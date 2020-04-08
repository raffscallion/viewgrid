library(shiny)
library(leaflet)
library(ncdf4)
library(raster)
library(tidyverse)
library(leaflet.extras)
library(plotly)
library(fs)
library(shinyWidgets)

data_path <- "C:/Users/sraffuse/Google Drive/Working/JVA/GOES-R Fire/R/HAQAST/data/CMAQ/fullgrid"

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("CMAQ View"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            fluidRow(
                column(6, selectInput("date", label = "Date",
                                      choices = seq.Date(from = as.Date("2017-10-08"),
                                                         to = as.Date("2017-10-16"),
                                                         by = "1 day"))
                       ),
                column(6, radioButtons("model", "Model Run:", inline = TRUE,
                                       choices = c("Baseline", "S1", "S2"))
                       )
            ),
            sliderInput("hour", "Hour:", min = 1, max = 24, value = 1),
            noUiSliderInput("height", "Height Layer:", min = 1, max = 28, value = 1,
                            orientation = "vertical", step = 1, height = "150px",
                            direction = "rtl"),
            numericInput("vlimit", "Vertical Limit", min = 500, max = 20000, value = 2000)
        ),

        # Show a plot of the generated distribution
        mainPanel(
            fluidRow(
                column(6, leafletOutput("map")),
                column(6, leafletOutput("integrated"))
            ),
            plotlyOutput("curtain")
        )
    )
)

server <- function(input, output) {

    g <- reactiveVal(NULL)
    proj_string <- reactiveVal(NULL)
    max_val <- reactiveVal(NULL)
    heights <- reactiveVal(NULL)
    
    get_nc <- reactive({

        date_string <- strftime(input$date, format = "%Y%m%d")
        filename <- paste(input$model, "out.combine", date_string, sep = "_")
        fullpath <- path(data_path, filename)
        nc <- nc_open(fullpath)
        
        g <- ncatt_get(nc, 0)
        proj_string <- glue::glue("+proj=lcc +lon_0={g$XCENT} +lat_1={g$P_ALP} ",
                                  "+lat_2={g$P_BET} +lat_0={g$YCENT} +units=m ",
                                  "+a=6370000.0 +b=6370000.0")
        
        # Get height from vertical levels in the global using the barometric formula
        # (https://en.wikipedia.org/wiki/Atmospheric_pressure)
        pressures <- g$VGLVLS
        heights <- log(1/pressures) / 0.000119

        # Set max height to 20000 (instead of Inf)
        heights[length(heights)] <- 20000
        height_table <- tibble(Layer = seq(1, length(heights), by = 1),
                                  Height = heights)

        # Store the global attributes, heights, and the proj string as global reactives
        g(g)
        proj_string(proj_string)
        heights(height_table)
        
        pm <- ncvar_get(nc, "ATOTIJ")

        # Also store the max value for use in color scales
        max_val(max(pm))
        
        # Permutate dimensions so that x and y are correct (will also need to flip after
        # converting to raster)
        pm <- aperm(pm, c(2, 1, 3, 4))

    })
    
    
    hourly_brick <- reactive({

        pm <- get_nc()
        stack <- pm[,,,input$hour]
        pmb <- brick(stack, xmn = g()$XORIG, ymn = g()$YORIG,
                            xmx = (g()$XORIG + g()$XCELL * g()$NROWS),
                            ymx = (g()$YORIG + g()$YCELL * g()$NCOLS),
                            crs = CRS(proj_string()))
        pmb <- flip(pmb, "y")
        
    })
    
    map_layer <- reactive({
        
        pmb <- hourly_brick()
        pmr <- pmb[[input$height]]
        
        # Squish and transform values before plotting so the scale looks just like
        # the curtain plot
        pmr[pmr < 1] <- 1
        pmr[pmr > 10000] <- 10000
        pmr <- log10(pmr)
        
    })
    
    integrated_layers <- reactive({
        
        # This is the sum of all the layers - the view from the satellite
        pmb <- hourly_brick()
        pmr <- sum(pmb)
        pmr[pmr < 1] <- 1
        pmr[pmr > 10000] <- 10000
        pmr <- log10(pmr)
        
    })
    
    output$map <- renderLeaflet({
        leaflet() %>%
            addProviderTiles(providers$Stamen.TonerLite,
                             providerTileOptions(noWrap = TRUE)) %>%
            addDrawToolbar(targetGroup = "draw",
                           polygonOptions = FALSE, circleOptions = FALSE,
                           rectangleOptions = FALSE, markerOptions = FALSE,
                           circleMarkerOptions = FALSE,
                           singleFeature = TRUE,
                           polylineOptions = drawPolylineOptions(
                               shapeOptions = 
                                   drawShapeOptions(color = "#000", weight = 2,
                                                    dashArray = "2 4"))) %>%
            setView(-120, 38, zoom = 6)
    })
    
    observe({
        leafletProxy("map") %>%
            clearImages() %>%
            addRasterImage(map_layer(), opacity = 0.8, colors = "viridis")
    })
    
    output$integrated <- renderLeaflet({
        leaflet() %>%
            addProviderTiles(providers$Stamen.TonerLite,
                             providerTileOptions(noWrap = TRUE)) %>%
            setView(-120, 38, zoom = 6)
    })
    
    observe({
        leafletProxy("integrated") %>%
            clearImages() %>%
            addRasterImage(integrated_layers(), opacity = 0.8, colors = "viridis")
    })
    
    curtain_line <- reactiveVal(NULL)
    
    observeEvent(input$map_draw_new_feature, {
        
        # Grab the line and make it spatial
        feature <- input$map_draw_new_feature
        coords <- unlist(feature$geometry$coordinates)
        mx <- matrix(coords, ncol = 2, byrow = TRUE)
        line <- spLines(mx, crs = '+proj=longlat +datum=WGS84')
        line <- spTransform(line, CRS(proj_string()))
        curtain_line(line)

    })
    
    curtain_data <- reactive({
        
        validate(need(!is.null(curtain_line()), "Draw a line"))

        # Extract the values from the brick along the line
        ex <- raster::extract(hourly_brick(), curtain_line(), nl = 28, along = TRUE,
                              cellnumbers = TRUE)[[1]]
        
        # Convert to a data frame for plotting as a curtain
        df <- as_tibble(ex)
        df$Order <- seq(1, nrow(df))
        
        df <- df %>%
            pivot_longer(starts_with("layer"), names_to = "Layer", values_to = "PM25") %>%
            separate(Layer, c("label", "y")) %>%
            mutate(Layer = as.numeric(y)) %>%
            select(-label, -y) %>%
            left_join(heights(), by = "Layer")

    })
    
    output$curtain <- renderPlotly({
        
        validate(need(!is.null(curtain_data()), "Draw a line"))
        
        height <- heights()$Height[heights()$Layer == input$height]
        
        gg <- ggplot(curtain_data(), aes(x = Order, y = Height, fill = PM25)) +
            geom_raster() +
            geom_hline(yintercept = height, linetype = "dotted") +
            scale_fill_viridis_c(limits =  c(1, 10000), trans = "log10",
                                  oob = scales::squish) +
            scale_y_continuous(limits = c(0, input$vlimit)) +
            theme_bw() +
            labs(x = "Transect",
                 y = "Height (m)",
                 caption = "Beginning of drawn line is on the left")
        ggplotly(gg, dynamicTicks = TRUE)
        
    })
    
}

# Run the application 
shinyApp(ui = ui, server = server)
