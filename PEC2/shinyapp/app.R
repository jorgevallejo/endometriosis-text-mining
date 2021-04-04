library(shiny)

# Starting value for data range
# Five years (in days) before current date
start_date <- Sys.Date()-(5*365.25)

# User interface
ui <- fluidPage(
  # Enter keywords
  textInput("keywords",
            label = "Palabras clave",
            value = "endometriosis",
            placeholder = "endometriosis"),
  # Enter date range
  dateRangeInput("fechas",
                 label = "Rango de fechas",
                 start = start_date,
                 format = "dd-mm-yyyy",
                 startview = "year",
                 weekstart = 1, # Monday
                 language = "es",
                 separator = "hasta"),
  textOutput("keyw")
)


# App behaviour
server <- function(input, output, session){
  # Generates the text for the query and displays
  output$keyw <- renderText({
    c(input$keywords, " AND ", format(input$fechas[1], "%Y/%m/%d"),
      ":",format(input$fechas[2], "%Y/%m/%d"), "[dp]")
    },
    sep = "")
  
  
}
  


# Execution
shinyApp(ui, server)