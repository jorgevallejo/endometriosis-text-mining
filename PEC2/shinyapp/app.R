library(shiny)


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
                 startview = "year",
                 weekstart = 1, # Monday
                 language = "es",
                 separator = "hasta"),
  textOutput("keyw")
)


# App behaviour
server <- function(input, output, session){
  output$keyw <- renderText({
    c(input$keywords, " desde ", format(input$fechas[1], "%d %b %Y"),
      " hasta ", format(input$fechas[2], "%d %b %Y"))
    })
}
  


# Execution
shinyApp(ui, server)