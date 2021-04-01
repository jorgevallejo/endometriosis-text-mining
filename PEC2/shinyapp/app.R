library(shiny)

animals <- c("dog", "cat", "mouse", "bird", "other", "I hate animals")

# User interface
ui <- fluidPage(
  "Collect free text",
  textInput("name", "What's your name?"),
  passwordInput("password", "What's your password?"),
  textAreaInput("story", "Tell me about yourself", rows = 3),
  "If you want to ensure that the text has certain properties you can use `validate()`",
  "\n",
  "\nCollect numeric inputs",
  numericInput("num", "Number one", value = 0, min = 0, max = 100),
  sliderInput("min", "Limit (minimum)", value = 50, min = 0, max = 100),
  sliderInput("rng", "Range", value = c(10, 20), min = 0, max = 100),
  "Collect dates",
  dateInput("dob", "When where you born?", language = "es", weekstart = 1),
  dateRangeInput("holiday", "When do you want to go on vacation next?", language = "es", weekstart = 1),
  "Collect limited choices",
  selectInput("state", "What's your favouirte state?", state.name
              #, multiple = TRUE # allows for selecting multiple options
              ),
  radioButtons("animal", "What's your favourite anima?", animals),
  radioButtons("rb", "Choose one:",
               choiceNames = list(
                 icon("angry"),
                 icon("smile"),
                 icon("sad-tear")
               ),
               choiceValues = list("angry", "happy", "sad")
               ),
  checkboxGroupInput("animal", "What animals do you like?", animals),
  checkboxInput("cleanup", "Clean up?", value = TRUE),
  checkboxInput("shutdown", "Shutdown?"),
  "File uploads",
  fileInput("upload", NULL),
  "Action buttons",
  actionButton("click", "Click me!"),
  actionButton("drink", "Drink me!", icon = icon("cocktail")),
  fluidRow(
    actionButton("click", "Click me!", class = "btn-danger"),
    actionButton("drink", "Drink me!", class = "btn-lg btn-success")
  ),
  fluidRow(
    actionButton("eat", "Eat me!", class = "btn-block")
  ),
  # OUTPUTS
  "Text",
  textOutput("text"),
  verbatimTextOutput("code"),
  "Tables",
  tableOutput("static"),
  dataTableOutput("dynamic"),
  "Plots",
  plotOutput("plot", width = "400px"),
  "Reactive programming example",
  textInput("name2", "What's your name?"),
  textOutput("greeting")
)

# App behaviour
server <- function(input, output, session){
  output$text <- renderText({
    "Hello friend!"
  })
  output$code <- renderPrint({
    summary(1:10)
  })
  output$static <- renderTable(head(mtcars))
  output$dynamic <- renderDataTable(mtcars, options = list(pageLength = 5))
  output$plot <- renderPlot(plot(1:5), res = 96)
  # Reactive programming example
  string <- reactive(paste0("Hello ", input$name2, "!"))
  output$greeting <- renderText(string())
  }
  


# Execution
shinyApp(ui, server)