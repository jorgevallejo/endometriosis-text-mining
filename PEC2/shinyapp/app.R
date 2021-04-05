library(shiny)
library(easyPubMed)
library(pubmed.mineR)

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
  textOutput("keyw"),
  # Search button
  actionButton("search", "Buscar en PubMed"),
  textOutput("n_archivos")
)


# App behaviour
server <- function(input, output, session){
  # # Generates the text for the query
  # query <- reactive(
  #   paste(c(input$keywords, " AND ", format(input$fechas[1], "%Y/%m/%d"),
  #         ":",format(input$fechas[2], "%Y/%m/%d"), "[dp]"), sep = "")
  # )
  
  # # Displays text of the query
  # output$keyw <- renderText(query())
  
  # Downloads search results
  # Generates the query text only at the press of search button
  # This duplicates the code in reactive expression "query"
    # there must be a way to simplify this
  pubmed_results <- eventReactive(input$search, {
    
    resultados_busqueda <- batch_pubmed_download(
      pubmed_query_string = "endometriosis AND 2021/03/31:2021/04/04[dp]",
      dest_file_prefix = "pubmed_",
      format = "abstract",
      batch_size = 5000)
    
    
    
    # ## Concatenate files
    # # List of files to be added together
    # files_list <- list.files(pattern = "pubmed_",
    #                          full.names = TRUE) # include path
    # # Create new file
    # out_file <- file(description = "todos_resultados.txt",
    #                  open = "w")
    # # Read each downloaded file and write into final file
    # for (i in files_list){
    #   x <- readLines(i)
    #   writeLines(x, out_file)
    # }
    # 
    # close(out_file)
    # # Generate object of class abstract
    # abstracts <- readabs("todos_resultados.txt")
    # 
    # return(abstracts)
  })
  # Muestra la cantidad de archivos descargados (no el numero de entradas en el archivo)
  output$n_archivos <- renderText(length(pubmed_results()))


  
}
  


# Execution
shinyApp(ui, server)