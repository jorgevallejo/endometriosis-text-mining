library(shiny)
library(easyPubMed)
library(pubmed.mineR)

# Starting value for data range
# Five years (in days) before current date
start_date <- Sys.Date()-5

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
  textOutput("n_archivos"),
  # Cites as a table
  tableOutput("titulos")
)


# App behaviour
server <- function(input, output, session){
  # Generates text for the query
  query <- reactive(
    paste(c(input$keywords, " AND " , format(input$fechas[1],"%Y/%m/%d"),":",
            format(input$fechas[2],"%Y/%m/%d"),"[dp]"), collapse="")
  )
  
  # # Displays text of the query
  # output$keyw <- renderText(query())
  
  # Downloads search results
  
  pubmed_results <- eventReactive(input$search, {
    
    resultados_busqueda <- batch_pubmed_download(
      pubmed_query_string = query(),
      dest_file_prefix = "pubmed_",
      format = "abstract",
      batch_size = 5000)
    
    ## Concatenate files
    # List of files to be added together
    files_list <- list.files(pattern = "pubmed_",
                             full.names = TRUE) # include path
    # Create new file
    out_file <- file(description = "todos_resultados.txt",
                     open = "w")
    # Read each downloaded file and write into final file
    for (i in files_list){
      x <- readLines(i)
      writeLines(x, out_file)
    }

    close(out_file)
    # Generate object of class abstract
    abstracts <- readabs("todos_resultados.txt")
    
    # Delete unnecesary text files
    files_to_delete <- list.files(pattern = "\\.txt$")
    file.remove(files_to_delete)

    abstracts
  })
  # Muestra la cantidad de citas recuperadas
  output$n_archivos <- renderText({
    paste0("Nº de citas recuperadas: ",
    length(pubmed_results()@PMID))
    })
  
  # Table of pmid plus title
  output$titulos <- renderTable({
    corpus <- pubmed_results()
    if (length(corpus@PMID) <10) {
      citas <- length(corpus@PMID)
    } else {
      citas <- 10
    }
    message(citas)
    tabla_titulos <- data.frame(corpus@PMID[1:citas], corpus@Journal[1:citas])
    colnames(tabla_titulos) <- c("PMID", "Publicaciones")
    tabla_titulos
  })


  
}
  


# Execution
shinyApp(ui, server)