library(shiny)
library(easyPubMed)
library(pubmed.mineR)
library(DT)

### Fixed variables ###

# Starting value for data range
# Ten years (in days) before current date
start_date <- Sys.Date() - (365.25 * 10)

### Custom functions ###

# Function for frequency barplots
freq_barplot <- function(varcat, varnum, main = ""){ # Categorical variable and numerical variable
  # Adjust width of left margin
  # https://stackoverflow.com/questions/10490763/automatic-adjustment-of-margins-in-horizontal-bar-chart
  par(mar=c(5.1, 
            max(4.1,max(nchar(as.character(varcat)))/1.5) ,
            4.1,
            2.1)
  )
  
  # The y object retrieves the coordinates of the categories
  # so they can be used for drawing text
  y <- barplot(varnum ~ varcat, 
               horiz = TRUE,
               las = 2,
               space = 0.1,
               main = main,
               ylab = "",
               xlab = "",
               xlim = c(0,max(varnum * 1.1)),
               axes = FALSE
  )
  # Writes the frequency of each gen at the end of the bar
  text(rev(varnum),
       y = y,
       labels = rev(varnum),
       adj = NULL,
       pos = 4,
       cex = 0.9
  )
}

# User interface
ui <- fluidPage(
  titlePanel("Endo-Mining",
             windowTitle = "Endo-Mining: minería de textos aplicada a la endometriosis"),
  navlistPanel(
    widths = c(2,10),
    tabPanel(title = "Buscar en PubMed",
  fluidRow(
    column(4,
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
  # textOutput("keyw"),
  # Search button
  actionButton("search", "Buscar en PubMed")
    ),
  column(8,
  textOutput("n_archivos"),
  # Cites as a table
  DT::dataTableOutput("titulos")
  )
  )
  ),
  tabPanel(title = "Frecuencia de palabras",
  fluidRow(
  # Table of words
    column(6,
           plotOutput("words_barplot"),
  tableOutput("palabras")
    ))),
  tabPanel(title = "Frecuencia de genes",
           fluidRow(
  column(6,
  # Barplot of genes
  plotOutput("genes_barplot"),
  # Table of genes
  tableOutput("genes_table")
    )
  )
  )
    )
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
    # Progress bar
    withProgress(message = "Descargando sumarios desde PubMed...",
                 detail = "Espere, por favor...",
                 value = 0, {
                   incProgress(7/15)
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
incProgress(15/15)
})
    abstracts
  })
  # Muestra la cantidad de citas recuperadas
  output$n_archivos <- renderText({
    paste0("Nº de citas recuperadas: ",
    length(pubmed_results()@PMID))
    })
  
  # Table of pmid plus title
  output$titulos <- DT::renderDataTable({
    corpus <- pubmed_results()
    # Table content
    tabla_titulos <- data.frame(corpus@PMID, corpus@Journal)
    colnames(tabla_titulos) <- c("PMID", "Publicaciones")
    datatable(tabla_titulos,
              selection = list(mode = 'single', selected = 1),
              options = list(language = list(url = 'spanish.json')))
    })
  
  
  ## Preprocesado del corpus primario
  # Word atomization
  words <- reactive({
    withProgress(message = "Recuperando palabras...",
                 value = 0, {
                   incProgress(1/2)
    words <- word_atomizations(pubmed_results())
    incProgress(2/2)
    words
                 })
    })
  # Table of words
  output$palabras <- renderTable({
    
    tabla_palabras <- data.frame(words()[1:10,])
    colnames(tabla_palabras) <- c("Palabra", "Frecuencia")
    tabla_palabras
  })
  
  # Barplot with frequency of words
  output$words_barplot <- renderPlot({
    tabla_frecuencias <- data.frame(words()[1:10,])
    tabla_frecuencias$words2 <- factor(tabla_frecuencias$words, 
                                 levels = rev(factor(tabla_frecuencias$words)))
    freq_barplot(varcat = tabla_frecuencias$words2,
                 varnum = tabla_frecuencias$Freq,
                 main = "Palabras más frecuentes")
  })
  
  # Gene atomization
  genes <- reactive({
    withProgress(message = 'Recuperando genes...',
                 detail = 'Suele tardar un rato...',
                 value = 0, {
                   incProgress(1/2)
    genes_data <- gene_atomization(pubmed_results())
    # Codify frequency of genes as numeric
    genes_table <- data.frame(genes_data,
                              stringsAsFactors = FALSE)
    colnames(genes_table) <- c("Symbol", "Nombre", "Frecuencia")
    genes_table$Frecuencia <- as.integer(genes_table$Frecuencia)
    incProgress(2/2)
                 })
    genes_table
    })

  # Table with frequency of genes
  output$genes_table <- renderTable({
    genes()
  })
  
  
  # Barplot with frequency of genes
  output$genes_barplot <- renderPlot({
    tabla_frecuencias <- genes()[1:10,]
    tabla_frecuencias$genes2 <- factor(tabla_frecuencias$Symbol,
                                       levels = rev(factor(tabla_frecuencias$Symbol)))
    freq_barplot(varcat = tabla_frecuencias$genes2,
                 varnum = tabla_frecuencias$Frecuencia,
                 main = "Genes más frecuentes")
  })


  
}
  


# Execution
shinyApp(ui, server)