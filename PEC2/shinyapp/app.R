library(shiny)
library(easyPubMed)
library(pubmed.mineR)
library(DT)
library(tokenizers)
library(org.Hs.eg.db) # GO over-representation test
library(clusterProfiler) # GO over-representation test

### Fixed variables ###

# Starting value for data range
# Five years (in days) before current date
start_date <- Sys.Date() - (5 * 365.25)

### Custom functions ###

# Function for frequency barplots
freq_barplot <- function(varcat, varnum, main = ""){ # Categorical variable and numerical variable
  # Adjust width of left margin
  # https://stackoverflow.com/questions/10490763/automatic-adjustment-of-margins-in-horizontal-bar-chart
  par(mar=c(5.1, 
            max(4.3,max(nchar(as.character(varcat)))/1.5) ,
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

  ## GO-over-representation test
  # Get human genes ID from pubtator contingency table
  universe_genes <- read.csv("human_geneID_universe.csv",
                             header = FALSE,
                             # Store as vector instead of dataframe
                             colClasses = "character")[,1]
  # Translate gene symbols to entrezId
  # Beware: the result is a dataframe
  # entrez <- select(org.Hs.eg.db,
  #                  keys = genes()$Symbol,
  #                  columns = c("SYMBOL", "ENTREZID"),
  #                  keytype = "SYMBOL") 
  ontology_aspect <- list("Función molecular" = "MF",
                          "Componente celular" = "CC",
                          "Proceso biológico" = "BP")
  adjust_methods <- list("Bonferroni" = "bonferroni",
                         "Holm" = "holm",
                         "Hommel" = "hommel",
                         "Benjamini & Hochberg" = "BH",
                         "Benjamini & Yekutieli" = "BY")
  # General GO overrepresentation function
  # ego_function <- enrichGO(gene = entrez()$ENTREZID,
  #                          universe = universe_genes(),
  #                          OrgDb = org.Hs.eg.db,
  #                          ont = ontology_aspect()$input$select_aspect,
  #                          pAdjustMethod =,
  #                          pvalueCutoff = 0.05,
  #                          qvalueCutoff = 0.5,
  #                          readable = FALSE)

### User interface ###
ui <- fluidPage(
  titlePanel("Endo-Mining",
             windowTitle = "Endo-Mining: minería de textos aplicada a la endometriosis"),
  navlistPanel(
    widths = c(2,10),
    tabPanel(title = "Buscar en PubMed",
             h1("Búsqueda en PubMed"),
  fluidRow(
    column(4,
  # Enter keywords
  textInput("keywords",
            label = "Palabras clave",
            value = "endometriosis",
            placeholder = "E.g., endometriosis"),
  # Enter date range
  dateRangeInput("fechas",
                 label = "Rango de fechas",
                 start = start_date,
                 format = "dd-mm-yyyy",
                 startview = "year",
                 weekstart = 1, # Monday
                 language = "es",
                 separator = "hasta"),
  # Do not select by date
  checkboxInput("check_all_dates",
                label = " Seleccionar máximo rango de fechas",
                value = FALSE),
  p(),
  p(strong("Texto consulta a PubMed")),
  # Query text
  verbatimTextOutput("keyw"),
  fluidRow(
    column(4,
           tabsetPanel(
             id = "SearchButton",
             type = "hidden",
             tabPanelBody(value = "button",
                          actionButton("search", "Buscar en PubMed")),
             tabPanelBody(value = "not_button",
                          "Búsqueda desactivada")
           )
    )
  )),
  column(8, # quiza deberia ser 6
  textOutput("n_archivos"),
  # Cites as a table
  DT::dataTableOutput("titulos"),
  # Abstract of selected cite
  htmlOutput("abstractText")
  )
  )),
  
  tabPanel(title = "Frecuencia de palabras",
           h1("Frecuencia de palabras"),
  fluidRow(
  # Table of words
    column(6,
           plotOutput("words_barplot"),
  tableOutput("palabras")
    ))),
  tabPanel(title = "Frecuencia de genes",
           h1("Frecuencia de genes"),
           fluidRow(
  column(6,
  # Barplot of genes
  plotOutput("genes_barplot"),
  # Table of genes
  tableOutput("genes_table")
    )
  )
  ),
  # Gráficas de frecuencia
  tabPanel(title = "Gráficas de frecuencia",
           h1("Placeholder")),
  # Caracterización de genes
  tabPanel(title = "Caracterización de genes",
           h1("Caracterización de genes por ontología génica"),
           column(4,
                  selectInput(inputId = "select_aspect",
                              "Aspecto funcional",
                              choices = c("Función molecular",
                                          "Componente celular",
                                          "Proceso biológico")),
                  numericInput(inputId = "go_categories",
                               "Categorías",
                               value = 20,
                               step = 1),
                  numericInput(inputId = "p_valor",
                               "Punto de corte: P-valor",
                               value = 0.01,
                               max = 1,
                               min = 0,
                               step = 0.005),
                  numericInput(inputId = "q_valor",
                               "Punto de corte: Q-valor",
                               value = 0.5,
                               max = 1,
                               min = 0,
                               step = 0.05),
                  selectInput(inputId = "metodo_ajuste",
                              "Método de ajuste del p-valor",
                              choices = c("Bonferroni",
                                          "Holm",
                                          "Hommel",
                                          "Benjamini & Hochberg",
                                          "Benjamini & Yekutieli"))
                  ),
           column(6,
                  renderPlot(barplot(height = ego,
                                     showCategory = 20,
                                     title = paste0("Términos GO enriquecidos \n(",
                                                    go_plot_titles[[ontology]],
                                                    ")")))
                  ))
    )
)




# App behaviour
server <- function(input, output, session){
  # Generates text for the query
  query <- reactive({
    validate(need(input$keywords != "", message = "POR FAVOR, INTRODUZCA LAS PALABRAS CLAVE DE SU INTERÉS" ),
             need(input$fechas[1] < input$fechas[2], message = "LA FECHA DE INICIO DEBE SER ANTERIOR A LA FECHA FINAL"))
    
    paste(c(input$keywords, " AND " , format(input$fechas[1],"%Y/%m/%d"),":",
            format(input$fechas[2],"%Y/%m/%d"),"[dp]"), collapse="")
  })
  
  # # Displays text of the query while being written
  output$keyw <-  renderText( query() )
  
  # Shows or hides search button
  # Hides when there are no keywords OR start date is bigger than finish date
  observe({
    if (input$keywords == "" || input$fechas[1] > input$fechas[2]) {
      updateTabsetPanel(inputId = "SearchButton",
                        selected = "not_button")
    }else{
      updateTabsetPanel(inputId = "SearchButton",
                        selected = "button")
    }
    }
  )
  
  # Updates date range when checkbox is ticked
  observe({
    if (input$check_all_dates == TRUE){
      updateDateRangeInput(inputId = "fechas",
                           start = "1800-01-01",
                           end = "3000-12-31")
    } else {
      updateDateRangeInput(inputId = "fechas",
                                 start = start_date,
                                 end = Sys.Date())}
               
  })
  
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
    # Display error message when input is wrong
    validate(need(input$SearchButton == "button", message = query() ))
    corpus <- pubmed_results()
    # Table content
    tabla_titulos <- data.frame(corpus@PMID, corpus@Journal)
    colnames(tabla_titulos) <- c("PMID", "Publicaciones")
    datatable(tabla_titulos,
              selection = list(mode = 'single', selected = 1),
              options = list(language = list(url = 'spanish.json')))
    })
  
  ## Abstract of selected pmid
  output$abstractText <- renderText({
    row_selected <- input$titulos_rows_selected
    abstracts <- pubmed_results()@Abstract[row_selected]
    abstractSentences <- tokenize_sentences(abstracts, simplify = TRUE)
    to_print <- paste('<p>', '<h4>', '<font_color = \"#4B04F6\"><b>', pubmed_results()@Journal[row_selected],
              '</b></font>', '</h4></p>', '\n')
    for (i in seq_along(abstractSentences)){
      if (i < 3) {
        to_print <- paste(to_print,
          '<p>', '<h4>', '<font_color = \"#4B04F6\"><b>', abstractSentences[i],
          '</b></font>', '</h4></p>', '\n')
      } else{
        to_print <- paste(to_print,
                          '<p><i>',abstractSentences[i],'</i></p>','\n')
      }
    }
    to_print <- paste(paste0('<p><a href="https://www.ncbi.nlm.nih.gov/pubmed/',pubmed_results()@PMID[row_selected],'" target=_blank>'
                            , 'Abrir publicación en PubMed', '</a></p>','\n'),
                      to_print)
    to_print
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
  
  ## GO-over-representation test
  ego_cc <- enrichGO(gene = entrez()$ENTREZID,
                           universe = universe_genes(),
                           OrgDb = org.Hs.eg.db,
                           ont = 'CC',
                           pAdjustMethod =,
                           pvalueCutoff = 0.05,
                           qvalueCutoff = 0.5,
                           readable = FALSE)
  
  output$ontologyCC <- renderPlot(height = ego_CC,
                                  showCategory = 20,
                                  title = "Título")
  
  # # Get human genes ID from pubtator contingency table
  # universe_genes <- read.csv("PEC2/shinyapp/human_geneID_universe.csv",
  #                            header = FALSE,
  #                            # Store as vector instead of dataframe
  #                            colClasses = "character")[,1]
  # # Translate gene symbols to entrezId
  # # Beware: the result is a dataframe
  # entrez <- select(org.Hs.eg.db,
  #                  keys = genes()$Symbol,
  #                  columns = c("SYMBOL", "ENTREZID"),
  #                  keytype = "SYMBOL") 


  
}
  


# Execution
shinyApp(ui, server)