# library(profvis) # Profiling the application

library(shiny)
library(easyPubMed)
library(pubmed.mineR)
library(DT)
library(tokenizers)
library(BiocManager) # Necessary for building org.Hs.eg.db into the app
options(repos = BiocManager::repositories()) # Necessary for building org.Hs.eg.db into the app
library(org.Hs.eg.db) # Obtaining EntrezID corresponding to gene symbol
library(enrichR) # GO over-representation test, interfaze for Enrichr webtool
library(wordcloud) # For wordcloud graphs

### Fixed variables ###

# Starting value for data range
# Five years (in days) before current date
end_date <- Sys.Date()
start_date <- end_date - (5 * 365.25)
# end_date <- "2021-05-20" # Temporal - only for test purposes
# start_date <- "2021-04-20" # Temporal - only for test purposes

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
               axes = FALSE,
               col = colorRampPalette(c("blue", "red"))(varcat)
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
  
  ontology_aspect <- list("Función molecular" = "GO_Molecular_Function_2018",  
                          "Componente celular" = "GO_Cellular_Component_2018",
                          "Proceso biológico" = "GO_Biological_Process_2018")
  
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
                 end = end_date,
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
           # h1("Frecuencia de palabras"),
           uiOutput("header_frecuencia_palabras"),
  fluidRow(
  # Table of words
    column(4,
           DT::dataTableOutput("palabras")
    ),
    column(6,
           DT::dataTableOutput("palabras_2ario")
           )),
  fluidRow(
    column(4,
    # Hyperlink to selected publication
    htmlOutput("HyperlinkPalabra")
    ),
    column(6,
           # Abstract of selected publication
           htmlOutput("abstractPalabra")
           )
  )),
  tabPanel(title = "Frecuencia de genes",
           h1("Frecuencia de genes"),
           fluidRow(
  column(5,
  # Table of genes
  DT::dataTableOutput("genes_table")
    ),
  column(5,
         # Secondary corpus of genes
         DT::dataTableOutput("genes_cites_table"))
  ),
          fluidRow(
            column(5,
                   htmlOutput("hyperlink_gene")),
            column(5,
                   # Abstract of selected publication by gene
                   htmlOutput("abstractGene"))
          )
  ),
  # Gráficas de frecuencia
  tabPanel(title = "Gráficas de frecuencia",
           h1("Gráficas de frecuencia"),
           fluidRow(
             column(5,
                    selectInput(inputId = "select_words_genes",
                                "Seleccionar resultados de",
                                choices = c("Palabras más frecuentes",
                                            "Genes más frecuentes")),
                    selectInput(inputId = "barplot_cloud",
                                "Típo de gráfica",
                                choices = c("Gráfico de barras",
                                            "Nube de palabras")),
                    # Optional UI with tabsets
                    # Will display results controls for barplot or wordcloud
                    # according to previous selection
                    ## What is achieved at the moment is just changing the default number
                    ## of categories/words displayed depending on barplot or wordcloud.
                    ## That could have been arranged in a more simple way using updateSliderInput,
                    ## but I have used a tabsetPanel because originally there where going
                    ## to be more and different controls for each kind of graph. I have not
                    ## had time to implement those, though.
                    tabsetPanel(
                      id = 'controles_barplot_wordcloud',
                      type = 'hidden',
                      tabPanelBody(
                        "Gráfico de barras",
                        sliderInput(inputId = 'genes_words_max',
                                    label = "Categorías en el gráfico",
                                    min = 1,
                                    max = 100,
                                    value = 20,
                                    step = 1)),
                      tabPanelBody(
                        "Nube de palabras",
                        sliderInput(inputId = 'words_cloud_max',
                                    label = "Límite de palabras",
                                    min = 1,
                                    max = 1000,
                                    value = 100,
                                    step = 1))
                        )),
             column(5,
                    # Optional UI with tabsets
                    # Will display results for words or genes
                    # according to user selection.
                    tabsetPanel(
                      id = 'graficas_frecuencia',
                      type = 'hidden',
                      tabPanelBody(
                        "Palabras más frecuentes - Gráfico de barras",
                        # Barplot de palabras
                        plotOutput("words_barplot")
                      ),
                      tabPanelBody(
                        "Genes más frecuentes - Gráfico de barras",
                        # Barplot of genes
                        plotOutput("genes_barplot")),
                      tabPanelBody(
                        "Palabras más frecuentes - Nube de palabras",
                        # Nube de palabras
                        plotOutput("words_wordcloud")),
                      tabPanelBody(
                        "Genes más frecuentes - Nube de palabras",
                        # Nube de genes
                        plotOutput("genes_wordcloud"))
                      )))),
  # Caracterización de genes
  tabPanel(title = "Caracterización de genes",
           h1("Caracterización de genes por ontología génica"),
           column(4,
                  selectInput(inputId = "select_display",
                              "Mostrar resultados como",
                              choices = c("Tabla",
                                          "Gráfico de barras")),
                  numericInput(inputId = "go_categories",
                               "Categorías mostradas",
                               value = 10,
                               step = 1),
                  selectInput(inputId = "select_aspect",
                              "Aspecto funcional",
                              choices = c("Componente celular",
                                          "Proceso biológico",
                                          "Función molecular")),
                  numericInput(inputId = "p_valor",
                               "Nivel de significatividad (p-valor ajustado)",
                               value = 0.05,
                               max = 1,
                               min = 0,
                               step = 0.005),
                  actionButton(inputId = "GO_button",
                               label = "Caracterizar"),
                  p(),
                  p("El método de ajuste del p-valor en esta aplicación 
                    (necesario para controlar la probabilidad de falsos positivos 
                    en las comparaciones múltiples) es el conocido como de", 
                    span("Benjamini & Hochberg", style = "font-style:italic"), "."), 
                    p("El botón de descarga proporciona un archivo CSV con todos los 
                    términos GO recuperados, junto con los p-valores originales;
                    permitiendo al usuario calcular los p-valores por su cuenta
                      si considera necesario usar un método diferente.")
                  ),
           column(6,
                  # Optional UI with tabsets
                  # Will display results as table or as barplot
                  # according to user selection.
                  tabsetPanel(
                    id = 'tabla_grafico',
                    type = 'hidden',
                    tabPanelBody(
                      "Tabla",
                      # GO terms as a table
                      DT::dataTableOutput("GOterms"),
                      # Download button
                      uiOutput("GO_download_ui"),
                      # Hyperlink to AmiGO website
                      htmlOutput("GO_link")
                    ),
                    tabPanelBody("Gráfico de barras",
                                 plotOutput(outputId = "GO_barplot"
                                            )
                                 
                                 )
                  )
                  
                  
                  
                  
                  )),
  tabPanel(title = "Acerca de",
           h1("Acerca de Endo-Mining"),
           p("Versión 0.6.0 (", a(href="version.txt"," visita el registro de cambios"), ")"),
           p("Endo-Mining es una aplicación web diseñada para llevar a cabo análisis exploratorios
             rápidos y ligeros de la información genética contenida en los sumarios de publicaciones
             biomédicas almacenados en la base de datos PubMed."),
           p("Diseñado por Jorge Vallejo Ortega como parte del Trabajo de Fin de Máster en el máster
             Bioinformática y Bioestadística de la Universitat Oberta de Catalunya."),
           p("Repositorio del proyecto en GitHub."),
           p(),
           p("Consultor: Romina Nebrij"),
           p("Responsable de área: Antoni Pérez Navarro")
    )
))




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
                                 end = end_date)}
               
  })
  
  # Downloads search results ## Temporal-comentado para tests en local

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
                            , 'Visitar página de la cita en PubMed', '</a></p>','\n'),
                      to_print)
    to_print
   })
  
  ## Preprocesado del corpus primario
  # Word atomization # Comentario temporal para tests
  words <- reactive({
    withProgress(message = "Recuperando palabras...",
                 value = 0, {
                   incProgress(1/2)
    words <- word_atomizations(pubmed_results())
    incProgress(2/2)
    words
                 })
    })
  
  ## Temporal for words in local
  # words <- reactive(readRDS("test_files/words.RDS"))
  ## Temporal for pubmed results in local
  # pubmed_results <- reactive(readRDS("test_files/pubmed_results_temporal.RDS"))
  
  # Header for frequency of words section
  output$header_frecuencia_palabras <- renderUI({
    query_keywords <- input$keywords
    if (is.null(query_keywords)) {h1("Frecuencia de palabras")}
    else {
      h1(
        HTML(
          paste0(
            "Frecuencia de palabras en",
            tags$br(),
            "publicaciones sobre ", query_keywords)))}
  })
  
  # Table of words
  output$palabras <- DT::renderDataTable({
    # Table content
    tabla_palabras <- data.frame(words())
    # tabla_palabras <- words() # Temporal mientras pruebo en local
    datatable(tabla_palabras,
              colnames = c("Palabra", "Frecuencia"),
              rownames = FALSE,
              caption = 'Haga click en las cabeceras de las columnas para cambiar el orden',
              selection = list(mode = 'single', selected = 1),
              options = list(language = list(url = 'spanish.json')))
  })
  
  # Secondary corpus based on selected word
  corpus_2ario <- reactive({
    withProgress(message = "Generando corpus secundario...",
                                 value = 0, {
    corpus <- pubmed_results()
    # corpus <- pubmed_results_temporal() # Temporal for testing in local
    setProgress(1/4)
    word_selected <- input$palabras_rows_selected
    setProgress(2/4)
    term <- words()[word_selected, 1] # Recover selected word from words dataframe
    setProgress(3/4)
    getabs(corpus, term, FALSE)
  })
  })
  

  
  
  # Table for secondary corpus on words
  output$palabras_2ario <- DT::renderDataTable({
    # Table content
    tabla_titulos_2ario <- data.frame(corpus_2ario()@PMID,
                                      corpus_2ario()@Journal)
    datatable(tabla_titulos_2ario,
              colnames = c("PMID", "Publicación"),
              rownames = FALSE,
              caption = "Citas que contienen la palabra seleccionada",
              selection = list(mode = 'single', selected = 1),
              options = list(language = list(url = 'spanish.json')))
  })
  
  ## Abstract of selected pmid for words
  ### This should be re-factored into a function because I am using
  ### the same code that in output$abstractText and output$abstractGene
  output$abstractPalabra <- renderText({
    row_selected <- input$palabras_2ario_rows_selected
    abstracts <- corpus_2ario()@Abstract[row_selected]
    abstractSentences <- tokenize_sentences(abstracts, simplify = TRUE)
    to_print <- paste('<p>', '<h4>', '<font_color = \"#4B04F6\"><b>', corpus_2ario()@Journal[row_selected],
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
    to_print <- paste(paste0('<p><a href="https://www.ncbi.nlm.nih.gov/pubmed/',corpus_2ario()@PMID[row_selected],'" target=_blank>'
                             , 'Visitar página de la cita en PubMed', '</a></p>','\n'),
                      to_print)
    to_print
  })
  
  # Display words or genes barplot
  observeEvent({input$select_words_genes
                input$barplot_cloud}, {
    updateTabsetPanel(
      inputId = "graficas_frecuencia",
      selected = paste0(input$select_words_genes,
                        " - ",
                        input$barplot_cloud))
  })
  
  # Display barplot or wordcloud controls
  observeEvent(input$barplot_cloud, {
      updateTabsetPanel(
        inputId = "controles_barplot_wordcloud",
        selected = input$barplot_cloud)
    })
  
  # Update slider of min and max represented words/genes
  observeEvent(input$select_words_genes,
               updateSliderInput(
                 inputId = 'genes_words_max',
                 max = min(100, nrow(genes()))
               ))
  
  
  # Barplot with frequency of words
  output$words_barplot <- renderPlot({
    tabla_frecuencias <- data.frame(words()[1:input$genes_words_max,])
    tabla_frecuencias$words2 <- factor(tabla_frecuencias$words, 
                                 levels = rev(factor(tabla_frecuencias$words)))
    freq_barplot(varcat = tabla_frecuencias$words2,
                 varnum = tabla_frecuencias$Freq,
                 main = "Palabras más frecuentes")
  },
  height = reactive(max(600, input$genes_words_max * 20)),
  res = 96,
  alt = 'Gráfica de barras de palabras más frecuentes')

  # Wordcloud with frequency of words
  output$words_wordcloud <- renderPlot({
    tabla_frecuencias <- data.frame(words()[1:input$words_cloud_max,])
    tabla_frecuencias$words2 <- factor(tabla_frecuencias$words, 
                                       levels = rev(factor(tabla_frecuencias$words)))
    wordcloud::wordcloud(words = tabla_frecuencias$words2,
                         freq = tabla_frecuencias$Freq,
                         random.order = FALSE,
                         colors = (colorRampPalette(c("blue", "red"))(100))) # Provisional
  },
  height = 600,
  res = 96,
  alt = 'Gráfica de barras de palabras más frecuentes')
  
  # Genes temporal
  # genes <- reactive({genes_data <- readRDS("test_files/genes.RDS")
  #                   genes_table <- data.frame(genes_data,
  #                                                                         stringsAsFactors = FALSE)
  #                                               colnames(genes_table) <- c("Símbolo", "Nombre", "Frecuencia")
  #                                               genes_table$Frecuencia <- as.integer(genes_table$Frecuencia)
  #                                               genes_table
  #                   })
    
  # Gene atomization ## Temporal - comentado para tests en local
  genes <- reactive({
    withProgress(message = 'Recuperando genes...',
                 detail = 'Suele tardar un rato...',
                 value = 0, {
                   incProgress(1/2)
    genes_data <- gene_atomization(pubmed_results())
    # Codify frequency of genes as numeric
    genes_table <- data.frame(genes_data,
                              stringsAsFactors = FALSE)
    colnames(genes_table) <- c("Símbolo", "Nombre", "Frecuencia")
    genes_table$Frecuencia <- as.integer(genes_table$Frecuencia)
    incProgress(2/2)
                 })
    genes_table
    })

  # Add EntrezID column into genes table
  genes_plus_entrez <- reactive({
    genes_table <- genes()
    keys <- genes_table[, "Símbolo"] # Char vector for looking up in database
    entrez <- mapIds(org.Hs.eg.db, # vector with correspondence symbol-entrezid
                     keys = keys,
                     column = "ENTREZID",
                     keytype = "SYMBOL",
                     multiVals = 'first'
    )
    genes_table$Entrez_ID <- entrez  # Add new column to genes dataframe
    genes_table <- genes_table[, c("Símbolo", "Entrez_ID", "Nombre", "Frecuencia")] # Rearrange columns
  })
  
  # Table with frequency of genes
  output$genes_table <- renderDataTable({
    req(genes_plus_entrez())
    datatable(genes_plus_entrez(),
              rownames = FALSE,
              caption = 'Haga click en las cabeceras de las columnas para cambiar el orden',
              selection = list(mode = 'single', selected = 1),
              options = list(language = list(url = 'spanish.json')))
  })
  
  # Secondary corpus based on selected gene symbol
  corpus_2ario_gene <- reactive({
    withProgress(message = "Generando corpus secundario...",
                 value = 0, {
                   corpus <- pubmed_results()
                   # corpus <- pubmed_results_temporal() # Temporal for testing in local
                   setProgress(1/4)
                   gene_selected <- input$genes_table_rows_selected
                   setProgress(2/4)
                   term <- genes()[gene_selected, 1] # Recover selected word from words dataframe
                   setProgress(3/4)
                   getabs(corpus, term, FALSE)
                 })
  })
  
  # Table with citations that include the gen
  output$genes_cites_table <- DT::renderDataTable({
    tabla_genes_2ario <- data.frame(corpus_2ario_gene()@PMID,
                                    corpus_2ario_gene()@Journal)
    datatable(tabla_genes_2ario,
              rownames = FALSE,
              colnames = c("PMID", "Publicación"),
              caption = 'Publicaciones que contienen el gen seleccionado',
              selection = list(mode = 'single', selected = 1),
              options = list(language = list(url = 'spanish.json')))
  })
  
  ## Abstract of selected pmid for gene
  ### This should be re-factored into a function because I am using
  ### the same code that in output$abstractText and output$abstractPalabra
  output$abstractGene <- renderText({
    row_selected <- input$genes_cites_table_rows_selected
    abstracts <- corpus_2ario_gene()@Abstract[row_selected]
    abstractSentences <- tokenize_sentences(abstracts, simplify = TRUE)
    to_print <- paste('<p>', '<h4>', '<font_color = \"#4B04F6\"><b>', corpus_2ario_gene()@Journal[row_selected],
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
    to_print <- paste(paste0('<p><a href="https://www.ncbi.nlm.nih.gov/pubmed/',corpus_2ario_gene()@PMID[row_selected],'" target=_blank>'
                             , 'Visitar página de la cita en PubMed', '</a></p>','\n'),
                      to_print)
    to_print
  })
  
  # Hyperlink for Entrez ID
  output$hyperlink_gene <- renderText({
    req(genes())
    row_selected <- input$genes_table_rows_selected
    # isolate Entrez ID for composing hyperlink
    gene_id <- genes_plus_entrez()[row_selected, c("Símbolo", "Entrez_ID")]
    # Build hyperlink
    paste0('<br /><br /><p><a href="https://www.ncbi.nlm.nih.gov/gene/', gene_id[1,2],'" target=_blank>',
           'Abrir enlace a la página de información del gen ', gene_id[1,1],
           ' en NCBI Gene', '</a></p>','\n')
  })
  
  # Barplot with frequency of genes
  output$genes_barplot <- renderPlot({
    tabla_frecuencias <- genes()[1:input$genes_words_max,]
    tabla_frecuencias$genes2 <- factor(tabla_frecuencias$Símbolo,
                                       levels = rev(factor(tabla_frecuencias$Símbolo)))
    freq_barplot(varcat = tabla_frecuencias$genes2,
                 varnum = tabla_frecuencias$Frecuencia,
                 main = "Genes más frecuentes")
  },
  height = reactive(max(600, input$genes_words_max * 20)),
  res = 96,
  alt = 'Gráfica de barras de genes más frecuentes')
  
  # Wordcloud with frequency of genes
  output$genes_wordcloud <- renderPlot({
    tabla_frecuencias <- genes()[1:input$words_cloud_max,]
    tabla_frecuencias$genes2 <- factor(tabla_frecuencias$Símbolo,
                                       levels = rev(factor(tabla_frecuencias$Símbolo)))
    wordcloud::wordcloud(words = tabla_frecuencias$genes2,
                         freq = tabla_frecuencias$Frecuencia,
                         random.order = FALSE,
                         colors = (colorRampPalette(c("blue", "red"))(100)) # Provisional
                         )
  },
  height = 600,
  res = 96,
  alt = 'Nube de palabras de los genes más frecuentes')
  
  # Display results of GO enrichment as table or barplot
  observeEvent(input$select_display, {
    updateTabsetPanel(
      inputId = "tabla_grafico",
      selected = input$select_display)
  })

  

# Compute enrichment of terms in gene set using enrichR
# as an interface for the web tool Enrichr
  ego_terms <- eventReactive(input$GO_button,{
    withProgress(message = "Calculando términos GO enriquecidos", {
    databases <- c("GO_Molecular_Function_2018", "GO_Cellular_Component_2018", "GO_Biological_Process_2018")
    genes_candidatos <- genes()[, "Símbolo"]
    incProgress(2/5)
    enriched <- enrichr(genes_candidatos,
                        databases = databases)
    })
  })
   
  # Ontology aspect that the user wants to explore
  ontology <- reactive(ontology_aspect[[input$select_aspect]])
  
  # GO terms dataframe
  go_dataframe <- reactive({
    # Create dataframe per ontology aspect
    dataframe <- ego_terms()[[ontology()]]
    # Adjusted P-value Cutoff
    dataframe[dataframe[,"Adjusted.P.value"] <= input$p_valor, ]
  })
   
  ### GO terms table
  output$GOterms <- DT::renderDataTable({
      datatable(go_dataframe()[, c("Term", "Adjusted.P.value", "Combined.Score", "Overlap")],
        #ego_terms()[[ontology()]][, c("Term", "Adjusted.P.value", "Combined.Score", "Overlap")],
                rownames = FALSE,
                colnames = c("Término GO", "p-valor ajustado", "Puntuación combinada", "Genes coincidentes"),
                selection = list(mode = 'single', selected = 1),
                options = list(language = list(url = 'spanish.json'),
                               # Number of rows in each page are determined by user 
                               pageLength = min(nrow(go_dataframe()),
                                 #nrow(ego_terms()[[ontology()]]),
                                                input$go_categories))) %>%
        formatSignif('Adjusted.P.value', 2) %>%  # Significative digits 
        formatRound('Combined.Score', 0) %>%  # Round Score to units
        formatStyle(columns = c("Adjusted.P.value", "Overlap"), `text-align` = 'left') # Center columns
    })
  
  # Hyperlink for GO term
  output$GO_link <- renderText({
    req(ego_terms())
    row_selected <- input$GOterms_rows_selected
    # isolate GO ID for composing hyperlink
    ego_term <- regmatches(ego_terms()[[ontology()]][row_selected,"Term"],  # Term selected in the table
                           regexec(pattern = 'GO:([[:digit:]]+)',  # Look for a substring of digits after 'GO:'
                                   ego_terms()[[ontology()]][row_selected,"Term"]))[[1]][1]  # Select the second term of the results vector
      
    # Build hyperlink
    paste0('<br /><br /><p><a href="http://amigo.geneontology.org/amigo/term/', ego_term,'" target=_blank>',
           'Abrir enlace a la página de información del término ', ego_terms()[[ontology()]][row_selected,"Term"],
           ' en AmiGO', '</a></p>','\n')
  })
  
  ## Prepare GO data for download
  
  output$GO_download_ui <- renderUI({
    req(ego_terms())
    downloadButton("GO_download",
                   label = "Descargar como archivo .csv")
  })
    
    output$GO_download <- downloadHandler(
    filename = function() {
      paste0(query(),'enrichedGOterms_',input$select_aspect,'.csv')
    },
    content = function(file) {
      # The table to download will not be cut off by p-value
      # so that the user will have access to all the info
      download_table <- ego_terms()[[ontology()]]
      # Subset columns (those that would make sense for the user)
      download_table <- download_table[, c("Term", "P.value", "Adjusted.P.value",
                                           "Combined.Score", "Overlap", "Genes")]
      # Change colnames to coincide with the ones in the web app
      colnames(download_table) <- c("Término GO", "p-valor", "p-valor ajustado",
                                    "Puntuación combinada", "Genes coincidentes", "Genes")
      # Generate the csv file
      write.csv(download_table,
                file = file,
                row.names = FALSE)
    },
    contentType = "text/csv"
  )
   
  # Plot enriched GO terms
  # y-axis is number of genes in each term
  # Order is by p-value
  output$GO_barplot <- renderPlot(
    plotEnrich(df = go_dataframe(), 
                 #ego_terms()[[ontology()]],
               # Number of bars in the plot is the minimum between actual
               # number of rows in the table or the number inputed by user
               showTerms = min(nrow(go_dataframe()),
                 #nrow(ego_terms()[[ontology()]]),
                               input$go_categories),
               numChar = 40, # Characters in x-axis labels
               xlab = paste0('Términos GO (',
                             min(nrow(go_dataframe()),input$go_categories),
                             ' de ', 
                             nrow(go_dataframe()), ' significativos)'),
               ylab = "Número de genes en la categoría",
               title = paste0("Términos GO enriquecidos \n(", input$select_aspect,")")),
    height = reactive(max(600, input$go_categories * 20)),
    res = 96,
    alt = 'Gráfica de barras de términos GO enriquecidos'
  )
}
  


# Execution
# profvis::profvis(runApp(shinyApp(ui, server)))

shinyApp(ui, server)

