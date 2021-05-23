# library(profvis)

library(shiny)
library(easyPubMed)
library(pubmed.mineR)
library(DT)
library(tokenizers)
# library(BiocManager) # Necessary for building clusterProfiler into the app
# options(repos = BiocManager::repositories()) # Necessary for building clusterProfiler into the app
# library(org.Hs.eg.db) # GO over-representation test
# library(clusterProfiler) # GO over-representation test
# library(ggplot2) # For putting xlabel in GO enrichment barplot
library(enrichR) # GO over-representation test, interfaze for Enrichr webtool

### Fixed variables ###

# Starting value for data range
# Five years (in days) before current date
# end_date <- Sys.Date()
# start_date <- end_date - (5 * 365.25)
end_date <- "2021-05-20" # Temporal - only for test purposes
start_date <- "2021-04-20" # Temporal - only for test purposes

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
  
  ontology_aspect <- list("Función molecular" = "GO_Molecular_Function_2018",  
                          "Componente celular" = "GO_Cellular_Component_2018",
                          "Proceso biológico" = "GO_Biological_Process_2018")
  
  # adjust_methods <- list("Bonferroni" = "bonferroni",
  #                        "Holm" = "holm",
  #                        "Hommel" = "hommel",
  #                        "Benjamini & Hochberg" = "BH",
  #                        "Benjamini & Yekutieli" = "BY")
  
  # # General GO overrepresentation function
  # ego_function <- function(genes, ontology, padjust, pvalue, qvalue) {
  #   enrichGO(gene = genes,
  #                          universe = universe_genes,
  #                          OrgDb = org.Hs.eg.db,
  #                          ont = ontology,
  #                          pAdjustMethod = padjust,
  #                          pvalueCutoff = pvalue,
  #                          qvalueCutoff = qvalue,
  #                          readable = FALSE)
  # }

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
                               "Punto de corte: P-valor",
                               value = 0.05,
                               max = 1,
                               min = 0,
                               step = 0.005),
                  # numericInput(inputId = "q_valor",
                  #              "Punto de corte: Q-valor",
                  #              value = 0.05,
                  #              max = 1,
                  #              min = 0,
                  #              step = 0.05),
                  selectInput(inputId = "metodo_ajuste",
                              "Método de ajuste del p-valor",
                              choices = "Benjamini & Hochberg"),
                  actionButton(inputId = "GO_button",
                               label = "GO test")
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
                                 end = end_date)}
               
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
    # genes_data <- gene_atomization(pubmed_results())
    # # Codify frequency of genes as numeric
    # genes_table <- data.frame(genes_data,
    #                           stringsAsFactors = FALSE)
    # colnames(genes_table) <- c("Símbolo", "Nombre", "Frecuencia")
    genes_table <- pubmed_results() %>% gene_atomization() %>%
      data.frame(stringsAsFactors = FALSE) %>%
      `colnames<-`(c("Símbolo", "Nombre", "Frecuencia"))
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
    tabla_frecuencias$genes2 <- factor(tabla_frecuencias$Símbolo,
                                       levels = rev(factor(tabla_frecuencias$Símbolo)))
    freq_barplot(varcat = tabla_frecuencias$genes2,
                 varnum = tabla_frecuencias$Frecuencia,
                 main = "Genes más frecuentes")
  })
  
  ## GO-over-representation test
  # GO enrichment analysis of the gene set
  # Configure sample genes as entrezID
  # keys <- reactive({
  #   genes()[, "Símbolo"]
  #   # load(file = "../intermediateData/genes.RData") # Temporal
  #   # genes[, "Gene_symbol"]
  # })

  #  entrezID <- reactive({
  #   entrez <- select(org.Hs.eg.db,
  #          keys = genes()[, "Símbolo"],
  #          columns = c("SYMBOL", "ENTREZID"),
  #          keytype = "SYMBOL")
  #   entrez[, "ENTREZID"]
  # })
  # # Compute enrichResult objects of each ontology aspect
  # ego_cc <- eventReactive(input$GO_button,{
  #   req(input$select_aspect == 'Componente celular')
  #   withProgress(message = "Computing enriched GO terms", {
  #     ego_function( genes = entrezID(),
  #                 ontology = 'CC',
  #                  padjust = adjust_methods[[input$metodo_ajuste]],
  #                  pvalue = input$p_valor,
  #                  qvalue = input$q_valor)
  #     })
  # })
  # 
  # ego_bp <- eventReactive(input$GO_button,{
  #   req(input$select_aspect == 'Proceso biológico')
  #   ego_function( genes = entrezID(),
  #                 ontology = 'BP',
  #                 padjust = adjust_methods[[input$metodo_ajuste]],
  #                 pvalue = input$p_valor,
  #                 qvalue = input$q_valor)
  # })
  # 
  # ego_mf <- eventReactive(input$GO_button,{
  #   req(input$select_aspect == 'Función molecular')
  #   ego_function( genes = entrezID(),
  #                 ontology = 'MF',
  #                 padjust = adjust_methods[[input$metodo_ajuste]],
  #                 pvalue = input$p_valor,
  #                 qvalue = input$q_valor)
  # })
  # 
  # # Selected ego results
  # ego_object <- reactive(
  #       if (input$select_aspect == 'Componente celular') {
  #         withProgress(message = "Calculando términos GO enriquecidos \npara componentes celulares",
  #         ego_cc())
  #       } else if (input$select_aspect == 'Proceso biológico') {
  #         withProgress(message = "Calculando términos GO enriquecidos \npara procesos biológicos",
  #       ego_bp())
  #       } else if (input$select_aspect == 'Función molecular') {
  #         withProgress(message = "Calculando términos GO enriquecidos \npara funciones moleculares",
  #       ego_mf())
  #       }
  #   )
  # 
  # # Compose data frame from eGO results
  # ego_table <- eventReactive(
  #   input$GO_button, {
      # withProgress(message = "Computing enriched GO terms", {
    #     if (input$select_aspect == 'Componente celular') {
    #       ego_object <- ego_cc()
    #     } else if (input$select_aspect == 'Proceso biológico') {
    #       ego_object <- ego_bp()
    #     } else if (input$select_aspect == 'Función molecular') {
    #       ego_object <- ego_mf()
    #     }

        # incProgress(1/2)


        # as.data.frame(ego_cc()[, c("ID", "Description", "GeneRatio", "BgRatio", "p.adjust")])
      # })
    # })

  # Display results in table or barplot
  observeEvent(input$select_display, {
    updateTabsetPanel(
      inputId = "tabla_grafico",
      selected = input$select_display)
  })

  # GO terms table
  # output$GOterms <- DT::renderDataTable({
  #   datatable(ego_table(),
  #             rownames = FALSE,
  #             colnames = c("GO_ID", "Descripción", "GeneRatio", "BgRatio", "p-valor ajustado"),
  #             selection = list(mode = 'single', selected = 1),
  #             options = list(language = list(url = 'spanish.json'))) %>%
  #     formatSignif('p.adjust', 2) %>%  # Significative digits for p.adjust column
  #     formatStyle(columns = c("GeneRatio", "BgRatio", "p.adjust"), `text-align` = 'center') # Center columns
  # })

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
  
    # output$GO_barplot <- renderPlot(barplot(height = ego_object(),
    #                      showCategory = input$go_categories,
    #                      title = paste0("Términos GO enriquecidos \n(",
    #                                     input$select_aspect,
    #                                     ")")) + labs(y = "Número de genes en cada categoría"),
    #                      # Plot size enough to display all chosen categories
    #                      height = reactive(
    #                        max(400, input$go_categories * 20)),
    #                      res = 96,
    #                      alt = 'Gráfica de barras')
  
  
  # 
  # output$ontologyCC <- renderPlot(height = ego_CC,
  #                                 showCategory = 20,
  #                                 title = "Título")
  
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
# profvis::profvis(runApp(shinyApp(ui, server)))

shinyApp(ui, server)

