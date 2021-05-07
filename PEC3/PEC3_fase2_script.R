### packages ###
library(easyPubMed)
library(pubmed.mineR)
library(clusterProfiler) # GO enrichment
library(org.Hs.eg.db)    # GO enrichment and
                         #    translation of gene ids
library(enrichplot)      # GO enrichment barplots
library(ggplot2)         # GO enrichment barplots

### Variables ###
# Default variables retrieve all results for endometriosis
keywords <- c("endometriosis")
# Dates must be in format YYYY/MM/DD
first_date <- "1800/12/31"
last_date <- format(Sys.Date() + 1, "%Y/%m/%d")

### Functions ###

# Generate query
query <- function (varkeywords=keywords, date1=first_date, date2=last_date) {
  paste(c(keywords, " AND " , date1, ":", date2,"[dp]"),
        collapse="")
}

# Horizontal barplot
freq_barplot <- function(varcat, varnum, main = ""){ # Categorical variable and numerical variable
  # Adjust width of left margin
  # https://stackoverflow.com/questions/10490763/automatic-adjustment-of-margins-in-horizontal-bar-chart
  par(mar=c(5.1, 
            max(4.3,max(nchar(as.character(varcat)))/1.8) ,
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
  text(rev(varnum),
       y = y,
       labels = rev(varnum),
       adj = NULL,
       pos = 4,
       cex = 0.9
  )
}

### Script main body

# Create directory structure

# 'data' contains raw source data.
# 'intermediateData' contains .RData objects with processed data.
# 'results' stores final report files.

directories <- c("script/data", "script/results", "script/intermediateData")

# Create directories
lapply(directories, function(x){
  if (!(dir.exists(x))){
    dir.create(x,
               recursive = TRUE)
  }
})

# Retrieve query results

batch_pubmed_download(query(),
                      dest_dir = "script/data/",
                      dest_file_prefix = "total_endometriosis",
                      format = "abstract",
                      batch_size = 5000
)

#concatenate text files
# List of files to be added together
files_list <- list.files(path = "script/data/",
                         pattern = "total",
                         full.names = TRUE) # include path
# Create new file
out_file <- file(description = "script/data/todos.txt",
                 open = "w")
# Read each downloaded file and write into final file
for (i in files_list){
  x <- readLines(i)
  writeLines(x, out_file)
}

close(out_file)

# Delete batch text files
file.remove(files_list)

# Generate S4 object of class 'Abstract' (corpus primario):
  
# Generate the object
abstracts <- readabs("script/data/todos.txt")
# Save object
save(abstracts, file = "script/intermediateData/abstracts.RData")

# Word atomization:
words <- word_atomizations(abstracts)
save(words, file = "script/intermediateData/words.RData")
# Move text file to results directory
file.rename(from = "word_table.txt",
            "script/results/words.txt")

# Barplot of word frequencies:

# Select the twelve most frequent
words2 <- words[1:12,]
# Reverse order factors
words2$words2 <- factor(words2$words, 
                             levels = rev(factor(words2$words)))
# Draw barplot and save as png
# Open png file
png(filename = "script/results/wordsbarplot.png")
# Create plot
freq_barplot(varcat = words2$words2,
             varnum = words2$Freq,
             main = "Palabras más frecuentes")
# Close file
dev.off()

# Gene extraction
genes <- gene_atomization(abstracts)
save(genes, file = "script/intermediateData/genes.RData")
# Move text file to results directory
file.rename(from = "table.txt",
            "script/results/genes.txt")

# Gene barplot:
  
# Select the twelve most frequent
genes2 <- data.frame(genes[1:12,],
                          stringsAsFactors = FALSE)
# Codify frequency as numeric
genes2$Freq <- as.numeric(genes2$Freq)
# Reverse order factors
genes2$genes2 <- factor(genes2$Gene_symbol, 
                             levels = rev(factor(genes2$Gene_symbol)))
# Draw barplot and save as png
# Open png file
png(filename = "script/results/genebarplot.png")
# Create plot
freq_barplot(varcat = genes2$genes2,
             varnum = genes2$Freq,
             main = "Genes más frecuentes")
# Close file
dev.off()

## Gene characterization (GO enrichment)

# Translate gene symbols to entrezId
# Beware: the result is a dataframe
keys <- genes[, "Gene_symbol"]

entrez <- select(org.Hs.eg.db,
                 keys = keys,
                 columns = c("SYMBOL", "ENTREZID"),
                 keytype = "SYMBOL")

# Get all genes ID from pubtator contingency table
load("PEC2/data/geneID_frequencies.RData")
# List of entrezID in org.Hs.eg.db
human_genes_entrezid <- keys(org.Hs.eg.db)
# Vector of all genes from pubtator list
universe <- names(geneID_frequencies)
# Filter by human genes vector
universe <- names(geneID_frequencies)[names(geneID_frequencies) %in% human_genes_entrezid]
# Enrichment test: biological processes
for (ontology in c("BP", "CC", "MF")) {
  ego <- enrichGO(gene = entrez$ENTREZID,
                   universe = universe,
                   OrgDb = org.Hs.eg.db,
                   ont = ontology,
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.5,
                   readable = TRUE)

  # Generate results as a csv file
  write.csv(ego,
            file= paste0("script/results/ego_", ontology, ".csv"))
  
  # Simplify redundant results
  ego_simplified <- simplify(ego, cutoff = 0.7, by = "p.adjust", select_fun = min)
  # Generate csv files from simplified results
  write.csv(ego,
            file= paste0("script/results/ego_", ontology, "_simplified.csv"))

  # Generate barplots of enrichment results
  # Terms for plot titles
  go_plot_titles <- c("BP" = "procesos biológicos",
                              "CC" = "componentes celulares",
                              "MF" = "funciones moleculares")
  
  # Create plot
  barplot(height = ego_simplified,
          showCategory = 20,
          title = paste0("Términos GO enriquecidos \n(",
                         go_plot_titles[[ontology]],
                         ")"))
  # Close png file
  ggsave(filename = paste0("script/results/go_", ontology, "barplot.png"))
  
}
