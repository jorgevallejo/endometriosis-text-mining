### packages ###
library(easyPubMed)
library(pubmed.mineR)

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
            max(4.1,max(nchar(as.character(varcat)))/1.8) ,
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

directories <- c("data", "results", "intermediateData")

# Create directories
lapply(directories, function(x){
  if (!(dir.exists(x))){
    dir.create(x)
  }
})

# Retrieve query results

batch_pubmed_download(query(),
                      dest_dir = "data/",
                      dest_file_prefix = "total_endometriosis", # Correct this #
                      format = "abstract",
                      batch_size = 5000
)

#concatenate text files
# List of files to be added together
files_list <- list.files(path = "data/",
                         pattern = "total",
                         full.names = TRUE) # include path
# Create new file
out_file <- file(description = "intermediateData/todos.txt",
                 open = "w")
# Read each downloaded file and write into final file
for (i in files_list){
  x <- readLines(i)
  writeLines(x, out_file)
}

close(out_file)

# Generate S4 object of class 'Abstract' (corpus primario):
  
# Generate the object
abstracts <- readabs("intermediateData/todos.txt")
# Save object
save(abstracts, file = "intermediateData/abstracts.RData")

# Word atomization:
words <- word_atomizations(abstracts)
save(words, file = "intermediateData/words.RData")
# Move text file to results directory
file.rename(from = "word_table.txt",
            "results/words.txt")

# Barplot of word frequencies:

# Select the twelve most frequent
words2 <- words[1:12,]
# Reverse order factors
words2$words2 <- factor(words2$words, 
                             levels = rev(factor(words2$words)))
# Draw barplot and save as png
# Open png file
png(filename = "results\wordsbarplot.png")
# Create plot
freq_barplot(varcat = words2$words2,
             varnum = words2$Freq,
             main = "Palabras más frecuentes")
# Close file
dev.off()

# Gene extraction
genes <- gene_atomization(abstracts)
save(genes, file = "intermediateData/genes.RData")
# Move text file to results directory
file.rename(from = "table.txt",
            "results/genes.txt")

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
png(filename = "results\genebarplot.png")
# Create plot
freq_barplot(varcat = genes2$genes2,
             varnum = genes2$Freq,
             main = "Genes más frecuentes")
# Close file
dev.off()