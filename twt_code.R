library("readxl")
library("plyr")
library("xlsx")

mut <- c("BRAF", "NRAS", "NF1", "KRAS", "GNAQ", "GNA11", "KIT")

s1 <- as.data.frame(read_excel("NIHMS362881-supplement-3.xlsx", sheet = "TABLE S1", skip = 6, col_names = T))
s1 <- s1[, c(1, 16, 17)]
s1[, 1] <- gsub("T", "", s1[,1])

s3 <- read.table("NIHMS393895-supplement-03.txt", header = T, sep = "\t", fill=T)
s3 <- s3[, c(1, 7, 8, 9, 10, 11, 14, 17, 24)]
s3[, c(3, 4, 5, 6, 8)] <- NULL
s3 <- unique(subset(s3, grepl("ME0", s3[, 3]) & s3[, 1] %in% mut))
names(s3)[3] <- "Sample"
s3[, 3] <- gsub("-.*", "", s3[, 3])
s3[, 4] <- gsub("p.", "", s3[, 4])

s4 <- as.data.frame(read_excel("NIHMS362881-supplement-3.xlsx", sheet = "TABLE S4", skip = 5, col_names = T))
s4 <- s4[, c(1, 2, 20)]
s4 <- unique(subset(s4, grepl("ME0", s4[, 1]) & s4[, 2] %in% mut))

s6 <- as.data.frame(read_excel("NIHMS362881-supplement-3.xlsx", sheet = "TABLE S6", skip = 5, col_names = T))
s6 <- s6[, c(1, 2, 19)]
s6 <- unique(subset(s6, grepl("ME0", s6[, 1]) & s6[, 2] %in% mut))
#s6 produce only KIT YIDPTQL570del, which is included in s3 so we can ignore s6

s10 <- as.data.frame(read_excel("NIHMS362881-supplement-3.xlsx", sheet = "TABLE S10", skip = 6, col_names = T))
s10 <- s10[, c(1, 6, 7)]
s10 <- subset(s10, grepl("ME0", s10[, 1]))
s10[, 1] <- gsub("T", "", s10[,1])

s4A <- read.table("NIHMS362881-supplement-2.txt", header = T, sep = "\t")
s4A <- s4A[, c(1, 7, 8)]
s4A <- unique(subset(s4A, grepl("ME0", s4A[, 1]) & s4A[, 3] %in% mut))
names(s4A)[1] <- "Sample"

s6_1 <- as.data.frame(read_excel("NIHMS393895-supplement-02.xlsx", sheet = "TABLE S6", skip = 5, col_names = T))
s6_1 <- s6_1[, c(1, 2, 6, 9)]
s6_1 <- subset(s6_1, grepl("ME0", s6_1[, 1]) & (s6_1[, 2] %in% mut))
s6_1[, 4] <- gsub("p.", "", s6_1[, 4])
names(s6_1)[1] <- "Sample"

#join s1 and s10
names(s1)[c(2, 3)] <- c("BRAF", "NRAS")
s1[, 2] <- gsub("wild type", "", s1[, 2])
s1[, 3] <- gsub("wild type", "", s1[, 3])
names(s10)[c(2, 3)] <- c("BRAF", "NRAS")
#keep NRAS with only Q61 ------> this will cause NRAS with other mutations considered as wild type
#s1 <- s1[!(!grepl("", s1$BRAF) & !grepl("Q61", s1$NRAS)),]
s110 <- join(s1, s10, type="full")

#join s4 and s6_1
s4[, 3] <- gsub("p.", "", s4[, 3])
names(s4)[1] <- "Sample"
#variant classification in s6_1 is not important so we can ignore this column
s6_1[, 3] <- NULL
names(s6_1)[2] <- "Hugo_Symbol"
names(s6_1)[3] <- "Protein_Change"
s461 <- join(s4, s6_1, type="full")
#replace NA values
s461$Protein_Change <- as.character(s461$Protein_Change)
s461$Protein_Change[is.na(s461$Protein_Change)] <- ""

#join s4A and s3
#remove rows with intron mutation
s3 <- s3[!(s3$Variant_Classification == "Intron"),]
s4A <- s4A[!(s4A$Variant_Classification == "Intron"),]
#maintain NRAS with only Q61
s3 <- s3[!(s3$Hugo_Symbol == "NRAS" & !grepl("Q61", s3$Protein_Change)),]
#remove rows in s4A contained in s3
s4A <- anti_join(s4A, s3, by=c("Sample", "Hugo_Symbol", "Variant_Classification"))
#add LoH to NF1 mut
s3$Variant_Classification <- as.character(s3$Variant_Classification)
s3$Variant_Classification[s3$Hugo_Symbol == "NF1"] <- paste0(s3[, 2], ", LoH")
s4A3 <- merge(s4A, s3, by=c("Sample", "Hugo_Symbol"), all=T)
#replace NA values by ""
s4A3$Variant_Classification.x <- as.character(s4A3$Variant_Classification.x)
s4A3$Variant_Classification.x[is.na(s4A3$Variant_Classification.x)] <- ""
s4A3$Variant_Classification.y <- as.character(s4A3$Variant_Classification.y)
s4A3$Variant_Classification.y[is.na(s4A3$Variant_Classification.y)] <- ""
#merge 2 Variant_Classification columns
s4A3$Variant_Classification <- paste(s4A3$Variant_Classification.x, s4A3$Variant_Classification.y)
#remove 2 Variant_Classification columns
s4A3$Variant_Classification.x <- NULL
s4A3$Variant_Classification.y <- NULL

#remove rows in s461 contained in s4A3
s461 <- anti_join(s461, s4A3, by=c("Sample", "Hugo_Symbol", "Protein_Change"))
#join s4A3 and s461
s4A3_461 <- merge(s4A3, s461, by=c("Sample", "Hugo_Symbol"), all=T)
#remove all na values
s4A3_461$Protein_Change.x <- as.character(s4A3_461$Protein_Change.x)
s4A3_461$Protein_Change.x[is.na(s4A3_461$Protein_Change.x)] <- ""
s4A3_461$Protein_Change.y <- as.character(s4A3_461$Protein_Change.y)
s4A3_461$Protein_Change.y[is.na(s4A3_461$Protein_Change.y)] <- ""
# merge 2 Protein Change columns
s4A3_461$Protein_Change <- paste(s4A3_461$Protein_Change.x, s4A3_461$Protein_Change.y)
#remove 2 Protein Change columns
s4A3_461$Protein_Change.x <- NULL
s4A3_461$Protein_Change.y <- NULL
#combine Variant Classification and Protein Change
s4A3_461$Specification <- paste(s4A3_461$Protein_Change, s4A3_461$Variant_Classification)
#remove columns
s4A3_461$Variant_Classification <- NULL
s4A3_461$Protein_Change <- NULL
#unmelt data
s4A3_461 <- reshape(s4A3_461, idvar = "Sample", timevar = "Hugo_Symbol", direction = "wide")
names(s4A3_461)[c(2:8)] <- gsub("Specification.", "", names(s4A3_461)[c(2:8)])

#remove rows of s110 contained in s4A3_416
result <- merge(s110, s4A3_461, by="Sample", all=T)
#delete BRAF.x and NRAS.x that are already in BRAF.y, NRAS.y
result$BRAF.x[grepl(result$BRAF.x, result$BRAF.y)] <- ""
result$NRAS.x[grepl(result$NRAS.x, result$NRAS.y)] <- ""

#remove NA values in BRAF and NRAS
result$BRAF.x[is.na(result$BRAF.x)] <- ""
result$BRAF.y[is.na(result$BRAF.y)] <- ""
result$NRAS.x[is.na(result$NRAS.x)] <- ""
result$NRAS.y[is.na(result$NRAS.y)] <- ""

#merge BRAF x y, NRAS x y
result$BRAF <- paste(result$BRAF.x, result$BRAF.y)
result$NRAS <- paste(result$NRAS.x, result$NRAS.y)

#remove BRAF, NRAS x y
result[, c("BRAF.x", "BRAF.y", "NRAS.x", "NRAS.y")] <- NULL

#remove all NAs
result$GNAQ[is.na(result$GNAQ)] <- ""
result$GNA11[is.na(result$GNA11)] <- ""
result$KRAS[is.na(result$KRAS)] <- ""
result$NF1[is.na(result$NF1)] <- ""
result$KIT[is.na(result$KIT)] <- ""

result <- result[, c(1, 7, 8, 4, 2, 3, 6, 5)]
write.xlsx(result, "TWT_specification.xlsx", sheetName="1", col.names=TRUE, row.names=TRUE, append=FALSE)

library("shiny")
ui <- fluidPage(
  titlePanel(strong("Identification of mutated samples")),
  
  sidebarLayout(
    sidebarPanel(
      h4("Genes and mutational status"),
      h6("*Check the box to view the mutational status (mutated or non-mutated) of the gene"),
      h6("*Uncheck the box for out-of-consideration gene"),
      tags$br(),
      h6("*BRAF mutations indicate only V600E ones"),
      lapply(mut, function(gene) {
        fluidRow(
          column(3, checkboxInput(inputId=toString(gene), label=toString(gene), TRUE)),
          column(9, selectInput(inputId=paste0("status", toString(gene)), label = NULL, choices = c("mutated", "non-mutated"))))
      })
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(align="center",
      
      # Output: Formatted text for caption ----
      h4(strong("Mutational status")),
      
      # Output: HTML table
      tableOutput("view")
      
    )
  )
)

server <- function(input, output) {
  datasetInput <- reactive({
    res <- result
    for (gene in mut) {
      if (input[[gene]]) {
        if (input[[paste0("status", gene)]] == "mutated") {
          temp <- subset(result, result[, toString(gene)] == TRUE)
        } 
        else {
          temp <- subset(result, result[, toString(gene)] == FALSE)
        }
        res <- join(res, temp, type="inner", match = "all")
      }
      else {
        # case of uncheck
        res[, gene] <- NULL
      }
    }
    # another case of uncheck
    #for (gene in mut) {
      # we do not care about this gene (mutated or not)
      #if (!input[[gene]]){
        #res[, gene] <- NULL
      #}
    #}
    res
  })
  
  output$view <- renderTable({
    datasetInput()
  })
}
shinyApp(ui = ui, server = server)
