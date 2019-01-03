library("readxl")
library("plyr")
library("xlsx")

mut <- c("BRAF", "NRAS", "NF1")

s6 <- as.data.frame(read_excel("NIHMS393895-supplement-02.xlsx", sheet = "TABLE S6", skip = 5, col_names = T))
s6 <- s6[, c(1, 2, 9)]
s6 <- subset(s6, grepl("ME0", s6[, 1]) & (s6[, 2] %in% mut))
s6[, c("BRAF", "NRAS", "NF1")] <- c(grepl("BRAF", s6[, 2]) & grepl("p.V600E", s6[, 3]), grepl("NRAS", s6[, 2]), grepl("NF1", s6[, 2]))
s6[, c(2:3)] <- NULL
s6 <- aggregate(s6, by=list(s6[, 1]), FUN = any)
s6[, 2] <- NULL
names(s6)[1] <- "Sample"

s4 <- as.data.frame(read_excel("NIHMS362881-supplement-3.xlsx", sheet = "TABLE S4", skip = 5, col_names = T))
s4 <- s4[, c(1, 2, 20)]
s4 <- subset(s4, grepl("ME0", s4[, 1]) & (s4[, 2] %in% mut))
s4[, c("BRAF", "NRAS", "NF1")] <- c(grepl("BRAF", s4[, 2]) & grepl("p.V600E", s4[, 3]), grepl("NRAS", s4[, 2]), grepl("NF1", s4[, 2]))
s4[, c(2:3)] <- NULL
s4 <- aggregate(s4, by=list(s4[, 1]), FUN = any)
s4[, 2] <- NULL
names(s4)[1] <- "Sample"

s10 <- as.data.frame(read_excel("NIHMS362881-supplement-3.xlsx", sheet = "TABLE S10", skip = 6, col_names = T))
s10 <- s10[, c(1, 6, 7)]
s10 <- subset(s10, grepl("ME0", s10[, 1]))
s10[, 1] <- gsub("T", "", s10[,1])
s10[, 2] <- grepl("V600E", s10[, 2])
s10[, 3] <- !is.na(s10[, 3])
names(s10)[c(2, 3)] <- c("BRAF", "NRAS")
s10$NF1 <- F

s1 <- as.data.frame(read_excel("NIHMS362881-supplement-3.xlsx", sheet = "TABLE S1", skip = 6, col_names = T))
s1 <- s1[, c(1, 16, 17)]
s1[, 1] <- gsub("T", "", s1[,1])
s1[, c(2, 3)] <- c(!grepl("wild type", s1[, 2]), !grepl("wild type", s1[, 3]))
names(s1)[c(2, 3)] <- c("BRAF", "NRAS")
s1$NF1 <- F

result <- aggregate(. ~ Sample, rbind.fill(s1, s4, s6, s10), FUN = any)
result$twt <- !(result$BRAF | result$NRAS | result$NF1)
write.xlsx(result, "TWT_res.xlsx", sheetName="1", col.names=TRUE, row.names=TRUE, append=FALSE)

#KIT mutation for twt
new_mut <- c("KIT")
twt <- subset(result, result$twt)

s4KIT <- as.data.frame(read_excel("NIHMS362881-supplement-3.xlsx", sheet = "TABLE S4", skip = 5, col_names = T))
s4KIT <- s4KIT[, c(1, 2)]
s4KIT <- subset(s4KIT, (s4KIT[, 1] %in% twt[, 1]) & (s4KIT[, 2] %in% new_mut))
#no rows selected

s6KIT <- as.data.frame(read_excel("NIHMS362881-supplement-3.xlsx", sheet = "TABLE S6", skip = 5, col_names = T))
s6KIT <- s6KIT[, c(1, 2)]
s6KIT <- subset(s6KIT, (s6KIT[, 1] %in% twt[, 1]) & (s6KIT[, 2] %in% new_mut))
s6KIT[, 2] <- T
names(s6KIT) <- c("Sample", "KIT")

s4AKIT <- read.table("NIHMS362881-supplement-2.txt", header = T, sep = "\t")
s4AKIT <- s4AKIT[, c(1, 8)]
s4AKIT <- unique(subset(s4AKIT, (s4AKIT[, 1] %in% twt[, 1]) & (s4AKIT[, 2] %in% new_mut)))
s4AKIT[, 2] <- T
names(s4AKIT) <- c("Sample", "KIT")

twt$KIT <- (twt[, 1] %in% s4KIT[, 1]) | (twt[, 1] %in% s6KIT[, 1]) | (twt[, 1] %in% s4AKIT[, 1])
write.xlsx(twt, "TWT_KIT.xlsx", sheetName="1", col.names=TRUE, row.names=TRUE, append=FALSE)

library("shiny")
ui <- fluidPage(
  titlePanel(strong("Identification of mutated samples")),
  
  sidebarLayout(
    sidebarPanel(
      h4("Genes and mutational status"),
      h6("*Check the box to view the mutational status (mutated or non-mutated) of the gene"),
      h6("*Uncheck the box for out-of-consideration gene"),
      tags$br(),
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