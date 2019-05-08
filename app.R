###########################################
# VCF Annotation Tool
###########################################

##############library######
library(shiny)
library(ggplot2)
library(plotly)
library(gridExtra)
library(DT)
library(GGally)

source("./scripts/functions.R")
#############################


###################ui##########
ui <- bootstrapPage(
    tags$head(tags$style(".shiny-progress {top: 50% !important;left: 50% !important;margin-top: -100px !important;margin-left: -250px !important; color: blue;font-size: 20px;font-style: italic;}")),
    titlePanel(div(" Xiu Huang's VCF Annotation Workshop"), windowTitle = "VCF Annotation"),
    fluidRow(
        sidebarPanel(
            width=3,
            wellPanel(
                p(" Please select vcf data"),
                fileInput("dataFile", "",
                          accept = c(".vcf")
                ),
                actionButton("load", "Start Annotation!")
            )
        ),
        mainPanel(
            width=9,
            tabsetPanel(
                tabPanel(
                    "Annotated Output", 
                    value = 1, 
                    dataTableOutput('data_adj'),
                    uiOutput("download")
                )
            )
        )
    )
)


server <- function(input, output, session) {
    
    data <- reactiveValues(List=NULL)
    
    observeEvent(input$load,{
      vcf_file <- input$dataFile$datapath
      data$vcf_file <- vcf_file
      data$df <- parse_vcf(vcf_file)
    })

    output$download <- renderUI({
        if(!is.null(data$df) & !is.null(data$vcf_file)) {
            downloadButton('dl_data_adj', 'Download Annotated Vcf')        }
    })
        
    output$data_adj <- renderDataTable({
        dat = data$df
        sketch = htmltools::withTags(table(
            class = 'display',
            thead(
                tr(
                    th('chr', title = 'Chromosome'),
                    th('start', title = 'Start position'),
                    th('end', title = 'End position'),
                    th('REF', title = 'Reference Allele'),
                    th('ALT', title = 'Alternative Allele'), 
                    th('TYPE', title = 'The type of allele, either snp, mnp, ins, del, or complex. '),
                    th('CIGAR', title = 'The extended CIGAR representation of each alternate allele'), 
                    th('WCONSEQUENCE', title = 'Type of variation (Substitution, Insertion, Silent, Intergenic, etc.) If there are multiple possibilities, annotate with the most deleterious possibility.'),
                    th('DP', title = 'Depth of sequence coverage at the site of variation.'), 
                    th('AO', title = 'Number of reads supporting the variant.'), 
                    th('AF', title = 'Percentage of reads supporting the variant versus those supporting reference reads.'), 
                    th('ExAC_var', title = 'Relevant variant IDs from ExAC database. '), 
                    th('allele_freq', title = 'Allele frequency of variant from Broad Institute ExAC Project API. '), 
                    th('category', title = 'Variant category from Broad Institute ExAC Project API. '), 
                    th('major_consequence', title = 'Major consequence from Broad Institute ExAC Project API. '), 
                    th('HGVSp', title = 'HGVSp from Broad Institute ExAC Project API. ')
                )
            )
        ))
        if (!is.null(dat)) dat[, ExAC_var := ifelse(!is.na(ExAC_var), sapply(ExAC_var, function(x) paste(sprintf('<a href="http://exac.broadinstitute.org/variant/%1$s"  target="_blank" >%1$s</a>', unlist(strsplit(unlist(strsplit(x, ",")), ";"))), collapse=", ")), "")]
        datatable(dat, filter = 'top', selection = 'single', options = list(lengthMenu = c(50, 100, 200), pageLength = 50, autoWidth = TRUE), rownames=FALSE, container = sketch, escape = FALSE)
    })
    
    output$dl_data_adj <- downloadHandler(
        filename = function() {
            "annotated_vcf.vcf"
        },
        content = function(file) {
            dat = data$df
            vcf_file = data$vcf_file
            vcf <- fread(file = vcf_file, sep = "\t", sep2 = "\t", header = T, skip = "#CHROM")
            header <- rbind(fread(vcf_file, header = F, sep = "\t")[grepl("^##", V1)], data.table(V1 = '##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations: \'Worst consequence | ExAc Var ID | ExAC_AF | ExAC_Category | ExAC_major_consequence | ExAC_HGVSp\' ">'))
            header[, V1 := ifelse(grepl("fileDate", V1), paste("##fileDate=", format(Sys.Date(), "%Y%m%d"), sep = ""), V1)]
            vcf[, INFO := paste(vcf$INFO, ";ANN=", apply(dat, 1, function(x){paste(x[c(8, 12:15)], collapse = "|")}), sep = "")]
            fwrite(header, file = file, col.names = F, append = F, eol = "\n", quote = F)
            fwrite(vcf, file = file, col.names = T, append = T, sep = "\t", quote = F, eol = "\n")
        }
    )
}

shinyApp(ui = ui, server = server)
