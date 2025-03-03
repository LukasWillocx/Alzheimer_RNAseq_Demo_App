# List of required packages
packages <- c("shiny", "shinythemes", "shinycssloaders", 
              "networkD3", "dplyr", "tidyr", "ggplot2",
              "grid", "htmlwidgets", "clusterProfiler", 
              "AnnotationHub", "DESeq2")

for(pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

source('functions_app2.R')

ui <- fluidPage(
  includeCSS("pf_styles.css"),
  tags$script(HTML("
    function toggleTheme() {
      const body = document.body;
      body.dataset.theme = body.dataset.theme === 'dark' ? 'light' : 'dark';
    }
  ")),
  
  titlePanel(
    div("Data Science Portfolio", 
        style = "display: flex; justify-content: space-between; align-items: center;",
        tags$button(
          id = "themeToggle",
          onclick = "toggleTheme()",
          "🌓",
          class = "theme-toggle-btn"
        )
    ), 
    windowTitle = "DS Portfolio"
  ),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Navigation"),
      selectInput("dataset", "Choose Dataset:",
                  choices = c("Iris", "mtcars", "diamonds")),
      hr(),
      tags$div(tags$div(style = "text-align: center",h4('Dataset Information'),
                        h5('Bulk RNA-seq of postmortem human prefrontal cortex and cerebellum brain tissues')),
               hr(),
               HTML(paste('<strong>Accession number:</strong> <a href="https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-14167" target="_blank">E-MTAB-14167</a>','<br>')),
               HTML(paste("<strong>Format:</strong>",'.txt','<br>')),
               HTML(paste("<strong>Assay technology:</strong>",'RNA sequencing','<br>')),
               HTML(paste("<strong>Sample Size:</strong> ",'91','<br>')),
               HTML(paste("<strong>Number of Observations:</strong> ", '62705','<br>')),
               HTML(paste('<strong>Reference Genome:</strong> <a href="https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/" target="_blank">GRCh38.p14</a>','<br>')),
               HTML(paste('<strong>Sample Attributes:</strong> <a href="https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-14167/sdrf" target="_blank">SDRF</a>','.txt','<br>')),
               HTML(paste("<strong>File Size:</strong>",'53.5 MB','<br>')),
               HTML(paste("<strong>Downloaded on:</strong>",'February 18th, 2025','<br>'))
      ),
      hr(),
      tags$div(style = "text-align: center",h5('R Libraries')),
      hr(),
      tags$ul(style = "list-style-type: none; padding-left: 0;", 
              lapply(packages, function(lib) {
                tags$li(style = "display: inline-block; margin: 5px; text-align: center;",
                        tags$span(
                          style = "border-radius: 20%; background-color: var(--background-color); padding: 5px;",
                          tags$b(lib)))})),
      hr(),
      tags$div(style = "text-align: center",h5('About the application')),
      hr(),
    ),
    
    mainPanel(
      width = 9,
      tabsetPanel(
        type = "tabs",
        tabPanel("Overview",
                 fluidRow(
                   column(12, div(class = "custom-panel",
                                 h3("Sample attributes"),
                                 withSpinner(sankeyNetworkOutput("sankey")))),
                 ),
                 fluidRow(
                   column(12, div(class = "custom-panel",
                                  h3("Reference genome for RNA seq mapping"),
                                  withSpinner(plotOutput("circos"))))
                 )
        ),
        tabPanel("Comparison ~sex",
                 div(class = "custom-panel",
                     h3("Volcano plot - female Vs male (cerebellum tissue)"),
                     withSpinner(plotOutput("DEA_sex_cerebellum"))),
                 div(class = "custom-panel",
                     h3("Volcano plot - female Vs male (prefrontal cortex tissue)"),
                     withSpinner(plotOutput("DEA_sex_pfc"))),
                 p(''),
                 p('The visualization in this panel serves as a sanity check. By considering sex as the grouping variable to perform the differential gene expression analysis, 
                   it is safe to assume the result would have a heavy emphasis on Y chromosome related gene expression. The statistically significant differential gene expressions
                   (an adjusted p-value <0.05) are labeled with their chromosome denomination. Besides the horizontally dotted line to identify the statistical significance threshold, vertically dotted lines
                   are added to denote biological significance. These are drawn at a log2 Fold Change of -2 and 2, to represent an absolute change in gene expression between the two groups of at least a factor four.'),
        ),
        tabPanel("Comparison ~genotype",
                 div(class = "custom-panel",
                     h3("Volcano plot - ApoE3,3 Vs ApoE4,4 (cerebellum tissue)"),
                     withSpinner(plotOutput("DEA_genotype_cerebellum"))),
                 div(class = "custom-panel",
                     h3("Volcano plot - ApoE3,3 Vs ApoE4,4 (prefrontal cortex tissue)"),
                     withSpinner(plotOutput("DEA_genotype_pfc")))
        ),
        tabPanel("Comparison ~disease",
                 div(class = "custom-panel",
                     h3("Volcano plot - Alzheimer's Vs none (cerebellum tissue)"),
                     withSpinner(plotOutput("DEA_disease_cerebellum"))),
                 div(class = "custom-panel",
                     h3("Volcano plot - Alzheimer's Vs none (prefrontal cortex tissue)"),
                     withSpinner(plotOutput("DEA_disease_pfc")))
        ),
      )
    )
  )
)

server <- function(input, output) {
  
  # loading prewritten data to reduce on the fly server calculations
  meta_data_sorted  <- read.csv('meta_data_sorted.csv')
  read_counts_clean <- read.csv('read_counts_clean.csv')
  
  DEA_sex_cerebellum_data<-read.csv('DEA_sex_results_in_cerebellum.csv')
  DEA_sex_pfc_data<-read.csv('DEA_sex_results_in_pfc.csv')
  
  DEA_genotype_cerebellum_data<-read.csv('DEA_genotype_results_in_cerebellum.csv')
  DEA_genotype_pfc_data<-read.csv('DEA_genotype_results_in_pfc.csv')
  
  DEA_disease_cerebellum_data<-read.csv('DEA_disease_results_in_cerebellum.csv')
  DEA_disease_pfc_data<-read.csv('DEA_disease_results_in_pfc.csv')
  
  
  output$sankey<-renderSankeyNetwork({
    sankey_design(meta_data_sorted)
  })
  #sex
  output$DEA_sex_cerebellum<-renderPlot({
    DEA_volcano_plotter(DEA_sex_cerebellum_data)+
      geom_label(aes(label = ifelse(stat_sign, Chr, NA),color=stat_sign),fill='#f0d09f')+
      scale_x_continuous(limits = c(-11, 11))
  },bg='transparent')
  
  output$DEA_sex_pfc<-renderPlot({
    DEA_volcano_plotter(DEA_sex_pfc_data)+
      geom_label(aes(label = ifelse(stat_sign, Chr, NA),color=stat_sign),fill='#f0d09f')+
      scale_x_continuous(limits = c(-11, 11))
  },bg='transparent')
  
  #genotype
  output$DEA_genotype_cerebellum<-renderPlot({
    DEA_volcano_plotter(DEA_genotype_cerebellum_data)+
      geom_label(aes(label = ifelse(stat_sign, Chr, NA),color=stat_sign),fill='#f0d09f')+
      scale_x_continuous(limits = c(-6, 6))
  },bg='transparent')
  
  output$DEA_genotype_pfc<-renderPlot({
    DEA_volcano_plotter(DEA_genotype_pfc_data)+
      geom_label(aes(label = ifelse(stat_sign, Chr, NA),color=stat_sign),fill='#f0d09f')+
      scale_x_continuous(limits = c(-6, 6))
  },bg='transparent')
  
  #disease
  output$DEA_disease_cerebellum<-renderPlot({
    DEA_volcano_plotter(DEA_disease_cerebellum_data)+
      geom_label(aes(label = ifelse(stat_sign, Chr, NA),color=stat_sign),fill='#f0d09f')+
      scale_x_continuous(limits = c(-5, 5))
  },bg='transparent')
  
  output$DEA_disease_pfc<-renderPlot({
    DEA_volcano_plotter(DEA_disease_pfc_data)+
      geom_label(aes(label = ifelse(stat_sign, Chr, NA),color=stat_sign),fill='#f0d09f')+
      scale_x_continuous(limits = c(-5, 5))
  },bg='transparent')
  
}

shinyApp(ui = ui, server = server)