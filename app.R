# List of required packages
packages <- c("shiny", "shinythemes", "shinycssloaders", 
              "networkD3", "dplyr", "tidyr", "ggplot2",
              "grid", "htmlwidgets", "clusterProfiler", 
              "AnnotationHub", "DESeq2","plotly")

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
                                  h3("Samples in principal component vector space, labeled according to their categorical attributes"),
                                  p('The principal component space entails the dimensional reduction of the 62705 distinct RNASeq reads in the data set. The first two
                                  principal components comprise nearly 97% of the variance observed in the data. Upon closer inspection it can be argued that the separation along principal component 2 quite neatly
                                  represents the biological difference between the two sampling tissues of cerebellum and prefrontal cortex. None of the other attributes show
                                  clear separation along principal component 1. 
                                    ')
                                  )),
                 ),
                 fluidRow(
                   column(3, div(class = "custom-panel",
                                  h5('Attribute: disease'),
                                  withSpinner(plotlyOutput("PCA_samples_plot_disease")))),
                   column(3, div(class = "custom-panel",
                                  h5('Attribute: genotype'),
                                  withSpinner(plotlyOutput("PCA_samples_plot_genotype")))),
                   column(3, div(class = "custom-panel",
                                  h5('Attribute: sex'),
                                  withSpinner(plotlyOutput("PCA_samples_plot_sex")))),
                   column(3, div(class = "custom-panel3",
                                  h5('Attribute: sampling location'),
                                  withSpinner(plotlyOutput("PCA_samples_plot_sampling_location")))),
                 )
        ),
        tabPanel("Cerebellum tissue analysis",
                 fluidRow(
                 column(4,div(class = "custom-panel",
                                h3("Volcano plot - female Vs male"),
                                withSpinner(plotOutput("DEA_sex_cerebellum")),
                              p('Excererbated log2 Fold Changes when comparing Female vs Male. The central lump of genes that are neither biologically or statistically significantly diferentially expressed gets squashed
                          toward the x axis. Insanely high -log10 adjusted p values are achieved by the unmistakable difference in some sex chromosome related gene titers.')),
                        ),
                 column(4,div(class = "custom-panel",
                                h3("Volcano plot - ApoE3,3 Vs ApoE4,4"),
                                withSpinner(plotOutput("DEA_genotype_cerebellum")))
                        ),
                 column(4,div(class = "custom-panel",
                                h3("Volcano plot - Alzheimer's Vs none"),
                                withSpinner(plotOutput("DEA_disease_cerebellum")))
                        ),
                 )
        ),
                 
        tabPanel("Prefrontal cortex tissue analysis",
                 fluidRow(
                   column(4,div(class = "custom-panel",
                                h3("Volcano plot - female Vs male"),
                                withSpinner(plotOutput("DEA_sex_pfc")))
                          ),
                   column(4,div(class = "custom-panel",
                                h3("Volcano plot - ApoE3,3 Vs ApoE4,4"),
                                withSpinner(plotOutput("DEA_genotype_pfc")))
                          ),
                   column(4,div(class = "custom-panel",
                                h3("Volcano plot - Alzheimer's Vs none"),
                                withSpinner(plotOutput("DEA_disease_pfc")))
                          ),
                 )
        ),
      )
    )
  )
)

server <- function(input, output) {
  
  # loading prewritten data to reduce on the fly server calculations
  meta_data_sorted  <- read.csv('meta_data_sorted.csv')
  rownames(meta_data_sorted)<-meta_data_sorted$samples
  
  read_counts_clean <- read.csv('read_counts_clean.csv')
  rownames(read_counts_clean) <- read_counts_clean$X
  read_counts_clean<-read_counts_clean %>%
    select(-X)
  rnaseq_data<-t(read_counts_clean)
  
  pca_result <<- prcomp(rnaseq_data)
  
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
  
  output$PCA_samples_plot_sex<-renderPlotly({
    p<-PCA_samples_plot(pca_result,meta_data_sorted,'sex')
    p
  })
  output$PCA_samples_plot_genotype<-renderPlotly({
    p<-PCA_samples_plot(pca_result,meta_data_sorted,'genotype')
    ggplotly(p)
  })
  output$PCA_samples_plot_disease<-renderPlotly({
    p<-PCA_samples_plot(pca_result,meta_data_sorted,'disease')
    ggplotly(p)
  })
  output$PCA_samples_plot_sampling_location<-renderPlotly({
    p<-PCA_samples_plot(pca_result,meta_data_sorted,'sampling_location')
    ggplotly(p)
  })
}

shinyApp(ui = ui, server = server)