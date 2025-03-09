# List of required packages
packages <- c("shiny", "shinythemes", "shinycssloaders", 
              "networkD3", "dplyr", "tidyr", "ggplot2",
              "grid", "htmlwidgets","AnnotationHub", "DESeq2","plotly", "DT",
              "clusterProfiler", "org.Hs.eg.db", "AnnotationDbi","biomaRt",
              "stringr")

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
    div("RNASeq analysis", 
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
      h4('Navigation - DGEA panel'),
      selectInput("attribute", "Choose the sample attribute for DGEA",
                  choices = c("Sex   (Female Vs Male)" = "sex", 
                              "Genotype   (ApoE3,3 Vs ApoE4,4)" = "genotype", 
                              "Disease   (Alzheimer's Vs none)" = "disease")),
       
      selectInput("tissue", "Choose tissue sampling location",
                  choices = c("Cerebellum" = "cerebellum", 
                              "Prefrontal cortex" = "prefrontal cortex")),
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
      p('This application concerns a demo, featuring a simple RNASeq analysis. The dataset consists of bulk RNASeq of postmortem human prefrontal cortex and cerebellum brain tissues.
        Differential gene expression is scrutinized, based on sample attributes. These are sex, genotype and disease status. The former is particularly interesting to evaluate whether the
        analysis is sound. In a comparison between male and female, a vast difference in expression at the level of the sex chromosomes should become apparent. In the other application it was
        outlined that the apolipoprotein E4 (ApoE4) allele was a risk factor for Alzheimer\'s disease. The comparison of a homozygous risk factor group vs a homozygous wild-type group is therefore
        an interesting consideration. It was opted to not pool samples from both tissues, as this could introduce biological tissue variance, which is outlined by the principal component analysis.'),
      hr(),
      tags$div(style = "text-align: center",h5('Known issues')),
      hr(),
      p('The Sankey diagram in the Overview panel currently doesn\'t scale appropriately to the full size of the box in a firefox browser.
        Chromium-based browsers should display everything as intended.')
    ),
    
    mainPanel(
      width = 9,
      tabsetPanel(
        type = "tabs",
        tabPanel("Overview",
                 fluidRow(
                   column(12, div(class = "custom-panel",
                                 h3("Sample attributes"),
                                 tags$div(style = "width: 100%; height: 100%;",
                                 withSpinner(sankeyNetworkOutput("sankey",width='100%',height='500px'))))),
                 ),
                 fluidRow(
                   column(12, div(class = "custom-panel",
                                  h3("Samples in principal component vector space, labeled according to their categorical attributes"),
                                  p('The principal component space entails the dimensional reduction of the 62705 distinct RNASeq reads in the data set. The first two
                                  principal components comprise nearly 97% of the variance observed in the data. Upon closer inspection it can be argued that the separation along principal component 2 quite neatly
                                  represents the biological difference between the two sampling tissues of cerebellum and prefrontal cortex. None of the other attributes show
                                  clear separation along either principal component. 
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
        tabPanel("Differential Gene Expression Analysis (DGEA)",
                 fluidRow(
                 column(12,div(class = "custom-panel",
                                h3("Volcano plot - identifying biologically (|log2 FoldChange| > 2) and statistically (p.adj < 0.05) relevant differential gene expression"),
                                withSpinner(plotOutput("volcano_plot")),
                        ),
                 )),
                 fluidRow(
                 column(12,div(class = "custom-panel",
                               h3("Significantly differentially expressed genes"),
                               withSpinner(dataTableOutput("signi_gene_table")),
                        ),
                 ),
                 column(12,div(class = "custom-panel",
                               h3("Gene ontology enrichment"),
                               withSpinner(uiOutput("GOE_plot_or_message")),
                 ),
                 ),
        ),
      ),
      tabPanel('EDA code',
               includeMarkdown('www/ApoE.Rmd'),
    ),
      tabPanel('Functions code',
               includeMarkdown('www/app_functions.Rmd'),
      )
  )
)
)
)

server <- function(input, output) {
  
  # loading pre-written data to reduce on the fly server calculations
  meta_data_sorted  <- read.csv('data/meta_data_sorted.csv')
  rownames(meta_data_sorted)<-meta_data_sorted$samples
  
  read_counts_clean <- read.csv('data/read_counts_clean.csv')
  rownames(read_counts_clean) <- read_counts_clean$X
  read_counts_clean<-read_counts_clean %>%
    dplyr::select(-X)
  rnaseq_data<-t(read_counts_clean)
  
  # caching principal component analysis
  pca_result <<- prcomp(rnaseq_data)
  
  # pre-written DEA results from the resource intensive DESeq2 function
  DEA_sex_cerebellum_data<-read.csv('data/DEA_sex_results_in_cerebellum.csv')
  DEA_sex_pfc_data<-read.csv('data/DEA_sex_results_in_pfc.csv')
  
  DEA_genotype_cerebellum_data<-read.csv('data/DEA_genotype_results_in_cerebellum.csv')
  DEA_genotype_pfc_data<-read.csv('data/DEA_genotype_results_in_pfc.csv')
  
  DEA_disease_cerebellum_data<-read.csv('data/DEA_disease_results_in_cerebellum.csv')
  DEA_disease_pfc_data<-read.csv('data/DEA_disease_results_in_pfc.csv')
  
  # pre-written GOE results from the resource intensive enrichGO function
  GOE_disease_cerebellum <- read.csv('data/GOE_disease_cerebellum.csv')
  GOE_disease_pfc <- read.csv('data/GOE_disease_pfc.csv')
  
  GOE_sex_cerebellum <- read.csv('data/GOE_sex_cerebellum.csv')
  GOE_sex_pfc <- read.csv('data/GOE_sex_pfc.csv')
  
  GOE_genotype_cerebellum <- read.csv('data/GOE_genotype_cerebellum.csv')
  
  # load the dataframe that converts ENSEMBL ids to Entrez ids (fetch takes 3 minutes otherwise)
  entrez_mapping<-read.csv('data/entrez_mapping.csv') 
  entrez_mapping<-rename(entrez_mapping,ENSEMBL=ensembl_gene_id)
  
  # Overview panel (top)
  output$sankey<-renderSankeyNetwork({
    sankey_design(meta_data_sorted)
  })
  
  # Overview panel (bottom)
  output$PCA_samples_plot_sex<-renderPlotly({
    PCA_samples_plot(pca_result,meta_data_sorted,'sex')
  })
  output$PCA_samples_plot_genotype<-renderPlotly({
    PCA_samples_plot(pca_result,meta_data_sorted,'genotype')
  })
  output$PCA_samples_plot_disease<-renderPlotly({
    PCA_samples_plot(pca_result,meta_data_sorted,'disease')
  })
  output$PCA_samples_plot_sampling_location<-renderPlotly({
    PCA_samples_plot(pca_result,meta_data_sorted,'sampling_location')
  })
  
  # DGEA panel top
  select_DEA_data <- reactive({
    
    attribute <- input$attribute
    tissue <- input$tissue
    
    if (attribute == "sex") {
      if (tissue == "cerebellum") {
        return(DEA_sex_cerebellum_data)
      } else if (tissue == "prefrontal cortex") {
        return(DEA_sex_pfc_data)
      }
    } else if (attribute == "genotype") {
      if (tissue == "cerebellum") {
        return(DEA_genotype_cerebellum_data)
      } else if (tissue == "prefrontal cortex") {
        return(DEA_genotype_pfc_data)
      }
    } else if (attribute == "disease") {
      if (tissue == "cerebellum") {
        return(DEA_disease_cerebellum_data)
      } else if (tissue == "prefrontal cortex") {
        return(DEA_disease_pfc_data)
      }
    }
  })
  
  select_GOE_data <- reactive({
    
    attribute <- input$attribute
    tissue <- input$tissue
    
    if (attribute == "sex") {
      if (tissue == "cerebellum") {
        return(GOE_sex_cerebellum)
      } else if (tissue == "prefrontal cortex") {
        return(GOE_sex_pfc)
      }
    } else if (attribute == "genotype") {
      if (tissue == "cerebellum") {
        return(GOE_genotype_cerebellum)
      } else if (tissue == "prefrontal cortex") {
        return(0)
      }
    } else if (attribute == "disease") {
      if (tissue == "cerebellum") {
        return(GOE_disease_cerebellum)
      } else if (tissue == "prefrontal cortex") {
        return(GOE_disease_pfc)
      }
    }
  })
  
  output$volcano_plot <- renderPlot({
    DEA_volcano_plotter(select_DEA_data())
  },bg='transparent')
  
  select_signi_genes <- reactive ({
    select_DEA_data()%>%
      filter(padj < 0.05 & abs(log2FoldChange)>2)%>%
      dplyr::select(c(ENSEMBL,Chr,log2FoldChange,padj))%>%
      arrange(padj) %>%
      left_join(entrez_mapping,by='ENSEMBL')
  })
  
  output$signi_gene_table <- renderDataTable({
    datatable(select_signi_genes()) %>% 
      formatStyle(columns = names(select_signi_genes()), color = "#7aa6a1") # Change 'blue' to your desired color
  })
  
  # Outputting the gene ontology bar graphs (conditional on there being enough significant genes)
  output$GOE_plot_or_message <- renderUI({
    if (nrow(select_signi_genes()) < 10) {
      textOutput("GOE_message")
    } else {
      plotlyOutput("GOE_plot")
    }
  })
  
  output$GOE_message <- renderText ({
    "Insufficient amount of significant genes present for Gene Ontology Enrichment.
    The minmal amount of siginificant genes is internally set at ten. \n Choose another combination of attribute or tissue.
    All cerebellum analyses should yield results."
  })
  
  output$GOE_plot <- renderPlotly({
    ggplotly(GOE_plotter(select_signi_genes(),select_GOE_data()),tooltip = c("ONTOLOGY", "y"))
  })
}

shinyApp(ui = ui, server = server)