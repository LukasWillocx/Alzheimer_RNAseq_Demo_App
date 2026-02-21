# Horizontal Navbar App - RNASeq Analysis

# List of packages
all_packages <- c(
  "shiny", "shinycssloaders", "networkD3", "dplyr", "tidyr", "ggplot2",
  "grid", "htmlwidgets", "AnnotationHub", "DESeq2", "bslib", "plotly", "DT",
  "clusterProfiler", "org.Hs.eg.db", "AnnotationDbi", "biomaRt",
  "stringr", "luwitemplate"
)

packages <- c(
  "shiny", "shinycssloaders", "networkD3", "dplyr", "tidyr", "ggplot2",
  "grid", "htmlwidgets", "AnnotationHub", "DESeq2", "bslib", "plotly", "DT",
  "clusterProfiler", "org.Hs.eg.db", "AnnotationDbi", "biomaRt",
  "stringr", "luwitemplate"
)

for(pkg in packages) {
  library(pkg, character.only = TRUE)
}

source('app_functions.R')

# UI
ui <- bslib::page_navbar(
  title = "RNASeq Analysis",
  theme = my_theme(),
  dark_mode_css(),
  window_title = "Alzheimer RNASeq Demo",
  tags$head(tags$style(HTML("
    #sankey svg { width: 100% !important; height: 100% !important; }
  "))),
  
  # Overview Tab
  nav_panel(
    "Overview",
    layout_sidebar(
      sidebar = sidebar(
        width = 590,
        tags$div(
          tags$div(style = "text-align: center", h4('Dataset Information')),
          hr(),
          HTML(paste('<strong>Accession number:</strong> <a href="https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-14167" target="_blank">E-MTAB-14167</a>', '<br>')),
          HTML(paste("<strong>Format:</strong>", '.txt', '<br>')),
          HTML(paste("<strong>Assay technology:</strong>", 'RNA sequencing', '<br>')),
          HTML(paste("<strong>Sample Size:</strong> ", '91', '<br>')),
          HTML(paste("<strong>Number of Observations:</strong> ", '62705', '<br>')),
          HTML(paste('<strong>Reference Genome:</strong> <a href="https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/" target="_blank">GRCh38.p14</a>', '<br>')),
          HTML(paste('<strong>Sample Attributes:</strong> <a href="https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-14167/sdrf" target="_blank">SDRF</a>', '.txt', '<br>')),
          HTML(paste("<strong>File Size:</strong>", '53.5 MB', '<br>')),
          HTML(paste("<strong>Downloaded on:</strong>", 'February 18th, 2025', '<br>'))
        ),
        tags$div(
          tags$div(style = "text-align: center", h5('R Libraries')),
          hr(),
          tags$ul(
            style = "list-style-type: none; padding-left: 0;",
            lapply(all_packages, function(lib) {
              if (lib == "luwitemplate") {
                doc_url <- "https://github.com/LukasWillocx/luwitemplate"
              } else if (lib %in% c("AnnotationHub", "DESeq2", "clusterProfiler",
                                    "org.Hs.eg.db", "AnnotationDbi", "biomaRt")) {
                doc_url <- paste0("https://bioconductor.org/packages/release/bioc/html/", lib, ".html")
              } else {
                doc_url <- paste0("https://CRAN.R-project.org/package=", lib)
              }
              
              tags$li(
                style = "display: inline-block; margin: 5px; text-align: center;",
                tags$a(
                  href = doc_url,
                  target = "_blank",
                  style = "text-decoration: none;",
                  tags$span(
                    style = "border-radius: 15px; background-color: rgba(128, 128, 128, 0.1);
                    border: 1px solid rgba(128, 128, 128, 0.3); padding: 5px 10px;
                    display: inline-block; cursor: pointer;
                    transition: background-color 0.2s, border-color 0.2s;",
                    onmouseover = "this.style.backgroundColor='rgba(128, 128, 128, 0.2)'; this.style.borderColor='rgba(128, 128, 128, 0.5)';",
                    onmouseout = "this.style.backgroundColor='rgba(128, 128, 128, 0.1)'; this.style.borderColor='rgba(128, 128, 128, 0.3)';",
                    tags$b(lib)
                  )
                )
              )
            })
          )
        ),
        tags$div(
          tags$div(style = "text-align: center", h5('About the Application')),
          hr(),
          p('This application features a simple RNASeq analysis. The dataset consists of bulk RNASeq of postmortem human prefrontal cortex and cerebellum brain
            tissues. Differential gene expression is shown, based on sample attributes.
            These are sex, genotype and disease status. Apolipoprotein E4
            (ApoE4) allele is a risk factor for Alzheimer\'s disease. The comparison of a
            homozygous risk factor group vs a homozygous wild-type group is therefore an
            interesting consideration. It was opted to not pool samples from both tissues, as
            this could introduce biological tissue variance, which is somewhat outlined by the principal
            component analysis.'),
        )
      ),
      layout_columns(
        col_widths = c(12, 12, 3, 3, 3, 3),
        fill = TRUE,
        card(
          full_screen = TRUE,
          card_header("Sample Attributes"),
          card_body(
            tags$div(style = "width: 100%; height: 100%;",
                     withSpinner(sankeyNetworkOutput("sankey", width = '100%', height = '360px')))
          )
        ),
        card(
          full_screen = TRUE,
          p("The samples in the 2D principal component vector space are labeled according to their
           categorical attributes. The first two principal components comprise nearly 97%
           of the variance observed in the data. Separation along PC2 seems to represent the
           biological difference between the two sampling tissues.")
        ),
        card(
          full_screen = TRUE,
          card_header("Attribute: disease"),
          card_body(withSpinner(plotlyOutput("PCA_samples_plot_disease", height = "300px")))
        ),
        card(
          full_screen = TRUE,
          card_header("Attribute: genotype"),
          card_body(withSpinner(plotlyOutput("PCA_samples_plot_genotype", height = "300px")))
        ),
        card(
          full_screen = TRUE,
          card_header("Attribute: sex"),
          card_body(withSpinner(plotlyOutput("PCA_samples_plot_sex", height = "300px")))
        ),
        card(
          full_screen = TRUE,
          card_header("Attribute: sampling location"),
          card_body(withSpinner(plotlyOutput("PCA_samples_plot_sampling_location", height = "300px")))
        )
      )
    )
  ),
  
  # DGEA Tab
  nav_panel(
    "DGEA",
    layout_sidebar(
      sidebar = sidebar(
        width = 250,
        h4("Navigation"),
        selectInput("attribute", "Choose the sample attribute for DGEA",
                    choices = c("Sex   (Female Vs Male)"            = "sex",
                                "Genotype   (ApoE3,3 Vs ApoE4,4)"  = "genotype",
                                "Disease   (Alzheimer's Vs none)"   = "disease"),
                    selected = 'genotype'
        ),
        selectInput("tissue", "Choose tissue sampling location",
                    choices = c("Cerebellum"        = "cerebellum",
                                "Prefrontal cortex" = "prefrontal cortex")),
        hr(),
      ),
      layout_columns(
        col_widths = c(7, 5, 12),
        fill = FALSE,
        card(
          full_screen = TRUE,
          card_header("Volcano plot — biologically (|log2 FC| > 2) and statistically (p.adj < 0.05) relevant differential gene expression"),
          card_body(withSpinner(plotlyOutput("volcano_plot", height = "400px")))
        ),
        card(
          full_screen = TRUE,
          card_header("Gene ontology enrichment"),
          card_body(withSpinner(uiOutput("GOE_plot_or_message", height= '400px')))
        ),
        card(
          full_screen = TRUE,
          card_header("Significantly differentially expressed genes"),
          card_body(withSpinner(dataTableOutput("signi_gene_table")))
        ),
        
      )
    )
  ),
  
  # EDA Code Tab
  nav_panel(
    "EDA Code",
    card(
      card_body(
        includeMarkdown('www/ApoE.Rmd')
      )
    )
  ),
  
  # Functions Code Tab
  nav_panel(
    "Functions Code",
    card(
      card_body(
        includeMarkdown('www/app_functions.Rmd')
      )
    )
  ),
  
  nav_spacer(),
  nav_item(input_dark_mode(id = "dark_mode")),
)

# Server
server <- function(input, output, session) {
  
  colors <- get_theme_colors()
  dm <- use_dark_mode(input, session)
  
  # Load pre-written data
  meta_data_sorted  <- read.csv('data/meta_data_sorted.csv')
  rownames(meta_data_sorted) <- meta_data_sorted$samples
  
  read_counts_clean <- read.csv('data/read_counts_clean.csv')
  rownames(read_counts_clean) <- read_counts_clean$X
  read_counts_clean <- read_counts_clean %>% dplyr::select(-X)
  rnaseq_data <- t(read_counts_clean)
  
  # Cache PCA
  pca_result <<- prcomp(rnaseq_data)
  
  # Pre-computed DEA results
  DEA_sex_cerebellum_data      <- read.csv('data/DEA_sex_results_in_cerebellum.csv')
  DEA_sex_pfc_data             <- read.csv('data/DEA_sex_results_in_pfc.csv')
  DEA_genotype_cerebellum_data <- read.csv('data/DEA_genotype_results_in_cerebellum.csv')
  DEA_genotype_pfc_data        <- read.csv('data/DEA_genotype_results_in_pfc.csv')
  DEA_disease_cerebellum_data  <- read.csv('data/DEA_disease_results_in_cerebellum.csv')
  DEA_disease_pfc_data         <- read.csv('data/DEA_disease_results_in_pfc.csv')
  
  # Pre-computed GOE results
  GOE_disease_cerebellum  <- read.csv('data/GOE_disease_cerebellum.csv')
  GOE_disease_pfc         <- read.csv('data/GOE_disease_pfc.csv')
  GOE_sex_cerebellum      <- read.csv('data/GOE_sex_cerebellum.csv')
  GOE_sex_pfc             <- read.csv('data/GOE_sex_pfc.csv')
  GOE_genotype_cerebellum <- read.csv('data/GOE_genotype_cerebellum.csv')
  
  # ENSEMBL -> Entrez mapping
  entrez_mapping <- read.csv('data/entrez_mapping.csv')
  entrez_mapping <- rename(entrez_mapping, ENSEMBL = ensembl_gene_id)
  
  # ── Overview panel ─────────────────────────────────────────────────────────
  output$sankey <- renderSankeyNetwork({
    sankey_design(meta_data_sorted, theme = dm$theme())
  })
  
  output$PCA_samples_plot_sex <- renderPlotly({
    PCA_samples_plot(pca_result, meta_data_sorted, 'sex', dm$theme())
  })
  output$PCA_samples_plot_genotype <- renderPlotly({
    PCA_samples_plot(pca_result, meta_data_sorted, 'genotype', dm$theme())
  })
  output$PCA_samples_plot_disease <- renderPlotly({
    PCA_samples_plot(pca_result, meta_data_sorted, 'disease', dm$theme())
  })
  output$PCA_samples_plot_sampling_location <- renderPlotly({
    PCA_samples_plot(pca_result, meta_data_sorted, 'sampling_location', dm$theme())
  })
  
  # ── DGEA panel ─────────────────────────────────────────────────────────────
  select_DEA_data <- reactive({
    attribute <- input$attribute
    tissue    <- input$tissue
    
    if (attribute == "sex") {
      if (tissue == "cerebellum")        return(DEA_sex_cerebellum_data)
      if (tissue == "prefrontal cortex") return(DEA_sex_pfc_data)
    } else if (attribute == "genotype") {
      if (tissue == "cerebellum")        return(DEA_genotype_cerebellum_data)
      if (tissue == "prefrontal cortex") return(DEA_genotype_pfc_data)
    } else if (attribute == "disease") {
      if (tissue == "cerebellum")        return(DEA_disease_cerebellum_data)
      if (tissue == "prefrontal cortex") return(DEA_disease_pfc_data)
    }
  })
  
  select_GOE_data <- reactive({
    attribute <- input$attribute
    tissue    <- input$tissue
    
    if (attribute == "sex") {
      if (tissue == "cerebellum")        return(GOE_sex_cerebellum)
      if (tissue == "prefrontal cortex") return(GOE_sex_pfc)
    } else if (attribute == "genotype") {
      if (tissue == "cerebellum")        return(GOE_genotype_cerebellum)
      if (tissue == "prefrontal cortex") return(0)
    } else if (attribute == "disease") {
      if (tissue == "cerebellum")        return(GOE_disease_cerebellum)
      if (tissue == "prefrontal cortex") return(GOE_disease_pfc)
    }
  })
  
  output$volcano_plot <- renderPlotly({
    p <- DEA_volcano_plotter(select_DEA_data(), dm$theme())
    fig <- luwi_ggplotly(p, theme = dm$theme(), tooltip = c("text", "x", "y")) %>%
      layout(showlegend = FALSE)
    # Strip 'hoveron' — not supported by scattergl
    for (i in seq_along(fig$x$data)) {
      fig$x$data[[i]]$hoveron <- NULL
    }
    fig %>% plotly::toWebGL()
  })
  
  select_signi_genes <- reactive({
    select_DEA_data() %>%
      filter(padj < 0.05 & abs(log2FoldChange) > 2) %>%
      dplyr::select(c(ENSEMBL, Chr, log2FoldChange, padj)) %>%
      arrange(padj) %>%
      left_join(entrez_mapping, by = 'ENSEMBL')
  })
  
  output$signi_gene_table <- renderDataTable({
    datatable(select_signi_genes())
  })
  
  output$GOE_plot_or_message <- renderUI({
    if (nrow(select_signi_genes()) < 10) {
      textOutput("GOE_message")
    } else {
      plotlyOutput("GOE_plot")
    }
  })
  
  output$GOE_message <- renderText({
    "Insufficient amount of significant genes present for Gene Ontology Enrichment.
     The minimal amount of significant genes is internally set at ten.
     Choose another combination of attribute or tissue.
     All cerebellum analyses should yield results."
  })
  
  output$GOE_plot <- renderPlotly({
    p <- GOE_plotter(select_signi_genes(), select_GOE_data(), dm$theme())
    luwi_ggplotly(p, theme = dm$theme(), tooltip = c("ONTOLOGY", "y"))
  })
}

# Run the application
shinyApp(ui = ui, server = server)