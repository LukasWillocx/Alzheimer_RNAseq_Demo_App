
sankey_design<-function(df){

  df <- df %>%
    mutate(age_group = cut(age,
                           breaks = c(60, 70, 80, 90, 100),
                           labels = c("60-69", "70-79", "80-89", "90+"),
                           include.lowest = TRUE))

  # Create ordered nodes (5 levels)
  nodes <- data.frame(
    name = c(
      as.character(unique(df$disease)),       # Level 1: Disease
      as.character(unique(df$genotype)),      # Level 2: Genotype
      as.character(unique(df$sex)),           # Level 3: Sex
      as.character(unique(df$age_group)),     # Level 4: Age Group
      as.character(unique(df$sampling_location))  # Level 5: Location
    )
  )
  # Create connection layers
  links <- bind_rows(
    # Disease → Genotype 
    df %>% 
      group_by(disease, genotype) %>% 
      summarise(value = n(), .groups = "drop") %>% 
      mutate(source = match(disease, nodes$name) - 1,
             target = match(genotype, nodes$name) - 1),
    
    # Genotype → Sex
    df %>% 
      group_by(genotype, sex) %>% 
      summarise(value = n(), .groups = "drop") %>% 
      mutate(source = match(genotype, nodes$name) - 1,
             target = match(sex, nodes$name) - 1),
    
    # Sex → Age Group
    df %>% 
      group_by(sex, age_group) %>% 
      summarise(value = n(), .groups = "drop") %>% 
      mutate(source = match(sex, nodes$name) - 1,
             target = match(age_group, nodes$name) - 1),
    
    # Age Group → Location
    df %>% 
      group_by(age_group, sampling_location) %>% 
      summarise(value = n(), .groups = "drop") %>% 
      mutate(source = match(age_group, nodes$name) - 1,
             target = match(sampling_location, nodes$name) - 1)
  ) %>% 
    dplyr::select(source, target, value)
  
  
  # Generate Sankey diagram
  sankeyPlot<-sankeyNetwork(
    Links = links,
    Nodes = nodes,
    Source = "source",
    Target = "target",
    Value = "value",
    NodeID = "name",
    units = "Samples",
    fontFamily = "arial",
    fontSize = 15,
    nodeWidth = 25,
    nodePadding = 10,
    sinksRight = T,
    iterations = 5
    )
  
  
  customSankey <- sankeyPlot %>%
    htmlwidgets::onRender("
    function(el) {
      // Define custom colors for bands
      const customColors = ['#d4b996'];

      // Select all band elements and apply custom colors
      d3.selectAll('path.link')
        .style('stroke', function(d, i) {
          return customColors[i % customColors.length];
        });
      // Add column headers below level blocks
      var cols_x = this.sankey.nodes().map(d => d.x).filter((v, i, a) => a.indexOf(v) === i).sort(function(a, b){return a - b});
      var labels = ['disease', 'genotype', 'sex', 'age', 'location'];
      var svgHeight = el.querySelector('svg').getBoundingClientRect().height;

      cols_x.forEach((d, i) => {
        if (i < labels.length) { // Ensure we don't go out of bounds
          d3.select(el).select('svg')
            .append('text')
            .attr('x', d + 31.8)
            .attr('y', svgHeight - 12) // Adjust the y position to be below the diagram
            .attr('dy', '.35em') // Adjust vertical alignment if needed
            .attr('text-anchor', 'middle') // Center the text horizontally
            .style('font-size', '16px') // Adjust font size if necessary
            .text(labels[i]);
        }
      });
      // Set the background color to transparent
      el.style.backgroundColor = 'transparent';
      // Change the text color to blue
      el.querySelectorAll('text').forEach(function(node) {
        node.style.fill = '#7aa6a1';
      });
      // Define custom colors for nodes
      const customNodeColors = 
      ['blue', 'red','lightgrey', 'grey','black', 'lightblue','pink', 'darkgray','gray', 'lightgrey','gainsboro', 'DarkSalmon','Teal'];
      // Select all node elements and apply custom colors
      d3.selectAll('g.node')
        .select('rect')
        .style('fill', function(d, i) {
          return customNodeColors[i % customNodeColors.length];
      });
      d3.selectAll('g.node')
        .select('rect')
        .style('stroke', function(d, i) {
          return customNodeColors[i % customNodeColors.length];
      });
      
    }"
    )
  
  # Display the custom Sankey diagram 
  htmltools::browsable(customSankey)
}

DEA_volcano_plotter <- function(DEA_res){
  ggplot(as.data.frame(DEA_res), aes(x=log2FoldChange, y=-log10(padj))) +
    geom_point(aes(color=stat_sign), alpha=0.7) +
    theme_minimal() +
    theme(legend.position = 'None',
          plot.background = element_rect(fill = "transparent", colour = NA),
          panel.background = element_rect(fill = "transparent", colour = NA),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_line(color = "#BEBEBE"),
          panel.grid.major.y = element_line(color = "#BEBEBE"),
          axis.title = element_text(color = "#7aa6a1"),
          axis.text = element_text(color = "#7aa6a1")) +
    labs(x="Log2 Fold Change", y="-Log10 Adjusted P-value") +
    geom_hline(yintercept = -log10(0.05), linetype="dashed", color="#7aa6a1") + # Statistical significance
    geom_vline(xintercept=c(-2, 2), linetype="dashed", color="#7aa6a1") + # Biological relevance lines
    scale_color_manual(values=c('TRUE'='#7aa6a1', 'FALSE'='black'))
}
