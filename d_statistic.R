d_statistic = function(source_graph, sample_graph, source_closeness, source_betweenness){
  #Degree-All Distribution
  source_graph_dd <- degree.distribution(source_graph, cumulative = T, mode="all")
  sample_graph_dd <- degree.distribution(sample_graph, cumulative = T, mode="all")
  test_degree_all <- ks.test(source_graph_dd, sample_graph_dd)
  
  #Degree-In Distribution
  source_dd_in <- degree.distribution(source_graph, cumulative = T, mode="in")
  sample_dd_in <- degree.distribution(sample_graph, cumulative = T, mode="in")
  test_degree_in <- ks.test(source_dd_in, sample_dd_in)
  
  #Degree-Out Distribution
  source_dd_out <- degree.distribution(source_graph, cumulative = T, mode="out")
  sample_dd_out <- degree.distribution(sample_graph, cumulative = T, mode="out")
  test_degree_out <- ks.test(source_dd_out, sample_dd_out)
  
  #Clustering Coefficient Distribution
  source_cd <- transitivity(source_graph, type="local")
  sample_cd <- transitivity(sample_graph, type="local")
  test_cd <- ks.test(source_cd, sample_cd)
  
  #Closeness
  sample_closeness <- closeness(sample_graph)
  test_closeness <- ks.test(scale(source_closeness), scale(sample_closeness))
  
  #Betweenness
  sample_betweenness <- betweenness(sample_graph)
  test_betweenness <- ks.test(scale(source_betweenness), scale(sample_betweenness))
  
  values = list("degree_in" = test_degree_in$statistic, "degree_out" = test_degree_out$statistic,
                "degree_all" = test_degree_all$statistic, "clustering_coeff" = test_cd$statistic,
                "closeness" = test_closeness$statistic, "betweenness" = test_betweenness$statistic)
  return(values)
}