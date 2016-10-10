evaluate_d = function(files, source_graph, source_closeness, source_betweenness){
  degree_in <- degree_out <- degree_all <- clustering_coeff <- closeness <- betweenness <- vector()
  for(file in files){
    sample_graph <- readRDS(file)
    statistic <- d_statistic(source_graph, sample_graph, source_closeness, source_betweenness)
    degree_in <- append(degree_in, as.numeric(statistic$degree_in))
    degree_out <- append(degree_out, as.numeric(statistic$degree_out))
    degree_all <- append(degree_all, as.numeric(statistic$degree_all))
    clustering_coeff <- append(clustering_coeff, as.numeric(statistic$clustering_coeff))
    closeness <- append(closeness, as.numeric(statistic$closeness))
    betweenness <- append(betweenness, as.numeric(statistic$betweenness))
  }
  values <- list("degree_in" = degree_in, "degree_in_mean" = round(mean(degree_in), digits=3),
                 "degree_out" = degree_out, "degree_out_mean" = round(mean(degree_out), digits=3),
                 "degree_all" = degree_all, "degree_all_mean" = round(mean(degree_all), digits=3),
                 "clustering_coeff" = clustering_coeff, "clustering_coeff_mean" = round(mean(clustering_coeff), digits=3),
                 "closeness" = closeness, "closeness_mean" = round(mean(closeness), digits=3),
                 "betweenness" = betweenness, "betweenness_mean" = round(mean(betweenness), digits=3))
  return(values)
}

evaluate_random_node_sample = function(graph, min_sample = 0.15, max_sample = 0.5, steps = 5){
  files <- vector()
  for(i in seq(min_sample, max_sample, length.out=steps)){
    sample_len <- i * length(V(graph))
    nodes_vector <- sample(V(graph)$name, replace=FALSE, size = sample_len)
    file_string = paste("RN_sample_",round(i*100,digits=2),".rds",sep="")
    files <- append(files, file_string)
    saveRDS(induced.subgraph(graph=graph, vids=nodes_vector), file=file_string)
    print(paste("Random Node",round(i*100, digits=2), "% done"))
  }
  return(files)
}

evaluate_random_page_rank_node_sample = function(graph, min_sample = 0.15, max_sample = 0.5, steps = 5){
  files <- vector()
  n <- length(V(graph))
  for(i in seq(min_sample, max_sample, length.out=steps)){
    sample_len <- i * n
    page_rank <- page.rank(graph, vids = V(graph))
    nodes_vector <- sample(V(graph)$name, prob = page_rank$vector, replace=FALSE, size = sample_len)
    file_string = paste("RPRN_sample_",round(i*100,digits=2),".rds",sep="")
    files <- append(files, file_string)
    saveRDS(induced.subgraph(graph=graph, vids=nodes_vector), file=file_string)
    print(paste("Random Page Rank Node",round(i*100, digits=2), "% done"))
  }
  return(files)
}

evaluate_random_degree_node_sample = function(graph, min_sample = 0.15, max_sample = 0.5, steps = 5){
  files <- vector()
  n <- length(V(graph))
  for(i in seq(min_sample, max_sample, length.out=steps)){
    sample_len <- i * n
    node_degrees <- degree(graph)
    node_degree_prob <- node_degrees / length(V(graph))
    nodes_vector <- sample(V(graph)$name, prob = node_degree_prob, replace=FALSE, size = sample_len)
    file_string = paste("RDN_sample_",round(i*100,digits=2),".rds",sep="")
    files <- append(files, file_string)
    saveRDS(induced.subgraph(graph=graph, vids=nodes_vector), file=file_string)
    print(paste("Random Degree Node",round(i*100, digits=2), "% done"))
  }
  return(files)
}

evaluate_random_edge_sample = function(graph, min_sample = 0.15, max_sample = 0.5, steps = 5) {
  files <- vector()
  n <- length(V(graph))
  current_graph <- graph.empty()
  m <- 0
  E(graph)$id <- seq_len(ecount(graph))
  for(i in seq(min_sample, max_sample, length.out=steps)){
    while(length(V(current_graph)) < n*i){
      m = m + 1
      random_edge <- sample(E(graph)$id, size=m)
      current_graph <- subgraph.edges(graph, random_edge)
    }
    file_string = paste("RE_sample_",round(i*100,digits=2),".rds",sep="")
    files <- append(files, file_string)
    saveRDS(current_graph, file=file_string)
    print(paste("Random Edge",round(i*100, digits=2), "% done"))
  }
  return(files)
}

evaluate_random_node_neighbour_sample = function(graph, min_sample = 0.15, max_sample = 0.5, steps = 5){
  n <- length(V(graph))
  list_to_create_graph <- files <- vector()
  for(i in seq(min_sample, max_sample, length.out=steps)){
    while(length(list_to_create_graph) < i * n){
      current_node <- sample(V(graph)$name, 1)
      node_neighbourhood <- neighbors(graph, v = current_node, mode="out")
      list_to_create_graph <- append(list_to_create_graph, current_node)
      list_to_create_graph <- append(list_to_create_graph, node_neighbourhood$name)
      list_to_create_graph <- unique(list_to_create_graph)
    }
    file_string = paste("RNN_sample_",round(i*100,digits=2),".rds",sep="")
    files <- append(files, file_string)
    saveRDS(induced_subgraph(graph, v = list_to_create_graph), file=file_string)
    print(paste("Random Node Neighbour",round(i*100, digits=2), "% done"))
  }
  return(files)
}

evaluate_random_walk_sample = function(graph, min_sample = 0.15, max_sample = 0.5, steps = 5) {
  n <- length(V(graph))
  list_to_create_graph <- files <- vector()
  source_node <- sample(V(graph), 1)
  step_count <- 0
  current_node <- source_node
  for(i in seq(min_sample, max_sample, length.out=steps)){
    while(length(list_to_create_graph) < n * i){
      step_count <- step_count + 1
      if(!(current_node$name %in% list_to_create_graph)){
        list_to_create_graph <- append(list_to_create_graph, current_node$name)
      }
      to_go_or_not_to_go <- sample(x = c(TRUE, FALSE), size = 1, prob = c(0.85, 0.15))
      if(to_go_or_not_to_go){
        node_ego_graph <- unique(neighbors(graph, v=current_node, mode="all"))
        if(length(node_ego_graph)>1){
          current_node <- sample(node_ego_graph, 1)
        }
      }
      else{
        current_node <- source_node
      }
      if (step_count > n){
        current_node <- sample(V(graph), 1)
        source_node <- current_node
        step_count <- 0
      }
    }
    
    file_string = paste("RW_sample_",round(i*100,digits=2),".rds",sep="")
    files <- append(files, file_string)
    sample_graph <- induced_subgraph(graph, v = list_to_create_graph)
    saveRDS(sample_graph, file=file_string)
    print(paste("Random Walk",round(i*100, digits=2), "% done"))
  }
  return(files)
}

evaluate_random_jump_sample = function(graph, min_sample = 0.15, max_sample = 0.5, steps = 5) {
  n <- length(V(graph))
  list_to_create_graph <- files <- vector()
  source_node <- sample(V(graph), 1)
  current_node <- source_node
  for(i in seq(min_sample, max_sample, length.out=steps)){
    while(length(list_to_create_graph) < n * i){
      if(!(current_node$name %in% list_to_create_graph)){
        list_to_create_graph <- append(list_to_create_graph, current_node$name)
      }
      to_jump_or_not_to_jump <- sample(x = c(TRUE, FALSE), size = 1, prob = c(0.15, 0.85))
      # print(to_jump_or_not_to_jump)
      if(!to_jump_or_not_to_jump){
        node_ego_graph <- unique(neighbors(graph, v=current_node, mode="all"))
        if(length(node_ego_graph)>1){
          current_node <- sample(node_ego_graph, 1)
        }
      }
      else{
        current_node <- sample(V(graph), 1)
      }
    }
    file_string = paste("RJ_sample_",round(i*100,digits=2),".rds",sep="")
    files <- append(files, file_string)
    saveRDS(induced_subgraph(graph, v = list_to_create_graph), file=file_string)
    print(paste("Random Jump",round(i*100, digits=2), "% done"))
  }
  return(files)
}

evaluate_forest_fire_sample = function(graph, min_sample = 0.15, max_sample = 0.5, steps = 5, pf=0.7){
  list_to_delete_from_graph <- files <- visited_list <- vector()
  n <- length(V(graph)) 
  vertices <- V(graph)$name
  geom_prob <- geom_distr(pf)
  for(i in rev(seq(min_sample, max_sample, length.out=steps))){
    to_cut <- n * i
    while(length(list_to_delete_from_graph) < n * (1 - i)){
      current_node <- sample(V(graph)[! V(graph) %in% list_to_delete_from_graph]$name, 1)
      list_to_delete_from_graph <- ff_spread(graph, current_node, visited_list = list_to_delete_from_graph, pf=pf, pb=pb, to_cut=to_cut, rec_count=0)
    }
    list_to_create_from_graph <- vertices[!vertices %in% list_to_delete_from_graph]
    file_string = paste("FF_sample_",round(i*100,digits=2),".rds",sep="")
    files <- append(files, file_string)
    saveRDS(induced.subgraph(graph, list_to_create_from_graph), file=file_string)
    print(paste("Forest Fire",round(i*100, digits=2), "% done"))
  }
  return(files)
}

ff_spread = function(graph, current_node, pf = 0.7, pb = 0, visited_list, to_cut, rec_count = 0, geom_prob = geom_distr(pf)){
  rec_count <- rec_count + 1
  if (!(current_node %in% visited_list)){
    visited_list <- append(visited_list, current_node)
    visited_list <- unique(visited_list)
    
  }
  if((length(visited_list) > to_cut) || (rec_count > 1000)){
    return(visited_list)
  }
  neighbours <- neighbors(graph, v=current_node)
  geom_prob <- geom_distr(pf)
  neighbours_vector <- (length(geom_prob)-1):0
  neighbours_num <- sample(x = neighbours_vector, prob = geom_prob, size = 1, replace = FALSE)
  no_visit_neighbours <- neighbours[!neighbours %in% visited_list]
  if ((length(no_visit_neighbours)>0) && (neighbours_num>0) && (length(no_visit_neighbours) > neighbours_num)){
    burnt_neighbours <- sample(x = no_visit_neighbours, size = neighbours_num)
    visited_list <- append(visited_list, burnt_neighbours)
  }
  else if((length(no_visit_neighbours>0)) && (neighbours_num>0) && (length(no_visit_neighbours) <= neighbours_num)){
    burnt_neighbours <- no_visit_neighbours
    visited_list <- append(visited_list, burnt_neighbours)
  }
  else{
    burnt_neighbours <- c()
  }
  
  for (neighbour in burnt_neighbours){
    visited_list <- append(visited_list, ff_spread(graph=graph, current_node = neighbour, visited_list = visited_list, to_cut=to_cut, pf=pf, pb=pb, rec_count=rec_count))
    visited_list <- unique(visited_list)
  }
  
  return(visited_list)
}

geom_distr = function(pf){
  mean = pf/(1-pf)
  y = 0:(2*round(mean)+1)
  py = dgeom(y,pf)
  return(py)
}
