random_node_sample = function(graph, sample_size = 0.15){
  sample_len <- sample_size * length(V(graph))
  nodes_vector <- sample(V(graph)$name, replace=FALSE, size = sample_len)
  sample_graph <- induced.subgraph(graph=graph, vids=nodes_vector)
  return(sample_graph)
}

random_page_rank_node_sample = function(graph, sample_size = 0.15){
  sample_len <- sample_size * length(V(graph))
  page_rank <- page.rank(graph, vids = V(graph))
  nodes_vector <- sample(V(graph)$name, prob = page_rank$vector, replace=FALSE, size = sample_len)
  sample_graph <- induced.subgraph(graph=graph, vids=nodes_vector)
  return(sample_graph)
}

random_degree_node_sample = function(graph, sample_size = 0.15){
  sample_len <- sample_size * length(V(graph))
  node_degrees <- degree(graph)
  node_degree_prob <- node_degrees / length(V(graph))
  nodes_vector <- sample(V(graph)$name, prob = node_degree_prob, replace=FALSE, size = sample_len)
  sample_graph <- induced.subgraph(graph=graph, vids=nodes_vector)
  return(sample_graph)
}

random_edge_sample = function(graph, sample_size = 0.15) {
  n <- length(V(graph))
  current_graph <- graph.empty()
  i <- 0
  while(length(V(current_graph)) < n*sample_size){
    i = i + 1
    random_edge <- sample(E(graph)$id, size=i)
    current_graph <- subgraph.edges(graph, random_edge)
  }
  return(current_graph)
}

random_node_neighbour_sample = function(graph, sample_size = 0.15){
  sample_len <- sample_size * length(V(graph))
  list_to_create_graph <- vector()
  while(length(list_to_create_graph) < sample_len){
    current_node <- sample(V(graph)$name, 1)
    node_neighbourhood <- neighbors(graph, v = current_node, mode="out")
    list_to_create_graph <- append(list_to_create_graph, current_node)
    list_to_create_graph <- append(list_to_create_graph, node_neighbourhood$name)
    list_to_create_graph <- unique(list_to_create_graph)
  }
  sample_graph <- induced.subgraph(graph, v = list_to_create_graph)
  return(sample_graph)
}

random_walk_sample = function(graph, sample_size = 0.15) {
  n <- length(V(graph))
  sample_size <- sample_size * n
  list_to_create_graph <- vector()
  source_node <- sample(V(graph), 1)
  step_count <- 0
  current_node <- source_node
  while(length(list_to_create_graph) < sample_size){
    if(length(list_to_create_graph) %% 1000 == 0){
      print(length(list_to_create_graph))
    }
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
  sample_graph <- induced_subgraph(graph, v = list_to_create_graph)
  return(sample_graph)
}

random_jump_sample = function(graph, sample_size = 0.15) {
  sample_size <- sample_size * length(V(graph))
  list_to_create_graph <- vector()
  source_node <- sample(V(graph), 1)
  current_node <- source_node
  while(length(list_to_create_graph) < sample_size){
    if(length(list_to_create_graph) %% 1000 == 0){
      print(length(list_to_create_graph))
    }
    if(!(current_node %in% list_to_create_graph)){
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
  sample_graph <- induced_subgraph(graph, v = list_to_create_graph)
  return(sample_graph)
}

forest_fire_sample = function(graph, sample_size = 0.15, pf = 0.7){
  list_to_delete_from_graph <- vector()
  visited_list <- vector() 
  to_cut <- length(V(graph)) * (1 - sample_size)
  while(length(list_to_delete_from_graph) < to_cut){
    current_node <- sample(V(graph)[! V(graph) %in% list_to_delete_from_graph]$name, 1)
    list_to_delete_from_graph <- ff_spread(graph, current_node, visited_list = list_to_delete_from_graph, pf=pf, pb=pb, to_cut=to_cut, rec_count=0)
  }
  list_to_create_from_graph <- V(graph)$name[!V(graph)$name %in% list_to_delete_from_graph]
  sample_graph <- induced.subgraph(graph, list_to_create_from_graph)
  return(sample_graph)
}

ff_spread = function(graph, current_node, pf = 0.7, pb = 0, visited_list, to_cut, rec_count = 0){
  rec_count <- rec_count + 1
  if(length(visited_list) %% 1000 == 0){
    print(length(visited_list))
  }
  if (!(current_node %in% visited_list)){
    visited_list <- append(visited_list, current_node)
    
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