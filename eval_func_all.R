evaluate_samples = function(graph, graph_name, source_closeness=NULL, source_betweenness=NULL,
                            runs = 3, min_sample=0.15, max_sample=0.5, steps=5){
  rw_stats <- rj_stats <- rnn_stats <- rdn_stats <- rprn_stats <- re_stats <- rn_stats <- ff_stats <- list()
  for(i in 1:runs){
    print(paste("RUN ", i))
    print(Sys.time())
    rw_eval <- evaluate_random_walk_sample(graph, min_sample, max_sample, steps)
    rj_eval <- evaluate_random_jump_sample(graph, min_sample, max_sample, steps)
    rnn_eval <- evaluate_random_node_neighbour_sample(graph, min_sample, max_sample, steps)
    rdn_eval <- evaluate_random_degree_node_sample(graph, min_sample, max_sample, steps)
    rprn_eval <- evaluate_random_page_rank_node_sample(graph, min_sample, max_sample, steps)
    re_eval <- evaluate_random_edge_sample(graph, min_sample, max_sample, steps)
    rn_eval <- evaluate_random_node_sample(graph, min_sample, max_sample, steps)
    ff_eval <- evaluate_forest_fire_sample(graph, 0.325, max_sample, steps)
    saveRDS(rw_eval, file=paste(graph_name, "_rw_eval_", i,".rds", sep=""))
    saveRDS(rj_eval, file=paste(graph_name, "_rj_eval_", i,".rds", sep=""))
    saveRDS(rnn_eval, file=paste(graph_name, "_rnn_eval_", i,".rds", sep=""))
    saveRDS(rdn_eval, file=paste(graph_name, "_rdn_eval_", i,".rds", sep=""))
    saveRDS(rprn_eval, file=paste(graph_name, "_rprn_eval_", i,".rds", sep=""))
    saveRDS(re_eval, file=paste(graph_name, "_re_eval_", i,".rds", sep=""))
    saveRDS(rn_eval, file=paste(graph_name, "_rn_eval_", i,".rds", sep=""))
    saveRDS(ff_eval, file=paste(graph_name, "_ff_eval_", i,".rds", sep=""))
    print(Sys.time())
    if(is.null(source_closeness)){
      source_closeness <- closeness(graph)
    }
    if(is.null(source_betweenness)){
      source_betweenness <- betweenness(graph)
    }
    
    rj_stat <- evaluate_d(rj_eval, graph, source_closeness, source_betweenness)
    rw_stat <- evaluate_d(rw_eval, graph, source_closeness, source_betweenness)
    rnn_stat <- evaluate_d(rnn_eval, graph, source_closeness, source_betweenness)
    rdn_stat <- evaluate_d(rdn_eval, graph, source_closeness, source_betweenness)
    rprn_stat <- evaluate_d(rprn_eval, graph, source_closeness, source_betweenness)
    re_stat <- evaluate_d(re_eval, graph, source_closeness, source_betweenness)
    rn_stat <- evaluate_d(rn_eval, graph, source_closeness, source_betweenness)
    ff_stat <- evaluate_d(ff_eval, graph, source_closeness, source_betweenness)

    rw_stats[[i]] <- rw_stat
    rj_stats[[i]] <- rj_stat
    rnn_stats[[i]] <- rnn_stat
    rdn_stats[[i]] <- rdn_stat
    rprn_stats[[i]] <- rprn_stat
    re_stats[[i]] <- re_stat
    rn_stats[[i]] <- rn_stat
    ff_stats[[i]] <- ff_stat
    print(Sys.time())
    gc()
  }
  df <- data.frame(
    algorithm=c("RandomNode","RandomPageRankNode","RandomDegreeNode","RandomEdge",
                "RandomNodeNeighbour","RandomJump","RandomWalk","ForestFire"),
    degree_in=c(em_di(rn_stats),em_di(rprn_stats),em_di(rdn_stats),em_di(re_stats),
                em_di(rnn_stats),em_di(rj_stats),em_di(rw_stats),em_di(ff_stats)),
    degree_out=c(em_do(rn_stats),em_do(rprn_stats),em_do(rdn_stats),em_do(re_stats),
                em_do(rnn_stats),em_do(rj_stats),em_do(rw_stats),em_do(ff_stats)),
    degree_all=c(em_da(rn_stats),em_da(rprn_stats),em_da(rdn_stats),em_da(re_stats),
                 em_da(rnn_stats),em_da(rj_stats),em_da(rw_stats),em_da(ff_stats)),
    clustering_coeff=c(em_cc(rn_stats),em_cc(rprn_stats),em_cc(rdn_stats),em_cc(re_stats),
                       em_cc(rnn_stats),em_cc(rj_stats),em_cc(rw_stats),em_cc(ff_stats)),
    closeness=c(em_cl(rn_stats),em_cl(rprn_stats),em_cl(rdn_stats),em_cl(re_stats),
                em_cl(rnn_stats),em_cl(rj_stats),em_cl(rw_stats),em_cl(ff_stats)),
    betweenness=c(em_be(rn_stats),em_be(rprn_stats),em_be(rdn_stats),em_be(re_stats),
                  em_be(rnn_stats),em_be(rj_stats),em_be(rw_stats),em_be(ff_stats))
  )
  return(df)
}

em_di = function(obj){
  sum <- 0
  for(i in obj){
    sum <- sum + i$degree_in_mean
  }
  return(sum/length(obj))
}
em_do = function(obj){
  sum <- 0
  for(i in obj){
    sum <- sum + i$degree_out_mean
  }
  return(sum/length(obj))
}
em_da = function(obj){
  sum <- 0
  for(i in obj){
    sum <- sum + i$degree_all_mean
  }
  return(sum/length(obj))
}
em_cc = function(obj){
  sum <- 0
  for(i in obj){
    sum <- sum + i$clustering_coeff_mean
  }
  return(sum/length(obj))
}
em_cl = function(obj){
  sum <- 0
  for(i in obj){
    sum <- sum + i$closeness_mean
  }
  return(sum/length(obj))
}
em_be = function(obj){
  sum <- 0
  for(i in obj){
    sum <- sum + i$betweenness_mean
  }
  return(sum/length(obj))
}