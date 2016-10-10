func.R -- All function used to sample graphs
Example usage:

sample_graph <- random_page_rank_node_sample(source_graph, sample_size=0.3)

d_statistic.R -- function used to test similarity of sampled graph using Kolmogorov-Smirnov two sample test
Example usage:

statistics_list <- d_statistic(source_graph, sample_graph, closeness(source_graph), betweenness(source_graph))

eval_func_all.R -- evaluate_samples function used to run all available sampling algorithms multiple times and return data.frame with mean D values
Example usage:

evaluate_samples(source_graph, output_graph_name, source_graph_closeness, source_graph_betweenness, runs=3, min_sample=0.15, max_sample=0.5, steps=5)

eval_func.R -- used by evaluate_samples function

