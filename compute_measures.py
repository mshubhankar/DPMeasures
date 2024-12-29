import constants
import networkx as nx
import numpy as np
import pickle
import os
from graph_algos import vertex_cover, nodedp_add_edge_nedges_lap, nodedp_add_edge_pdnodes_lap
import util
import glob
import sys
import time

# take as argument database name and noise type
database = sys.argv[1] # database input
noise_type = sys.argv[2] # 'rnoise' or 'conoise'
algo_version = sys.argv[3] # either 'bound_hier', 'hier_expomech', 'expomech', 'baseline_maxdeg', 'baseline_truedeg'

print("database: ", database)

measures = [
            "no_of_edges", #I_MI measure 
            "positive_degree_nodes", #I_P measure
            "vertex_cover" #I_R measure
]

epsilons = [0.1, 0.2, 0.5, 1.0, 2.0, 3.0, 5.0]
repeats = 10
data_directory = 'datasets/' + database + '/conflict_graph/'
n_rows = 10000
optimization_eps_percentage = 0.4 # percentage of epsilon to be used for optimizations (exponential mechanism + upper bound calculation)
bound_eps_percentage = 0.25 # percentage of epsilon to be used for upper bound calculation from optimization
storing_interval = 10 # graphs stored interval


results_folder = f'Results/{database}'

if not os.path.exists(results_folder):
    os.makedirs(results_folder)


for epsilon in epsilons:
    if 'baseline' in algo_version:
        measure_epsilon = epsilon
    if 'expomech' in algo_version:
        expo_eps = epsilon * optimization_eps_percentage
        measure_epsilon = epsilon - expo_eps
    if 'bound_hier' in algo_version:
        
        expo_eps = epsilon * optimization_eps_percentage
        measure_epsilon = epsilon - expo_eps
        bound_eps = bound_eps_percentage * expo_eps
        expo_eps = expo_eps - bound_eps
       
    

    for measure in measures:

        #start writing in two csv files for random and rnoise results
        result_csv = open(f'{results_folder}/{algo_version}_{measure}_{n_rows}_{noise_type}_eps_{epsilon}.csv', 'w')
        
        result_csv.write('iter,conflicts,error_nodes,lin_R,measure_value,privacy_noise,error,std,time\n')
        
        print(f"Computing {measure} measure for {database} with {noise_type} noise and epsilon {epsilon}")
               
        file_regex = f'graph_{n_rows}_{noise_type}_*_0.pkl'
        # count files that match file_regex and iterate over them
        file_count = len(glob.glob(data_directory + file_regex))
        

        for iter in range(file_count):
            # load graph
            G = pickle.load(open(data_directory + f'graph_{n_rows}_{noise_type}_{iter}_0.pkl', 'rb'))
            num_nodes = len(G.nodes())
            max_deg = num_nodes - 1
            true_histogram = util.get_deg_his(G, max_deg) 

            conflicts = G.edges()
            non_private_lin_r, _ = util.I_lin_R(conflicts) # non-private lin_R
            print('Iteration:', iter, 'Conflicts:', len(conflicts))
            temp_errors = []
            temp_results = []
            temp_measures = []
            temp_pnoise = []
            temp_time = []
            
            if measure == "vertex_cover": #I_R measure has same algorithm for all versions
                
                for i in range(repeats):
                    start_time = time.time()
                    true_size, noisy_size, vertex_noise = vertex_cover(G,epsilon)
                    end_time = time.time()
                    if true_size == 0:
                        error = 0
                    else:
                        error = abs(noisy_size-(non_private_lin_r*2))/(non_private_lin_r*2)
                    noisy_result = noisy_size            
                    print('error: ', error)
                    measure_value = true_size
                    privacy_noise = vertex_noise

                    temp_errors.append(error)
                    temp_results.append(noisy_result)
                    temp_measures.append(measure_value)
                    temp_pnoise.append(privacy_noise)
                    temp_time.append(end_time-start_time)

            else: #I_MI and I_P measures
                if 'baseline' not in algo_version:
                    #Initializing candidate set for theta
                    cand_interval = n_rows//10
                    nearest_thousand = int(n_rows/cand_interval)*cand_interval
                    theta_candidates = [i for i in range(cand_interval, nearest_thousand, cand_interval)] #theta=n_rows is added later in expomech computation in graph_algos.py
                    small_candidates = [1, 5, 10, 100, 500]
                    theta_candidates = small_candidates + theta_candidates
                elif algo_version == 'baseline_maxdeg':
                    theta_candidates = [n_rows]
                elif algo_version == 'baseline_truedeg':
                    theta_candidates = [max([i for i, val in enumerate(true_histogram) if val != 0])] # true max degree (non-private)
            
            
                temp_errors = []
                temp_results = []
                temp_measures = []
                temp_pnoise = []


                temp_errors = []
                temp_results = []
                temp_measures = []
                temp_pnoise = []

                if algo_version == 'bound_hier':
                    
                    upper_bound = util.get_upperbound(f'datasets/{database}', data_directory + f'graph_{n_rows}_{noise_type}_0_0.csv' , bound_eps)
                    
                    # post_processing as upper bound can't be less than 0 or greater than n_rows
                    if upper_bound <= 0:
                        upper_bound = 1
                
                    if upper_bound > n_rows:
                        upper_bound = n_rows

                    # remove all elements in theta_candidates that are larger than the sum_theta
                    theta_candidates = [i for i in theta_candidates if i <= upper_bound]
                    # add sum_theta to theta_candidates if it is not already in there
                    if upper_bound not in theta_candidates:
                        theta_candidates.append(upper_bound)

                for i in range(repeats):
                
                    if measure == "no_of_edges":
                        start_time = time.time()
                        if 'baseline' in algo_version:
                            noisy_edges, true_theta_edges = nodedp_add_edge_nedges_lap(G, num_nodes, measure_epsilon, theta_candidates, algo_version=algo_version)
                        else:
                            noisy_edges, true_theta_edges = nodedp_add_edge_nedges_lap(G, num_nodes, epsilon, theta_candidates, algo_version=algo_version, expo_eps=expo_eps)
                        end_time = time.time()
                        true_edges = G.number_of_edges()
                        noisy_result = noisy_edges
                        if noisy_result < 0: #post-processing as number of edges can't be less than 0 
                            noisy_result = 0
                        error = abs(noisy_result-true_edges)/true_edges
                        measure_value = true_theta_edges
                        privacy_noise = true_theta_edges - noisy_result
                        print('error: ', error)
                    
                    elif measure == "positive_degree_nodes":
                        start_time = time.time()
                        if 'baseline' in algo_version:
                            noisy_pdnodes, truncated_pdnodes = nodedp_add_edge_pdnodes_lap(G, num_nodes, measure_epsilon, theta_candidates, algo_version=algo_version)
                        else:          
                            noisy_pdnodes, truncated_pdnodes = nodedp_add_edge_pdnodes_lap(G, num_nodes, epsilon, theta_candidates, algo_version=algo_version, expo_eps=expo_eps)
                        end_time = time.time()
                        true_pos_nodes = num_nodes - true_histogram[0]
                        
                        error = abs(noisy_pdnodes-true_pos_nodes)/true_pos_nodes
                        noisy_result = noisy_pdnodes
                        if noisy_result < 0: #post-processing as number of positive degree nodes can't be less than 0 
                            noisy_result = 0
                        measure_value = truncated_pdnodes
                        privacy_noise = truncated_pdnodes - noisy_pdnodes
                        print('error: ', error)

                    else:
                        print("no valid algo")
                    temp_errors.append(error)
                    temp_results.append(noisy_result)
                    temp_measures.append(measure_value)
                    temp_pnoise.append(privacy_noise)
                    temp_time.append(end_time-start_time)
    
        result_csv.write(f'{iter}, {len(conflicts)},{num_nodes - true_histogram[0]},{non_private_lin_r},{np.mean(temp_measures)},{np.mean(temp_pnoise)},{np.mean(temp_errors)},{np.std(temp_errors)}, {np.mean(temp_time)}\n')

    result_csv.close()
        