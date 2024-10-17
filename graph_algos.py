import numpy as np
import util
from sklearn.linear_model import LinearRegression 


def learn_theta_nedges(G, num_nodes, epsilon, expo_eps, theta_list, n_hiers=1):
    
    expo_eps = expo_eps/n_hiers

    for i in range(n_hiers):
        prob_list = []
        max_theta = max(theta_list)
        sensitivity = max_theta
        deg_list_max_theta = deg_list_by_edge_addition(G, max_theta)
        max_edges = sum(deg_list_max_theta)/2
        scores = []
        for theta in theta_list:
            deg_list_theta = deg_list_by_edge_addition(G, theta)
            theta_edges = sum(deg_list_theta)/2
            score = - (max_edges - theta_edges) - np.sqrt(2) * theta/ epsilon  #e_bias + lap error
            scores.append(score)
            prob = np.exp(expo_eps * score /(2.0 * sensitivity))
            prob_list.append(prob)
        if i == 0:
            # adding final candidate = num nodes 
            theta_cand_list = theta_list + [num_nodes]
            score = - np.sqrt(2) * num_nodes/ epsilon # no e_bias
            scores.append(score)
            prob_list.append(np.exp(expo_eps * score /(2.0 * sensitivity)))
        
        selected_theta = theta_cand_list[util.sample_prob_list(prob_list)]
        theta_cand_list = [theta for theta in theta_cand_list if selected_theta >= theta] # preparing for next expo mech    
        
    return selected_theta


def learn_theta_pdnodes(G, num_nodes, epsilon, expo_eps, theta_list, n_hiers=1):
    
    expo_eps = expo_eps/n_hiers
    
    for i in range(n_hiers):
        prob_list = []
        max_theta = max(theta_list)
        sensitivity = max_theta
        deg_list_max_theta = deg_list_by_edge_addition(G, max_theta)
        deg_his_max_theta = util.deg_seq_to_deg_his(deg_list_max_theta, max_theta)
        pd_max_theta = len(G.nodes()) - deg_his_max_theta[0]
        scores = []
        for theta in theta_list:
            deg_list_theta = deg_list_by_edge_addition(G, theta)
            deg_his = util.deg_seq_to_deg_his(deg_list_theta, theta)
            pd_theta = len(G.nodes()) - deg_his[0]
            score = - (pd_max_theta - pd_theta) - np.sqrt(2) * theta/ epsilon  #n_bias + lap error
            prob = np.exp(expo_eps * score /(2.0 * sensitivity))
            prob_list.append(prob)
            scores.append(score)
    if i == 0:
        # adding final candidate = num nodes 
        num_nodes = len(G.nodes())
        theta_cand_list = theta_list + [num_nodes]
        score = - np.sqrt(num_nodes) * (num_nodes)/ epsilon # no n_bias
        prob_list.append(np.exp(expo_eps * score /(2.0 * sensitivity)))

    selected_theta = theta_cand_list[util.sample_prob_list(prob_list)]
    theta_cand_list = [theta for theta in theta_cand_list if selected_theta >= theta] # preparing for next expo mech 
     
    return selected_theta


# Day et al. SIGMOD'16 , algo 1
# projection by edge-addition, variant: output degSeq instead of graph, part of Algo 4
# Edge addition algorithm from empty graph till no edges can be added; keep degree bounded by theta
def deg_list_by_edge_addition(G, theta):    
    
    num_nodes = len(G.nodes())
    nodes = list(G.nodes()) # the node ids are not strictly from 0 to |nodesNum|-1
    inv_nodes_list = {} # map node id to index
    for id in range(num_nodes):
        v = nodes[id]
        inv_nodes_list[v]=id
    
    nodes_random_indices = np.random.permutation(num_nodes)
    deg_list_temp = np.zeros(num_nodes)
    
    for vid in nodes_random_indices: # edge addition
        v = nodes[vid]
        for u in G.neighbors(v):
            uid = inv_nodes_list[u]
            if uid < vid and deg_list_temp[uid]<theta and deg_list_temp[vid]<theta:
                deg_list_temp[uid] = deg_list_temp[uid]+1
                deg_list_temp[vid] = deg_list_temp[vid]+1
                
    deg_list_temp=sorted(deg_list_temp, reverse=False) # small to large degrees
    
    return deg_list_temp # return the degree sequence and the number of edges   


def nodedp_add_edge_nedges_lap(G, n_nodes, epsilon, theta_list, algo_version, expo_eps = 0.4):
    n_nodes = len(G.nodes())
    # Learning theta 
    if len(theta_list) >1: #many choices
        n_hiers = 1
        if 'hier' in algo_version:
            n_hiers = 2
        theta = learn_theta_nedges(G, n_nodes, epsilon, expo_eps, theta_list, algo_version=algo_version, n_hiers=n_hiers)
        print('chosen theta: ', theta)
    else:
        theta = theta_list[0]
    
    degree_list = deg_list_by_edge_addition(G,theta)
    n_edges = sum(degree_list)/2
    sens = theta
    noisy_n_edges = n_edges + np.random.laplace(0, sens/epsilon)
    
    return noisy_n_edges, n_edges


def nodedp_add_edge_pdnodes_lap(G, n_nodes, epsilon, theta_list, algo_version, expo_eps = 0.4):
   
    # Learning theta 
    if len(theta_list) >1: #many choices
        n_hiers = 1
        if 'hier' in algo_version:
            n_hiers = 2
        theta = learn_theta_pdnodes(G, n_nodes, epsilon, expo_eps, theta_list, algo_version=algo_version, n_hiers=n_hiers)
        print('chosen theta: ', theta)
    else:
        theta = theta_list[0]   

    #Edge addition algorithm from empty graph till no edges can be added keep degree bounded by theta
    degree_list = deg_list_by_edge_addition(G,theta)
    deg_his = util.deg_seq_to_deg_his(degree_list, theta) # get degree histogram
    pd_nodes = n_nodes - deg_his[0] # number of nodes with positive degree = total nodes - number of nodes with degree 0
    sens = theta 
    noisy_pd_nodes = pd_nodes + np.random.laplace(0, sens/epsilon)
    
    return noisy_pd_nodes, pd_nodes
    


def vertex_cover(G, epsilon):
    G = G.copy()
    num_nodes = len(G.nodes())
    nodes = list(G.nodes()) # the node ids are not strictly from 0 to |nodesNum|-1

    c = 0 # vertex cover size

    # get all edges of the graph
    edges = list(G.edges())

    random_perm_edges = list(np.random.permutation(edges))
    
    for (u,v) in random_perm_edges:
        #check if (u,v) edge exists in G
        if G.has_edge(u,v):
            if (u == v):
                print('self loop', u)
            c = c + 2
            G.remove_node(u) 
            G.remove_node(v)
    privacy_noise = np.random.laplace(0, 2/epsilon)
    noisy_c = c + privacy_noise
    return c, noisy_c, privacy_noise 
        
            
