import pandas as pd
from util import build_dynamic_queries, constraints_check
import networkx as nx
import pickle
import random
import re
import os
import COnoiseViolationsAlgorithm as cvio
import RnoiseViolationsAlgorithm as rvio
random.seed(0)

datasets = [
            'Adult',
            'Stock',
            'Flight',
            'Hospital',
            'Tax', 
            ]


repeat = 1 # number of repeats with different seeds
conoise_iter = 200 # number of iterations for conoise
storing_interval = 10 # interval to store the graph
rnoise_alpha = 0.01 # percentage of cells to be violated with rnoise
rnoise_beta = 0 #skew of the Zipfiand distribution used to select values from the active domain
rnoise_typo_prob = 0.5 #probability of a typo or random value 
type_noise = 'rnoise' # 'rnoise' or 'conoise'
n_rows = 10000

# pick random seeds for each repeat
random_seeds = [random.randint(0, 1000) for _ in range(repeat)]

for testDirectoryPath in datasets:
    df_full = pd.read_csv('datasets/'+ testDirectoryPath + '/' + "inputDB2.csv", keep_default_na=False, na_values=['-1.#IND', '1.#QNAN', '1.#IND', '-1.#QNAN', '#N/A N/A', '#N/A', 'N/A', 'n/a','', '#NA', 'NULL','null', 'NaN', '-NaN', 'nan', '-nan', ''] , header=0)
    constraints_raw = open('datasets/'+ testDirectoryPath+'/dcs.txt', 'r')
    constraints = [line.strip() for line in constraints_raw.readlines()]
    constraints = [x.replace(' ', '_') for x in constraints] #in case the columns names include spaces

    print('Dataset:', testDirectoryPath, 'DCs:', len(constraints))

    for rep in range(repeat):
        
        df = df_full.sample(n=n_rows, random_state=random_seeds[rep])
        
        if type_noise == 'conoise':
        
            for iter in range(conoise_iter):
                global t1,t2
                sample = df.sample(n=2)
                t1 = sample.iloc[0]
                t2 = sample.iloc[1]

                # clean constraint from excessive chars
                constraintSetRaw = random.choice(constraints)
                constraintSet = constraintSetRaw[4:-1].split('&')
                constraintSet = [re.split('(!=|>=|<=|>|<|=)', i) for i in constraintSet]

                # in case the constraint refers to a single tuple
                if "t2" not in constraintSet:
                    t2 = t1
                
                # generate violations using the fittingViolationAlgorithm in ViolationsAlgorithm.py
                t = cvio.fittingViolationAlgorithm(constraintSet,df,t1,t2)
                cvio.updateTable(df,t[0],t[1],sample)

                if iter % storing_interval == 0:
                    allColumns = {}
                    orig_columns = df.columns
                    for col in df.columns: 
                        allColumns[col] = col.replace(' ','_')
                    df = df.rename(columns=allColumns)

                    allConstraints = build_dynamic_queries(constraints,df)


                    # violating_pairs, violating_tuples, time1, time2 = constraints_check(df,constraints, allColumns, allConstraints[0], allConstraints[1])
                    violating_pairs, time1 = constraints_check(df,constraints, allColumns, allConstraints[0], allConstraints[1])

                    G = nx.Graph()

                    for i in range(len(df)):
                            G.add_node(i)
                    
                    # iterate all rows of dataframe violating_pairs and add edges to graph
                    for i, row in violating_pairs.iterrows():
                        if row['id1'] != row['id2']: #don't add if self loop
                            G.add_edge(row['id1'], row['id2'])
                    
                    saving_dir = 'datasets/'+ testDirectoryPath + '/conflict_graph' + '/'

                    #create directory if not exists
                    if not os.path.exists(saving_dir):
                        os.makedirs(saving_dir)

                    # save the graph
                    with open(saving_dir + f'graph_{n_rows}_conoise_{iter}_{rep}.pkl', 'wb') as f:
                        pickle.dump(G, f)

                    df.to_csv(saving_dir + f'graph_{n_rows}_conoise_{iter}_{rep}.csv', index=False, columns=orig_columns)
                    conflicts = G.edges()
                    print('Iteration:', iter, 'Conflicts:', len(conflicts))
                    
        elif type_noise == 'rnoise':

            # calculate all propabilities
            listToStr = ' '.join([str(elem) for elem in constraints])
            colomnsInConstraints = rvio.col_in_constraints(listToStr,df)
            all_probs = rvio.calculate_all_probs(df,colomnsInConstraints,rnoise_beta)

            cells_count = len(df.columns) * df.shape[0]
            iterations = int(rnoise_alpha * cells_count) # number of cells to be violated
            print('Number of cells to be violated:', iterations)
            
            for iter in range(1, iterations):
                rvio.rand_vio_algorithm(df,colomnsInConstraints,all_probs,rnoise_typo_prob)
                
                if iter % storing_interval == 0:
                    allColumns = {}
                    orig_columns = df.columns
                    for col in df.columns: 
                        allColumns[col] = col.replace(' ','_')
                    df = df.rename(columns=allColumns)

                    allConstraints = build_dynamic_queries(constraints,df)


                    violating_pairs, time1 = constraints_check(df,constraints, allColumns, allConstraints[0], allConstraints[1])
                    G = nx.Graph()

                    for i in range(len(df)):
                            G.add_node(i)

                    # iterate all rows of dataframe violating_pairs and add edges to graph
                    for i, row in violating_pairs.iterrows():
                        if row['id1'] != row['id2']: #don't add if self loop
                            G.add_edge(row['id1'], row['id2'])
                    
                    saving_dir = 'datasets/'+ testDirectoryPath + '/conflict_graph' + '/'

                    #create directory if not exists
                    if not os.path.exists(saving_dir):
                        os.makedirs(saving_dir)

                    # save the graph
                    with open(saving_dir + f'graph_{n_rows}_rnoise_{iter}_{rep}.pkl', 'wb') as f:
                        pickle.dump(G, f)

                    df.to_csv(saving_dir + f'graph_{n_rows}_rnoise_{iter}_{rep}.csv', index=False, columns=orig_columns)

                    conflicts = G.edges()
                    print('Iteration:', iter, 'Conflicts:', len(conflicts))
            
    print('Graphs saved for dataset:', testDirectoryPath)
    