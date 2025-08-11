import pickle
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
font = {'size'   : 18}

plt.rc('font', **font)


databases = [
    'Stock',
    'Hospital',
    'Flight',
    # 'Adult',
    
    # 'Tax',
    
]

ntypes = [
    # 'conoise',
    'rnoise'
]
queries = [
            'vertex_cover', 
           'no_of_edges',
           'positive_degree_nodes'
           ]

size = 10000
epsilons = [0.1, 0.2, 0.5, 1.0, 2.0, 3.0, 5.0]
epsilon = 1.0

def repair_df(df):
    df['error'] = df['error'].replace(' inf',0.0)
    df['std'] = df['std'].replace(' nan',0.0)
    # convert df['error'] and df['std'] to float if in string
    df['error'] = df['error'].astype(float)
    df['std'] = df['std'].astype(float)
    return df


ntype = ntypes[0]
for query in queries:

    
    plt.figure(figsize=(8,6))
    
        
    for database in databases:
        results_folder = f'Results/{database}'

        values = []
        errors = []
        for epsilon in epsilons:
                # if query != 'vertex_cover':
            result_file = f'{results_folder}/thetasum_hier_{query}_samegraph_{size}_{ntype}_eps_{epsilon}.csv'
            # else:
                # result_file = f'{results_folder}/{query}_samegraph_{size}_{ntype}_eps_{epsilon}.csv'
            
            df = pd.read_csv(result_file)
            
            # convert all inf values to 0 in df['error'] and df['std']
            df = repair_df(df)

            s = df.shape[0]
            values.append(df.iloc[0]['error'])
            errors.append(df.iloc[0]['std'])

        
        
        # plot the df with x axis from col 'iter' and y axis from 'error' with std deviation 'std'
        x_eps = [str(eps) for eps in epsilons]
        # plt.errorbar(df['iter'], df['error'], yerr=df['std'], label=f'{ntype}', fmt='-o')
        if query == 'no_of_edges':
            
            eb = plt.errorbar(x_eps, values,yerr=errors, label=database, capsize=5, capthick=2)
            eb[-1][0].set_linestyle('--')
        elif query == 'positive_degree_nodes':
            eb =plt.errorbar(x_eps, values, yerr=errors, label=database,  capsize=5, capthick=2)
            eb[-1][0].set_linestyle('--')    
        elif query == 'vertex_cover':
            # plt.errorbar(x_eps, values, yerr=errors, label=r"$I_{R}$")
            eb =plt.errorbar(x_eps, values, yerr=errors, label=database,  capsize=5, capthick=2)
            eb[-1][0].set_linestyle('--')
        
        print(database)
        print(errors)

        plt.xlabel('Privacy Budget')
       
        plt.xticks(x_eps)
        # cutoff lower bound for y axis at 0
        # plt.ylim(bottom=0)
        if query == 'vertex_cover':
            plt.yscale('log')
            plt.ylabel('Error (log scale)')
        else:
            plt.ylabel('Error')
            
    # plt.legend()
    
    if query == 'no_of_edges':
        plt.title(r"$I_{MI}$")
    elif query == 'positive_degree_nodes':
        plt.title(r"$I_{P}$")
    elif query == 'vertex_cover':
        plt.title(r"$I_{R}$")
    # plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.5), fancybox=True, ncol=3)
    plt.savefig(f'Plots/epsilon/{query}.jpg', dpi=300, bbox_inches='tight')
    plt.show()
            