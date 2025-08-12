import pickle
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
font = {'size'   : 30}

plt.rc('font', **font)
plt.rcParams.update({'figure.autolayout': True})

databases = [
    # 'Flight',
    'Adult',
    # 'Stock',
    # 'Hospital',
    # 'Tax'
    
]

ntypes = [
    'conoise',
    'rnoise'
]
queries = [
            'vertex_cover', 
           'no_of_edges',
           'positive_degree_nodes'
           ]

size = 100
epsilon = 1.0

def repair_df(df):
    df['error'] = df['error'].replace(' inf',0.0)
    df['std'] = df['std'].replace(' nan',0.0)
    # convert df['error'] and df['std'] to float if in string
    df['error'] = df['error'].astype(float)
    df['std'] = df['std'].astype(float)
    return df

query_titles = {
            'positive_degree_nodes': r'$\mathcal{I}_P$',
            'no_of_edges': r'$\mathcal{I}_{MI}$',
            'vertex_cover': r'$\mathcal{I}_R$'
        }
for database in databases:
    results_folder = f'Results/{database}'

    
    for ntype in ntypes:
        fig, axs = plt.subplots(1, len(queries), figsize=(24, 6))
        # fig.suptitle(f'{database}', fontsize=16)
        
        for i, query in enumerate(queries):
            ax = axs[i]

            baselineV_file = f'{results_folder}/baseline_truedeg_{query}_{size}_{ntype}_eps_{epsilon}.csv'
            baselineV_df = pd.read_csv(baselineV_file)
            baselineV_df = repair_df(baselineV_df)

            if query != 'vertex_cover':
                baselinetruedeg_file = f'{results_folder}/baseline_maxdeg_{query}_{size}_{ntype}_eps_{epsilon}.csv'
                thetaexpomech_file = f'{results_folder}/expomech_{query}_{size}_{ntype}_eps_{epsilon}.csv'
                thetaexpomechhier_file = f'{results_folder}/hier_expomech_{query}_{size}_{ntype}_eps_{epsilon}.csv'
                thetaexpomechsumhier_file = f'{results_folder}/bound_hier_{query}_{size}_{ntype}_eps_{epsilon}.csv'

                baselinetruedeg_df = pd.read_csv(baselinetruedeg_file)
                thetaexpomech_df = pd.read_csv(thetaexpomech_file)
                thetaexpomechhier_df = pd.read_csv(thetaexpomechhier_file)
                thetaexpomechsumhier_df = pd.read_csv(thetaexpomechsumhier_file)

                # Clean the data
                baselineV_df = repair_df(baselineV_df)
                baselinetruedeg_df = repair_df(baselinetruedeg_df)
                thetaexpomech_df = repair_df(thetaexpomech_df)
                thetaexpomechhier_df = repair_df(thetaexpomechhier_df)
                thetaexpomechsumhier_df = repair_df(thetaexpomechsumhier_df)

                # Print differences
                diff = abs(thetaexpomechhier_df['error'] - baselinetruedeg_df['error']) / baselinetruedeg_df['error']
                print(f"Mean diff between hierarchical and baseline ({database}, {query}): {diff.mean():.4f}")
                diff = abs(baselineV_df['error'] - baselinetruedeg_df['error']) / baselinetruedeg_df['error']
                print(f"Mean diff between |V| and baseline ({database}, {query}): {diff.mean():.4f}")
            
            # Plot
            if query == 'no_of_edges':
                ax.plot(baselineV_df['conflicts'], baselineV_df['error'], label=r"$\theta = |V|$")
                ax.plot(baselinetruedeg_df['conflicts'], baselinetruedeg_df['error'], label=r"$\theta = \max \deg(v)$", linestyle='--')
                ax.plot(thetaexpomech_df['conflicts'], thetaexpomech_df['error'], label=r"Exponential mech.")
                ax.plot(thetaexpomechhier_df['conflicts'], thetaexpomechhier_df['error'], label=r"Hierarchical exp. mech.")
                ax.plot(thetaexpomechsumhier_df['conflicts'], thetaexpomechsumhier_df['error'], label=r"Upper bound + Hierarchical")

            elif query == 'positive_degree_nodes':
                ax.plot(baselineV_df['error_nodes'], baselineV_df['error'], label=r"$\theta = |V|$")
                ax.plot(baselinetruedeg_df['error_nodes'], baselinetruedeg_df['error'], label=r"$\theta = \max \deg(v)$", linestyle='--')
                ax.plot(thetaexpomech_df['error_nodes'], thetaexpomech_df['error'], label=r"Exponential mech.")
                ax.plot(thetaexpomechhier_df['error_nodes'], thetaexpomechhier_df['error'], label=r"Hierarchical exp. mech.")
                ax.plot(thetaexpomechsumhier_df['error_nodes'], thetaexpomechsumhier_df['error'], label=r"Upper bound + Hierarchical")

            elif query == 'vertex_cover':
                ax.plot(2 * baselineV_df['lin_R'], baselineV_df['error'], label='True vertex cover')

            ax.set_xlabel('True value of measure')
            ax.set_ylabel('Error')
            ax.set_yscale('log')
            ax.set_title(query_titles[query])

            if database == 'Adult' and query == 'positive_degree_nodes':
                ax.set_ylim(0.1, 10)
            
            ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x:,.0f}' if x < 1000 else f'{x/1000:.0f}k'))

            if i == len(queries) - 1:
                ax.legend(loc='center left', bbox_to_anchor=(1.05, 0.5))
        
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        plt.savefig(f'Plots/combined_comparingtheta_{database}_samegraph_{size}_{ntype}_eps_{epsilon}.jpg', dpi=300, bbox_inches='tight')
        plt.show()
