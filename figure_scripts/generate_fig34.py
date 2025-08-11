import pickle
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
font = {'size'   : 30}

plt.rc('font', **font)
plt.rcParams.update({'figure.autolayout': True})

databases = [
    'Flight',
    'Adult',
    'Stock',
    'Hospital',
    'Tax'
    
]

ntypes = [
    'conoise',
    'rnoise'
]
queries = [
             
           'no_of_edges',
           'positive_degree_nodes',
           'vertex_cover',
           ]

def repair_df(df):
    df['error'] = df['error'].replace(' inf',0.0)
    df['std'] = df['std'].replace(' nan',0.0)
    # convert df['error'] and df['std'] to float if in string
    df['error'] = df['error'].astype(float)
    df['std'] = df['std'].astype(float)
    return df

size = 10000
epsilon = 1.0

stats_dict = {}

for database in databases:
    results_folder = f'Results/{database}'
    
    stats_dict[database] = {}

    for ntype in ntypes:
        from matplotlib.gridspec import GridSpec

        fig = plt.figure(figsize=(24, 6))
        gs = GridSpec(1, 3, figure=fig, width_ratios=[1, 1, 1])

        axs = [fig.add_subplot(gs[0, i]) for i in range(3)]
 

        query_titles = {
            'positive_degree_nodes': r'$\mathcal{I}_P$',
            'no_of_edges': r'$\mathcal{I}_{MI}$',
            'vertex_cover': r'$\mathcal{I}_R$'
        }

        for i, query in enumerate(queries):
            ax = axs[i]
            stats_dict[database][query] = {}

            baselineV_file = f'{results_folder}/baseline|V|_{query}_samegraph_{size}_{ntype}_eps_{epsilon}.csv'
            baselineV_df = repair_df(pd.read_csv(baselineV_file))

            if query != 'vertex_cover':
                baselinetruedeg_file = f'{results_folder}/baselinetruemaxdeg_{query}_samegraph_{size}_{ntype}_eps_{epsilon}.csv'
                thetaexpomech_file = f'{results_folder}/expomech_{query}_samegraph_{size}_{ntype}_eps_{epsilon}.csv'
                thetaexpomechhier_file = f'{results_folder}/expomech_hier_{query}_samegraph_{size}_{ntype}_eps_{epsilon}.csv'
                thetaexpomechsumhier_file = f'{results_folder}/thetasum_hier_{query}_samegraph_{size}_{ntype}_eps_{epsilon}.csv'

                baselinetruedeg_df = repair_df(pd.read_csv(baselinetruedeg_file))
                thetaexpomech_df = repair_df(pd.read_csv(thetaexpomech_file))
                thetaexpomechhier_df = repair_df(pd.read_csv(thetaexpomechhier_file))
                thetaexpomechsumhier_df = repair_df(pd.read_csv(thetaexpomechsumhier_file))

                thetaexpomechsumhier_df['result'] = thetaexpomechsumhier_df['measure_value'] + thetaexpomechsumhier_df['privacy_noise']
                thetaexpomechsumhier_df['std_dev'] = thetaexpomechsumhier_df['measure_value'] * thetaexpomechsumhier_df['std']
                thetaexpomechsumhier_df['lower_bound'] = (thetaexpomechsumhier_df['result'] - thetaexpomechsumhier_df['std_dev']).clip(lower=0)
                thetaexpomechsumhier_df['upper_bound'] = thetaexpomechsumhier_df['result'] + thetaexpomechsumhier_df['std_dev']

                stats_dict[database][query]['average_error'] = thetaexpomechsumhier_df['error'].mean()
                stats_dict[database][query]['max_error'] = thetaexpomechsumhier_df['error'].max()

            if query == 'no_of_edges':
                r2t_baseline_file = f'Results_improved/r2t/{database}/nedges_samegraph_{size}_{ntype}_eps_{epsilon}.csv'
                r2t_baseline_df = pd.read_csv(r2t_baseline_file)
                r2t_baseline_df['result'] = r2t_baseline_df['private_result'] / 2
                r2t_baseline_df['upper_bound'] = r2t_baseline_df['result'] + r2t_baseline_df['privacy_result_std'] / 2
                r2t_baseline_df['lower_bound'] = (r2t_baseline_df['result'] - r2t_baseline_df['privacy_result_std'] / 2).clip(lower=0)

                stats_dict[database][query]['average_error_r2t'] = r2t_baseline_df['error'].mean()
                stats_dict[database][query]['max_error_r2t'] = r2t_baseline_df['error'].max()

                ax.plot(thetaexpomechsumhier_df['iter'], thetaexpomechsumhier_df['result'], label="Our approach")
                ax.plot(thetaexpomechsumhier_df['iter'], thetaexpomechsumhier_df[' conflicts'], label="True value")
                ax.plot(r2t_baseline_df['iter'], r2t_baseline_df['result'], label="R2T")
                ax.fill_between(thetaexpomechsumhier_df['iter'],
                                thetaexpomechsumhier_df['lower_bound'],
                                thetaexpomechsumhier_df['upper_bound'], alpha=0.2)
                ax.fill_between(r2t_baseline_df['iter'],
                                r2t_baseline_df['lower_bound'],
                                r2t_baseline_df['upper_bound'], alpha=0.2, color='green')

            elif query == 'positive_degree_nodes':
                ax.plot(thetaexpomechsumhier_df['iter'], thetaexpomechsumhier_df['result'], label="Our approach")
                ax.plot(thetaexpomechsumhier_df['iter'], thetaexpomechsumhier_df['error_nodes'], label="True value")
                ax.fill_between(thetaexpomechsumhier_df['iter'],
                                thetaexpomechsumhier_df['lower_bound'],
                                thetaexpomechsumhier_df['upper_bound'], alpha=0.2)

            else:
                baselineV_df['result'] = baselineV_df['measure_value'] + baselineV_df['privacy_noise']
                baselineV_df['std_dev'] = baselineV_df['measure_value'] * baselineV_df['std']
                baselineV_df['lower_bound'] = (baselineV_df['result'] - baselineV_df['std_dev']).clip(lower=0)
                baselineV_df['upper_bound'] = baselineV_df['result'] + baselineV_df['std_dev']

                stats_dict[database][query]['average_error'] = baselineV_df['error'].mean()
                stats_dict[database][query]['max_error'] = baselineV_df['error'].max()

                ax.plot(baselineV_df['iter'], baselineV_df['result'], label="Our approach")
                ax.plot(baselineV_df['iter'], 2 * baselineV_df['lin_R'], label="True value")
                ax.fill_between(baselineV_df['iter'],
                                baselineV_df['lower_bound'],
                                baselineV_df['upper_bound'], alpha=0.2)

            ax.set_title(query_titles[query])
            ax.set_xlabel("Iterations")
            if i == 0:
                ax.set_ylabel("Measure value")
            else:
                ax.set_ylabel("")

            ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x:,.0f}' if x < 1000 else f'{x/1000:.0f}k'))
            ax.set_ylim(bottom=0)

            # Custom limit for Stock dataset on no_of_edges
            if database == 'Stock' and query == 'no_of_edges':
                if ntype == 'rnoise':
                    ax.set_ylim(top=20)
                elif ntype == 'conoise':
                    ax.set_ylim(top=200)


        # fig.suptitle(f'{database}')
        handles, labels = ax.get_legend_handles_labels()
        # fig.legend(handles, labels, loc='lower center', ncol=3, bbox_to_anchor=(0.5, -0.02))
        fig.tight_layout(rect=[0, 0.05, 1, 0.95])
        # fig.subplots_adjust(wspace=0.3)
        fig.savefig(f'Plots/combined_{database}_samegraph_{size}_{ntype}_eps_{epsilon}.jpg', dpi=300, bbox_inches='tight')
        plt.show()