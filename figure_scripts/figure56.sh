# !/bin/sh
# This script is used to generate Figure 5 and 6 in the paper.
# Written only for adult dataset. Repeat for other datasets if needed.
# Usage: ./figure56.sh
# Make sure to run this script from the root directory of the project.


#For conoise
# python create_graph.py  --repeat 1 --conoise_iter 200 --storing_interval 10  --type_noise conoise --n_rows 100
# python compute_measures.py --database Adult --noise_type conoise --algo_version bound_hier --n_rows 100 #for our most improved method -- upperbound + hierarchical exponential mechanism
# python compute_measures.py --database Adult --noise_type conoise --algo_version hier_expomech --n_rows 100 #for our method -- hierarchical exponential mechanism
# python compute_measures.py --database Adult --noise_type conoise --algo_version expomech --n_rows 100 #for our method -- exponential mechanism
# python compute_measures.py --database Adult --noise_type conoise --algo_version baseline_maxdeg --n_rows 100 #for baseline max degree
# python compute_measures.py --database Adult --noise_type conoise --algo_version baseline_truedeg --n_rows 100 #for baseline true degree (nonprivate)

#For rnoise
# python create_graph.py  --repeat 1 --rnoise_alpha 0.01 --rnoise_beta 0 --rnoise_typo_prob 0.5 --type_noise rnoise --n_rows 100
# python compute_measures.py --database Adult --noise_type rnoise --algo_version bound_hier --n_rows 100 #for our most improved method -- upperbound + hierarchical exponential mechanism
# python compute_measures.py --database Adult --noise_type rnoise --algo_version hier_expomech --n_rows 100 #for our method -- hierarchical exponential mechanism   
# python compute_measures.py --database Adult --noise_type rnoise --algo_version expomech --n_rows 100 #for our method -- exponential mechanism
# python compute_measures.py --database Adult --noise_type rnoise --algo_version baseline_maxdeg --n_rows 100 #for baseline max degree
# python compute_measures.py --database Adult --noise_type rnoise --algo_version baseline_truedeg --n_rows 100 #for baseline true degree (nonprivate)

# python figure_scripts/generate_fig5.py
python figure_scripts/generate_fig6.py