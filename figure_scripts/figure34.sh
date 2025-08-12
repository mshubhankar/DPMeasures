#!/bin/sh
# This script is used to generate Figure 3 in the paper.
# Written only for adult dataset. Repeat for other datasets if needed.
# Usage: ./figure3.sh
# Make sure to run this script from the root directory of the project.

#Adult dataset

#For conoise
# python create_graph.py --storing_interval 10 --type_noise conoise --n_rows 100
# python compute_measures.py --database Adult --noise_type conoise --algo_version bound_hier --n_rows 100 #for our most improved method
# python compute_measures.py --database Adult --noise_type conoise --algo_version baseline_truedeg --n_rows 100 #for our most improved method
#For rnoise
# python create_graph.py --storing_interval 10 --type_noise rnoise --n_rows 100
# python compute_measures.py --database Adult --noise_type rnoise --algo_version bound_hier --n_rows 100 #for our most improved method
# python compute_measures.py --database Adult --noise_type rnoise --algo_version baseline_truedeg --n_rows 100 #for our most improved method

python figure_scripts/generate_fig34.py