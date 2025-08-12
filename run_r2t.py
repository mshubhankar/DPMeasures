import os
import glob
import pandas
import sqlalchemy
from sqlalchemy import text
from ExtractInfo import extract_run
from R2T_multiple import r2t_run
from tqdm import tqdm
import argparse

parser = argparse.ArgumentParser(description="Run R2T experiments on specified databases.")
parser.add_argument('--databases', nargs='+', default=['Adult'], help='List of databases to process')
parser.add_argument('--noise_type', type=str, default='conoise', help='Type of noise to use')
parser.add_argument('--n_rows', type=int, default=10000, help='Number of rows to consider in the dataset')
parser.add_argument('--epsilon', type=float, default=1.0, help='Privacy budget epsilon')
parser.add_argument('--runs', type=int, default=10, help='Number of runs for averaging results')

args = parser.parse_args()

databases = args.databases
noise_type = args.noise_type
n_rows = args.n_rows
epsilon = args.epsilon
runs = args.runs

algo = 'no_of_edges' # R2T doesn't support vertex_cover and positive_degree_nodes yet.


for database in databases:

    data_directory = "../datasets/"+database+"/"+"conflict_graph/"
    results_folder = f"../Results/{database}"
    if not os.path.exists(results_folder):
        os.makedirs(results_folder)


    print('database: ', database)
    
    result_csv_path = f'{results_folder}/r2t_{algo}_{n_rows}_{noise_type}_eps_{epsilon}.csv'
    result_csv = open(f'{results_folder}/r2t_{algo}_samegraph_{n_rows}_{noise_type}_eps_{epsilon}.csv', 'w')
    result_csv.write('iter,queried_result,private_result,privacy_result_std,error,error_std\n')
    print("algo: ", algo)
    
    file_regex = f'graph10k_{noise_type}_*_0.csv'
    # count files that match file_regex and iterate over them
    file_count = len(glob.glob(data_directory + file_regex))

    
    
    for iter in tqdm(range(10, file_count * 10, 10), desc="Processing files"):
        database_name = database.lower()
        con = sqlalchemy.create_engine(f'postgresql:///r2t')

        csv_name = data_directory + f'graph10k_{noise_type}_{iter}_0.csv'
        
        df = pandas.read_csv(csv_name)
        # Convert header to lower case for csv file
        df.columns = map(str.lower, d.columns)
        # Save csv back to file
        df.to_csv(csv_name, index=False)

        # Delete the table if it already exists
        df.to_sql(database_name, con, if_exists='replace')

        output_file = f'../results/{database}/output.txt'
        extract_run('r2t', f'../Query/{database_name}_{algo}.txt', database_name, f'../Query/id.txt', output_file)
        real_result, private_result, private_std, error, error_std, _, _ = r2t_run(output_file, epsilon, 0.1, n_rows, runs)
        result_csv.write(f'{iter},{real_result},{private_result},{private_std},{error},{error_std}\n')