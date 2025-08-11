# DPMeasures
The codebase contains code for the paper "Computing Inconsistency Measures Under Differential Privacy"

## Hardware requirements

This codebase has been tested on a MacBook Pro with an Apple M1 chip. The main hardware requirements are:

- Apple Silicon (M1 or later) or equivalent CPU
- At least 16 GB RAM recommended for large datasets
- Sufficient disk space for storing generated graphs and results

**Note:** The most time-consuming hardware requirement is during the graph generation step (`create_graph.py`), especially when working with larger datasets and a higher number of constraints. Ensure your system has adequate resources to handle these workloads efficiently.

## Installation

To set up the Python environment using conda and install the required packages, follow these steps:

1. Create a new conda environment:
    ```sh
    conda create --name dpmeasures python=3.11
    ```

2. Activate the environment:
    ```sh
    conda activate dpmeasures
    ```

3. Install the required packages from the `requirements.txt` file:
    ```sh
    pip install -r requirements.txt
    ```
 4. Install Gurobi and gurobipy:
     ```sh
     conda install -c gurobi gurobi
     ```

     Note: A licensed version of Gurobi is required for datasets with larger rows and constraints.

## Usage

To run the code, follow these steps:


1. Create the graphs using `create_graph.py`:
   ```sh
   python create_graph.py
   ```

   The `create_graph.py` script defines several variables that you can adjust to customize the graph creation process:

  - `repeat`: Number of repeats with different seeds (default: 1)
  - `conoise_iter`: Number of iterations for conoise (default: 200)
  - `storing_interval`: Interval to store the graph (default: 10)
  - `rnoise_alpha`: Percentage of cells to be violated with rnoise (default: 0.01)
  - `rnoise_beta`: Skew of the Zipfian distribution used to select values from the active domain (default: 0)
  - `rnoise_typo_prob`: Probability of a typo or random value (default: 0.5)
  - `type_noise`: Type of noise, either 'rnoise' or 'conoise' (default: 'rnoise')
  - `n_rows`: Number of rows (default: 10000)
   You can modify these variables directly in the script or pass them as command-line arguments.

2. Compute the measures using `compute_measures.py`:
  ```sh
  python compute_measures.py <database> <noise_type> <algo_version>
  ```

  The `compute_measures.py` script requires the following arguments:

  - `<database>`: Database name
  - `<noise_type>`: Type of noise, either 'rnoise' or 'conoise'.
  - `<algo_version>`: Algorithm version to use, either 'bound_hier', 'hier_expomech', 'expomech', 'baseline_maxdeg', or 'baseline_truedeg'.

## Generate figures using scripts in the `figure_scripts` directory:
   ```sh
   python figure_scripts/figureX.sh
   ```
    Replace `figureX.sh` with the appropriate script for the figure you want to generate (e.g., `figure1.sh`, `figure2.sh`, etc.). Figure 3-4 and Figure 5-6 have same script. 

 **Note:** The code introduces privacy noise by sampling from noisy distributions. As a result, the figures you generate may differ from those presented in the paper due to the inherent randomness in the noise sampling process.