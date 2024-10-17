# DPMeasures
The codebase contains code for the paper "Computing Inconsistency Measures Under Differential Privacy"
## Installation

To set up the Python environment using conda and install the required packages, follow these steps:

1. Create a new conda environment:
    ```sh
    conda create --name dpmeasures python=3.8
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

