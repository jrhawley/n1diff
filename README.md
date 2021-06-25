# Hedging uncertainty in differential gene expression analyses with James-Stein estimators

Code, data, and analysis for a chapter of my PhD thesis.
The code for calculating the James-Stein fold change estimates can be found in [`code/jse-shrinkage/jse.R`](code/jse-shrinkage/jse.R).

## Usage

This repository is structured as follows:

```shell
.
├── data/            # directory where all non-analysis data is stored
└── code/
    ├── result1/     # analysis scripts for `result1`
    ├── result2/     # analysis scripts for `result2`
    └── ...
└── results/
    ├── result1/     # results and plots, generated from `/code/result1/`
    ├── result2/     # results and plots, generated from `/code/result2/`
    └── ...
├── README.md        # this file
└── environment.yaml # Anaconda environment YAML file for the entire project
```

To re-run any of the analyses in the `code/` folders:

1. Build and activate the `conda` environment stored in `environment.yaml`
    ```shell
    conda create --file environment.yaml -n <ENV_NAME>
    conda activate <ENV_NAME>
    ```
2. Navigate to the result directory of interest
3. Run `snakemake`

That will regenerate the entire set of results for that specific folder.
You can preview that needs to be run by running `snakemake -n` before running the analyses.

