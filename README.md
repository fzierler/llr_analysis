# A note on algorithms and observables for finite temperature Yang-Mills theories on the lattice - Analysis workflow

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)

The workflow in this repository performs the analyses presented in the paper
[A note on algorithms and observables for finite temperature Yang-Mills theories on the lattice [XXXX.YYYYY]](https://arxiv.org/abs/XXXX.YYYYY).

## Requirements

- Conda, for example, installed from [Miniforge][miniforge]
- [Snakemake][snakemake], which may be installed using Conda

## Setup

1. Install the dependencies above.
2. Clone this repository including submodules
   (or download its Zenodo release and `unzip` it) 
   and `cd` into it:

   ```shellsession
   git clone --recurse-submodules https://github.com/fzierler/llr_analysis
   cd llr_analysis
   git checkout devel
   ```

3. The raw data and metadata can be downloaded from Zenodo
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX).
Download `metadata.tar` and `raw_data.tar` from Zenodo and decompress the 
archives in this directory. On slow and/or unstable connections consider using 
`wget -c` to download large files. 

## Running the workflow

The workflow is run using Snakemake with a separate snakefile per gauge group:

``` shellsession
snakemake --use-conda --cores 1 --snakefile workflow/Snakefile_sp4.smk
snakemake --use-conda --cores 1 --snakefile workflow/Snakefile_su3.smk
snakemake --use-conda --cores 1 --snakefile workflow/Snakefile_su4.smk
```

where the number `1` may be replaced by the number of CPU cores you wish to 
allocate to the computation.

Snakemake will automatically download and install all required Python packages.
This requires an Internet connection; if you are running in an HPC environment 
where you would need to run the workflow without Internet access, details on how
to preinstall the environment can be found in the [Snakemake documentation][snakemake-conda].

Using all 16 cores on an AMD 5950x CPU the analysis takes roughly XX minutes.

``` shellsession
time snakemake --use-conda --cores 16 --snakefile workflow/Snakefile_sp4.smk --forceall
real    XXmYY.ZZZs
user    XXmYY.ZZZs
sys     XXmYY.ZZZs

time snakemake --use-conda --cores 16 --snakefile workflow/Snakefile_su3.smk --forceall
real    XXmYY.ZZZs
user    XXmYY.ZZZs
sys     XXmYY.ZZZs
```

## Output

Output plots, tables, and definitions
are placed in the `assets/{gauge_group}/plots`, `assets/{gauge_group}/tables`, 
and `assets/{gauge_group}/definitions` directories, where {gauge_group} 
specifies the gauge group. 


Output data assets are placed into the `data_assets` directory.
Intermediary data are placed in the `intermediary_data` directory.

## Reusability

This workflow is relatively tailored to the data which it was originally written
to analyse. Additional ensembles may be added to the analysis by adding relevant
files to the `raw_data` directory, and adding corresponding entries to the files
in the `metadata` directory. However, extending the analysis in this way has not
been as fully tested as the rest of the workflow, and is not guaranteed to be 
trivial for someone not already familiar with the code.

[miniforge]: https://github.com/conda-forge/miniforge
[snakemake]: https://snakemake.github.io
[snakemake-conda]: https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html
