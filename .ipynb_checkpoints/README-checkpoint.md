# Project documentation

Data analysis of the GFP fluctuation project is structured in different scripts that are detailed here.

## Environment

### R 

R 4.1.2-foss-2021a (14 cores, 90 Gb)
All packages and corresponding versions are specified in /R_packages.csv
All packages are in /scicore/home/nimwegen/rocasu25/R/x86_64-pc-linux-gnu-library/4.1

### Python

Python 3.8.13 (GCC 7.5.0), with packages os, pandas, numpy, matplotlib, re, shutil, subprocess, seaborn, skimage, glob, tifffile, zarr and scipy.
Environment is stored in /scicore/home/nimwegen/rocasu25/Documents/Projects/biozentrum/.conda_environment.

### gfp_gaussian_process

Should be installed within this git directory.
https://github.com/bjks/gfp_gaussian_process
Commit 5be92d0.

### ggp_notebooks

Should be installed within this git directory.
https://github.com/bjks/ggp_notebooks
Commit 82ce0a2.

## Project structure

The essential scripts and files are listed below:

```
gfp_fluctuations_project
│   README.md
│   all_experiments.csv: list all experiments I have performed so far, with corresponding experimental details
│
└───analysis: contains Rmd files were data are analyzed to adress specific questions, these Rmd files are also used to produce figures.
|      └────figures
│      do_concentration_distributions_depend_on_gl_coordinates.Rmd
│      is_cv_constant_over_the_course_of_the_experiments.Rmd
│      quantifying_noise_contributions.Rmd
│      comparing_fullnoise_and_experimental_distributions.Rmd
|      comparing_distributions_between_parameter_regimes.Rmd
│      growth_rate_production_heatmap.Rmd
|      are_cv_gr_q_constant_over_the_course_of_the_experiments.Rmd
|      inspect_raw_data_distributions.Rmd
|      
│      
└───mother_machine_experiments_toolbox
│      preprocessing_curating_tyding.ipynb
│      load_functions_and_packages.R
│      read_MM_data.R
│      denoise_data.Rmd
│      import_denoised_data.R: denoised_data_files_dir has to be specified
│      load_denoised_data.R: denoised_data_dir has to be specified
|      load_denoised_data_complete_cycles.R: denoised_data_dir has to be specified
|      load_simulation_data.R: 'simulation_data_folder' has to be specified (path to forward integration folder)
│      import_simulated_data.R
│      mm.properties
│
└───import_scripts: R scripts to import, curate and export raw data
│      promoter_condition.Rmd: import, curated and save in ../raw_curated_data_complete_cycles, relevant data. NB: uncomplete cell cycles can also be saved.
|          uncomplete data are saved in ../raw_curated_data
│      datalist_promoter_condition.csv: information used by promoter_condition.Rmd to import data
│      promoter_import_all_data.Rmd: run all promoter_condition.Rmd files
│      file_list.txt: list of all files that were curated using deepMoma 
│      make_yaml_files_for_new_curation.Rmd: temporary, follow efforts to recurate datasets
│
└───raw_curated_data: contains raw data that can be denoised (**complete and uncomplete** cell cycles), as well as plots of cell traces (log(length) versus time)
│      e.g: hi1_glycerol040_rawdata.csv and hi1_glycerol040.pdf
|
└───raw_curated_data_complete_cycles: contains raw data that can be denoised (**only complete cell cycles**), as well as plots of cell traces (log(length) versus time)
│      e.g: hi1_glycerol040_rawdata.csv and hi1_glycerol040.pdf
|
└───denoising_raw_date: contains the output of the denoising procedure, experimental (data + logs, each denoising attempt is specified by a timestamp: e.g denoising_20220428095255)
│      └───parameters: contains input parameters for the denoising procedure
│      │     promoter_condition_parameters.txt
│      │     e.g: hi1_acetate_parameters.txt
│      │
│      └───denoising_timestamp: denoising procedure output
       denoising_readme.md: content of the denoising procedure output folder
│
└───denoised_data: contains denoised data, e.g: glucose020_hi1_denoised.csv
│
└───denoised_data_complete_cycles: contains denoised data, e.g: glucose020_hi1_denoised.csv
│   
└───gfp_gaussian_process: repository containing the code necessary to run denoising procedure 
│   
└───ggp_notebooks: repository containing the code to run forward integration
    
```

# Notes about autocorrelation functions

From a terminal, in gfp_gaussian_process  conda activated, use:

```
python correlation_from_joint.py -d /scicore/home/nimwegen/rocasu25/Documents/Projects/biozentrum/gfp_fluctuations_project/denoising_raw_data/denoising_20230201134400/ -k acetate005 acetate020 glycerol040 glucose020 glucoseaa020 -dt 12 12 6 3 1.5 -delimiter _
```

# TEMP

/scicore/home/nimwegen/GROUP/Moma/MM_Analysis/builds/unstable_20230131__FOR_PREPROCESSING_ONLY__mmpreprocesspy_v0.2.0-98d609a7/Modules


