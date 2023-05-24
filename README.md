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

## Dataset handling

### Raw data

The list of raw data is available in the following table: `./mother_machine_experiments_toolbox/all_mmexperiments_datalist.csv`.
Columns names are the following:

```
date,description,f_start,f_end,condition,t_interval,data_path,promoter,vector,cell_detection_offset,experimental_setup,time_threshold_h,replica
```

Where:
- `f_start`, `f_end` are frame numbers of start and end of the experiment
- `t_interval` is the time step
- `data_path` the folder where curated data (using deepMoma) can be found
- `cell_detection_offset` is the offset used (in pixels), to filter out cells too close to the top of the growth lane exit
- `experimental_setup` is the setup that was used to generate the data. 1 stands for `Stacy` (Ti-Eclipse epi-fluo microscope), 3 stands for `Momo` (Ti-2 epi-fluo microscope)
- `time threshold` is an arbitrary number of hours before which cells are discarded, to ensure MM experiment runs in a "permanent regime".
- `replica` is an arbitrary number to facilitate data analysis.

Other columns are pretty clear.

### Import, tidy, filter raw data

To import raw data, run `read_all_raw_MM_Data.R`. Necessary packages will be loaded, and `import_tidy_filter_raw_MMData.R` will be sourced.
A tidied and filtered of the raw data (one file per condition, promoter, date) will saved on the disk in `./raw_curated_data_complete_cycles`.

NB: `cell_detection_offset` and `time_threshold` are used in `import_tidy_filter_raw_MMData.R` to filter out cell-cycles that occur at the beginning of the experiment, or that get too close to the edge of the growth-lane.

### Denoise raw data

Once raw data are saved in `./mother_machine_experiments_toolbox/raw_curated_data_complete_cycles`, one can run the following script to denoise the data:
`./mother_machine_experiments_toolbox/denoise_raw_data.R`. This script will call the `gfp_gaussian` software (`https://github.com/bjks/gfp_gaussian_process`, Bjoern Kscheschinski) to denoised the data, using scicore nodes (1 node per job).

Data will be stored in `./mother_machine_experiments_toolbox/denoising_raw_data`, in a folder that can be identified with a timestamp (e.g `denoising_20230522144637`). Importantly, joints are also saved in the same folder.

### Tidy and import denoised data

Run `./mother_machine_experiments_toolbox/import_denoised_data.R`. Data will be tidied with all necessary information contained in raw data and then saved in `./denoised_data_complete_cycles`.
NB: the same data will also be saved in the folder: `denoised_data_complete_cycles_corrected` together with inferred parameters, with q corrected for autofluorescence volumic production.

### Proceed with forward integration

Currently, I run Bjoern's jupyter notebook to perform forward integration:

```
./ggp_notebooks/simulations/forward_integration.ipynb 
```

From (https://github.com/bjks/ggp_notebooks).
NB: the following piece of code is currently being used to call correct files and run the integration.

```{python}
def get_input_files(directory, keyword=None, ext=".csv"):
    entries = os.listdir(directory)
    final_files = []
    if keyword == None:
        for e in entries:
            if e.endswith(ext):
                final_files.append(os.path.join(directory,e))
    else:
        for e in entries:
            if e.endswith(ext) and keyword in e:
                final_files.append(os.path.join(directory,e))   
    return sorted(final_files)   

### get data and set output dir

data_files = get_input_files('/scicore/home/nimwegen/rocasu25/Documents/Projects/biozentrum/gfp_fluctuations_project/denoised_data_complete_cycles_corrected/', keyword='complete.csv')

len(data_files)
data_files
```

And then:

```
import warnings;
warnings.filterwarnings('ignore');

"""
HOW TO SET THE SETTING PARAMETERS
----------------------------------

- growth_noise and q_noise can be either True or False, 
where True will take the inferred lambda/q traces and False takes the mean growth/production rate

-var_dx and var_dg can be either a float, None or "data":
        - float: will enable the binomial sampling noise in gfp and gaussian sampling in log_length with therespective parameters
        - None: will turn off the noise source
        - "data": will 'measure' the ratio in length and gfp of mother and daughter cell in the data and transfer it to the simulation

"""
for data_file in sorted(data_files):
    sample = '_'.join(data_file.split('/')[-1].split('_')[:3])
    print("integrating:", sample)
    cells = ggp_df2cells(pd.read_csv(data_file, skiprows=header_lines(data_file)))
    cells_raw = ggp_df2cells(pd.read_csv(data_file, skiprows=header_lines(data_file)), log_length="log_length", gfp="fp")
    params = read_header(data_file)
    params = read_header(data_file.replace("prediction", "final"))

    out_dir = mk_missing_dir(os.path.join(*data_file.split('/')[:-1], 'integration', sample), depth=1)


    #### inferred var_dx, var_dg, beta
    var_dx = params["var_dx"][1]
    var_dg = params["var_dg"][1]
    beta   = params["beta"][1]
    
    

    # keys gives a name (for labels, file names) to the integration variant
    settings = {"allnoisesources": {"growth_noise": True, 
                                      "q_noise": True, 
                                      "var_dx": "data", 
                                      "var_dg": "data"},
               "nogrowthnoise": {"growth_noise": False, 
                                     "q_noise": True, 
                                      "var_dx": "data", 
                                      "var_dg": "data"},
               "noprodnoise": {"growth_noise": True, 
                                      "q_noise": False, 
                                      "var_dx": "data", 
                                      "var_dg": "data"},
               "nodivisionnoise": {"growth_noise": True, 
                                      "q_noise": True, 
                                      "var_dx": None, 
                                      "var_dg": None},
                "nogrowthdivnoise": {"growth_noise": False, 
                                     "q_noise": True, 
                                      "var_dx": None, 
                                      "var_dg": None},
                "nonoise": {"growth_noise": False, 
                                      "q_noise": False, 
                                      "var_dx": None, 
                                      "var_dg": None}}


    # r
    for s in settings.keys():
        
        out_dir_s = mk_missing_dir(os.path.join(out_dir, s), depth=1)
        cells_integrated = forward_integration(cells, params, 
                                               growth_noise=settings[s]["growth_noise"], 
                                               q_noise=settings[s]["q_noise"], 
                                               var_dx=settings[s]["var_dx"], 
                                               var_dg=settings[s]["var_dg"], 
                                               beta=beta)

        outfile_base = os.path.join(out_dir_s, sample)+'_'+s
        outfile_base = '_'.join(outfile_base.split()) # remove whitespaces

        save_cells(cells_integrated, outfile_base+'.csv')
        # plots
        plot_concentration(cells_integrated, cells,  get_longest_path(cells_integrated), cells_raw=None,
                           label_i=s,plot_file = outfile_base+'.pdf')
        plot_concentration_hist(cells_integrated, cells,
                           label_i=s,
                          plot_file = outfile_base+'_hist.pdf')
        plot_x_g(cells_integrated, cells,  get_longest_path(cells_integrated), cells_raw=None,
                           label_i=s,
                          plot_file = outfile_base+'_x_g.pdf')
        
#     break
```
Output files are stored in `/scicore/home/nimwegen/rocasu25/Documents/Projects/biozentrum/gfp_fluctuations_project/denoised_data_complete_cycles_corrected/integration`.

### Compute autocorrelation functions

We make use of the following python script (again from `gfp_gaussian_process`): `https://github.com/bjks/ggp_notebooks`.

```
./gfp_gaussian_process/python_src/correlation_from_joint.py
```
To run all computation (1 node per dataset), on the cluster, simply run: `./mother_machine_experiments_toolbox/computing_correlation_functions.R`.

### Load all data
- Denoised data together with parameters
- Forward integration data to assess contribution to overall noise
- Correlations

Run `./mother_machine_experiments_toolbox/load_simulated_and_denoised_data.R`

### Make figures

All scripts to make figures are kept in `./figures`. I try to keep it organized by naming script with the same name as their outputs.

## Project structure

The essential scripts and files are listed below:

```
gfp_fluctuations_project
│   README.md
│
└───analysis: contains Rmd files were data are analyzed to adress specific questions, these Rmd files are also used to produce figures.
|      └────figures      
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

