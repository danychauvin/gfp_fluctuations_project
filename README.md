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

### Denoising using gfp_gaussian_process

Should be installed within this git directory.
https://github.com/bjks/gfp_gaussian_process
Commit 5be92d0.

### Forward integration using ggp_notebooks

Should be installed within this git directory.
https://github.com/bjks/ggp_notebooks
Commit 82ce0a2.

### Microscopy data curation using DeepMoMA

I used the following:

Preprocessing:
`module use /scicore/home/nimwegen/GROUP/Moma/MM_Analysis/builds/unstable_202304\05__FOR_PREPROCESSING_ONLY__mmpreprocesspy_v0.4.0-5b062368/Modules`

Curation:
`module use /scicore/home/nimwegen/GROUP/Moma/MM_Analysis/builds/unstable_20230112__moma_v0.6.0-beta14-02af6055/Modules`

## Dataset handling

Data were treated the following way:

### Curating and tyding raw data 

Preprocessing and curation were conducted using my "curation manager" jupyter notebook (`./mother_machine_experiments_toolbox/mother_machine_curation_manager.ipynb`) and the above mentionned version of deepMoma (preprocessing and curation). Once relevant data (i.e CellStats csv files) were copied to a location specified in `./mother_machine_experiments_toolbox/all_mmexperiments_datalist.csv`.

The list of raw data is available in the following table: `./mother_machine_experiments_toolbox/all_mmexperiments_datalist.csv`. This table was written manually, as I accumulated data.

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

### Filtering raw data for denoising

Done running: `./mother_machine_experiments_toolbox/read_all_raw_MM_Data.R`. Necessary packages will be loaded, and `import_tidy_filter_raw_MMData.R` will be sourced.
A tidied and filtered version of each data (one file per condition, promoter, date) will saved on disk in `./raw_curated_data_complete_cycles`.

NB: `cell_detection_offset` and `time_threshold` are used in `import_tidy_filter_raw_MMData.R` to filter out cell-cycles that occur at the beginning of the experiment, or that get too close to the edge of the growth-lane.

### Denoising raw data

Once filtered raw data are saved in `./mother_machine_experiments_toolbox/raw_curated_data_complete_cycles`, one can run the following script to denoise the data:
`./mother_machine_experiments_toolbox/denoise_raw_data.R`.

This script will call the `gfp_gaussian` software (`https://github.com/bjks/gfp_gaussian_process`, Bjoern Kscheschinski) to denoise the data, using scicore nodes (1 node per job).

Output data will be stored in `./mother_machine_experiments_toolbox/denoising_raw_data`, in a folder that can be identified with a timestamp (e.g `denoising_20230522144637`). Importantly, joints necessary to compute correlation functions are also saved in the same folder.

### Saving denoised data together with raw data on disk

Run `./mother_machine_experiments_toolbox/import_denoised_data.R`. Data will be tidied with all necessary information contained in raw data and then saved in `./denoised_data_complete_cycles`.
In these csv, `gfp_nb_corrected` and `q_corrected` are corrected for the autofluorescence.

NB: the same data will also be saved in the folder: `denoised_data_complete_cycles_corrected` to proceed with forward integration, together with inferred parameters. Beware, in these files, inferred `mean_q` has been corrected for autofluorescence volumic production. It is smaller than inferred `mean_q`, which accounts for all sources of fluorescent molecules. Also, `gfp_nb` and `q` ARE corrected for the autofluorescence.

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

I use of the following python script (again from `gfp_gaussian_process`): `https://github.com/bjks/ggp_notebooks`.

```
./gfp_gaussian_process/python_src/correlation_from_joint.py
```
To run all computation (1 node per dataset), on the cluster, simply run: `./mother_machine_experiments_toolbox/computing_correlation_functions.R`.

## Load data before producing figures

### Load only raw data (filtered)

Run `./mother_machine_experiments_toolbox/load_raw_data.R` to load `raw_data` dataframe in memory.

### Load all data

Run `./mother_machine_experiments_toolbox/load_simulated_and_denoised_data.R` to load the following in memory:
- Denoised data together with parameters
- Forward integration data to assess contribution to overall noise
- Correlations data

## Produce figures

All scripts to make figures are kept in `./figures/scripts`. I try to keep it organized by naming script with the same name as their outputs, which should be in `./figures/pdfs`.
