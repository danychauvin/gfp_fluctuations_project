# constitexprImportdata

constitexprImportdata is a repository to keep my data import routines.

## List of scripts and description

### loadFunctionsAndPackages.R

Calls the packages and necessary functions to import and transform the data in this project.
Bear in mind: everything should be compatible with R 4.1.0.

### readMMData.R

Warnings: There are currently three versions:
- readMMData_202008.R is compatible with the old deepMoma output (202008).
- readMMData_202111.R is compatible with the newest deepMoma output (202111).
- readMMData_20220113.R is a development version aimed at testing the rotating box (Michael).


#### Input

Imports the output of deepMoma (output folders whose paths are specified in a summary csv file, whose path has to be specified in the variable 'path_to_MM_data_summary').

The summary csv file should have the following format:

```
date,description,f_start,f_end,condition,t_interval,orientation,data_path,promoter,vector
20200728,glucoseaa,0,Inf,glucoseaa,1.5,b,/scicore/home/nimwegen/rocasu25/MM_Data/Dany/20200728/20200728_bottom_chr_hi1_curated,hi1,ch
```

Where f_start and f_end refer to the first and last frame that are to be considered, t_interval is in minutes, and data_path is the absolute path of the folder containing deepMoma output.

#### Description and warnings

Depending on the deepMoma version, readMMData.R had distinct versions. The current version is compatible with deepMoma as in August 2020. This version enforces some filtering on the data: cells not fulfilling at least one of these criteria are discarded:

- Did not touch the top (as defined in the variable vertical_cutoff)
- More than 4 data points
- Division is the cell ending event ('div')
- Parents are known (id != -1)
- All frames are occuring more than 2hours after the beginning or the switching of environmental conditions within the chip.
- Log(length_um) vs time Pearson R2 value >0.95

#### Outputs

readMMData.R outputs two tidied dataframes: 
- myframes_to_mycells:
	+ Each row describes a cell at a frame f/time t. The number of rows is the sum of f_c over c, where f_c is the number of frames accounting for cell c full cell cycle, from birth to division.
- mycells:
	+ Each row describes a single cell and full cell cycles observables. There are as many rows as cells.

myframes_to_mycells and mycells, are representing the exact same set of cells.

These dataframes are then fed to the next script, transformMMData.R

### transformMMData.R

Modifies myframes_to_mycells and mycells, to compute total of number of gfp molecules from raw fluorescence, computes full cell cycles variables such as exponential growth rate, and populational quantities such as pseudo growth-rate and pseudo-concentration.

## How to use these scripts

In your RMD file/R script, simply use:

```
source("./constitexprImportdata/loadFunctionsAndPackages.R")
path_to_MM_data_summary <- "./dataLists/constitexprData.csv"
source("./constitexprImportdata/readMMData.R")
source("./constitexprImportdata/transformMMData.R")
```







