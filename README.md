# xpert-ultra

## The model will read data from the following files: 
- ultra_params.csv: a list of all the parameters to be sampled once per simulation (or with fixed values), with a modal estimate and high and low value for the sampled range of each parameter
- lowinc_params.csv: like ultra_params.csv, but with TB specificity values for Ultra calculated only based on results from study sites in countries with TB incidence <100/100,000/year. 
- highinc_params.csv: like ultra_params.csv, but with TB specificity values for Ultra calculated only based on results from study sites in countries with TB incidence >200/100,000/year. (Note that for both lowinc_params.csv and highinc_params.csv, the values for specificity without trace and with repeated trace call are not accurate; the values with trace are merely copied over, as these entries were not required for the results reported in the manuscript.)
-  agesexdistrib.csv: A calcuation of the age and sex distribution of notified adult TB cases for each country setting: This takes the number of notified cases in each decade between ages 15-64, for each sex, and calculates what fraction of total adult cases are in each age-and-sex category, binning age in 5-year increments by splitting each decade in half and applying the same age distribution for cases age >65 as is reported for the country's overall population age >65. 
- agesexmort.csv: supplies individual background mortality rates based on age, sex, and country (using south african data on HIV+ individuals for the South Africa cohort; those south-africa data will be applied to the HIV+ individuals in the other countries also.)


## To run the model: 

Edit first three lines of Xpert_Ultra_run_model_plosmed.R to specificy the setting, any fixed parameters ("NA" if none), and a 'tag' for the version of the parameters or model to run (see list below). 

Can also edit 'date' (line 7) to change the date label on output files, and cohortsize and nsims (lines 9-10) to change the number of patients modeled in the cohort and the number of simulations to be run of the cohort.

Output will be saved to file 'paste0(savelocation,"markovoutput_",runname, ".csv")', ncolumns=length(header), append = TRUE, sep=",")'



## Interpreting results:

The output file, which contains a large number of different results for each simulation of the cohort, can then be processed using the function in 'Xpert Ultra results function plosmed.R', to generate the main results reported in the paper (the number of TB deaths and number of unnecessary treatments for Xpert and for Ultra, and their incremental values [difference between xpert and ultra] and ratio).



## Possible 'tags' (see above) that can be supplied when running the model include: 
- "overall" for the primary analysis
- "highinc" or "lowinc" to use the specificity estimates calculated from study sites in only high- or low-incidence countries, respectively
- "pessimisticrr" for a sensitivity analysis assuming worse treatment outcomes for rifampin-resistant TB
- "chinaempiric" to add some empiric treatment in the China setting
- "changebehavior" for the sensitivity analyses in which clinicians use a second test with Ultra to try to improve specificity, with characteristics similar to xray
- "changebehavior2" to use a highly accurate second test with Ultra to improve specificity with little cost to sensitivity

## Underlying data
The data from the diagnostic accuracy (Dorman, Schumacher, et al 2017) that directly informed this model can be found in the file, Ultra study data for model.pptx.

