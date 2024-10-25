# Changes in the brain after psilocybin therapy for depression

This repo contains code for the statistical analysis and figure generation associated with the results described in: 
https://www.nature.com/articles/s41591-022-01744-z
<br>

BLUF: The evidence from two independent studies suggests that psilocybin therapy can alleviate depression via an increase in the brain's functional network integration.
<br>

![Energy landscape](figures/psilo_energy_schematic-01.jpg)

## Aims and usage

1. To extract measures of brain network modularity (Q) from preprocessed resting-state fMRI. 
2. Compare the Q of patients with depression before and after psilocybin therapy.
3. Compare individual changes in Q to changes in depression severity. 


## Usage
Add `third_party` to the Matlab path for the code in `functions` to run.

## 1. Modularity estimation
Timecourses are extracted from a set of regions to measure functional connectivity as estimate normalized brain network modularity. eg:

```
% define imaging matrix
imaging_volume = 4D_brain_imaging_brain

% define atlas matrix
atlas_volume = 3D_regional_brain_atlas

% Estimate norm. modulatrity
Q_norm = calculate_static_modularity(imaging_volume, atlas_volume)
```
NB:- Imaging and atlas data need to be in the same space (e.g., MNI)


The modularity estimation was heavily inspired by Karolina Finc's Nature Communications paper on working memory: https://www.nature.com/articles/s41467-020-15631-z. To normalize Q, the maximum true Q is divided by the mean of a distribution of shuffled functional connectivity matrices. 



## 2. Stats & plotting
The extracted Q and Beck Depression Inventory (BDI) measures are in `.mat` files here:
```
data/
|-- psilodep1/
|   |-- dat_1.mat
|
|-- psilodep2/
|   |-- dat_2.mat
```

### BDI 
To recreate the BDI figure that visualizes the reductions in depression severity following psilocybin therapy, run:
```
>>> BDI_figure
```
<img src="figures/BDI_figure.jpg" width="350" />

### Imaging figure
To recreate the brain imaging figures that visualises reductions in modularity (Q) following psilocybin therapy, run:
```
>>> Imaging_figure
```
<img src="figures/psilodep1_Q_figure.jpg" width="350" />
<img src="figures/psilodep2_Q_figure.jpg" width="350" />

### Stats
The additional stats are calculated in:
```
functions/
|-- psilodep1.m
|-- psilodep2.m
```

## This code is no longer supported
After leaving academia, I no longer have access to Matlab licenses! 

### Get in touch if you:
- find a bug 
- would like to translate this repo into an open source language
- want help with applying this sort of analysis to your data.






