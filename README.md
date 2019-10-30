# SDS
Accurate prediction of prognosis is important for cancer patients. Integrating clinical data with molecular data helps increase the accuracy. We propose novel methods for feature selection and feature reduction to lead to improved prognostic performance in case-only analysis of time-to-event endpoint studies.



## Prerequisites
R packages:
* survival
* glmnet
* survcomp
* matrixStats

Data sets (order of subjects have to match in all 3 datasets)
* survival data (data.frame) in form of: OVERALL.SURVIVAL, overall.survival.indicator
* gene expression data (matrix) already filtered and normalized
* clinical data (data.frame) with dummy variables for categorical variables

## Example workflow
An example how to use the package can be found in Example.pdf

## Citation
The paper can be found at this location: 

## Acknowledgments

* KUMC 'omics working group
* This work is supported by the National Cancer Institute (NCI) Cancer Center Support Grant P30 CA168524
