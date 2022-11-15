# tBNMA

This repository contains the code to reproduce the numerical simulations and data analysis from LeBlanc and Banks (2022). 

# Directories

There are four main directories in this repository: Data, Figures, Old Scripts, and jags.  

## Data

The Data directory contains the dataset used in the paper. This dataset was agglomerated from the datasets of multiple review papers: Thom et al. (2015); Liu et al. (2016); Guest et al. (2017); Mccool et al. (2017); Li and Xu (2018); Zhang et al. (2019); Lan et al. (2019); Brown et al. (2021); Feng et al. (2021).  The file is labeled `data_tbnma.csv`, and is a long datastructure.  Each row corresponds to one arm of one trial.  There are a number of covariates included:

* `Study`: The study this treatment arm originates from, in the style of the in-line citations from the review paper from which it was sourced.  
* `Year`: The year this treatment arm took place in.
* `Treatment`: The treatment tested in this treatment arm.
* `Mean Age`: The mean age of participants in this treatment arm.  Not available for all studies.
* `Male`: The proportion of participants that were male in this treatment arm.  Not available for all studies.
* `Success`: The number of successes in this treatment arm.
* `Trials`: The number of trials in this treatment arm. 
* `Source`: The review paper from which the study containing this treatment arm was sourced.  Note that most studies are present in more than one review paper --- when this occurs, we only highlight a single source.  

## Figures

This directory contains all of the figures appearing in the paper, as well as some additional figures which do not.

## Old Scripts

This directory contains scripts which were part of the development process and used in earlier versions of tBNMA models, but which have no direct use now.  They are not guaranteed to run without modifying, e.g., file paths.

## jags

All of the jags scripts for various versions of BNMA models are included.  The four which are relevant are

* `BNMA_Like_Bin_Trial_Multi_Arm.bug`: an implemenation of standard BNMA for binomial likelihoods and multiple treatment arms.
* `BNMA_Like_Bin_Trial_Multi_Arm_Time_Z.bug`: an implemenation of Meta-BNMA for binomial likelihoods, multiple treatment arms, and with time effects only on those treatments with at least $5$ appearences in the dataset.
* `BNMA_Like_Bin_Trial_Multi_Arm_Time_Sigmoidal_Z.bug`: an implemenation of Sig-BNMA for binomial likelihoods, multiple treatment arms, and with time effects only on those treatments with at least $5$ appearences in the dataset.
* `BNMA_Like_Bin_Trial_Multi_Arm_Time_GP_Z.bug`: an implemenation of GP-BNMA for binomial likelihoods, multiple treatment arms, and with time effects only on those treatments with at least $5$ appearences in the dataset.

# Scripts

In addition to the directories, there are three scripts located in the main repository which reproduce the results in the paper.

## `Simulate_BNMA_Bin_MultiArm.R`

This script uses the networks, timepoints, comparisons, and study sizes from the agglomerate dataset to simulate a standard BNMA dataset with no time-varying effects.  The four methods (BNMA, Meta-BNMA, Sig-BNMA, and GP-BNMA) are run on this dataset.  Plots of results are produced.

## `Simulate_BNMA_Bin_MultiArm_Sigmoidal.R`

This script uses the networks, timepoints, comparisons, and study sizes from the agglomerate dataset to simulate a Sig-BNMA dataset with time-varying effects on linezolid.  The four methods (BNMA, Meta-BNMA, Sig-BNMA, and GP-BNMA) are run on this dataset.  Plots of results are produced.

## `Data_Analysis.R`

The four methods (BNMA, Meta-BNMA, Sig-BNMA, and GP-BNMA) are run on the agglomerated dataset.  Plots of results are produced.

# References

* Brown, N., Goodman, A., Horner, C., Jenkins, A. and Brown, E. (2021), ‘Treatment of methicillin-resistant staphylococcus aureus (mrsa): updated guidelines from the uk’, JAC- Antimicrobial Resistance 3
* Diekema, D. J., Pfaller, M. A., Shortridge, D., Zervos, M. and Jones, R. N. (2019), ‘Twenty-year trends in antimicrobial susceptibilities among staphylococcus aureus from the sentry antimicrobial surveillance program’, Open Forum Infectious Diseases 6(Supplement 1), S47–S53.
* Guest, J., Esteban, J., Manganelli, A., Novelli, A., Rizzardini, G. and Serra-Burriel, M. (2017), 20 ‘Comparative efficacy and safety of antibiotics used to treat acute bacterial skin and skin structure infections: Results of a network meta-analysis’, PLOS ONE 12, e0187792
* Feng, J., Xiang, F., Cheng, J., Gou, Y. and Li, J. (2021), ‘Comparative efficacy and safety of vancomycin, linezolid, tedizolid, and daptomycin in treating patients with suspected or proven complicated skin and soft tissue infections: An updated network meta-analysis’, Infectious Diseases and Therapy 10.
* Lan, S.-H., Lin, W.-T., Chang, S.-P., Lu, L.-C., Chao, C.-M., Lai, C.-C. and Wang, J.-H.(2019), ‘Tedizolid versus linezolid for the treatment of acute bacterial skin and skin structure infection: A systematic review and meta-analysis’, Antibiotics 8, 137
* LeBlanc and Banks.  Time-Varying Bayesian Meta-Analysis.  2022. 
*Li, Y. and Xu, W. (2018), ‘Efficacy and safety of linezolid compared with other treatments for skin and soft tissue infections: A meta-analysis’, Bioscience Reports 38, BSR20171125
* Liu, C., Mao, Z., Yang, M., Kang, H., Liu, H., Pan, L., Hu, J., Luo, J. and Zhou, F. (2016), ‘Efficacy and safety of daptomycin for skin and soft tissue infections: A systematic review with trial sequential analysis’, Therapeutics and Clinical Risk Management Volume 12, 1455–1466.
* Mccool, R., Eales, J., Barata, T., Arber, M., Cikalo, M., Fleetwood, K., Glanville, J., Gould,I. and Kauf, T. (2017), ‘Pin17. systematic review and network meta-analysis of tedizolid for the treatment of acute bacterial skin and skin structure infection (absssi) due to methicillin-resistant staphylococcus aureus (mrsa)’, Value in Health 18, A231
* Thom, H., Thompson, J., Scott, D., Halfpenny, N., Sulham, K. and Corey, G. (2015), ‘Comparative efficacy of antibiotics for the treatment of acute bacterial skin and skin structure infections (absssi): A systematic review and network meta-analysis’, Current medical research and opinion 31, 1–34.
* Zhang, Y., Wang, Y., van Driel, M., Mcguire, T., Zhang, T., Dong, Y., Liu, Y., Liu, L., Hao,
24 R., Cao, L., Xing, J. and Dong, Y. (2019), ‘Network meta-analysis and pharmacoeconomic evaluation of antibiotics for the treatment of patients infected with complicated skin and soft structure infection and hospital-acquired or ventilator-associated penumonia’, Antimicrobial Resistance Infection Control 8.
