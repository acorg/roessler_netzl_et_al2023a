# roessler_netzl_et_al2023a
This repository holds the code for the species comparison between human and hamster SARS-CoV-2 neutralisation data, published by Rössler, Netzl et al. Please cite the original publication if any data or code from this repository is used. 

The repository's DOI was created with Zenodo (https://docs.github.com/en/repositories/archiving-a-github-repository/referencing-and-citing-content). Please cite the publication and repository when using code from this repository.

Raw data can be found in the `data` directory. The code for the analyses shown in the main manuscript can be found in the `code` directory. To obtain a titer table for antigenic map construction, execute the `excel_to_titertable.R` script. To construct maps execute the `make_map.R` script. The STAN model is stored in the `model` directory and executing the `code/model_map_magnitude_distribution.R` script saves the sampled and optimized serum and magnitude effects shown in the main manuscript and SOM.

SOM figures can be found in the `som` directory. 

The `function` directory contains utility functions.

All analyses were performed in R version 4.2.2 (2022-10-31).
R Core Team (2022). R: A language and environment for statistical
  computing. R Foundation for Statistical Computing, Vienna,
  Austria. URL https://www.R-project.org/.
  
Antigenic maps were constructed using the Racmacs package, Version 1.1.35:
Wilks S (2022). _Racmacs: R Antigenic Cartography Macros_. https://acorg.github.io/Racmacs,
  https://github.com/acorg/Racmacs.
  
The bayesian modelling was performed with cmdstanr: 
Gabry J, Češnovar R (2022). _cmdstanr: R Interface to 'CmdStan'_.
  https://mc-stan.org/cmdstanr/, https://discourse.mc-stan.org.
