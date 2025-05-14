# Errors of attention adaptively warp spatial cognition
This repository contains analysis code supporting Brissenden, Yin, Vesia, and Lee (2025) Nature Human Behaviour (https://doi.org/10.1038/s41562-025-02109-5)

## Usage ##
### Data ###
The data supporting the findings of this study are available at https://osf.io/egskw/.

### Analysis ###
Install data from the link above into a directory named `./data`. Create an environmental variable `LOCAL` in a `.Renviron` file that points to the directory where the repository is saved. The analyses and figures presented in the manuscript can be recreated by executing the respective R scripts (e.g. `exp1_analysis.R`) in `./code`.

### Required Packages ###
`brms` (2.21.0)  

`rstan` (2.32.6) 

`BayesFactor` (0.9.12.4.7)  

`tidyverse` (2.0.0)  
  - `ggplot` (3.5.1)  
  - `dplyr` (1.1.4)  
  - `readr` (2.1.5)

 
`tidybayes` (3.0.6)   

`boot` (1.3.30)  

`gazerjb` (0.0.1.2; https://github.com/brissend/gazerjb)  



