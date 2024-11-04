# load libraries
library(rstan)
library(brms)
library(tidybayes)
library(tidyverse)
library(viridis)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# set environment
LOCAL = Sys.getenv("LOCAL")
setwd(LOCAL)

# load data
if (!file.exists('WM_adapt/data/exp1/exp1_group_fixed.csv')) {
  source('WM_adapt/code/exp1_analysis.R') # subject recall computed and saved by this script
}
grpdf = read.csv('WM_adapt/data/exp1/exp1_group_fixed.csv')

# exponential decay function prior
singleexpprior = prior(normal(0,0.2), nlpar = 'amp') +
  prior(normal(0,1000),nlpar = "rate",lb=0) + 
  prior(normal(0.5,0.2),nlpar = "asymptote") 

# Individual subject model fits
fitdf = data.frame()
for (ii in 1:length(unique(grpdf$subject))) {
  

  # single exponential model best fit group data
  set.seed(1111)
  fit_sub_single = brm(bf(y ~ (amp * (2^(-x/rate))) + asymptote, amp + rate + asymptote ~ 1, nl = TRUE),
                       data = filter(grpdf,subject == unique(grpdf$subject)[ii]),
                       prior = singleexpprior,
                       iter = 6000,
                       warmup = 2000,
                       chains = 4,
                       control = list(adapt_delta = 0.8))
  
  
  for (x in seq(min(grpdf$x),max(grpdf$x),by = 20)) {
    
    tmp = epred_draws(fit_sub_single,newdata = data.frame(x = x),ndraws=5000)
    
    fitdf = rbind(fitdf,
                  data.frame(subject = unique(grpdf$subject)[ii],
                             x = x, 
                             y = mean(tmp$.epred),
                             y_low = bayestestR::hdi(tmp$.epred,ci = 0.95)$CI_low,
                             y_hi = bayestestR::hdi(tmp$.epred,ci = 0.95)$CI_high))
  }
  
}

# plot
grpdf$block = factor(grpdf$block)
ggplot(grpdf,aes(x = x, y = (y - 0.5) / 0.17 * 100 )) + 
  geom_point(aes(color = block)) + 
  geom_ribbon(data = fitdf, alpha = 0.25,color = NA,aes(ymin = (y_low - 0.5)/0.17 * 100, ymax = (y_hi - 0.5)/0.17  *100)) + 
  geom_line(data = fitdf,size = 1,aes(x = x, y = (y - 0.5)/0.17 * 100)) + 
  geom_hline(yintercept = 0,linetype = 'dashed') + 
  scale_color_viridis(discrete = T ,begin = 0.05,end = 0.95) + 
  labs(color = 'Block',y = 'Recall Location (% of Backstep)', x = 'Trial') +
  facet_wrap(vars(subject),nrow = 8,ncol = 5,scales = 'free_y') + 
  theme_classic() +
  theme(legend.position = 'none')

