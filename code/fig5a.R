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
if (!file.exists('WM_adapt/data/exp4/exp4_group_mean_fixed.csv')) {
  source('WM_adapt/code/exp4_analysis.R') # group average recall computed and saved by this script
}
mndf = read.csv('WM_adapt/data/exp4/exp4_group_mean_fixed.csv')

# load model fit
fit_adapt = readRDS('WM_adapt/model_fits/Exp4_dblexpdecay_fit_group_mean.rds')

# fit linear model to pre-adapt and post-adapt trials
if (!file.exists('WM_adapt/model_fits/Exp4_preadapt_linear_fit_group_mean.rds')) {
  
  fit_preadapt = brm(memResponse ~ trial,
                     data = filter(mndf,phase == 'pre-adapt'),
                     prior = set_prior("student_t(1,0,4)",class = "b"),
                     iter = 6000,
                     warmup = 2000,
                     chains = 4,
                     control = list(adapt_delta = 0.8))
  saveRDS(fit_preadapt,file = 'WM_adapt/model_fits/Exp4_preadapt_linear_fit_group_mean.rds')
  
} else {
  fit_preadapt = readRDS('WM_adapt/model_fits/Exp4_preadapt_linear_fit_group_mean.rds')
}

mndf$trial_post = mndf$trial - 880 
if (!file.exists('WM_adapt/model_fits/Exp4_postadapt_linear_fit_group_mean.rds')) {
  
  fit_postadapt = brm(y ~ trial_post,
                      data = filter(mndf,phase == 'post-adapt'),
                      prior = set_prior("student_t(1,0,4)",class = "b"),
                      iter = 6000,
                      warmup = 2000,
                      chains = 4,
                      control = list(adapt_delta = 0.8))
  saveRDS(fit_postadapt,file = 'WM_adapt/model_fits/Exp4_postadapt_linear_fit_group_mean.rds')
  
} else {
  fit_postadapt = readRDS('WM_adapt/model_fits/Exp4_postadapt_linear_fit_group_mean.rds')
}

# posterior expected distribution (posterior predictive distribution sans residual)
# adaptation blocks
fitdf = data.frame()
minval = 4 
maxval = max(filter(mndf,phase == 'adapt')$trial_adapt)
for (x in seq(minval,maxval,length.out = 100)) {
  
  tmp = epred_draws(fit_adapt,newdata = data.frame(trial_adapt = x), ndraws = 5000)
  
  fitdf = rbind(fitdf,
                data.frame(trial = x + 100,  
                           memResponse = mean(tmp$.epred),
                           y_low = bayestestR::hdi(tmp$.epred,ci = 0.95)$CI_low,
                           y_hi = bayestestR::hdi(tmp$.epred,ci = 0.95)$CI_high,
                           phase = 'Adapt'))
}

# pre-adapt
minval = min(filter(mndf,phase == 'pre-adapt')$trial)
maxval = max(filter(mndf,phase == 'pre-adapt')$trial)
for (x in seq(minval,maxval,length.out = 20)) {
  
  tmp = epred_draws(fit_preadapt,newdata = data.frame(trial = x), ndraws = 5000)
  
  fitdf = rbind(fitdf,
                data.frame(trial = x, 
                           memResponse = mean(tmp$.epred),
                           y_low = bayestestR::hdi(tmp$.epred,ci = 0.95)$CI_low,
                           y_hi = bayestestR::hdi(tmp$.epred,ci = 0.95)$CI_high,
                           phase = 'Pre-Adapt'))
}

# post-adapt
minval = min(filter(mndf,phase == 'post-adapt')$trial) - 880
maxval = max(filter(mndf,phase == 'post-adapt')$trial) - 880
for (x in seq(minval,maxval,length.out = 20)) {
  
  tmp = epred_draws(fit_postadapt,newdata = data.frame(trial_post = x), ndraws = 5000)
  
  fitdf = rbind(fitdf,
                data.frame(trial = x + 880, 
                           memResponse = mean(tmp$.epred),
                           y_low = bayestestR::hdi(tmp$.epred,ci = 0.95)$CI_low,
                           y_hi = bayestestR::hdi(tmp$.epred,ci = 0.95)$CI_high,
                           phase = 'Post-Adapt'))
}


mndf$phase = factor(mndf$phase, levels = c('pre-adapt','adapt','post-adapt'))
mndf$block = factor(mndf$block)
ggplot(mndf,aes(x = trial, y = memResponse,group = phase )) + 
  geom_point(aes(color = block,shape = phase)) + 
  geom_ribbon(data = fitdf, alpha = 0.25,color = NA,aes(ymin = y_low, ymax = y_hi)) + 
  geom_line(data = fitdf,linewidth = 1,aes(x = trial, y = memResponse)) + 
  geom_hline(yintercept = 9,linetype = 'dashed') + 
  scale_color_viridis(discrete = T ,begin = 0.05,end = 0.95) + 
  scale_shape_discrete(name = 'Phase') + 
  scale_x_continuous(breaks= seq(0,1000,by=250)) + 
  scale_y_continuous(breaks = seq(9.6,6.6,by = -0.6),labels = paste0(seq(20,-80,by = -20),'%'),limits = c(6.6,9.6)) + 
  labs(color = 'Block',y = 'Recall Location (% of Backstep)', x = 'Trial') +
  theme_classic() +
  theme(legend.position = 'none')
