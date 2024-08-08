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
if (!file.exists('WM_adapt/data/exp2/exp2_group_fixed.csv')) {
  source('WM_adapt/code/exp2_analysis.R') # subject recall computed and saved by this script
}
grpdf = read.csv('WM_adapt/data/exp2/exp2_group_fixed.csv')

# exponential decay function prior
singleexpprior = prior(normal(0,0.2), nlpar = 'amp') +
  prior(normal(0,1000),nlpar = "rate",lb=0) + 
  prior(normal(0.5,0.2),nlpar = "asymptote") 

# Individual subject model fits
fitdf = data.frame()
for (ii in 1:length(unique(grpdf$subject))) {
  
  subdf = filter(grpdf, subject == unique(grpdf$subject)[ii])
  subdf$adaptfit = c(rep(1,85),rep(0,25))
  subdf$trial = subdf$x - 100 
  subdf$trial[1:25] = 0 # treat all pre-adapt trials as trial 0
  
  # pre-adapt fit
  set.seed(1111)
  fit_preadapt = brm(y ~ x,
                     data = filter(subdf,phase == 'pre-adapt'),
                     prior = set_prior("student_t(1,0,0.2)",class = "b"),
                     iter = 6000,
                     warmup = 2000,
                     chains = 4,
                     control = list(adapt_delta = 0.8))
  
  # single exponential model best fit group data
  set.seed(1111)
  fit_adapt = brm(bf(y ~ (amp * (2^(-trial/rate))) + asymptote, amp + rate + asymptote ~ 1, nl = TRUE),
                       data = filter(subdf, adaptfit == 1),
                       prior = singleexpprior,
                       iter = 6000,
                       warmup = 2000,
                       chains = 4,
                       control = list(adapt_delta = 0.8))
  
  # post-adapt
  set.seed(1111)
  subdf$trial_post = subdf$x - 880 
  fit_postadapt = brm(y ~ trial_post,
                      data = filter(subdf,phase == 'post-adapt'),
                      prior = set_prior("student_t(1,0,0.2)",class = "b"),
                      iter = 6000,
                      warmup = 2000,
                      chains = 4,
                      control = list(adapt_delta = 0.8))
  
  
  minval = min(filter(subdf,phase == 'pre-adapt')$x)
  maxval = max(filter(subdf,phase == 'pre-adapt')$x)
  for (x in seq(minval,maxval,length.out = 20)) {
    
    tmp = epred_draws(fit_preadapt,newdata = data.frame(x = x), ndraws = 5000)
    
    fitdf = rbind(fitdf,
                  data.frame(subject = unique(grpdf$subject)[ii],
                             x = x, 
                             y = mean(tmp$.epred),
                             y_low = bayestestR::hdi(tmp$.epred,ci = 0.95)$CI_low,
                             y_hi = bayestestR::hdi(tmp$.epred,ci = 0.95)$CI_high,
                             phase = 'Pre-Adapt'))
  }
  
  minval = 4 
  maxval = max(filter(subdf,phase == 'adapt')$trial)
  for (x in seq(minval,maxval,length.out = 100)) {
    
    tmp = epred_draws(fit_adapt,newdata = data.frame(trial = x), ndraws = 5000)
    
    fitdf = rbind(fitdf,
                  data.frame(subject = unique(grpdf$subject)[ii],
                             x = x + 100,  
                             y = mean(tmp$.epred),
                             y_low = bayestestR::hdi(tmp$.epred,ci = 0.95)$CI_low,
                             y_hi = bayestestR::hdi(tmp$.epred,ci = 0.95)$CI_high,
                             phase = 'Adapt'))
  }
  
  minval = min(filter(subdf,phase == 'post-adapt')$x) - 880
  maxval = max(filter(subdf,phase == 'post-adapt')$x) - 880
  for (x in seq(minval,maxval,length.out = 20)) {
    
    tmp = epred_draws(fit_postadapt,newdata = data.frame(trial_post = x), ndraws = 5000)
    
    fitdf = rbind(fitdf,
                  data.frame(subject = unique(grpdf$subject)[ii],
                             x = x + 880, 
                             y = mean(tmp$.epred),
                             y_low = bayestestR::hdi(tmp$.epred,ci = 0.95)$CI_low,
                             y_hi = bayestestR::hdi(tmp$.epred,ci = 0.95)$CI_high,
                             phase = 'Post-Adapt'))
  }
  
  
}

# plot
grpdf$block = factor(grpdf$block)
grpdf$phase = factor(grpdf$phase, levels = c('pre-adapt','adapt','post-adapt'))
ggplot(grpdf,aes(x = x, y = (y - 0.5) / 0.17 * 100 ),group = phase) + 
  geom_point(aes(color = block,shape = phase)) + 
  geom_ribbon(data = fitdf, alpha = 0.25,color = NA,aes(ymin = (y_low - 0.5)/0.17 * 100, ymax = (y_hi - 0.5)/0.17  *100)) + 
  geom_line(data = fitdf,size = 1,aes(x = x, y = (y - 0.5)/0.17 * 100, group = phase)) + 
  geom_hline(yintercept = 0,linetype = 'dashed') + 
  scale_color_viridis(discrete = T ,begin = 0.05,end = 0.95) + 
  labs(color = 'Block',y = 'Reported Location (% of Backstep)', x = 'Trial') +
  facet_wrap(vars(subject),nrow = 8,ncol = 5,scales = 'free_y') + 
  theme_classic() +
  theme(legend.position = 'none')



