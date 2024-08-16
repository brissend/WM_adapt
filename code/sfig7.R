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
if (!file.exists('WM_adapt/data/exp4/exp4_nosaccade_stim_group_mean_fixed.csv')) {
  source('WM_adapt/code/exp4_eyetrack_analysis.R')
}
mndf_nosaccade_stim = read_csv('WM_adapt/data/exp4/exp4_nosaccade_stim_group_mean_fixed.csv')
mndf_nosaccade_delay = read_csv('WM_adapt/data/exp4/exp4_nosaccade_delay_group_mean_fixed.csv')

# ============================================
# Fit models (stimulus presentation + 100 ms)
# ============================================
# linear fit
if (!file.exists('WM_adapt/model_fits/Exp4_nosaccade_stim_linear_fit_group_mean.rds')) {
  fit_linear = brm(memResponse ~ trial_adapt,
                   data = filter(mndf_nosaccade_stim,adaptfit == 1),
                   prior = set_prior("student_t(1,0,0.2)",class = "b"),
                   iter = 6000,
                   warmup = 2000,
                   chains = 4,
                   control = list(adapt_delta = 0.8))
  saveRDS(fit_linear,file = 'WM_adapt/model_fits/Exp4_nosaccade_stim_linear_fit_group_mean.rds')
} else {
  fit_linear = readRDS('WM_adapt/model_fits/Exp4_nosaccade_stim_linear_fit_group_mean.rds')
}

#  single exponential decay fit (Robinson, Soetedjo & Noto, 2006)
if (!file.exists('WM_adapt/model_fits/Exp4_nosaccade_stim_expdecay_fit_group_mean.rds')) {
  singleexpprior = prior(normal(0,5), nlpar = 'amp') +
    prior(normal(0,1000),nlpar = "rate",lb=0) +
    prior(normal(9,5),nlpar = "asymptote")
  
  fit_single = brm(bf(memResponse ~ (amp * (2^(-trial_adapt/rate))) + asymptote, amp + rate + asymptote ~ 1, nl = TRUE),
                   data = filter(mndf_nosaccade_stim,adaptfit == 1),
                   prior = singleexpprior,
                   iter = 6000,
                   warmup = 2000,
                   chains = 4,
                   control = list(adapt_delta = 0.8))
  saveRDS(fit_single,file = 'WM_adapt/model_fits/Exp4_nosaccade_stim_expdecay_fit_group_mean.rds')
} else {
  fit_single = readRDS('WM_adapt/model_fits/Exp4_nosaccade_stim_expdecay_fit_group_mean.rds')
}

# double exponential
if (!file.exists('WM_adapt/model_fits/Exp4_nosaccade_stim_dblexpdecay_fit_group_mean.rds')) {
  dbl_exp_prior = prior(normal(0,5), nlpar = 'amp1') +
    prior(normal(0,5), nlpar = 'amp2') + 
    prior(normal(0,50),nlpar = "rate1",lb = 0) +
    prior(normal(0,1000),nlpar = "rate2",lb = 0) + 
    prior(normal(9,5),nlpar = "plateau")
  
  fit_dbl = brm(bf(memResponse ~ (amp1 * (2^(-trial_adapt/rate1))) + (amp2*(2^(-trial_adapt/rate2))) + plateau , amp1 + amp2 + rate1 + rate2 + plateau ~ 1, nl = TRUE),
                data = filter(mndf_nosaccade_stim,adaptfit == 1),
                prior = dbl_exp_prior,
                iter = 6000,
                warmup = 2000,
                chains = 4,
                control = list(adapt_delta = 0.8))
  saveRDS(fit_dbl,file = 'WM_adapt/model_fits/Exp4_nosaccade_stim_dblexpdecay_fit_group_mean.rds')
} else {
  fit_dbl = readRDS('WM_adapt/model_fits/Exp4_nosaccade_stim_dblexpdecay_fit_group_mean.rds')
}

# model comparison
exp4_nosaccade_stim_loo = loo(fit_linear,fit_single,fit_dbl)
exp4_nosaccade_stim_loo
fit_adapt = fit_dbl # double exponential fits best


# fit linear model to pre-adapt and post-adapt trials
if (!file.exists('WM_adapt/model_fits/Exp4_nosaccade_stim_preadapt_linear_fit_group_mean.rds')) {
  
  fit_preadapt = brm(memResponse ~ trial,
                     data = filter(mndf_nosaccade_stim,phase == 'pre-adapt'),
                     prior = set_prior("student_t(1,0,4)",class = "b"),
                     iter = 6000,
                     warmup = 2000,
                     chains = 4,
                     control = list(adapt_delta = 0.8))
  saveRDS(fit_preadapt,file = 'WM_adapt/model_fits/Exp4_nosaccade_stim_preadapt_linear_fit_group_mean.rds')
  
} else {
  fit_preadapt = readRDS('WM_adapt/model_fits/Exp4_nosaccade_stim_preadapt_linear_fit_group_mean.rds')
}

mndf_nosaccade_stim$trial_post = mndf_nosaccade_stim$trial - 880 
if (!file.exists('WM_adapt/model_fits/Exp4_nosaccade_stim_postadapt_linear_fit_group_mean.rds')) {
  
  fit_postadapt = brm(memResponse ~ trial_post,
                      data = filter(mndf_nosaccade_stim,phase == 'post-adapt'),
                      prior = set_prior("student_t(1,0,4)",class = "b"),
                      iter = 6000,
                      warmup = 2000,
                      chains = 4,
                      control = list(adapt_delta = 0.8))
  saveRDS(fit_postadapt,file = 'WM_adapt/model_fits/Exp4_nosaccade_stim_postadapt_linear_fit_group_mean.rds')
  
} else {
  fit_postadapt = readRDS('WM_adapt/model_fits/Exp4_nosaccade_stim_postadapt_linear_fit_group_mean.rds')
}

# posterior expected distribution (posterior predictive distribution sans residual)
# adaptation blocks
fitdf = data.frame()
minval = 4 
maxval = max(filter(mndf_nosaccade_stim,phase == 'adapt')$trial_adapt)
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
minval = min(filter(mndf_nosaccade_stim,phase == 'pre-adapt')$trial)
maxval = max(filter(mndf_nosaccade_stim,phase == 'pre-adapt')$trial)
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
minval = min(filter(mndf_nosaccade_stim,phase == 'post-adapt')$trial) - 880
maxval = max(filter(mndf_nosaccade_stim,phase == 'post-adapt')$trial) - 880
for (x in seq(minval,maxval,length.out = 20)) {
  
  tmp = epred_draws(fit_postadapt,newdata = data.frame(trial_post = x), ndraws = 5000)
  
  fitdf = rbind(fitdf,
                data.frame(trial = x + 880, 
                           memResponse = mean(tmp$.epred),
                           y_low = bayestestR::hdi(tmp$.epred,ci = 0.95)$CI_low,
                           y_hi = bayestestR::hdi(tmp$.epred,ci = 0.95)$CI_high,
                           phase = 'Post-Adapt'))
}


mndf_nosaccade_stim$phase = factor(mndf_nosaccade_stim$phase, levels = c('pre-adapt','adapt','post-adapt'))
mndf_nosaccade_stim$block = factor(mndf_nosaccade_stim$block)
ggplot(mndf_nosaccade_stim,aes(x = trial, y = memResponse,group = phase )) + 
  geom_point(aes(color = block,shape = phase)) + 
  geom_ribbon(data = fitdf, alpha = 0.25,color = NA,aes(ymin = y_low, ymax = y_hi)) + 
  geom_line(data = fitdf,linewidth = 1,aes(x = trial, y = memResponse)) + 
  geom_hline(yintercept = 9,linetype = 'dashed') + 
  scale_color_viridis(discrete = T ,begin = 0.05,end = 0.95) + 
  scale_shape_discrete(name = 'Phase') + 
  scale_x_continuous(breaks= seq(0,1000,by=250)) + 
  scale_y_continuous(breaks = seq(9.6,6.0,by = -0.6),labels = paste0(seq(20,-100,by = -20),'%'),limits = c(6.0,9.9)) + 
  labs(color = 'Block',y = 'Recall Location (% of Backstep)', x = 'Trial') +
  theme_classic() +
  theme(legend.position = 'none')


# ============================================
# Fit models (stimulus presentation + delay period)
# ============================================
# linear fit
if (!file.exists('WM_adapt/model_fits/Exp4_nosaccade_delay_linear_fit_group_mean.rds')) {
  fit_linear = brm(memResponse ~ trial_adapt,
                   data = filter(mndf_nosaccade_delay,adaptfit == 1),
                   prior = set_prior("student_t(1,0,0.2)",class = "b"),
                   iter = 6000,
                   warmup = 2000,
                   chains = 4,
                   control = list(adapt_delta = 0.8))
  saveRDS(fit_linear,file = 'WM_adapt/model_fits/Exp4_nosaccade_delay_linear_fit_group_mean.rds')
} else {
  fit_linear = readRDS('WM_adapt/model_fits/Exp4_nosaccade_delay_linear_fit_group_mean.rds')
}

#  single exponential decay fit (Robinson, Soetedjo & Noto, 2006)
if (!file.exists('WM_adapt/model_fits/Exp4_nosaccade_delay_expdecay_fit_group_mean.rds')) {
  singleexpprior = prior(normal(0,5), nlpar = 'amp') +
    prior(normal(0,1000),nlpar = "rate",lb=0) +
    prior(normal(9,5),nlpar = "asymptote")
  
  fit_single = brm(bf(memResponse ~ (amp * (2^(-trial_adapt/rate))) + asymptote, amp + rate + asymptote ~ 1, nl = TRUE),
                   data = filter(mndf_nosaccade_delay,adaptfit == 1),
                   prior = singleexpprior,
                   iter = 6000,
                   warmup = 2000,
                   chains = 4,
                   control = list(adapt_delta = 0.8))
  saveRDS(fit_single,file = 'WM_adapt/model_fits/Exp4_nosaccade_delay_expdecay_fit_group_mean.rds')
} else {
  fit_single = readRDS('WM_adapt/model_fits/Exp4_nosaccade_delay_expdecay_fit_group_mean.rds')
}

# double exponential
if (!file.exists('WM_adapt/model_fits/Exp4_nosaccade_delay_dblexpdecay_fit_group_mean.rds')) {
  dbl_exp_prior = prior(normal(0,5), nlpar = 'amp1') +
    prior(normal(0,5), nlpar = 'amp2') + 
    prior(normal(0,50),nlpar = "rate1",lb = 0) +
    prior(normal(0,1000),nlpar = "rate2",lb = 0) + 
    prior(normal(9,5),nlpar = "plateau")
  
  fit_dbl = brm(bf(memResponse ~ (amp1 * (2^(-trial_adapt/rate1))) + (amp2*(2^(-trial_adapt/rate2))) + plateau , amp1 + amp2 + rate1 + rate2 + plateau ~ 1, nl = TRUE),
                data = filter(mndf_nosaccade_delay,adaptfit == 1),
                prior = dbl_exp_prior,
                iter = 6000,
                warmup = 2000,
                chains = 4,
                control = list(adapt_delta = 0.8))
  saveRDS(fit_dbl,file = 'WM_adapt/model_fits/Exp4_nosaccade_delay_dblexpdecay_fit_group_mean.rds')
} else {
  fit_dbl = readRDS('WM_adapt/model_fits/Exp4_nosaccade_delay_dblexpdecay_fit_group_mean.rds')
}

# model comparison
exp4_nosaccade_delay_loo = loo(fit_linear,fit_single,fit_dbl)
exp4_nosaccade_delay_loo
fit_adapt = fit_single # exponential models indistinguishable (using simpler single exponential model)


# fit linear model to pre-adapt and post-adapt trials
if (!file.exists('WM_adapt/model_fits/Exp4_nosaccade_delay_preadapt_linear_fit_group_mean.rds')) {
  
  fit_preadapt = brm(memResponse ~ trial,
                     data = filter(mndf_nosaccade_delay,phase == 'pre-adapt'),
                     prior = set_prior("student_t(1,0,4)",class = "b"),
                     iter = 6000,
                     warmup = 2000,
                     chains = 4,
                     control = list(adapt_delta = 0.8))
  saveRDS(fit_preadapt,file = 'WM_adapt/model_fits/Exp4_nosaccade_delay_preadapt_linear_fit_group_mean.rds')
  
} else {
  fit_preadapt = readRDS('WM_adapt/model_fits/Exp4_nosaccade_delay_preadapt_linear_fit_group_mean.rds')
}

mndf_nosaccade_delay$trial_post = mndf_nosaccade_delay$trial - 880 
if (!file.exists('WM_adapt/model_fits/Exp4_nosaccade_delay_postadapt_linear_fit_group_mean.rds')) {
  
  fit_postadapt = brm(memResponse ~ trial_post,
                      data = filter(mndf_nosaccade_delay,phase == 'post-adapt'),
                      prior = set_prior("student_t(1,0,4)",class = "b"),
                      iter = 6000,
                      warmup = 2000,
                      chains = 4,
                      control = list(adapt_delta = 0.8))
  saveRDS(fit_postadapt,file = 'WM_adapt/model_fits/Exp4_nosaccade_delay_postadapt_linear_fit_group_mean.rds')
  
} else {
  fit_postadapt = readRDS('WM_adapt/model_fits/Exp4_nosaccade_delay_postadapt_linear_fit_group_mean.rds')
}

# posterior expected distribution (posterior predictive distribution sans residual)
# adaptation blocks
fitdf = data.frame()
minval = 4 
maxval = max(filter(mndf_nosaccade_delay,phase == 'adapt')$trial_adapt)
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
minval = min(filter(mndf_nosaccade_delay,phase == 'pre-adapt')$trial)
maxval = max(filter(mndf_nosaccade_delay,phase == 'pre-adapt')$trial)
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
minval = min(filter(mndf_nosaccade_delay,phase == 'post-adapt')$trial) - 880
maxval = max(filter(mndf_nosaccade_delay,phase == 'post-adapt')$trial) - 880
for (x in seq(minval,maxval,length.out = 20)) {
  
  tmp = epred_draws(fit_postadapt,newdata = data.frame(trial_post = x), ndraws = 5000)
  
  fitdf = rbind(fitdf,
                data.frame(trial = x + 880, 
                           memResponse = mean(tmp$.epred),
                           y_low = bayestestR::hdi(tmp$.epred,ci = 0.95)$CI_low,
                           y_hi = bayestestR::hdi(tmp$.epred,ci = 0.95)$CI_high,
                           phase = 'Post-Adapt'))
}


mndf_nosaccade_delay$phase = factor(mndf_nosaccade_delay$phase, levels = c('pre-adapt','adapt','post-adapt'))
mndf_nosaccade_delay$block = factor(mndf_nosaccade_delay$block)
ggplot(mndf_nosaccade_delay,aes(x = trial, y = memResponse,group = phase )) + 
  geom_point(aes(color = block,shape = phase)) + 
  geom_ribbon(data = fitdf, alpha = 0.25,color = NA,aes(ymin = y_low, ymax = y_hi)) + 
  geom_line(data = fitdf,linewidth = 1,aes(x = trial, y = memResponse)) + 
  geom_hline(yintercept = 9,linetype = 'dashed') + 
  scale_color_viridis(discrete = T ,begin = 0.05,end = 0.95) + 
  scale_shape_discrete(name = 'Phase') + 
  scale_x_continuous(breaks= seq(0,1000,by=250)) + 
  scale_y_continuous(breaks = seq(9.6,6.0,by = -0.6),labels = paste0(seq(20,-100,by = -20),'%'),limits = c(6.0,9.9)) + 
  labs(color = 'Block',y = 'Recall Location (% of Backstep)', x = 'Trial') +
  theme_classic() +
  theme(legend.position = 'none')

