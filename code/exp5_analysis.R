# load libraries
library(rstan)
library(brms)
library(BayesFactor)
library(tidybayes)
library(tidyverse)
library(boot)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# set environment
LOCAL = Sys.getenv("LOCAL")
setwd(LOCAL)

# load data
subs = list.files(path = "WM_adapt/data/exp5/", include.dirs = TRUE, pattern = '[0-9]{6}[A-Za-z]{2}$')

grpdf = list()
grprandf = list()
adaptpct = rep(NA,length(subs))
washoutpct = rep(NA,length(subs))
for (s in seq_along(subs)) {
  subdf = read.csv(file = file.path('WM_adapt/data/exp5',subs[s],paste0(subs[s],'.csv')),na.strings = 'nan')
  subdf$blockNum = factor(subdf$blockNum)
  subdf$trial = 1:nrow(subdf)
  subdf$subject = rep(str_extract(subs[s],pattern = '[A-Z]{2}$'),nrow(subdf))
  
  # WM-fixed
  sub_df_fixed = filter(subdf,trialType == 'memory')
  grpdf[[s]] = sub_df_fixed
  
  # WM-random
  sub_df_random = filter(subdf, trialType == 'memory_random')
  sub_df_random$error = sub_df_random$memResponse - sub_df_random$memTargetPos
  grprandf[[s]] = sub_df_random
  
  # compute adaptation and washout magnitude
  tmp = grpdf[[s]] %>% group_by(blockNum) %>% dplyr::summarise(mn = mean(memResponse,na.rm=T)) %>% select(mn) %>% unlist()
  block1 = tmp[1]
  block4 = tmp[4]
  block5 = tmp[5]
  adaptpct[s] = (block4 - block1)/3 # back step size = 3° of visual angle
  washoutpct[s] = (block5 - block4)/3
}

# concatenate subjects into single data frame
grpdf = do.call(rbind, grpdf) 
grprandf = do.call(rbind,grprandf)

# save WM-fixed and WM-random group data frame
if (!file.exists('WM_adapt/data/exp5/exp5_group_fixed.csv')) {
  write.csv(grpdf,
            file = 'WM_adapt/data/exp5/exp5_group_fixed.csv',
            row.names = F)
}
if (!file.exists('WM_adapt/data/exp5/exp5_group_random.csv')) {
  write.csv(grprandf,
            file = 'WM_adapt/data/exp5/exp5_group_random.csv',
            row.names = F)
}

# compute group average WM-fixed recall
mndf = grpdf %>% 
  group_by(trial) %>% 
  dplyr::summarize(memResponse = mean(memResponse,na.rm=T),n = n()) %>% 
  filter(n > 1)
mndf$block = factor(c(rep(1,25),rep(2:4,each=20),rep(5,25)))
mndf$phase = factor(c(rep('pre-adapt',25),rep('adapt',60),rep('post-adapt',25)))

# adaptation magnitude (bootstrap SE)
bootmu = function(sample,i) mean(sample[i])
set.seed(1111)
boot(adaptpct,bootmu,100)

# washout magnitude
set.seed(1111)
boot(washoutpct,bootmu,100)

# correlation Bayes factor between adaptation and washout magnitude
set.seed(1111)
cor.test(adaptpct,washoutpct)
correlationBF(adaptpct,washoutpct)

# Bayes Factor analysis comparing pre-adapt block and last adaptation block
set.seed(1111)
df_pre_vs_adapt = filter(grpdf,blockNum %in% c(1,4)) %>% droplevels()
df_pre_vs_adapt$blockNum = factor(df_pre_vs_adapt$blockNum)
df_pre_vs_adapt$subject = factor(df_pre_vs_adapt$subject)

bf_full = lmBF(memResponse ~ blockNum + subject, 
               whichRandom = "subject",
               rscaleRandom = 'nuisance',
               data = df_pre_vs_adapt %>% drop_na(memResponse))
bf_int = lmBF(memResponse ~ subject, whichRandom = "subject",
              rscaleRandom = 'nuisance',
              data = df_pre_vs_adapt %>% drop_na(memResponse))

bf_full/bf_int

# Bayes Factor analysis comparing last adaptation block and post-adaptation block
set.seed(1111)
df_adapt_vs_post = filter(grpdf,blockNum %in% c(4,5)) %>% droplevels()
df_adapt_vs_post$block = factor(df_adapt_vs_post$blockNum)
df_adapt_vs_post$subject = factor(df_adapt_vs_post$subject)

bf_full = lmBF(memResponse ~ blockNum + subject, 
               whichRandom = "subject",
               rscaleRandom = 'nuisance',
               data = df_adapt_vs_post %>% drop_na(memResponse))
bf_int = lmBF(memResponse ~ subject, whichRandom = "subject",
              rscaleRandom = 'nuisance',
              data = df_adapt_vs_post %>% drop_na(memResponse))

bf_full/bf_int

# Bayes factor analysis examining linear effect of trial in pre-adapt phase
bf_full = lmBF(memResponse ~ trial + subject,
               whichRandom = "subject",
               data = filter(grpdf,phase == 'pre-adapt') %>% drop_na(memResponse))
bf_int = lmBF(memResponse ~ subject,
              whichRandom = "subject",
              data = filter(grpdf,phase == 'pre-adapt') %>% drop_na(memResponse))

bf_full/bf_int

# timecourse model fits
mndf$adaptfit = c(rep(1,85),rep(0,25))
mndf$trial_adapt = mndf$trial - 100 # renumber so that trial 1 corresponds w/ onset of attentional errors
mndf$trial_adapt[1:25] = 0 # treat all pre-adapt trials as trial 0

# save group average WM-fixed recall
if (!file.exists('WM_adapt/data/exp5/exp5_group_mean_fixed.csv')) {
  write.csv(mndf,
            file = 'WM_adapt/data/exp5/exp5_group_mean_fixed.csv',
            row.names = F)
}

# linear fit
if (!file.exists('WM_adapt/model_fits/Exp5_linear_fit_group_mean.rds')) {
  fit_linear = brm(memResponse ~ trial_adapt,
                   data = filter(mndf,adaptfit == 1),
                   prior = set_prior("student_t(1,0,0.2)",class = "b"),
                   iter = 6000,
                   warmup = 2000,
                   chains = 4,
                   control = list(adapt_delta = 0.8))
  saveRDS(fit_linear,file = 'WM_adapt/model_fits/Exp5_linear_fit_group_mean.rds')
} else {
  fit_linear = readRDS('WM_adapt/model_fits/Exp5_linear_fit_group_mean.rds')
}

#  single exponential decay fit (Robinson, Soetedjo & Noto, 2006)
if (!file.exists('WM_adapt/model_fits/Exp5_expdecay_fit_group_mean.rds')) {
  singleexpprior = prior(normal(0,5), nlpar = 'amp') +
    prior(normal(0,1000),nlpar = "rate",lb=0) +
    prior(normal(-9,5),nlpar = "asymptote")
  
  fit_single = brm(bf(memResponse ~ (amp * (2^(-trial_adapt/rate))) + asymptote, amp + rate + asymptote ~ 1, nl = TRUE),
                   data = filter(mndf,adaptfit == 1),
                   prior = singleexpprior,
                   iter = 6000,
                   warmup = 2000,
                   chains = 4,
                   control = list(adapt_delta = 0.8))
  saveRDS(fit_single,file = 'WM_adapt/model_fits/Exp5_expdecay_fit_group_mean.rds')
} else {
  fit_single = readRDS('WM_adapt/model_fits/Exp5_expdecay_fit_group_mean.rds')
}

# double exponential
if (!file.exists('WM_adapt/model_fits/Exp5_dblexpdecay_fit_group_mean.rds')) {
  dbl_exp_prior = prior(normal(0,5), nlpar = 'amp1') +
    prior(normal(0,5), nlpar = 'amp2') + 
    prior(normal(0,50),nlpar = "rate1",lb = 0) +
    prior(normal(0,1000),nlpar = "rate2",lb = 0) + 
    prior(normal(-9,5),nlpar = "plateau")
  
  fit_dbl = brm(bf(memResponse ~ (amp1 * (2^(-trial_adapt/rate1))) + (amp2*(2^(-trial_adapt/rate2))) + plateau , amp1 + amp2 + rate1 + rate2 + plateau ~ 1, nl = TRUE),
                data = filter(mndf,adaptfit == 1),
                prior = dbl_exp_prior,
                iter = 6000,
                warmup = 2000,
                chains = 4,
                control = list(adapt_delta = 0.8))
  saveRDS(fit_dbl,file = 'WM_adapt/model_fits/Exp5_dblexpdecay_fit_group_mean.rds')
} else {
  fit_dbl = readRDS('WM_adapt/model_fits/Exp5_dblexpdecay_fit_group_mean.rds')
}

# model comparison
exp5_loo = loo(fit_linear,fit_single,fit_dbl)
exp5_loo

# WM-random bin analysis (adaptation field)
ranbindf = grprandf %>% mutate(bin = cut(memTargetPos,breaks = c(-Inf,seq(-10.5,-4.5,length.out = 7),Inf))) %>%
  dplyr::group_by(subject,blockNum,bin) %>% dplyr::summarise(error = mean(error,na.rm = T))
blockdiffdf = ranbindf %>% filter(blockNum %in% c('1','4')) %>%
  ungroup %>% spread(blockNum, error) %>% 
  group_by(subject,bin) %>% 
  mutate(memdiff = `4` - `1`)
blockdiffdf$bin_index = blockdiffdf$bin
levels(blockdiffdf$bin_index) = c('+2','+1','0','-1','-2','-3','-4','-5')

# bin distance regression 
blockdiffdf$abs_bin_pos = rep(NA,nrow(blockdiffdf))
blockdiffdf$abs_bin_pos[ blockdiffdf$bin_index == '-5'] = 5
blockdiffdf$abs_bin_pos[ blockdiffdf$bin_index == '-4'] = 4
blockdiffdf$abs_bin_pos[ blockdiffdf$bin_index == '-3'] = 3
blockdiffdf$abs_bin_pos[ blockdiffdf$bin_index == '-2'] = 2
blockdiffdf$abs_bin_pos[ blockdiffdf$bin_index == '-1'] = 1
blockdiffdf$abs_bin_pos[ blockdiffdf$bin_index == '0'] = 0
blockdiffdf$abs_bin_pos[ blockdiffdf$bin_index == '+1'] = 1
blockdiffdf$abs_bin_pos[ blockdiffdf$bin_index == '+2'] = 2

set.seed(1111)
blockdiffdf$subject = factor(blockdiffdf$subject)
bf_full = lmBF(memdiff ~ abs_bin_pos + subject,
               whichRandom = "subject",
               data = blockdiffdf %>% drop_na(memdiff))
bf_int = lmBF(memdiff ~ subject,
              whichRandom = "subject",
              data = blockdiffdf %>% drop_na(memdiff))
bf_full/bf_int

post_samples = posterior(bf_full,iterations = 10000)
mean(post_samples[,"abs_bin_pos"]) * 1 / 3 * 100 # slope * bin_distance / backstep_size * 100
hdi(as.vector(post_samples[,"abs_bin_pos"])) * 1 / 3 * 100

# exp 4 and 5 combined analysis
grpdf4 = read.csv('WM_adapt/data/exp4/exp4_group_fixed.csv')
grpdf5 = read.csv('WM_adapt/data/exp5/exp5_group_fixed.csv')
grpdf4$hemifield = 'Right'
grpdf5$hemifield = 'Left'
grpdf5$memResponse = grpdf5$memResponse * -1 # flip sign
grpdf = rbind(grpdf4,grpdf5)

set.seed(1111)
df_pre_vs_adapt = filter(grpdf,blockNum %in% c(1,4)) %>% droplevels()
df_pre_vs_adapt$blockNum = factor(df_pre_vs_adapt$blockNum)
df_pre_vs_adapt$subject = factor(df_pre_vs_adapt$subject)
df_pre_vs_adapt$hemifield = factor(df_pre_vs_adapt$hemifield)

bf_full = lmBF(memResponse ~ blockNum*hemifield + subject, 
               whichRandom = "subject",
               rscaleRandom = 'nuisance',
               data = df_pre_vs_adapt %>% drop_na(memResponse))
bf_main = lmBF(memResponse ~ blockNum + hemifield + subject, whichRandom = "subject",
              rscaleRandom = 'nuisance',
              data = df_pre_vs_adapt %>% drop_na(memResponse))
bf_full/bf_main

post_samples = posterior(bf_full,iterations = 10000)
summary(post_samples)
rightdiff = post_samples[,"blockNum:hemifield-4.&.Right"] - post_samples[,"blockNum:hemifield-1.&.Right"]
leftdiff = post_samples[,"blockNum:hemifield-4.&.Left"] - post_samples[,"blockNum:hemifield-1.&.Left"]
mean(rightdiff - leftdiff) / 3 * 100
hdi(as.vector(rightdiff - leftdiff)) / 3 * 100
