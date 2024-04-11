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
subs = paste0('WM_adapt/data/exp1/',list.files('WM_adapt/data/exp1/',pattern = '^data'))

# row numbers for WM-fixed and WM-random trials for pavlovia output files (difficult to parse via other means)
fixed_trial_nums =   c( 13, 26, 39, 52, 65, 78, 91, 104, 117, 130, 143, 156, 169, 182, 195, 208, 221, 
                        234, 247, 260, 273, 286, 299, 312, 325, 338, 351, 364, 377, 390, 403, 416,
                        429, 442, 455, 468, 481, 494, 507, 520, 533, 546, 559, 572, 585,  598,  611,  
                        624,  637,  650,  663,  676,  689,  702,  715,  728,  741,  754,  767,  780,
                        793,  806,  819,  832,  845,  858, 871,  884,  897,  910,  923,  936,  949,
                        962,  975,  988, 1001, 1014, 1027, 1040, 1053, 1066, 1079, 1092, 1105, 1118,
                        1131, 1144, 1157, 1170, 1183, 1196, 1209, 1222, 1235, 1248, 1261, 1274, 1287, 1300)
random_trial_nums = c(6,   19,   32,   45,   58,   71,   84,   97,  110,  123,  136,  149,  162,  175,  188,
                      201,  214,  227,  240,  253,  266,  279, 292,  305,  318,  331,  344,  357,  370,  383,  
                      396,  409,  422,  435,  448,  461, 474,  487,  500,  513,  526,  539,  552,  565,
                      578,  591,  604,  617,  630,  643,  656,  669,  682,  695,  708,  721,  734,  747,
                      760,  773,  786,  799,  812,  825,  838,  851, 864,  877,  890,  903,  916,  929,  
                      942,  955, 968,  981,  994, 1007, 1020, 1033, 1046, 1059, 1072, 1085, 1098, 1111, 
                      1124, 1137, 1150, 1163, 1176, 1189, 1202, 1215, 1228, 1241, 1254, 1267, 1280, 1293)

grpdf = list()
grprandf = list()
adaptpct = rep(NA,length(subs))
for (s in seq_along(subs)) {
  subdf = read.csv(subs[s]) # load subject data file
  
  # compute accuracy for attention trials immediately prior to fixed WM trial
  attaccuracy = rep(NA, length(fixed_trial_nums))
  attrt = rep(NA,length(fixed_trial_nums))
  for (t in seq_along(fixed_trial_nums)) {
    if (t == 1) {
      attaccuracy[t] = mean(subdf$key_resp.corr[subdf$trialCounter < fixed_trial_nums[t]],na.rm = T)
      attrt[t] = mean(subdf$key_resp.rt[subdf$trialCounter < fixed_trial_nums[t]],na.rm = T)
    } else {
      attaccuracy[t] = mean(subdf$key_resp.corr[(subdf$trialCounter < fixed_trial_nums[t]) & (subdf$trialCounter > fixed_trial_nums[t-1])],
                            na.rm=T)
      attrt[t] = mean(subdf$key_resp.rt[(subdf$trialCounter < fixed_trial_nums[t]) & (subdf$trialCounter > fixed_trial_nums[t-1])],
                      na.rm=T)
    }
  }
  
  # extract WM-fixed trials
  subfixeddf = data.frame(y = subdf$probe_slider_fixed.response[subdf$trialCounter %in% fixed_trial_nums],
                     x = fixed_trial_nums,
                     attacc = attaccuracy,
                     mnrt = attrt,
                     block = factor(c(rep(1:5,each=20))))
  
  # extract WM-random trials
  subrandf = data.frame(y = subdf$probe_slider_random.response[subdf$trialCounter %in% random_trial_nums] -
                          subdf$memTargetPos[subdf$trialCounter %in% random_trial_nums],
                        x = random_trial_nums,
                        block = factor(c(rep(1:5,each=20))),
                        memTargetPos = subdf$memTargetPos[subdf$trialCounter %in% random_trial_nums])
  
  # data quality check
  print(subs[s])
  print(sum(!is.na(subdf$key_resp.rt))/1100 > (2/3)) # Att-error response rate
  print(sum(!is.na(subfixeddf$y))/100 > (2/3)) # WM-fixed response rate
  print(sum(!is.na(subrandf$y))/100 > (2/3)) # WM-random response rate
  print(mean(subdf$key_resp.corr,na.rm=T) > (2/3)) # Att-error accuracy (including non-response as incorrect)
  print(mean(abs(subrandf$y),na.rm=T) < 0.15) # WM-random MAE
  print(all(subdf$key_resp_block_break.rt[!is.na(subdf$key_resp_block_break.rt)]/60 < 10)) # inter-block break time
  print(subdf$quiz_slider_1.response[1411] < 3) # head movement question
  print(subdf$quiz_slider_2.response[1411] >= 8) # eye movement question
  print('')
  
  # add to group list
  grpdf[[s]] = subfixeddf
  grpdf[[s]]$subject = rep(s,nrow(grpdf[[s]]))
  grprandf[[s]] = subrandf
  grprandf[[s]]$subject = rep(s,nrow(grprandf[[s]]))
  
  # compute adaptation magnitude
  tmp = grpdf[[s]] %>% group_by(block) %>% dplyr::summarise(mn = mean(y,na.rm=T)) %>% select(mn) %>% unlist()
  block1 = tmp[1]
  block5 = tmp[5]
  adaptpct[s] = (block5 - block1)/0.17

}

# concatenate subjects into single data frame
grpdf = do.call(rbind, grpdf) 
grprandf = do.call(rbind,grprandf)

# save WM-fixed and WM-random group data frame
if (!file.exists('WM_adapt/data/exp1/exp1_group_fixed.csv')) {
  write.csv(grpdf,
            file = 'WM_adapt/data/exp1/exp1_group_fixed.csv',
            row.names = F)
}
if (!file.exists('WM_adapt/data/exp1/exp1_group_random.csv')) {
  write.csv(grprandf,
            file = 'WM_adapt/data/exp1/exp1_group_random.csv',
            row.names = F)
}

# compute group average WM-fixed recall
mndf = grpdf %>% group_by(x) %>% dplyr::summarize(y = mean(y,na.rm=T),rt = mean(mnrt,na.rm = T),attacc = mean(attacc,na.rm = T),n = n()) %>% filter(n > 1)
mndf$block = factor(rep(seq(1,5),each = 20))

# save group average WM-fixed recall
if (!file.exists('WM_adapt/data/exp1/exp1_group_mean_fixed.csv')) {
  write.csv(mndf,
            file = 'WM_adapt/data/exp1/exp1_group_mean_fixed.csv',
            row.names = F,
            col.names = T)
}

# adaptation magnitude (bootstrap SE)
set.seed(1111)
bootmu = function(sample,i) mean(sample[i])
boot(adaptpct,bootmu,1000)

# Bayes Factor analysis comparing 1st and last block
set.seed(1111)
df_pre_vs_adapt = filter(grpdf,block %in% c(1,5)) %>% droplevels()
df_pre_vs_adapt$block = factor(df_pre_vs_adapt$block)
df_pre_vs_adapt$subject = factor(df_pre_vs_adapt$subject)

bf_full = lmBF(y ~ block + subject, 
               whichRandom = "subject",
               rscaleRandom = 'nuisance',
               data = df_pre_vs_adapt %>% drop_na(y))
bf_int = lmBF(y ~ subject, whichRandom = "subject",
              rscaleRandom = 'nuisance',
              data = df_pre_vs_adapt %>% drop_na(y))

bf_full/bf_int

# timecourse model fits
# linear fit
if (!file.exists('WM_adapt/model_fits/Exp1_linear_fit_group_mean.rds')) {
  fit_linear = brm(y ~ x,
                   data = mndf,
                   prior = set_prior("student_t(1,0,0.2)",class = "b"),
                   iter = 6000,
                   warmup = 2000,
                   chains = 4,
                   control = list(adapt_delta = 0.8))
  saveRDS(fit_linear,file = 'WM_adapt/model_fits/Exp1_linear_fit_group_mean.rds')
} else {
  fit_linear = readRDS('WM_adapt/model_fits/Exp1_linear_fit_group_mean.rds')
}

#  single exponential decay fit (Robinson, Soetedjo & Noto, 2006)
if (!file.exists('WM_adapt/model_fits/Exp1_expdecay_fit_group_mean.rds')) {
  singleexpprior = prior(normal(0,0.2), nlpar = 'amp') +
    prior(normal(0,1000),nlpar = "rate",lb=0) + 
    prior(normal(0.5,0.2),nlpar = "asymptote") 
  
  fit_single = brm(bf(y ~ (amp * (2^(-x/rate))) + asymptote, amp + rate + asymptote ~ 1, nl = TRUE),
                   data = mndf,
                   prior = singleexpprior,
                   iter = 6000,
                   warmup = 2000,
                   chains = 4,
                   control = list(adapt_delta = 0.8))
  saveRDS(fit_single,file = 'WM_adapt/model_fits/Exp1_expdecay_fit_group_mean.rds')
} else {
  fit_single = readRDS('WM_adapt/model_fits/Exp1_expdecay_fit_group_mean.rds')
}

# double exponential
if (!file.exists('WM_adapt/model_fits/Exp1_dblexpdecay_fit_group_mean.rds')) {
  dbl_exp_prior = prior(normal(0,0.2), nlpar = 'amp1',lb = 0) +
    prior(normal(0,0.2), nlpar = 'amp2',lb = 0) + 
    prior(normal(0,50),nlpar = "rate1",lb = 0) +
    prior(normal(0,500),nlpar = "rate2",lb = 0) + 
    prior(normal(0.5,0.2),nlpar = "plateau")
  
  fit_dbl = brm(bf(y ~ (amp1 * (2^(-x/rate1))) + (amp2*(2^(-x/rate2))) + plateau , amp1 + amp2 + rate1 + rate2 + plateau ~ 1, nl = TRUE),
                data = mndf,
                prior = dbl_exp_prior,
                iter = 6000,
                warmup = 2000,
                chains = 4,
                control = list(adapt_delta = 0.8))
  saveRDS(fit_dbl,file = 'WM_adapt/model_fits/Exp1_dblexpdecay_fit_group_mean.rds')
} else {
  fit_dbl = readRDS('WM_adapt/model_fits/Exp1_dblexpdecay_fit_group_mean.rds')
}

# model comparison
exp1_loo = loo(fit_linear,fit_single,fit_dbl)
exp1_loo

# WM-random bin analysis (adaptation field)
ranbindf = grprandf %>% mutate(bin = cut(memTargetPos,breaks = c(-Inf,seq(0.15,0.65,0.1),Inf))) %>%
  group_by(subject,block,bin) %>% dplyr::summarise(y = mean(y,na.rm = T))
blockdiffdf = ranbindf %>% filter(block %in% c('1','5')) %>%
  ungroup %>% spread(block, y) %>% 
  group_by(subject,bin) %>% 
  mutate(ydiff = `5` - `1`)
blockdiffdf$bin_index = blockdiffdf$bin
levels(blockdiffdf$bin_index) = c('-4','-3','-2','-1','0','+1','+2')

# bin distance regression 
blockdiffdf$abs_bin_pos = rep(NA,nrow(blockdiffdf))
blockdiffdf$abs_bin_pos[ blockdiffdf$bin_index == '-3'] = 0.3
blockdiffdf$abs_bin_pos[ blockdiffdf$bin_index == '-2'] = 0.2
blockdiffdf$abs_bin_pos[ blockdiffdf$bin_index == '-1'] = 0.1
blockdiffdf$abs_bin_pos[ blockdiffdf$bin_index == '0'] = 0
blockdiffdf$abs_bin_pos[ blockdiffdf$bin_index == '+1'] = 0.1
blockdiffdf$abs_bin_pos[ blockdiffdf$bin_index == '+2'] = 0.2

set.seed(1111)
blockdiffdf$subject = factor(blockdiffdf$subject)
bf_full = lmBF(ydiff ~ abs_bin_pos + subject,
               whichRandom = "subject",
               data = blockdiffdf %>% drop_na(ydiff))
bf_int = lmBF(ydiff ~ subject,
              whichRandom = "subject",
              data = blockdiffdf %>% drop_na(ydiff))
bf_full/bf_int

post_samples = posterior(bf_full,iterations = 10000)
mean(post_samples[,"abs_bin_pos"]) * 0.1 / 0.17 * 100 # slope * bin_distance / backstep_size * 100
hdi(as.vector(post_samples[,"abs_bin_pos"])) * 0.1 / 0.17 * 100

