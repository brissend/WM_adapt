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
subs = paste0('WM_adapt/data/exp3/',list.files('WM_adapt/data/exp3/'))

# trial numbers for WM-fixed and WM-random trials for pavlovia output files (difficult to parse via other means)
fixed_trial_nums =  c(4,   8,  12,  16,  20,  24,  28,  32,  36,  40,  44,  48,  52,  56,  60,  64,  68,  72,  76,  80,  84,  88,  92,  96,
                      100, 113, 126, 139, 152, 165, 178, 191, 204, 217, 230, 243, 256, 269, 282, 295, 308, 321, 334, 347, 360, 373, 386, 399,
                      412, 425, 438, 451, 464, 477, 490, 503, 516, 529, 542, 555, 568, 581, 594, 607, 620, 633, 646, 659, 672, 685, 698, 711,
                      724, 737, 750, 763, 776, 789, 802, 815, 828, 841, 854, 867, 880, 884, 888, 892, 896, 900, 904, 908, 912, 916, 920, 924,
                      928, 932, 936, 940, 944, 948, 952, 956, 960, 964, 968, 972, 976, 980)
random_trial_nums = c(1,   2,   3,   5,   6,  7,   9,  10,  11,  13,  14,  15,  17,  18,  19,  21,  22,  23,  25,  26,  27,  29,  30,  31,
                      33,  34,  35,  37,  38,  39,  41,  42,  43,  45,  46,  47,  49,  50,  51,  53,  54,  55,  57,  58,  59,  61,  62,  63,
                      65,  66,  67,  69,  70,  71,  73,  74,  75,  77,  78,  79,  81,  82,  83,  85,  86,  87,  89,  90,  91,  93,  94,  95,
                      97,  98,  99, 106, 119, 132, 145, 158, 171, 184, 197, 210, 223, 236, 249, 262, 275, 288, 301, 314, 327, 340, 353, 366,
                      379, 392, 405, 418, 431, 444, 457, 470, 483, 496, 509, 522, 535, 548, 561, 574, 587, 600, 613, 626, 639, 652, 665, 678,
                      691, 704, 717, 730, 743, 756, 769, 782, 795, 808, 821, 834, 847, 860, 873, 881, 882, 883, 885, 886, 887, 889, 890, 891,
                      893, 894, 895, 897, 898, 899, 901, 902, 903, 905, 906, 907, 909, 910, 911, 913, 914, 915, 917, 918, 919, 921, 922, 923,
                      925, 926, 927, 929, 930, 931, 933, 934, 935, 937, 938, 939, 941, 942, 943, 945, 946, 947, 949, 950, 951, 953, 954, 955,
                      957, 958, 959, 961, 962, 963, 965, 966, 967, 969, 970, 971, 973, 974, 975, 977, 978, 979)

grpdf = list()
grprandf = list()
for (s in seq_along(subs)) {
  
  subdf = read.csv(subs[s]) 
  
  subfixeddf = data.frame(y = subdf$probe_slider_fixed.response[subdf$trialCounter %in% fixed_trial_nums],# * ((0.75-0.17)/(0.7-0.17)),
                     x = fixed_trial_nums,
                     block = factor(c(rep(1,25),rep(2:4,each=20),rep(5,25))),
                     adapt = factor(c(rep('pre-adapt',25),rep('adapt',60),rep('post-adapt',25))))
  
  subrandf = data.frame(y = subdf$probe_slider_random.response[subdf$trialCounter %in% random_trial_nums] -
                          subdf$memTargetPos[subdf$trialCounter %in% random_trial_nums],
                        x = random_trial_nums,
                        block = factor(c(rep(1,75),rep(2:4,each=20),rep(5,75))),
                        adapt = factor(c(rep('pre-adapt',75),rep('adapt',60),rep('post-adapt',75))),
                        memTargetPos = subdf$memTargetPos[subdf$trialCounter %in% random_trial_nums])
  
  # data quality check
  print(subs[s])
  print(sum(!is.na(subdf$key_resp.rt))/660 > 0.75) # Att-adapt response rate
  print(sum(!is.na(subfixeddf$y))/100 > 0.75) # WM-fixed response rate
  print(sum(!is.na(subrandf$y))/100 > 0.75) # WM-random response rate
  print(mean(df$key_resp.corr,na.rm=T) > (2/3)) # Att-adapt accuracy (including non-response as incorrect)
  print(mean(abs(subrandf$y),na.rm=T) < 0.15) # WM-random MAE
  print(all(df$key_resp_block_break.rt[!is.na(df$key_resp_block_break.rt)]/60 < 10)) # inter-block break time
  print(df$quiz_slider_1.response[1104] <= 3) # head movement question
  print(df$quiz_slider_2.response[1104] >= 8) # eye movement question
  print('')
  
  # add to group list
  grpdf[[s]] = subfixeddf
  grpdf[[s]]$subject = rep(s,nrow(grpdf[[s]]))
  grprandf[[s]] = subrandf
  grprandf[[s]]$subject = rep(s,nrow(grprandf[[s]]))
  
  # compute adaptation and washout magnitude
  tmp = grpdf[[s]] %>% group_by(block) %>% dplyr::summarise(mn = mean(y,na.rm=T)) %>% select(mn) %>% unlist()
  block1 = tmp[1]
  block4 = tmp[4]
  block5 = tmp[5]
  adaptpct[s] = (block4 - block1)/0.17
  washoutpct[s] = (block5 - block4)/0.17
}

# concatenate subjects into single data frame
grpdf = do.call(rbind, grpdf) 
grprandf = do.call(rbind,grprandf)

# adaptation magnitude (bootstrap SE)
bootmu = function(sample,i) mean(sample[i])
set.seed(1111)
boot(adaptpct,bootmu,100)

# Bayes factor analysis comparing 1st and last block
set.seed(1111)
df_pre_vs_adapt = filter(grpdf,block %in% c(1,4)) %>% droplevels()
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
