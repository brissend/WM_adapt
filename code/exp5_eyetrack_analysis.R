# load libraries
library(tidyverse)
library(boot)
library(BayesFactor)

# set environment
LOCAL = Sys.getenv("LOCAL")
setwd(LOCAL)

# get subject list
datadir = 'WM_adapt/data/exp5/'
subs = list.files(path = datadir, include.dirs = TRUE, pattern = '[0-9]{6}[A-Za-z]{2}$')

# load saccade data
saccade_adapt = read_csv(file.path(datadir,'exp5_saccade_att.csv'))
saccade_fixed = read_csv(file.path(datadir,'exp5_saccade_fixed.csv'))
saccade_random = read_csv(file.path(datadir,'exp5_saccade_random.csv'))

### Att - error ###
# count number of trials w/ at least 1 saccade during stimulus presentation period + 100 ms
adapt_pct_saccade_df = saccade_adapt %>% 
  filter(onset >= 600, # stimulus onset
         onset <= 917, # exogenous cue (116.67 ms) + target (100 ms) + 100 ms)
         diffr >= 1, # amplitude > 1 dva
         (difftheta > 90) & (difftheta < 270)) %>% # left hemifield
  group_by(ID,trial) %>% 
  dplyr::summarise(n = n()) %>% 
  ungroup() %>% 
  group_by(ID) %>% 
  dplyr::summarise(ntrial = length(unique(trial))) %>%
  mutate(pct_nosaccade = ((660 - ntrial) / 660) * 100) %>%
  mutate(pct_saccade = 100 - pct_nosaccade)

mean(adapt_pct_saccade_df$pct_nosaccade)
bootmu = function(sample,i) mean(sample[i])
set.seed(1111)
boot(adapt_pct_saccade_df$pct_nosaccade,bootmu,100)

### WM - fixed ###
# all adaptation blocks
wm_fixed_pct_saccade_df = saccade_fixed %>%
  filter(onset >= 600, # stimulus onset
         onset <= 817, # memory sample (116.67 ms) + 100 ms
         trial > 100, trial <= 776, # adaptation blocks
         diffr >= 1, # amplitude > 1 dva
         (difftheta > 90) & (difftheta < 270)) %>% # left hemifield
  group_by(ID,trial) %>% 
  dplyr::summarise(n = n()) %>% 
  ungroup() %>% 
  group_by(ID) %>% 
  dplyr::summarise(ntrial = length(unique(trial))) %>%
  mutate(pct_nosaccade = ((60 - ntrial) / 60) * 100) %>% # 60 total WM-fixed trials across all adaptation blocks
  mutate(pct_saccade = 100 - pct_nosaccade)

pct_vec = c(wm_fixed_pct_saccade_df$pct_nosaccade,rep(100,11 - length(wm_fixed_pct_saccade_df$pct_nosaccade))) # add subjects who made no saccades during critical period (0 trials w/ saccade)
mean(pct_vec)
set.seed(1111)
boot(pct_vec,bootmu,100)

# last adaptation block
wm_fixed_pct_saccade_df = saccade_fixed %>% 
  filter(onset >= 600,
         onset <= 817,
         trial >= 529, trial <= 776, # last adaptation block
         diffr >= 1,
         (difftheta > 90) & (difftheta <  270)) %>%
  group_by(ID,trial) %>% 
  dplyr::summarise(n = n()) %>% 
  ungroup() %>% 
  group_by(ID) %>% 
  dplyr::summarise(ntrial = length(unique(trial))) %>%
  mutate(pct_nosaccade = ((20 - ntrial) / 20) * 100) %>% # 20 WM-fixed trials in last adaptation block
  mutate(pct_saccade = 100 - pct_nosaccade)

pct_vec = c(wm_fixed_pct_saccade_df$pct_nosaccade,rep(100,11 - length(wm_fixed_pct_saccade_df$pct_nosaccade)))
mean(pct_vec)
set.seed(1111)
boot(pct_vec,bootmu,100)

### Adaptation magnitude excluding trials with saccades ###
# load behavioral data
# load data
if (!file.exists('WM_adapt/data/exp5/exp5_group_mean_fixed.csv')) {
  source('WM_adapt/code/exp5_analysis.R') # group average recall computed and saved by this script
}
mndf = read.csv('WM_adapt/data/exp5/exp5_group_mean_fixed.csv')

if (!file.exists('WM_adapt/data/exp5/exp5_group_fixed.csv')) {
  source('WM_adapt/code/exp5_analysis.R') # group recall computed and saved by this script
}
grpdf = read.csv('WM_adapt/data/exp5/exp5_group_fixed.csv')

# determine if any saccade was made during stimulus presentation and delay period
saccade_fixed_stim = saccade_fixed %>%
  filter(onset >= 600, # stimulus onset
         onset <= 817, # memory sample (116.67 ms) + 100 ms
         diffr >= 1, # amplitude > 1 dva
         (difftheta > 90) & (difftheta <  270)) # left hemifield

saccade_fixed_delay = saccade_fixed %>%
  filter(onset >= 600, # stimulus onset
         onset <= 1217, # memory sample (116.67 ms) + delay period (500 ms)
         diffr >= 1, # amplitude > 1 dva
         (difftheta > 90) & (difftheta <  270)) # left hemifield

grpdf$saccade_stim_bool = rep(F, nrow(grpdf))
grpdf$saccade_delay_bool = rep(F,nrow(grpdf))
for (s in unique(grpdf$subject)) {
  for (t in unique(grpdf$trial)) {
    
    # save whether saccade was made during stimulus presentation (+100 ms stimulus offset)
    grpdf$saccade_stim_bool[which((grpdf$subject == s) & (grpdf$trial == t))] = 
      nrow(filter(saccade_fixed_stim,ID == s,trial == t)) > 0
    
    # save whether saccade was made during stimulus presentation + delay period
    grpdf$saccade_delay_bool[which((grpdf$subject == s) & (grpdf$trial == t))] = 
      nrow(filter(saccade_fixed_delay,ID == s,trial == t)) > 0
    
  }
}

# compute average timecourse excluding trials with a saccade during either stimulus presentation or delay period
# stimulus presentation + 100 ms fixation
mndf_nosaccade_stim = grpdf %>% 
  filter(!saccade_stim_bool) %>%
  group_by(trial) %>% 
  dplyr::summarize(memResponse = mean(memResponse,na.rm=T),n = n()) %>% 
  filter(n > 1)
mndf_nosaccade_stim = left_join(mndf_nosaccade_stim,as_tibble(mndf), # get block and phase info from original mean timecourse
                                by = 'trial',
                                suffix=c("",".y")) %>%
  select(-ends_with(".y")) 

# save 
if (!file.exists('WM_adapt/data/exp5/exp5_nosaccade_stim_group_mean_fixed.csv')) {
  write.csv(mndf_nosaccade_stim,
            file = 'WM_adapt/data/exp5/exp5_nosaccade_stim_group_mean_fixed.csv',
            row.names = F)
}

# stimulus presentation + delay period
mndf_nosaccade_delay = grpdf %>% 
  filter(!saccade_delay_bool) %>%
  group_by(trial) %>% 
  dplyr::summarize(memResponse = mean(memResponse,na.rm=T),n = n()) %>% 
  filter(n > 1)
mndf_nosaccade_delay = left_join(mndf_nosaccade_delay,as_tibble(mndf), # get block and phase info from original mean timecourse
                                 by = 'trial',
                                 suffix=c("",".y")) %>% 
  select(-ends_with(".y")) 

# save 
if (!file.exists('WM_adapt/data/exp5/exp5_nosaccade_delay_group_mean_fixed.csv')) {
  write.csv(mndf_nosaccade_delay,
            file = 'WM_adapt/data/exp5/exp5_nosaccade_delay_group_mean_fixed.csv',
            row.names = F)
}

subs = unique(grpdf$subject)
adaptpct_nosaccade_stim = rep(NA,length(subs))
washoutpct_nosaccade_stim = rep(NA,length(subs))
adaptpct_nosaccade_delay = rep(NA,length(subs))
washoutpct_nosaccade_delay = rep(NA,length(subs))
for (s in seq_along(subs)) {
  
  # compute adaptation and washout magnitude
  tmp = grpdf %>% filter(subject == subs[s],!saccade_stim_bool) %>% group_by(blockNum) %>% dplyr::summarise(mn = mean(memResponse,na.rm=T)) %>% select(mn) %>% unlist()
  block1 = tmp[1]
  block4 = tmp[4]
  block5 = tmp[5]
  adaptpct_nosaccade_stim[s] = (block4 - block1)/3 # back step size = 3° of visual angle
  washoutpct_nosaccade_stim[s] = (block5 - block4)/3
  
  tmp = grpdf %>% filter(subject == subs[s],!saccade_delay_bool) %>% group_by(blockNum) %>% dplyr::summarise(mn = mean(memResponse,na.rm=T)) %>% select(mn) %>% unlist()
  block1 = tmp[1]
  block4 = tmp[4]
  block5 = tmp[5]
  adaptpct_nosaccade_delay[s] = (block4 - block1)/3 # back step size = 3° of visual angle
  washoutpct_nosaccade_delay[s] = (block5 - block4)/3
}

# adaptation magnitude (bootstrap SE)
bootmu = function(sample,i) mean(sample[i])
set.seed(1111)
boot(adaptpct_nosaccade_stim,bootmu,100)

set.seed(1111)
boot(adaptpct_nosaccade_delay,bootmu,100)

# washout magnitude
set.seed(1111)
boot(washoutpct_nosaccade_stim,bootmu,100)

set.seed(1111)
boot(washoutpct_nosaccade_delay,bootmu,100)

# correlation Bayes factor between adaptation and washout magnitude
set.seed(1111)
cor.test(adaptpct_nosaccade_stim,washoutpct_nosaccade_stim)
correlationBF(adaptpct_nosaccade_stim,washoutpct_nosaccade_stim)

set.seed(1111)
cor.test(adaptpct_nosaccade_delay,washoutpct_nosaccade_delay)
correlationBF(adaptpct_nosaccade_delay,washoutpct_nosaccade_delay)

# Bayes Factor analysis comparing pre-adapt block and last adaptation block
set.seed(1111)
df_pre_vs_adapt = filter(grpdf,blockNum %in% c(1,4),!saccade_stim_bool) %>% droplevels()
df_pre_vs_adapt$block = factor(df_pre_vs_adapt$blockNum)
df_pre_vs_adapt$subject = factor(df_pre_vs_adapt$subject)

bf_full = lmBF(memResponse ~ blockNum + subject, 
               whichRandom = "subject",
               rscaleRandom = 'nuisance',
               data = df_pre_vs_adapt %>% drop_na(memResponse))
bf_int = lmBF(memResponse ~ subject, whichRandom = "subject",
              rscaleRandom = 'nuisance',
              data = df_pre_vs_adapt %>% drop_na(memResponse))

bf_full/bf_int

set.seed(1111)
df_pre_vs_adapt = filter(grpdf,blockNum %in% c(1,4),!saccade_delay_bool) %>% droplevels()
df_pre_vs_adapt$block = factor(df_pre_vs_adapt$blockNum)
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
df_adapt_vs_post = filter(grpdf,blockNum %in% c(4,5),!saccade_stim_bool) %>% droplevels()
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

set.seed(1111)
df_adapt_vs_post = filter(grpdf,blockNum %in% c(4,5),!saccade_delay_bool) %>% droplevels()
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
set.seed(1111)
bf_full = lmBF(memResponse ~ trial + subject,
               whichRandom = "subject",
               data = filter(grpdf,phase == 'pre-adapt',!saccade_stim_bool) %>% drop_na(memResponse))
bf_int = lmBF(memResponse ~ subject,
              whichRandom = "subject",
              data = filter(grpdf,phase == 'pre-adapt',!saccade_stim_bool) %>% drop_na(memResponse))

bf_full/bf_int

set.seed(1111)
bf_full = lmBF(memResponse ~ trial + subject,
               whichRandom = "subject",
               data = filter(grpdf,phase == 'pre-adapt',!saccade_delay_bool) %>% drop_na(memResponse))
bf_int = lmBF(memResponse ~ subject,
              whichRandom = "subject",
              data = filter(grpdf,phase == 'pre-adapt',!saccade_delay_bool) %>% drop_na(memResponse))

bf_full/bf_int


