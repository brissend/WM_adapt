# load libraries
library(tidyverse)
library(boot)

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
