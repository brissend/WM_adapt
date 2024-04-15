# load libraries
library(gazerjb) # forked and edited version of gazer package
library(tidyverse)

# set environment
LOCAL = Sys.getenv("LOCAL")
setwd(LOCAL)

# get subject list
datadir = 'WM_adapt/data/exp4/'
subs = list.files(path = datadir, include.dirs = TRUE, pattern = '[0-9]{6}[A-Za-z]{2}$')

# load saccade data
saccade_adapt = read_csv(file.path(datadir,'exp4_saccade_att.csv'))
saccade_fixed = read_csv(file.path(datadir,'exp4_saccade_fixed.csv'))
saccade_random = read_csv(file.path(datadir,'exp4_saccade_random.csv'))

# % of trials without a saccade during and immediately after stimulus presentation
saccade_adapt_filter = filter(saccade_adapt,
                              onset >= 600, # stimulus onset
                              onset <= 917, # exogenous cue (116.67 ms) + target (100 ms) + 100 ms)
                              diffr >= 1, # amplitude > 1 dva
                             (difftheta < 90) | (difftheta >  270)) # right hemifield

