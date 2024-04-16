# load libraries
library(gazerjb) # forked and edited version of gazer package
library(zoo)
library(tidyverse)

# set environment
LOCAL = Sys.getenv("LOCAL")
setwd(LOCAL)

# get subject list
datadir = 'WM_adapt/data/exp4/'
subs = list.files(path = datadir, include.dirs = TRUE, pattern = '[0-9]{6}[A-Za-z]{2}$')

# parse eye tracking ascii files
for (s in seq_along(subs)) {
  
  # extract samples
  if (!file.exists(file.path(datadir,subs[s],paste0(subs[s],'_saccade.csv')))) {
    # original gazer function only extracted gaze data from ascii file
    # this version extracts gaze, saccade, and fixation data and saves to separate files
    parse_asc(subs[s], homeDir = datadir, overwriteBlinks = F, cutPreview = 0) 
  }
  
  # extract condition
  if (!file.exists(file.path(datadir,subs[s],paste0(subs[s],'_messages.csv')))) {
    find_messages_asc_fix(subs[s], 
                          homeDir = datadir,
                          vars2extract = 'Condition')
  }
  
  # merging eye data and behavioral condition files
  # gaze
  if (!file.exists(file.path(datadir,subs[s],paste0(subs[s],'_combined.csv')))) {
    merge_asc_files(subs[s],datadir)
  }
  
  # saccade
  if (!file.exists(file.path(datadir,subs[s],paste0(subs[s],'_saccade_combined.csv')))) {
    sdf = read_csv(file.path(datadir,subs[s],paste0(subs[s],'_saccade.csv')))
    mdf = read_csv(file.path(datadir,subs[s],paste0(subs[s],'_messages.csv')))
    mdf$ID = rep(sdf$ID[1],length.out = nrow(mdf))
    smdf = full_join(sdf,mdf,by = c('trial','ID')) %>% arrange(.,trial)
    write_csv(smdf,
              file = file.path(datadir,subs[s],paste0(subs[s],'_saccade_combined.csv')))
  }
  
  # fixation
  if (!file.exists(file.path(datadir,subs[s],paste0(subs[s],'_fixation_combined.csv')))) {
    fdf = read_csv(file.path(datadir,subs[s],paste0(subs[s],'_fixation.csv')))
    mdf = read_csv(file.path(datadir,subs[s],paste0(subs[s],'_messages.csv')))
    mdf$ID = rep(fdf$ID[1],length.out = nrow(mdf))
    fmdf = full_join(fdf,mdf,by = c('trial','ID')) %>% arrange(.,trial)
    write_csv(fmdf,file = file.path(datadir,subs[s],paste0(subs[s],'_fixation_combined.csv')))
  }
}
  