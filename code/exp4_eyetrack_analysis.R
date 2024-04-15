# load libraries
library(gazerjb) # forked and edited version of gazer package (https://github.com/brissend/gazerjb)
library(zoo)
library(tidyverse)

# set environment
LOCAL = Sys.getenv("LOCAL")
setwd(LOCAL)

# get subject list
datadir = 'WM_adapt/data/exp4/'
subs = list.files(path = datadir, include.dirs = TRUE, pattern = '[0-9]{6}[A-Za-z]{2}$')

# screen parameters
res = c(1920,1200)
screen_dist = 85 # cm
screen_width = 51.84 # cm
ppd = res[1] / (2 * atan(screen_width/(2 * screen_dist)) * (180/pi)) # pixels per degree of visual angle
center = res/2

# load eye tracking data for each subject
gaze = list()
saccade = list()
fixation = list()
for (s in seq_along(subs)) {
  gaze[[s]] = fread(file.path(datadir,subs[s],paste0(subs[s],'_combined.csv'))) %>% as_tibble()
  saccade[[s]] = read_csv(file.path(datadir,subs[s],paste0(subs[s],'_saccade_combined.csv')),
                          na = c(".","","NA"))
  fixation[[s]] = read_csv(file.path(datadir,subs[s],paste0(subs[s],'_fixation_combined.csv')),
                           na = c(".","","NA"))
  
}
gaze = do.call(rbind,gaze)
saccade = do.call(rbind,saccade)
fixation = do.call(rbind,fixation)

# identify and extend blinks based on the velocity of pupil size changes
gaze$pupil_na = gaze$pupil
gaze$pupil_na[gaze$pupil_na == 0] = NA
gaze$extendpupil = extend_blinks(gaze$pupil_na,fillback = 100,fillforward = 100, hz = 1000) # van Ede, Chekroud, Nobre, 2019 Nature Human Behavior

# set x and y coords corresponding with blinks to NA
gaze$x[is.na(gaze$extendpupil)] = NA
gaze$y[is.na(gaze$extendpupil)] = NA

# also set x and y coords recorded as missing by parse_asc back to NA
gaze$x[gaze$x == 1e08] = NA
gaze$y[gaze$y == 1e08] = NA

# interpolate across blinks and missing data (linear much better than cubic)
# gaze = gaze %>%
#   dplyr::mutate(x_interp = zoo::na.spline(x, na.rm = FALSE),
#                 y_interp = zoo::na.spline(y, na.rm = FALSE))
gaze = gaze %>%
  dplyr::mutate(x_interp = zoo::na.approx(x, na.rm = FALSE, rule = 2),
                y_interp = zoo::na.approx(y, na.rm = FALSE, rule = 2))

# separate gaze for each condition
gaze_adapt = gaze %>%
  filter(Condition == 'adapt')
gaze_fixed = gaze %>%
  filter(Condition == 'memory')
gaze_random = gaze %>% 
  filter(Condition == 'memory_random')

# drift correction each trial 
# (perform separately for each condition due to different timing)
# compute median gaze position for first 200 ms of trial (100 ms trial start cue + 100 ms fixation)
gaze_adapt = gaze_adapt %>%
  filter(time > 0,time < 200) %>%
  group_by(ID, trial) %>%
  dplyr::summarise(xbaseline = median(x_interp,na.rm = T),
                   ybaseline = median(y_interp, na.rm = T)) %>%
  dplyr::full_join(., gaze_adapt) %>% 
  ungroup()

# recenter gaze coordinates for rest of trial based on median gaze position at start of trial
gaze_adapt = gaze_adapt %>%
  dplyr::mutate(x_interp_driftcorr = x_interp - xbaseline + center[1],
                y_interp_driftcorr = y_interp - ybaseline + center[2]) %>%
  dplyr::arrange(subject, trial, time)

# wm-fixed; gaze data for wm trials included 1000 ms ITI prior to the trial start cue
gaze_fixed = gaze_fixed %>%
  filter(time > 1000,time < 1200) %>%
  group_by(ID, trial) %>%
  dplyr::summarise(xbaseline = median(x_interp,na.rm = T),
                   ybaseline = median(y_interp, na.rm = T)) %>%
  dplyr::full_join(., gaze_fixed) %>% 
  ungroup()

# recenter
gaze_fixed = gaze_fixed %>%
  dplyr::mutate(x_interp_driftcorr = x_interp - xbaseline + center[1],
                y_interp_driftcorr = y_interp - ybaseline + center[2]) %>%
  dplyr::arrange(subject, trial, time)

# wm-random
gaze_random = gaze_random %>%
  filter(time > 1000,time < 1200) %>%
  group_by(ID, trial) %>%
  dplyr::summarise(xbaseline = median(x_interp,na.rm = T),
                   ybaseline = median(y_interp, na.rm = T)) %>%
  dplyr::full_join(., gaze_random) %>% 
  ungroup()

# recenter
gaze_random = gaze_random %>%
  dplyr::mutate(x_interp_driftcorr = x_interp - xbaseline + center[1],
                y_interp_driftcorr = y_interp - ybaseline + center[2]) %>%
  dplyr::arrange(subject, trial, time)

# convert x and y coords from pixel units to degrees of visual angle
gaze_adapt$xdva_interp = (gaze_adapt$x_interp - center[1]) / ppd
gaze_adapt$ydva_interp = ((gaze_adapt$y_interp - center[2]) / ppd) * -1
gaze_adapt$x_interp_driftcorr_dva = (gaze_adapt$x_interp_driftcorr - center[1]) / ppd
gaze_adapt$y_interp_driftcorr_dva = ((gaze_adapt$y_interp_driftcorr - center[2]) / ppd) * -1

gaze_fixed$xdva_interp = (gaze_fixed$x_interp - center[1]) / ppd
gaze_fixed$ydva_interp = ((gaze_fixed$y_interp - center[2]) / ppd) * -1
gaze_fixed$x_interp_driftcorr_dva = (gaze_fixed$x_interp_driftcorr - center[1]) / ppd
gaze_fixed$y_interp_driftcorr_dva = ((gaze_fixed$y_interp_driftcorr - center[2]) / ppd) * -1

gaze_random$xdva_interp = (gaze_random$x_interp - center[1]) / ppd
gaze_random$ydva_interp = ((gaze_random$y_interp - center[2]) / ppd) * -1
gaze_random$x_interp_driftcorr_dva = (gaze_random$x_interp_driftcorr - center[1]) / ppd
gaze_random$y_interp_driftcorr_dva = ((gaze_random$y_interp_driftcorr - center[2]) / ppd) * -1