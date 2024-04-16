# load libraries
library(gazerjb) # forked and edited version of gazer package (https://github.com/brissend/gazerjb)
library(zoo)
library(tidyverse)

# set environment
LOCAL = Sys.getenv("LOCAL")
setwd(LOCAL)

# get subject list
datadir = 'WM_adapt/data/exp5/'
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

# ===================
# gaze preprocessing
# ===================

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

# remove ITI from beginning of memory trials (gaze data for wm trials included 1000 ms ITI prior to the trial start cue)
gaze_fixed = filter(gaze_fixed, time >= 1000) %>% mutate(time = time - 1000)
gaze_random = filter(gaze_random, time >= 1000) %>% mutate(time = time - 1000)

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

# wm-fixed median gaze
gaze_fixed = gaze_fixed %>%
  filter(time > 0,time < 200) %>%
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

# wm-random median gaze
gaze_random = gaze_random %>%
  filter(time > 0,time < 200) %>%
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

# convert cartesian coordinates to polar coordinates
gaze_adapt = gaze_adapt %>%
  mutate(r = useful::cart2pol(x_interp_driftcorr_dva,y_interp_driftcorr_dva,degrees = T)$r,
         theta = useful::cart2pol(x_interp_driftcorr_dva,y_interp_driftcorr_dva,degrees = T)$theta %% 360)

gaze_fixed = gaze_fixed %>%
  mutate(r = useful::cart2pol(x_interp_driftcorr_dva,y_interp_driftcorr_dva,degrees = T)$r,
         theta = useful::cart2pol(x_interp_driftcorr_dva,y_interp_driftcorr_dva,degrees = T)$theta %% 360)

gaze_random = gaze_random %>%
  mutate(r = useful::cart2pol(x_interp_driftcorr_dva,y_interp_driftcorr_dva,degrees = T)$r,
         theta = useful::cart2pol(x_interp_driftcorr_dva,y_interp_driftcorr_dva,degrees = T)$theta %% 360)

# save preprocessed gaze data
write_csv(gaze_adapt,file = file.path(datadir,'exp5_gaze_att.csv'))
write_csv(gaze_fixed,file = file.path(datadir,'exp5_gaze_fixed.csv'))
write_csv(gaze_random,file = file.path(datadir,'exp5_gaze_random.csv'))

# ======================
# saccade preprocessing
# ======================
saccade_adapt = saccade %>%
  filter(Condition == 'adapt')
saccade_fixed = saccade %>%
  filter(Condition == 'memory')
saccade_random = saccade %>% 
  filter(Condition == 'memory_random')

# remove ITI from beginning of memory trials
saccade_fixed = filter(saccade_fixed, onset >= 1000) %>% mutate(onset = onset - 1000,
                                                                offset = offset - 1000)
saccade_random = filter(saccade_random, onset >= 1000) %>% mutate(onset = onset - 1000,
                                                                  offset = offset - 1000)

# drift correction
# att-error
saccade_adapt_baseline = gaze_adapt %>%
  filter(time > 0,time < 200) %>%
  group_by(ID, trial) %>%
  dplyr::summarise(xbaseline = median(x_interp,na.rm = T),
                   ybaseline = median(y_interp, na.rm = T))

saccade_adapt = left_join(saccade_adapt, saccade_adapt_baseline,by = c('ID','trial')) %>%
  dplyr::mutate(startx_driftcorr = startx - xbaseline + center[1],
                starty_driftcorr = starty - ybaseline + center[2],
                endx_driftcorr = endx - xbaseline + center[1],
                endy_driftcorr = endy - ybaseline + center[2]) %>%
  dplyr::arrange(ID, trial, onset)

# wm-fixed
saccade_fixed_baseline = gaze_fixed %>%
  filter(time > 0,time < 200) %>%
  group_by(ID, trial) %>%
  dplyr::summarise(xbaseline = median(x_interp,na.rm = T),
                   ybaseline = median(y_interp, na.rm = T))

saccade_fixed = left_join(saccade_fixed, saccade_fixed_baseline,by = c('ID','trial')) %>%
  dplyr::mutate(startx_driftcorr = startx - xbaseline + center[1],
                starty_driftcorr = starty - ybaseline + center[2],
                endx_driftcorr = endx - xbaseline + center[1],
                endy_driftcorr = endy - ybaseline + center[2]) %>%
  dplyr::arrange(ID, trial, onset)

# wm-random
saccade_random_baseline = gaze_random %>%
  filter(time > 0,time < 200) %>%
  group_by(ID, trial) %>%
  dplyr::summarise(xbaseline = median(x_interp,na.rm = T),
                   ybaseline = median(y_interp, na.rm = T))

saccade_random = left_join(saccade_random, saccade_random_baseline,by = c('ID','trial')) %>%
  dplyr::mutate(startx_driftcorr = startx - xbaseline + center[1],
                starty_driftcorr = starty - ybaseline + center[2],
                endx_driftcorr = endx - xbaseline + center[1],
                endy_driftcorr = endy - ybaseline + center[2]) %>%
  dplyr::arrange(ID, trial, onset)


# convert x and y coords to degrees of visual angle
# att-error
saccade_adapt$startx_driftcorr_dva = (saccade_adapt$startx_driftcorr - center[1]) / ppd
saccade_adapt$starty_driftcorr_dva = ((saccade_adapt$starty_driftcorr - center[2]) / ppd) * -1
saccade_adapt$endx_driftcorr_dva = (saccade_adapt$endx_driftcorr - center[1]) / ppd
saccade_adapt$endy_driftcorr_dva = ((saccade_adapt$endy_driftcorr - center[2]) / ppd) * -1
saccade_adapt$diffx_driftcorr_dva = saccade_adapt$endx_driftcorr_dva - saccade_adapt$startx_driftcorr_dva
saccade_adapt$diffy_driftcorr_dva = saccade_adapt$endy_driftcorr_dva - saccade_adapt$starty_driftcorr_dva

# wm-fixed
saccade_fixed$startx_driftcorr_dva = (saccade_fixed$startx_driftcorr - center[1]) / ppd
saccade_fixed$starty_driftcorr_dva = ((saccade_fixed$starty_driftcorr - center[2]) / ppd) * -1
saccade_fixed$endx_driftcorr_dva = (saccade_fixed$endx_driftcorr - center[1]) / ppd
saccade_fixed$endy_driftcorr_dva = ((saccade_fixed$endy_driftcorr - center[2]) / ppd) * -1
saccade_fixed$diffx_driftcorr_dva = saccade_fixed$endx_driftcorr_dva - saccade_fixed$startx_driftcorr_dva
saccade_fixed$diffy_driftcorr_dva = saccade_fixed$endy_driftcorr_dva - saccade_fixed$starty_driftcorr_dva

# wm-random
saccade_random$startx_driftcorr_dva = (saccade_random$startx_driftcorr - center[1]) / ppd
saccade_random$starty_driftcorr_dva = ((saccade_random$starty_driftcorr - center[2]) / ppd) * -1
saccade_random$endx_driftcorr_dva = (saccade_random$endx_driftcorr - center[1]) / ppd
saccade_random$endy_driftcorr_dva = ((saccade_random$endy_driftcorr - center[2]) / ppd) * -1
saccade_random$diffx_driftcorr_dva = saccade_random$endx_driftcorr_dva - saccade_random$startx_driftcorr_dva
saccade_random$diffy_driftcorr_dva = saccade_random$endy_driftcorr_dva - saccade_random$starty_driftcorr_dva

# convert start and end points to polar coordinates
# att-error
saccade_adapt = saccade_adapt %>%
  mutate(startr = useful::cart2pol(startx_driftcorr_dva,starty_driftcorr_dva,degrees = T)$r,
         starttheta = useful::cart2pol(startx_driftcorr_dva,starty_driftcorr_dva,degrees = T)$theta %% 360,
         endr = useful::cart2pol(endx_driftcorr_dva,endy_driftcorr_dva,degrees = T)$r,
         endtheta = useful::cart2pol(endx_driftcorr_dva,endy_driftcorr_dva,degrees = T)$theta %% 360,
         diffr = useful::cart2pol(diffx_driftcorr_dva,diffy_driftcorr_dva,degrees = T)$r,
         difftheta = useful::cart2pol(diffx_driftcorr_dva,diffy_driftcorr_dva,degrees = T)$theta %% 360)

saccade_fixed = saccade_fixed %>%
  mutate(startr = useful::cart2pol(startx_driftcorr_dva,starty_driftcorr_dva,degrees = T)$r,
         starttheta = useful::cart2pol(startx_driftcorr_dva,starty_driftcorr_dva,degrees = T)$theta %% 360,
         endr = useful::cart2pol(endx_driftcorr_dva,endy_driftcorr_dva,degrees = T)$r,
         endtheta = useful::cart2pol(endx_driftcorr_dva,endy_driftcorr_dva,degrees = T)$theta %% 360,
         diffr = useful::cart2pol(diffx_driftcorr_dva,diffy_driftcorr_dva,degrees = T)$r,
         difftheta = useful::cart2pol(diffx_driftcorr_dva,diffy_driftcorr_dva,degrees = T)$theta %% 360)

saccade_random = saccade_random %>%
  mutate(startr = useful::cart2pol(startx_driftcorr_dva,starty_driftcorr_dva,degrees = T)$r,
         starttheta = useful::cart2pol(startx_driftcorr_dva,starty_driftcorr_dva,degrees = T)$theta %% 360,
         endr = useful::cart2pol(endx_driftcorr_dva,endy_driftcorr_dva,degrees = T)$r,
         endtheta = useful::cart2pol(endx_driftcorr_dva,endy_driftcorr_dva,degrees = T)$theta %% 360,
         diffr = useful::cart2pol(diffx_driftcorr_dva,diffy_driftcorr_dva,degrees = T)$r,
         difftheta = useful::cart2pol(diffx_driftcorr_dva,diffy_driftcorr_dva,degrees = T)$theta %% 360)

# save preprocessed saccade data
write_csv(saccade_adapt,file = file.path(datadir,'exp5_saccade_att.csv'))
write_csv(saccade_fixed,file = file.path(datadir,'exp5_saccade_fixed.csv'))
write_csv(saccade_random,file = file.path(datadir,'exp5_saccade_random.csv'))
