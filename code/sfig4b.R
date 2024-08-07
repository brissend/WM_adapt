# load libraries
library(tidyverse)
library(viridis)

# set environment
LOCAL = Sys.getenv("LOCAL")
setwd(LOCAL)

# load gaze data
gaze_fixed = read_csv('WM_adapt/data/exp4/exp4_gaze_fixed.csv')

# gaze density plot
gaze_fixed_agg = gaze_fixed %>% 
  filter(time >=600,time <=817,trial >= 529, trial <= 776) %>% 
  group_by(ID,trial) %>%
  dplyr::summarise(x_interp_driftcorr_dva = mean(x_interp_driftcorr_dva,na.rm = T),
                   y_interp_driftcorr_dva = mean(y_interp_driftcorr_dva,na.rm = T))

# create 1 dva radius ring around fixation
angles = seq(0,359,1)
polardf = data.frame(x = 1 * cos((angles) * (pi/180)),
                     y = 1 * sin((angles) * (pi/180)),
                     angle = angles)

ggplot(NULL) + 
  coord_fixed() + 
  xlim(-1.1,10) + 
  ylim(-3,3) + 
  stat_density_2d(data = gaze_fixed_agg,
                  aes(x = x_interp_driftcorr_dva, y = y_interp_driftcorr_dva,
                      fill = after_stat(ndensity)),
                  geom = 'raster',contour = FALSE,n = 500) +
  scale_fill_viridis_c() +
  geom_point(data = polardf,aes(x = x,y = y),color = 'white',size = 0.05,alpha = 0.1,shape = 16) + 
  geom_vline(xintercept = 0,linetype = 'dashed',color = 'white', size = 0.5,alpha = 0.5) + 
  geom_hline(yintercept = 0,linetype = 'dashed',color = 'white',size = 0.5,alpha = 0.5) +
  labs(x = 'Visual Angle (°)',y = 'Visual Angle (°)') +
  facet_wrap(~ID,nrow = 6,ncol = 2) + 
  theme_classic() + 
  theme(legend.position = 'none')

