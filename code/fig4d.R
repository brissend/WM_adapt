# load libraries
library(tidyverse)

# set environment
LOCAL = Sys.getenv("LOCAL")
setwd(LOCAL)

# load data
if (!file.exists('WM_adapt/data/exp5/exp5_group_random.csv')) {
  source('WM_adapt/code/exp5_analysis.R') # group aggregated WM-random data frame saved by this script
}
grprandf = read.csv('WM_adapt/data/exp5/exp5_group_random.csv')

# bin WM-random location
ranbindf = grprandf %>% 
  mutate(bin = cut(memTargetPos,breaks = c(-Inf,seq(-10.5,-4.5,length.out = 7),Inf))) %>%
  group_by(subject,blockNum,bin) %>%
  dplyr::summarise(error = mean(error,na.rm = T))

blockdiffdf = ranbindf %>% 
  filter(blockNum %in% c('1','4')) %>%
  ungroup %>% 
  spread(blockNum, error) %>% 
  group_by(subject,bin) %>% 
  mutate(ydiff = `4` - `1`)
levels(blockdiffdf$bin) =  c('+2','+1' ,'0','-1','-2','-3','-4','-5')
blockdiffdf$ydiffpct = blockdiffdf$ydiff/3 * 100

# plot
blockdiffdf_mean_cl_boot = blockdiffdf %>% group_by(bin) %>% summarise(mean_cl_boot(ydiffpct))
ggplot(blockdiffdf_mean_cl_boot,aes(x = bin,y = y,color = bin,fill = bin)) + 
  geom_hline(yintercept = 0,linetype = 'dashed') + 
  geom_linerange(aes(ymin = ymin,ymax = ymax)) + 
  geom_point(size = 5) + 
  labs( y = 'Recall Error (Block 4 - Block 1)',x = 'WM Random Location Bin Distance (°)') + 
  scale_y_continuous(breaks =  seq(-20,60,by = 20),labels = paste0(seq(20,-60,by = -20),'%'),limits = c(-30,70)) + 
  theme_classic() + 
  theme(legend.position = 'none')



