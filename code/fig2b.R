# load libraries
library(tidyverse)

# set environment
LOCAL = Sys.getenv("LOCAL")
setwd(LOCAL)

# load data
if (!file.exists('WM_adapt/data/exp1/exp1_group_random.csv')) {
  source('WM_adapt/code/exp1_analysis.R') # group aggregated WM-random data frame saved by this script
}
grprandf = read.csv('WM_adapt/data/exp1/exp1_group_random.csv')

# Bin WM-random location
ranbindf = grprandf %>% 
  mutate(bin = cut(memTargetPos,breaks = c(-Inf,seq(0.15,0.65,0.1),Inf))) %>%
  group_by(subject,block,bin) %>% 
  dplyr::summarise(y = mean(y,na.rm = T))

blockdiffdf = ranbindf %>% filter(block %in% c('1','5')) %>%
  ungroup %>% spread(block, y) %>% 
  group_by(subject,bin) %>% 
  mutate(ydiff = `5` - `1`)
levels(blockdiffdf$bin) = c('-0.4','-0.3','-0.2','-0.1','0','+0.1','+0.2')

ggplot(blockdiffdf,aes(x = bin,y = ydiff,color = bin)) + 
  geom_hline(yintercept = 0,linetype = 'dashed') + 
  stat_summary(fun=mean, geom="point", shape=23, size=2, color="red", fill="red") + 
  stat_summary(fun.data = "mean_cl_boot", size = 1) + 
  labs( y = 'Recall Error (Block 5 - Block 1)',x = 'WM Random Location Bin Distance') + 
  scale_y_continuous(breaks = seq(0.017,-0.051,by = -0.017),labels = paste0(seq(10,-30,by=-10),'%')) +
  coord_cartesian(ylim = c(-0.053,0.017)) + 
  theme_classic() + 
  theme(legend.position = 'none')

