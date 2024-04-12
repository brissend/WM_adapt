# load libraries
library(tidybayes)
library(tidyverse)
library(viridis)

# set environment
LOCAL = Sys.getenv("LOCAL")
setwd(LOCAL)

# load data
if (!file.exists('WM_adapt/data/exp1/exp1_group_mean_fixed.csv')) {
  source('WM_adapt/code/exp1_analysis.R') # group average recall computed and saved by this script
}
mndf = read.csv('WM_adapt/data/exp1/exp1_group_mean_fixed.csv')

# load model fit
fit = readRDS('WM_adapt/model_fits/Exp1_expdecay_fit_group_mean.rds')

# posterior expected distribution (posterior predictive distribution sans residual)
fitdf = data.frame()
for (x in seq(min(mndf$x),max(mndf$x),by = 20)) {
  
  tmp = epred_draws(fit,newdata = data.frame(x = x), ndraws = 5000)
  
  fitdf = rbind(fitdf,
                data.frame(x = x, 
                           y = mean(tmp$.epred),
                           y_low = bayestestR::hdi(tmp$.epred,ci = 0.95)$CI_low,
                           y_hi = bayestestR::hdi(tmp$.epred,ci = 0.95)$CI_high))
}

# plot
mndf$block = factor(mndf$block)
ggplot(mndf,aes(x = x, y = y )) + 
  geom_point(aes(color = block)) + 
  geom_ribbon(data = fitdf, alpha = 0.25,color = NA,aes(ymin = y_low, ymax = y_hi)) + 
  geom_line(data = fitdf,size = 1,aes(x = x, y = y)) + 
  geom_hline(yintercept = 0.5,linetype = 'dashed') + 
  scale_color_viridis(discrete = T ,begin = 0.05,end = 0.95) + 
  scale_x_continuous(breaks = c(0,250,500,750,1000,1250)) + 
  scale_y_continuous(breaks = seq(0.517,0.449,by = -0.017),labels = paste0(seq(10,-30,by=-10),'%')) + 
  coord_cartesian(ylim = c(0.447,0.517)) + 
  labs(color = 'Block',y = 'Recall Location (% of Backstep)', x = 'Trial') +
  theme_classic() +
  theme(legend.position = 'none')

