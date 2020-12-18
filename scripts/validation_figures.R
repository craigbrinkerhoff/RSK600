#Fall 2020
#Craig Brinkerhoff
#Produce figures for pepsi validation for BIGER algorithm

library(tidyr)
library(tidyverse)
library(cowplot)
library(RColorBrewer)
library(colorspace)
library(ggtext)
theme_set(theme_cowplot())

# Sum of squares function------------------------------
sumsq <- function(x) sum(x^2)

#read in results----------------------------------
full_output <- read.csv('inputs//validation//results_SWOT_11day_slopeK.csv')
full_output <- filter(full_output, river != 'guez_etal_2020/ArialKhan')
r2 <- round(summary(lm(log(full_output$kest_mean)~log(full_output$kobs)))$r.squared, 2)
full_output <- drop_na(full_output)
rmse <- round(exp(Metrics::rmse(log(full_output$kest_mean), log(full_output$kobs))), 2)

#Model Validation across all rivers and timesteps--------------------------------------------------------
valPlot <- ggplot(full_output, aes(x=(kobs), y=(kest_mean))) +
  geom_abline(size=2, linetype='dashed', color='grey') +
 # geom_point(size=5, pch=21, color='black', fill='#1b9e77') +
  geom_pointrange(aes(ymin = kest_low, ymax = kest_high), fatten=10, fill='#1b9e77', pch=21, color='black') +
  geom_smooth(size=2, color='black', method='lm', se=F)+
  xlab('Scaling Model k [m/dy]') +
  ylab('Remotely Sensed k [m/dy]') +
    scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x)))+
   scale_x_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  coord_cartesian(xlim=c(10^-0.75, 10^0.75), ylim=c(10^-0.75,10^0.75))+
  scale_color_discrete_qualitative(palette = 'Harmonic') +
  theme(legend.position = "none",
        axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold')) +
  geom_richtext(aes(x=10^-0.5, y=10^0.6), color='black', label=paste0('RMSE: ', rmse, ' m/dy')) +
 geom_richtext(aes(x=10^-0.5, y=10^0.5), color='black', label=paste0('r<sup>2</sup>: ', r2), size=5)
ggsave('outputs//validation//validation.jpg', valPlot, width=8, height=8)

#by river metrics---------------------------------------------
stats_by_reach <- group_by(full_output, river) %>%
  summarise(kge = hydroGOF::KGE(kest_mean, kobs),
            nse =   hydroGOF::NSE(kest_mean, kobs),
            nrmse = sqrt(mean((kobs - kest_mean)^2)) / mean(kobs, na.rm=T),
            rBIAS =   mean(kest_mean - kobs) / mean(kobs, na.rm=T),
            rrmse =   sqrt(mean((kest_mean - kobs)^2 / kobs^2)),
            CVobs = round(sd(kobs)/mean(kobs, na.rm=T), 2))

plot_stats <- gather(stats_by_reach, key=key, value=value, c('nse', 'nrmse', 'rBIAS', 'rrmse', 'kge'))
plot_stats <- filter(plot_stats, key %in% c('nrmse', 'rrmse', 'rBIAS', 'kge'))
plotSWOTreaches <- ggplot(plot_stats, aes(x=key, y=value, fill=key)) +
  geom_boxplot(size=1) +
  geom_hline(yintercept=1, linetype='dashed') +
  geom_hline(yintercept=0, linetype='dashed') +
  geom_hline(yintercept=-0.41, linetype='dashed') +
  xlab('Metric') +
  ylab('Value') +
  coord_cartesian(ylim = c(-1,1))+
  scale_fill_brewer(palette = 'Dark2', name='') +
  theme(legend.position = "none",
        axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'))
 
# example timeseries plot per river------------------------------------------------------------------------------------
set.seed(700)

kge_bins <- quantile(stats_by_reach$kge, c(0.33, 0.66), na.rm = T)
kge_bad <- filter(stats_by_reach, kge <= kge_bins[1]) %>% select(river)
kge_good <- filter(stats_by_reach, kge >= kge_bins[2]) %>% select(river)
kge_eh <- filter(stats_by_reach, kge >= kge_bins[1] & kge <= kge_bins[2]) %>% select(river)

badRiver <- filter(full_output, river == sample(kge_bad$river, 1)) %>%
  gather(key=key, value=value, c(kobs, kest_mean))
badRiverPlot <- ggplot(badRiver, aes(x=time, y=value, color=key)) + #model
  geom_pointrange(aes(ymin = kest_low, ymax = kest_high), fatten=10) +
  geom_line(size=1) + 
  ylab('k600 [m/dy]') +
  xlab('Timestep') +
  scale_color_brewer(palette='Set2') +
  theme(legend.position = "none")

ehRiver <- filter(full_output, river == sample(kge_eh$river, 1)) %>%
  gather(key=key, value=value, c(kobs, kest_mean))
ehRiverPlot <- ggplot(ehRiver, aes(x=time, y=value, color=key)) + #model
  geom_pointrange(aes(ymin = kest_low, ymax = kest_high), fatten=10) +
  geom_line(size=1) +
  ylab('k600 [m/dy]') +
  xlab('Timestep') +
  scale_color_brewer(palette='Set2') +
  theme(legend.position = "none")

goodRiver <- filter(full_output, river == sample(kge_good$river, 1)) %>%
  gather(key=key, value=value, c(kobs, kest_mean))
goodRiverPlot <- ggplot(goodRiver, aes(x=time, y=value, color=key)) + #model
  geom_pointrange(aes(ymin = kest_low, ymax = kest_high), fatten=10) +
  geom_line(size=1) +
  ylab('k600 [m/dy]') +
  xlab('Timestep') +
  scale_color_brewer(palette='Set2', name='k600 [m/dy]', labels=c('Modeled', 'Observed'))

# extract the legend from one of the plots
legend <- get_legend(goodRiverPlot + theme(legend.box.margin = margin(0, 0, 0, 100)))

#bring it alllllll together via cowplot
timeseriesPlot <- plot_grid(goodRiverPlot + theme(legend.position = 'none'), ehRiverPlot, badRiverPlot, legend, ncol=2, labels=c('b','c','d',NA), label_size = 18)
plot2 <- plot_grid(plotSWOTreaches, timeseriesPlot, ncol=2, labels=c('a', NA), label_size = 18)
ggsave('outputs//validation//validation_by_river.jpg', plot2, width=14, height=8)


#Compare BIKER uncertainties against MC uncertanties-------------------------------------------------------
prior_k600_uncertainity <- 1.029
mean_modelSD_by_riv <- group_by(full_output, river) %>%
  summarise(meanSD = mean(kest_sd, na.rm=T))
uncertainity_comparison_alltimesteps <- ggplot(full_output, aes(x=kest_sd)) +
  geom_histogram(size=1, color='black', fill='darkgreen', bins=30) +
  geom_vline(xintercept = prior_k600_uncertainity, linetype='dashed', size=1.3, color='darkblue') +
  geom_vline(xintercept = mean(full_output$kest_sd), size=1.3, color='darkred') +
  xlab('BIKER k600 ln(sd)') +
  ylab('Count')
uncertainity_comparison_byRiv <- ggplot(mean_modelSD_by_riv, aes(x=meanSD)) +
  geom_histogram(size=1, color='black', fill='darkgreen', bins=30) +
  geom_vline(xintercept = prior_k600_uncertainity, linetype='dashed', size=1.3, color='darkblue') +
  geom_vline(xintercept = mean(mean_modelSD_by_riv$meanSD), size=1.3, color='darkred') +
  xlab('BIKER k600 ln(sd)') +
  ylab('Count')

RS_uncertainity <- plot_grid(uncertainity_comparison_alltimesteps, uncertainity_comparison_byRiv, ncol=2, labels = 'auto')
ggsave('outputs//validation/validation_uncertainity.jpg', RS_uncertainity)
