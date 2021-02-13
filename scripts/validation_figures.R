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
results <- read.csv('outputs//validation//BIKER_validation_results.csv')
full_output <- filter(results, errFlag == 0)

lm_kfit <- lm(log10(full_output$kest_mean)~log10(full_output$kobs))
r2 <- round(summary(lm_kfit)$r.squared, 2)
full_output <- drop_na(full_output)
rmse <- round((Metrics::rmse(log10(full_output$kest_mean), log10(full_output$kobs))), 2)

predInts <- predict(lm_kfit, interval='prediction')
full_output <- cbind(full_output, predInts)

#Model Validation across all rivers and timesteps--------------------------------------------------------
valPlot <- ggplot(full_output, aes(x=(kobs), y=(kest_mean))) +
  geom_abline(size=2, linetype='dashed', color='black') +
  geom_pointrange(aes(ymin = kest_low, ymax = kest_high), fatten=10, fill='#1b9e77', pch=21, color='black') +
  geom_smooth(size=2, color='darkgrey', method='lm', se=F)+
  geom_line(aes(y=10^(lwr)), color='darkgrey', linetype='dashed', size=1.75) +
  geom_line(aes(y=10^(upr)), color='darkgrey', linetype='dashed', size=1.75) +
  xlab('k600 via observed velocity [m/dy]') +
  ylab('BIKER k600 [m/dy]') +
  scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  scale_x_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  scale_color_discrete_qualitative(palette = 'Harmonic') +
  theme(legend.position = "none",
        axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold')) +
  geom_richtext(aes(x=10^-0.25, y=10^1), color='black', label=paste0('RMSE: 10^', rmse, ' m/dy'), size=6) +
  geom_richtext(aes(x=10^-0.25, y=10^0.8), color='black', label=paste0('r<sup>2</sup>: ', r2), size=6)

#cdfs---------------------
t <- gather(full_output, key=key, value=value, c(kobs, kest_mean, kest_high, kest_low))
t$flag <- ifelse(t$key == 'kobs', 1, 0)
k600_cdfs <- ggplot(t, aes(x=value, color=key, linetype=factor(flag))) +
  stat_ecdf(size=1.25) +
  scale_color_manual(values=c('#67a9cf','#67a9cf', '#66c2a5', 'black')) +
  xlab('k600 [m/dy]') +
  ylab('Percentile') +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  geom_hline(yintercept = 0.50, size=1.2) +
  theme(legend.position = "none",
        axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'))

valPlot <- plot_grid(valPlot, k600_cdfs, ncol=2, labels=c('a', 'b'))
ggsave('outputs//validation//validation.jpg', valPlot, width=12, height=7)

#by river metrics---------------------------------------------
#errors <- full_output_errors$kest_mean
#full_output2 <- cbind(full_output, errors)
#full_output2 <- gather(results, key=key, value=value, c('kest_mean', 'errFlag'))
stats_by_reach <- group_by(results, river, errFlag) %>%
  summarise(kge = hydroGOF::KGE(kest_mean, kobs),
            nse =   hydroGOF::NSE(kest_mean, kobs),
            nrmse = sqrt(mean((kobs - kest_mean)^2)) / mean(kobs, na.rm=T),
            rBIAS =   mean(kest_mean- kobs) / mean(kobs, na.rm=T),
            rrmse =   sqrt(mean((kest_mean- kobs)^2 / kobs^2)))

plot_stats <- gather(stats_by_reach, key=key, value=value, c('nse', 'nrmse', 'rBIAS', 'rrmse', 'kge'))
plot_stats <- filter(plot_stats, key %in% c('nrmse', 'rrmse', 'rBIAS', 'kge'))
print(plot_stats)
plotSWOTreaches <- ggplot(plot_stats, aes(x=key, y=value, fill=factor(errFlag))) +
  geom_boxplot(size=1) +
  geom_hline(yintercept=1, linetype='dashed') +
  geom_hline(yintercept=0, linetype='dashed') +
  geom_hline(yintercept=-0.41, linetype='dashed') +
  xlab('Metric') +
  ylab('Value') +
  coord_cartesian(ylim = c(-1,1))+
  scale_fill_brewer(palette = 'Accent', name='', labels=c('No Error', 'Internal + Layover Errors')) +
  theme(legend.position = "bottom",
        axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'))
 
# example timeseries plot per river------------------------------------------------------------------------------------
set.seed(700)

kge_bins <- quantile(stats_by_reach[stats_by_reach$errFlag ==0,]$kge, c(0.33, 0.66), na.rm = T)
kge_bad <- filter(stats_by_reach[stats_by_reach$errFlag ==0,], kge <= kge_bins[1]) %>% select(river)
kge_good <- filter(stats_by_reach[stats_by_reach$errFlag ==0,], kge >= kge_bins[2]) %>% select(river)
kge_eh <- filter(stats_by_reach[stats_by_reach$errFlag ==0,], kge >= kge_bins[1] & kge <= kge_bins[2]) %>% select(river)

badRiver <- filter(full_output, river == sample(kge_bad$river, 1)) %>%
  gather(key=key, value=value, c(kobs, kest_mean))
riv <- substr(as.character(badRiver[1,]$river), 16, nchar(as.character(badRiver[1,]$river)))
badRiverPlot <- ggplot(badRiver, aes(x=time, y=value, color=key)) + #model
  geom_pointrange(aes(ymin = kest_low, ymax = kest_high), fatten=10) +
  geom_line(size=1) + 
  ylab('k600 [m/dy]') +
  xlab('Timestep') +
  scale_color_brewer(palette='Set2') +
  theme(legend.position = "none") +
  ggtitle(riv)

ehRiver <- filter(full_output, river == sample(kge_eh$river, 1)) %>%
  gather(key=key, value=value, c(kobs, kest_mean))
riv <- substr(as.character(ehRiver[1,]$river), 16, nchar(as.character(ehRiver[1,]$river)))
ehRiverPlot <- ggplot(ehRiver, aes(x=time, y=value, color=key)) + #model
  geom_pointrange(aes(ymin = kest_low, ymax = kest_high), fatten=10) +
  geom_line(size=1) +
  ylab('k600 [m/dy]') +
  xlab('Timestep') +
  scale_color_brewer(palette='Set2') +
  theme(legend.position = "none") +
  ggtitle(riv)

goodRiver <- filter(full_output, river == sample(kge_good$river, 1)) %>%
  gather(key=key, value=value, c(kobs, kest_mean))
riv <- substr(as.character(goodRiver[1,]$river), 16, nchar(as.character(goodRiver[1,]$river)))
goodRiverPlot <- ggplot(goodRiver, aes(x=time, y=value, color=key)) + #model
  geom_pointrange(aes(ymin = kest_low, ymax = kest_high), fatten=10) +
  geom_line(size=1) +
  ylab('k600 [m/dy]') +
  xlab('Timestep') +
  scale_color_brewer(palette='Set2', name='k600 [m/dy]', labels=c('BIKER', 'Model Using Observed \nVelocity')) +
  ggtitle(riv)

# extract the legend from one of the plots
legend <- get_legend(goodRiverPlot + theme(legend.box.margin = margin(0, 0, 0, 100)))

#bring it alllllll together via cowplot
timeseriesPlot <- plot_grid(goodRiverPlot + theme(legend.position = 'none'), ehRiverPlot, badRiverPlot, legend, ncol=2, labels=c('b','c','d',NA), label_size = 18)
plot2 <- plot_grid(plotSWOTreaches, timeseriesPlot, ncol=2, labels=c('a', NA), label_size = 18)
ggsave('outputs//validation//validation_by_river.jpg', plot2, width=14, height=8)

#save results for ms---------------------------------------------------
results_3_2 <- data.frame('rmse'=rmse, 'r2'=r2)
write.csv(results_3_2, 'outputs//validation//results_all_riv.csv')
write.csv(stats_by_reach, 'outputs//validation//results_by_riv.csv')



#for finesst
#slope versus temporal variation of k600
# temp <- gather(stats_by_reach, key=key, value=value, c('CVobs', 'CVest'))
# ggplot(temp, aes(x=meanS, y=value*100, color=key)) + 
#   geom_point(size=5) +
#   geom_smooth(se=F, method='lm', size=1.5)+
#   scale_color_brewer(palette='Dark2', name='', labels=c('BIKER', 'Observed'))+
#   ylim(0,100)+
#   scale_x_log10(
#     breaks = scales::trans_breaks("log10", function(x) 10^x),
#     labels = scales::trans_format("log10", scales::math_format(10^.x))) +
#   ylab('CV K600 [%]') +
#   xlab('Mean SWOT slope') +
#   theme(legend.position = 'bottom')
