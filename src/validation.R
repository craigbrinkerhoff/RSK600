########################
#Description: Produce figures for BIKER runs on SWOT-simulated rivers
#Creator: Craig Brinkerhoff
#Date: Fall 2020
######################

print('validating BIKER...')

################
##READ IN RESULTS----------------------------------
###############
results <- read.csv('cache/validation/BIKER_validation_results.csv')
full_output <- filter(results, errFlag == 0) #remove results with SWOT measurement error

lm_kfit <- lm(log10(full_output$kest_mean)~log10(full_output$kobs))
r2 <- round(summary(lm_kfit)$r.squared, 2)
full_output <- drop_na(full_output)
rmse <- round((Metrics::rmse((full_output$kest_mean), (full_output$kobs))), 2)

predInts <- predict(lm_kfit, interval='prediction')
full_output <- cbind(full_output, predInts)

####################
##MODEL VALIDATION ACROSS ALL RIVERS AND TIMESTEPS--------------------------------------------------------
####################
valPlot <- ggplot(full_output, aes(x=(kobs), y=(kest_mean))) +
  geom_abline(size=2, linetype='dashed', color='black') +
  geom_pointrange(aes(ymin = kest_low, ymax = kest_high), fatten=10, fill='#1b9e77', pch=21, color='black') +
  geom_smooth(size=2, color='darkgrey', method='lm', se=F)+
  geom_line(aes(y=10^(lwr)), color='darkgrey', linetype='dashed', size=1.75) +
  geom_line(aes(y=10^(upr)), color='darkgrey', linetype='dashed', size=1.75) +
  xlab('ko2 via observed \nshear velocity [m/dy]') +
  ylab('BIKER ko2 [m/dy]') +
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
  annotate("text", label = paste0('RMSE: ', rmse, ' m/dy'), x = 10^0, y = 10^1, size = 5, colour = "black")+
  annotate("text", label = paste0('r2: ', r2), x = 10^0, y = 10^0.8, size = 5, colour = "black")

###################
##CDFS OF BIKER VS OBSERVED---------------------
##################
t <- gather(full_output, key=key, value=value, c(kobs, kest_mean, kest_high, kest_low))
t$flag <- ifelse(t$key == 'kobs', 1, 0)
k_cdfs <- ggplot(t, aes(x=value, color=key, linetype=factor(flag))) +
  stat_ecdf(size=1.25) +
  scale_color_manual(values=c('#67a9cf','#67a9cf', '#66c2a5', 'black')) +
  xlab('ko2 [m/dy]') +
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

valPlot <- plot_grid(valPlot, k_cdfs, ncol=2, labels=c('a', 'b'), label_size=18)
ggsave('cache/validation/validation.jpg', valPlot, width=12, height=7)

########################
##CALCULATE BY-RIVER ERROR METRICS---------------------------------------------
########################
stats_by_reach <- group_by(results, river, errFlag) %>%
  summarise(kge = hydroGOF::KGE(kest_mean, kobs),
            nrmse = sqrt(mean((kobs - kest_mean)^2, na.rm=T)) / mean(kobs, na.rm=T),
            rBIAS =   mean(kest_mean- kobs, na.rm=T) / mean(kobs, na.rm=T),
            rrmse =   sqrt(mean((kobs- kest_mean)^2 / kobs^2, na.rm=T)),
            meanKobs = mean(kobs, na.rm=T),
            meanWobs = mean(Wobs, na.rm=T),
            meanSobs = mean(Sobs, na.rm=T),
            RE = hydroGOF::rmse(kest_mean, kobs, na.rm=T)/mean(kobs, na.rm=T))

plot_stats <- gather(stats_by_reach, key=key, value=value, c('nrmse', 'rBIAS', 'rrmse', 'kge'))
plot_stats <- filter(plot_stats, key %in% c('nrmse', 'rrmse', 'rBIAS', 'kge'))

########################
##PLOT BY-RIVER ERROR METRICS---------------------------------------
########################
plotSWOTreaches <- ggplot(plot_stats, aes(x=key, y=value, fill=factor(errFlag))) +
  geom_boxplot(size=1) +
  geom_hline(yintercept=1, linetype='dashed') +
  geom_hline(yintercept=0, linetype='dashed') +
  geom_hline(yintercept=-0.41, linetype='dashed') +
  xlab('Metric') +
  ylab('Value') +
#  coord_cartesian(ylim = c(-1.5,1))+
  scale_fill_brewer(palette = 'Accent', name='', labels=c('No Error', 'SWOT Measurement Error')) +
  theme(legend.position = "bottom",
        axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'))

#################
##PLOT EXAMPLE TIMESERIES PLOTS FOR 3 RIVERS------------------------------------------------------------------------------------
#################
set.seed(5)

#Randomly select 3 rivers representative of good, eh, and bad KGE scores
kge_bins <- quantile(stats_by_reach[stats_by_reach$errFlag ==0,]$kge, c(0.33, 0.66), na.rm = T)
kge_bad <- filter(stats_by_reach[stats_by_reach$errFlag ==0,], kge <= kge_bins[1]) %>% select(river)
kge_good <- filter(stats_by_reach[stats_by_reach$errFlag ==0,], kge >= kge_bins[2]) %>% select(river)
kge_eh <- filter(stats_by_reach[stats_by_reach$errFlag ==0,], kge >= kge_bins[1] & kge <= kge_bins[2]) %>% select(river)

#bad river plot
badRiver <- filter(full_output, river == sample(kge_bad$river, 1)) %>%
  gather(key=key, value=value, c(kobs, kest_mean))
badRiver$kest_low <- ifelse(badRiver$key == 'kobs', NA, badRiver$kest_low)
badRiver$kest_high <- ifelse(badRiver$key == 'kobs', NA, badRiver$kest_high)
riv <- as.character(badRiver[1,]$river)
badRiverPlot <- ggplot(badRiver, aes(x=time, y=value, color=key)) + #model
  geom_point(size=3) +
  geom_ribbon(aes(ymin = kest_low, ymax = kest_high), alpha=0.75, fill='grey')+
  geom_line(size=1) +
  ylab('ko2 [m/dy]') +
  xlab('Timestep') +
  scale_color_brewer(palette='Set2') +
  theme(legend.position = "none") +
  ggtitle(riv)

#eh river plot
ehRiver <- filter(full_output, river == sample(kge_eh$river, 1)) %>%
  gather(key=key, value=value, c(kobs, kest_mean))
ehRiver$kest_low <- ifelse(ehRiver$key == 'kobs', NA, ehRiver$kest_low)
ehRiver$kest_high <- ifelse(ehRiver$key == 'kobs', NA, ehRiver$kest_high)
riv <- as.character(ehRiver[1,]$river)
ehRiverPlot <- ggplot(ehRiver, aes(x=time, y=value, color=key)) + #model
  geom_point(size=3) +
  geom_ribbon(aes(ymin = kest_low, ymax = kest_high), alpha=0.75, fill='grey')+
  geom_line(size=1) +
  ylab('ko2 [m/dy]') +
  xlab('Timestep') +
  scale_color_brewer(palette='Set2') +
  theme(legend.position = "none") +
  ggtitle(riv)

#good river plot
goodRiver <- filter(full_output, river == sample(kge_good$river, 1)) %>%
  gather(key=key, value=value, c(kobs, kest_mean))
goodRiver$kest_low <- ifelse(goodRiver$key == 'kobs', NA, goodRiver$kest_low)
goodRiver$kest_high <- ifelse(goodRiver$key == 'kobs', NA, goodRiver$kest_high)
riv <- as.character(goodRiver[1,]$river)
goodRiverPlot <- ggplot(goodRiver, aes(x=time, y=value, color=key)) + #model
  geom_point(size=3) +
  geom_ribbon(aes(ymin = kest_low, ymax = kest_high), alpha=0.75, fill='grey')+
  geom_line(size=1) +
  ylab('ko2 [m/dy]') +
  xlab('Timestep') +
  scale_color_brewer(palette='Set2', name='ko2 [m/dy]', labels=c('BIKER', 'Model Using \nObserved Velocity')) +
  ggtitle(riv)

# extract the legend from one of the plots
legend <- get_legend(goodRiverPlot + theme(legend.box.margin = margin(0, 0, 0, 100)))

#bring it alllllll together via cowplot
timeseriesPlot <- plot_grid(goodRiverPlot + theme(legend.position = 'none'), ehRiverPlot, badRiverPlot, legend, ncol=2, labels=c('b','c','d',NA), label_size = 18)
plot2 <- plot_grid(plotSWOTreaches, timeseriesPlot, ncol=2, labels=c('a', NA), label_size = 18)
ggsave('cache/validation/validation_by_river.jpg', plot2, width=14, height=8)

####################
##SAVE RESULTS FOR MS---------------------------------------------------
###################
results_3_2 <- data.frame('rmse'=rmse, 'r2'=r2)
write.csv(results_3_2, 'cache/validation/results_all_riv.csv')
write.csv(stats_by_reach, 'cache/validation/results_by_riv.csv')

################
##RBIAS VS HYDRAULIC PROPERTIES---------------------------------------
################
riverPropertiesPlot_k <- ggplot(stats_by_reach, aes(y=rBIAS, x=meanKobs)) +
  geom_point(size=3) +
  geom_hline(yintercept=0, linetype='dashed') +
  xlab('Mean observed k [m/dy]') +
  ylab('rBIAS') +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'))
riverPropertiesPlot_W <- ggplot(stats_by_reach, aes(y=rBIAS, x=meanWobs)) +
  geom_point(size=3) +
  geom_hline(yintercept=0, linetype='dashed') +
  xlab('Mean observed width [m]') +
  ylab('rBIAS') +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  theme(legend.position = "none",
        axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'))
riverPropertiesPlot_S <- ggplot(stats_by_reach, aes(y=rBIAS, x=meanSobs)) +
  geom_point(size=3) +
  geom_hline(yintercept=0, linetype='dashed') +
  xlab('Mean observed slope') +
  ylab('rBIAS') +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  theme(legend.position = "none",
        axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'))

# extract the legend from one of the plots
legend <- get_legend(riverPropertiesPlot_k + theme(legend.box.margin = margin(0, 0, 0, 100)))

#bring it alllllll together via cowplot
riverPropertiesPlot <- plot_grid(riverPropertiesPlot_k + theme(legend.position = 'none'), riverPropertiesPlot_W, riverPropertiesPlot_S, legend, ncol=2, labels=c('b','c','d',NA), label_size = 18)
ggsave('cache/validation/riverProperties_rbias.jpg', riverPropertiesPlot, width=14, height=8)

#########################
##KGE VS HYDRAULIC properties-------------------------------------------------------
#########################
riverPropertiesPlot_k <- ggplot(stats_by_reach, aes(y=kge, x=meanKobs)) +
  geom_point(size=3) +
  xlab('Mean observed k [m/dy]') +
  ylab('KGE') +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'))
riverPropertiesPlot_W <- ggplot(stats_by_reach, aes(y=kge, x=meanWobs)) +
  geom_point(size=3) +
  xlab('Mean observed width [m]') +
  ylab('KGE') +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  theme(legend.position = "none",
        axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'))
riverPropertiesPlot_S <- ggplot(stats_by_reach, aes(y=kge, x=meanSobs)) +
  geom_point(size=3) +
  xlab('Mean observed slope') +
  ylab('KGE') +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  theme(legend.position = "none",
        axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'))

# extract the legend from one of the plots
legend <- get_legend(riverPropertiesPlot_k + theme(legend.box.margin = margin(0, 0, 0, 100)))

#bring it alllllll together via cowplot
riverPropertiesPlot <- plot_grid(riverPropertiesPlot_k + theme(legend.position = 'none'), riverPropertiesPlot_W, riverPropertiesPlot_S, legend, ncol=2, labels=c('b','c','d',NA), label_size = 18)
ggsave('cache/validation/riverProperties_kge.jpg', riverPropertiesPlot, width=14, height=8)

#########################
##RRMSE VS HYDRAULIC properties-------------------------------------------------------
#########################
riverPropertiesPlot_k <- ggplot(stats_by_reach, aes(y=rrmse, x=meanKobs)) +
  geom_point(size=3) +
  xlab('Mean observed k [m/dy]') +
  ylab('RRMSE') +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'))
riverPropertiesPlot_W <- ggplot(stats_by_reach, aes(y=rrmse, x=meanWobs)) +
  geom_point(size=3) +
  xlab('Mean observed width [m]') +
  ylab('RRMSE') +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  theme(legend.position = "none",
        axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'))
riverPropertiesPlot_S <- ggplot(stats_by_reach, aes(y=rrmse, x=meanSobs)) +
  geom_point(size=3) +
  xlab('Mean observed slope') +
  ylab('RRMSE') +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  theme(legend.position = "none",
        axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'))

# extract the legend from one of the plots
legend <- get_legend(riverPropertiesPlot_k + theme(legend.box.margin = margin(0, 0, 0, 100)))

#bring it alllllll together via cowplot
riverPropertiesPlot <- plot_grid(riverPropertiesPlot_k + theme(legend.position = 'none'), riverPropertiesPlot_W, riverPropertiesPlot_S, legend, ncol=2, labels=c('b','c','d',NA), label_size = 18)
ggsave('cache/validation/riverProperties_rrmse.jpg', riverPropertiesPlot, width=14, height=8)

#########################
##NRMSE VS HYDRAULIC properties-------------------------------------------------------
#########################
riverPropertiesPlot_k <- ggplot(stats_by_reach, aes(y=nrmse, x=meanKobs)) +
  geom_point(size=3) +
  xlab('Mean observed k [m/dy]') +
  ylab('NRMSE') +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'))
riverPropertiesPlot_W <- ggplot(stats_by_reach, aes(y=nrmse, x=meanWobs)) +
  geom_point(size=3) +
  xlab('Mean observed width [m]') +
  ylab('NRMSE') +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  theme(legend.position = "none",
        axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'))
riverPropertiesPlot_S <- ggplot(stats_by_reach, aes(y=nrmse, x=meanSobs)) +
  geom_point(size=3) +
  xlab('Mean observed slope') +
  ylab('NRMSE') +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  theme(legend.position = "none",
        axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'))

# extract the legend from one of the plots
legend <- get_legend(riverPropertiesPlot_k + theme(legend.box.margin = margin(0, 0, 0, 100)))

#bring it alllllll together via cowplot
riverPropertiesPlot <- plot_grid(riverPropertiesPlot_k + theme(legend.position = 'none'), riverPropertiesPlot_W, riverPropertiesPlot_S, legend, ncol=2, labels=c('b','c','d',NA), label_size = 18)
ggsave('cache/validation/riverProperties_nrmse.jpg', riverPropertiesPlot, width=14, height=8)

###################################
##NRMSE/% RELATIVE ERROR VERSUS WANG ETAL 2021 RELATIVE ERROR--------------------------------------------
###################################
relativeError_plot <- ggplot(stats_by_reach[stats_by_reach$errFlag == 0,], aes(y=RE*100, x=river)) +
  geom_col(color='black') +
  coord_flip()+
  geom_hline(yintercept=mean(c(57.31, 57.55)), linetype='dashed', size=2) +
  xlab('River') +
  ylab('% Relative Error') +
  theme(legend.position = "none",
        axis.text.x=element_text(size=20),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'))
ggsave('cache/validation/relativeErr.jpg', relativeError_plot, width=14, height=8)
print(stats_by_reach[stats_by_reach$errFlag == 0,])
