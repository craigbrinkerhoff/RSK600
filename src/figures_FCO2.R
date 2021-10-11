########################
#Description: Produce figures for BIKER FCO2 runs on SWOT-simulated rivers
#Creator: Craig Brinkerhoff
#Date: Fall 2021
######################

#for local runs
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())
#setwd('C:\\Users\\craig\\Documents\\OneDrive - University of Massachusetts\\Ongoing Projects\\RSK600')

print('validating BIKER + CO2 data...')

`%notin%` <- Negate(`%in%`)

###########
##READ IN RESULTS----------------------------------
###########
output <- read.csv('cache/FCO2/CO2_results.csv')
output$river <- as.character(output$river)

#Calculate total wetted surface area over all rivers
total_sa <- output[!duplicated(output$sa), ] #[m2]
total_sa <- sum(total_sa$sa)

#CO2 data is only for 29 samples, so ignore k values beyond that (applies to a few rivers)
output <- drop_na(output)

#Calculate total mass fluxes of CO2 from rivers------------------------------------------------
massFluxes <- gather(output, key=key, value=value, c(FCO2_obs, FCO2_BIKER, FCO2_BIKER_low, FCO2_BIKER_high, FCO2_raymond2012, FCO2_raymond2013, FCO2_brinkerhoff2019)) %>%
  group_by(key) %>%
  summarise(sumFCO2 = mean(value * 1e-12 * 365, na.rm=T)*total_sa, #tG-C/yr from all rivers
            medianFCO2 = median(value * 1e-12 * 365, na.rm=T)*total_sa) #tG-C/yr from all rivers
write.csv(massFluxes, 'cache/FCO2/massFluxes.csv')

##########################
##PLOT MEDIAN/MEAN EFFLUX AND BY-RIVER METRICS PLOT----------------------------------------------------------------------
###########################
#bar plots of mass fluxes----------------------------------------
lowCI_sum <- massFluxes[massFluxes$key == 'FCO2_BIKER_low',]$sumFCO2
lowCI_med <- massFluxes[massFluxes$key == 'FCO2_BIKER_low',]$medianFCO2
highCI_sum <- massFluxes[massFluxes$key == 'FCO2_BIKER_high',]$sumFCO2
highCI_med <- massFluxes[massFluxes$key == 'FCO2_BIKER_high',]$medianFCO2

massFluxes <- massFluxes[massFluxes$key != 'FCO2_BIKER_low' & massFluxes$key != 'FCO2_BIKER_high',]

massFluxes$lowCI_sum <- c(lowCI_sum, NA, NA, NA, NA)
massFluxes$lowCI_med <- c(lowCI_med, NA, NA, NA, NA)
massFluxes$highCI_sum <- c(highCI_sum, NA, NA, NA, NA)
massFluxes$highCI_med <- c(highCI_med, NA, NA, NA, NA)
massFluxes$key <- factor(massFluxes$key,levels = c("FCO2_BIKER", "FCO2_brinkerhoff2019", "FCO2_raymond2012", "FCO2_raymond2013", "FCO2_obs"))

barPlot_mean <- ggplot(massFluxes, aes(y=sumFCO2, x=key, fill=key)) +
  geom_bar(stat='identity', color='black', size=1.2) +
  geom_errorbar(aes(x=key, ymin=lowCI_sum, ymax=highCI_sum), width=0.4, colour="black", size=1.3)+
  geom_hline(yintercept = massFluxes[massFluxes$key=='FCO2_obs', ]$sumFCO2, size=1.2, linetype='dashed', color='#ff7f00') +
  ylim(0,7.25)+
  ylab('Mean [tG-C/yr]') +
  scale_fill_manual(values=c('#e31a1c', '#3f007d', '#6a51a3', '#9e9ac8', '#ff7f00'), name='', labels=c('BIKER \n ', 'Brinkerhoff \n2019', 'Raymond \n2012', 'Raymond \n2013', 'Observed \n ')) +
#  scale_fill_brewer(palette = 'Set1', name='', labels=c('BIKER \n ', 'Brinkerhoff \n2019', 'Raymond \n2012', 'Raymond \n2013', 'Observed \n ')) +
  theme(legend.position = "right",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_text(size=20),
        axis.title.y=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'))

barPlot_median <- ggplot(massFluxes, aes(y=medianFCO2, x=key, fill=key)) +
  geom_bar(stat='identity', color='black', size=1.2) +
  geom_errorbar(aes(x=key, ymin=lowCI_med, ymax=highCI_med), width=0.4, colour="black", size=1.3)+
  geom_hline(yintercept = massFluxes[massFluxes$key=='FCO2_obs', ]$medianFCO2, size=1.2, linetype='dashed', color='#ff7f00') +
  ylim(0,7.25)+
  ylab('Median [tG-C/yr]') +
  scale_fill_manual(values=c('#e31a1c', '#3f007d', '#6a51a3', '#9e9ac8', '#ff7f00'), name='', labels=c('BIKER \n ', 'Brinkerhoff \n2019', 'Raymond \n2012', 'Raymond \n2013', 'Observed \n ')) +
  #scale_fill_brewer(palette = 'Set1', name='', labels=c('BIKER \n ', 'Brinkerhoff \n2019', 'Raymond \n2012', 'Raymond \n2013', 'Observed \n ')) +
  theme(legend.position = "right",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_text(size=20),
        axis.title.y=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'))

barPlots <- plot_grid(barPlot_mean, barPlot_median, ncol=1, labels=c('b', 'c'), label_size=18)

#################
#Scatterplot of differnt HG model results-------------------
################
#massFluxes_byRiver <- group_by(output, river) %>%
#  filter(n() >= 3) %>% #only keep rivers with at least 3 measurements to get by-river medians
#  gather(key=key, value=value, c(FCO2_BIKER, FCO2_raymond2012, FCO2_raymond2013, FCO2_brinkerhoff2019)) %>%
#  group_by(river, key) %>%
#  summarise(sumFCO2 = mean(value * 1e-12 * 365, na.rm=T)*total_sa, #tG-C/yr from river
#            medianFCO2 = median(value * 1e-12 * 365, na.rm=T)*total_sa,  #tG-C/yr from river
#            obsFCO2 = median(FCO2_obs * 1e-12 * 365, na.rm=T)*total_sa) #tG-C/yr from river

#scatterPlot <- ggplot(massFluxes_byRiver, aes(x= obsFCO2, y=medianFCO2, color=key)) +
#  geom_point(size=5, alpha=0.50)+
#  geom_smooth(method='lm', se=F, size=5) +
#  geom_abline(size=2, linetype='dashed') +
#  scale_color_brewer(palette = 'Set1', name='', labels=c('BIKER \n ', 'Brinkerhoff \n2019', 'Raymond \n2012', 'Raymond \n2013')) +
#  xlab('Median river efflux via\nobserved hydraulics [tG-C/yr]') +
#  ylab('Median predicted\nriver efflux [tG-C/yr]')+
#  theme(legend.position = "none",
#        axis.text=element_text(size=20),
#        axis.title=element_text(size=24,face="bold"),
#        legend.text = element_text(size=17),
#        legend.title = element_text(size=17, face='bold'))

stats <- gather(output, key=key, value=value, c(FCO2_BIKER, FCO2_raymond2012, FCO2_raymond2013, FCO2_brinkerhoff2019)) %>%
  group_by(river) %>%
  filter(n() >= 3) %>%
  group_by(river, key) %>%
  summarise(r2 = summary(lm(value~FCO2_obs))$r.squared,
            NRMSE = sqrt(mean((FCO2_obs - value)^2, na.rm=T)) / mean(FCO2_obs, na.rm=T),
            rBIAS =   mean(value- FCO2_obs, na.rm=T) / mean(FCO2_obs, na.rm=T),
            KGE =   KGE(value, FCO2_obs), #sqrt(mean((kobs- kest_mean)^2 / kobs^2, na.rm=T)),
            n_data=n())

plot_stats <- gather(stats, key=key2, value=value2, c('NRMSE', 'rBIAS', 'KGE', 'r2'))
plot_stats <- filter(plot_stats, key2 %in% c('NRMSE', 'KGE', 'rBIAS', 'r2'))
lowKGE <- nrow(plot_stats[plot_stats$key2 == 'KGE' & plot_stats$value2 < -2,])

boxPlots <- ggplot(plot_stats, aes(x=key2, y=value2, fill=key))+
  geom_boxplot(size=1, alpha=0.75) +
  geom_hline(yintercept=1, linetype='dashed') +
  geom_hline(yintercept=0, linetype='dashed') +
  geom_hline(yintercept=0.50, linetype='dashed') +
  xlab('Metric') +
  ylab('Value') +
  scale_fill_manual(values=c('#e31a1c', '#3f007d', '#6a51a3', '#9e9ac8'), name='', labels=c('BIKER \n ', 'Brinkerhoff \n2019', 'Raymond \n2012', 'Raymond \n2013')) +
  annotate("text", x = 'KGE', y = -1.8, label = paste0('# rivers \n< -2: ', lowKGE), size=5)+
  coord_cartesian(ylim=c(-2,2))+
  theme(legend.position = "none",
        axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'))

#######################
#BRING VARIOUS FCO2 FIGURES TOGETHER INTO SINGLE FILE AND WRITE TO DISK-----------------------------
#######################
FCO2_models <- plot_grid(boxPlots, barPlots, ncol=2, labels=c('a', NA), label_size=18)
ggsave('cache/FCO2/FCO2_models.jpg', FCO2_models, width=15, height=10)

##################
##SAVE STATS TO FILES FOR MS MARKDOWN FILE--------------------------------------
#################
#results_all_rivs <- data.frame('rmse'=rmse, 'r2'=lmr2)
#write.csv(results_all_rivs, 'cache/FCO2/fco2_stats_all.csv')


#SCRAPS-----------------
#######################
##PLOT FCO2 VALIDATION ACROSS ALL RIVERS AND TIMESTEPS-------------------------
#######################
#get error stats
#lm <- lm(log10(FCO2_BIKER)~log10(FCO2_obs), data=output)
#lmr2 <- round(summary(lm)$r.squared,2)
#rmse <- round((Metrics::rmse((output$FCO2_obs), (output$FCO2_BIKER))), 2) #g/m2/dy

#get prediction intervals
#predInts <- predict(lm, interval='prediction')
#output <- cbind(output, predInts)

#overall FCO2 validation plot
#flux_plot<- ggplot(output, aes(x=(FCO2_obs), y=(FCO2_BIKER))) +
#  geom_pointrange(aes(ymin = FCO2_BIKER_low, ymax = FCO2_BIKER_high, fill = river), fatten=10, pch=21, color='black') +
#  geom_smooth(size=2, color='grey', method='lm', se=F)+
#  geom_line(aes(y=10^(lwr)), color='grey', linetype='dashed', size=1.75) +
#  geom_line(aes(y=10^(upr)), color='grey', linetype='dashed', size=1.75) +
#  geom_abline(size=2, linetype='dashed', color='darkgrey') +
#  xlab('FCO2 via observed \ndepth [g-C/m2*dy]') +
#  ylab('FCO2 via BIKER [g-C/m2*dy]') +
#  scale_y_log10(
#    limits = c(10^-2, 10^2),
#    breaks = scales::trans_breaks("log10", function(x) 10^x),
#    labels = scales::trans_format("log10", scales::math_format(10^.x)))+
#  scale_x_log10(
#    limits = c(10^-2, 10^2),
#    breaks = scales::trans_breaks("log10", function(x) 10^x),
#    labels = scales::trans_format("log10", scales::math_format(10^.x)))+
#  scale_fill_discrete_qualitative(palette = 'Dark3') +
#  annotate("text", x = 10^-0.4, y = 10^1.5, label = '(Colors correspond \nto rivers)', size=6)+
#  theme(legend.position = "none",
#        axis.text=element_text(size=20),
#        axis.title=element_text(size=24,face="bold"),
#        legend.text = element_text(size=17),
#        legend.title = element_text(size=17, face='bold')) +
  #annotate("text", label = paste0('RMSE: ', rmse, ' g-C/m2*dy'), x = 10^-1, y = 10^1.5, size = 7, colour = "black")+
  #annotate("text", label = paste0('r2: ', lmr2), x = 10^-1, y = 10^1.3, size = 7, colour = "black")
