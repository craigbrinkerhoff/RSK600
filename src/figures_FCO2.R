########################
#Description: Produce figures for BIKER FCO2 runs on SWOT-simulated rivers
#Creator: Craig Brinkerhoff
#Date: Winter 2022
######################

print('validating BIKER + CO2 data...')

`%notin%` <- Negate(`%in%`)

theme_set(theme_classic())

###########
##READ IN RESULTS----------------------------------
###########
output <- read.csv('cache/FCO2/CO2_results.csv')
output$river <- as.character(output$river)

#Calculate total wetted surface area over all rivers
total_sa <- output[!duplicated(output$sa), ] #[m2]
total_sa <-sum(total_sa$sa)

#CO2 data is only for 29 samples, so ignore k values beyond that (applies to a few rivers)
output <- drop_na(output)

#Calculate total mass fluxes of CO2 from rivers------------------------------------------------
massFluxes <- gather(output, key=key, value=value, c(FCO2_obs, FCO2_BIKER, FCO2_BIKER_low, FCO2_BIKER_high, FCO2_raymond2012, FCO2_raymond2013, FCO2_brinkerhoff2019, FCO2_brinkerhoff2019_low, FCO2_brinkerhoff2019_high)) %>%
  group_by(river, key) %>%
  summarise(yearlyCFux_byriv = mean(value*1e-12*sa_m2, na.rm=T)) %>% #tG-C/yr from all rivers
  group_by(key) %>%
  summarise(yearlyCFux = sum(yearlyCFux_byriv, na.rm=T))
write.csv(massFluxes, 'cache/FCO2/massFluxes.csv')

##########################
##PLOT MEDIAN/MEAN EFFLUX AND BY-RIVER METRICS PLOT----------------------------------------------------------------------
###########################
#bar plots of mass fluxes----------------------------------------
lowCI_sum <- massFluxes[massFluxes$key == 'FCO2_BIKER_low',]$yearlyCFux
highCI_sum <- massFluxes[massFluxes$key == 'FCO2_BIKER_high',]$yearlyCFux
massFluxes <- massFluxes[massFluxes$key != 'FCO2_BIKER_low' & massFluxes$key != 'FCO2_BIKER_high',]

lowCI_sum_brink <- massFluxes[massFluxes$key == 'FCO2_brinkerhoff2019_low',]$yearlyCFux
highCI_sum_brink <- massFluxes[massFluxes$key == 'FCO2_brinkerhoff2019_high',]$yearlyCFux
massFluxes <- massFluxes[massFluxes$key != 'FCO2_brinkerhoff2019_low' & massFluxes$key != 'FCO2_brinkerhoff2019_high',]

#MAKE ANOTHER COLUMN FOR BRINK 2019 ci AND THENADD IT
massFluxes$lowCI_sum <- c(lowCI_sum, lowCI_sum_brink, NA, NA, NA)
massFluxes$highCI_sum <- c(highCI_sum, highCI_sum_brink, NA, NA, NA)


#massFluxes$lowCI_sum_brink <- c(NA, NA, NA, lowCI_sum_brink, NA, NA, NA)
#massFluxes$highCI_sum_brink <- c(NA, NA, NA, highCI_sum_brink, NA, NA, NA)
#print(massFluxes)
massFluxes$key <- factor(massFluxes$key,levels = c("FCO2_BIKER", "FCO2_brinkerhoff2019", "FCO2_raymond2012", "FCO2_raymond2013", "FCO2_obs"))

barPlot <- ggplot(massFluxes, aes(y=yearlyCFux, x=key, fill=key)) +
  geom_bar(stat='identity', color='black', size=1.2) +
  geom_errorbar(aes(x=key, ymin=lowCI_sum, ymax=highCI_sum), width=0.4, colour="black", size=1.3)+
  geom_hline(yintercept = massFluxes[massFluxes$key=='FCO2_obs', ]$yearlyCFux, size=1.2, linetype='dashed', color='#ff7f00') +
  ylab('CO2 Flux [Tg-C/yr]') +
  scale_fill_manual(values=c('#e31a1c', '#3f007d', '#6a51a3', '#9e9ac8', '#ff7f00'), name='', labels=c('BIKER \n ', 'Brinkerhoff \n2019', 'Raymond \n2012', 'Raymond \n2013', 'Observed \n ')) +
  theme(legend.position = "right",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_text(size=20),
        axis.title.y=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'))

# stats <- gather(output, key=key, value=value, c(FCO2_BIKER, FCO2_raymond2012, FCO2_raymond2013, FCO2_brinkerhoff2019)) %>%
#   group_by(river) %>%
#   filter(n() >= 3) %>%
#   group_by(river, key) %>%
#   summarise(r = sqrt(summary(lm(value~FCO2_obs))$r.squared),
#             NRMSE = sqrt(mean((FCO2_obs - value)^2, na.rm=T)) / mean(FCO2_obs, na.rm=T),
#             rBIAS =   mean(value- FCO2_obs, na.rm=T) / mean(FCO2_obs, na.rm=T),
#             KGE =   KGE(value, FCO2_obs),
#             n_data=n())
# 
# plot_stats <- gather(stats, key=key2, value=value2, c('NRMSE', 'rBIAS', 'KGE', 'r'))
# plot_stats <- filter(plot_stats, key2 %in% c('NRMSE', 'KGE', 'rBIAS', 'r'))
# lowKGE <- nrow(plot_stats[plot_stats$key2 == 'KGE' & plot_stats$value2 < -2,])
# 
# boxPlots <- ggplot(plot_stats, aes(x=key2, y=value2, fill=key))+
#   geom_boxplot(size=1, alpha=0.75) +
#   geom_hline(yintercept=1, linetype='dashed') +
#   geom_hline(yintercept=0, linetype='dashed') +
#   geom_hline(yintercept=0.50, linetype='dashed') +
#   xlab('Metric') +
#   ylab('Value') +
#   scale_fill_manual(values=c('#e31a1c', '#3f007d', '#6a51a3', '#9e9ac8'), name='', labels=c('BIKER \n ', 'Brinkerhoff \n2019', 'Raymond \n2012', 'Raymond \n2013')) +
#   annotate("text", x = 'KGE', y = -1.8, label = paste0('# rivers \n< -2: ', lowKGE), size=5)+
#   coord_cartesian(ylim=c(-2,2))+
#   theme(legend.position = "none",
#         axis.text=element_text(size=20),
#         axis.title=element_text(size=24,face="bold"),
#         legend.text = element_text(size=17),
#         legend.title = element_text(size=17, face='bold'))

#######################
#BRING VARIOUS FCO2 FIGURES TOGETHER INTO SINGLE FILE AND WRITE TO DISK-----------------------------
#######################
#FCO2_models <- plot_grid(boxPlots, barPlot, ncol=2, labels=c('a', 'b'), label_size=18)
ggsave('cache/FCO2/fig8.jpg', barPlot, width=8, height=8)
#write.csv(stats, 'cache/FCO2/results_by_riv.csv')
