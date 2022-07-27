########################
#Description: Produce figures for BIKER FCO2 runs on SWOT-simulated rivers
#Creator: Craig Brinkerhoff
#Date: Winter 2022
######################

print('validating BIKER + CO2 data...')

`%notin%` <- Negate(`%in%`)

theme_set(theme_classic())

###################
##READ IN RESULTS----------------------------------
####################
output <- read.csv('cache/FCO2/CO2_results.csv')
output$river <- as.character(output$river)

#Calculate total wetted surface area over all rivers
total_sa <- output[!duplicated(output$sa), ] #[m2]
total_sa <-sum(total_sa$sa)

#CO2 data is only for 29 samples, so ignore k values beyond that (applies to some rivers)
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

#MAKE ANOTHER COLUMN FOR CIs AND THEN ADD THEM
massFluxes$lowCI_sum <- c(lowCI_sum, lowCI_sum_brink, NA, NA, NA)
massFluxes$highCI_sum <- c(highCI_sum, highCI_sum_brink, NA, NA, NA)

massFluxes$key <- factor(massFluxes$key,levels = c("FCO2_BIKER", "FCO2_brinkerhoff2019", "FCO2_raymond2012", "FCO2_raymond2013", "FCO2_obs"))

barPlot <- ggplot(massFluxes, aes(y=yearlyCFux, x=key, fill=key)) +
  geom_bar(stat='identity', color='black', size=1.2) +
  geom_errorbar(aes(x=key, ymin=lowCI_sum, ymax=highCI_sum), width=0.4, colour="black", size=1.3)+
  geom_hline(yintercept = massFluxes[massFluxes$key=='FCO2_obs', ]$yearlyCFux, size=1.2, linetype='dashed', color='#ff7f00') +
  ylab('CO2 Flux [Tg-C/yr]') +
  scale_fill_manual(values=c('#e31a1c', '#3f007d', '#6a51a3', '#9e9ac8', '#ff7f00'), name='', labels=c('BIKER \n ', 'Brinkerhoff \n2019', 'Raymond \n2012', 'Raymond \n2013', '"Observed" \n ')) +
  theme(legend.position = "right",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_text(size=20),
        axis.title.y=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'))

#######################
#BRING VARIOUS FCO2 FIGURES TOGETHER INTO SINGLE FILE AND WRITE TO DISK-----------------------------
#######################
ggsave('cache/FCO2/fig7.jpg', barPlot, width=8, height=8)
