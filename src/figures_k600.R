########################
#Description: Produce figures for BIKER runs on SWOT-simulated rivers
#Creator: Craig Brinkerhoff
#Date: Fall 2021
######################

print('validating BIKER...')

`%notin%` <- Negate(`%in%`)

################
##READ IN RESULTS----------------------------------
###############
results <- read.csv('cache/validation/BIKER_validation_results.csv')
full_output <- filter(results, errFlag == 0) #remove results with SWOT measurement error
full_output <- group_by(full_output, river) %>%
  filter(Wobs >= 100) %>%
  ungroup()

#full_output <- group_by(full_output, river) %>%
#  slice(seq(1, n(), samplingRate))

#Calculate r2 and rmse metrics
#r2 <- round(summary(lm_kfit)$r.squared, 2)
#full_output <- drop_na(full_output)
#rmse <- round((Metrics::rmse((full_output$kest_mean), (full_output$kobs))), 2)

#Calculate prediction intervals for plot
#predInts <- predict(lm_kfit, interval='prediction')
#full_output <- cbind(full_output, predInts)

####################
##MODEL VALIDATION ACROSS ALL RIVERS AND TIMESTEPS--------------------------------------------------------
####################
valPlot <- ggplot(full_output, aes(x=(kobs), y=(kest_mean))) +
  geom_pointrange(aes(ymin = kest_low, ymax = kest_high, fill=river), fatten=10, pch=21, color='black') +
#  geom_smooth(size=2, color='darkgrey', method='lm', se=F)+
#  geom_line(aes(y=10^(lwr)), color='darkgrey', linetype='dashed', size=1.75) +
#  geom_line(aes(y=10^(upr)), color='darkgrey', linetype='dashed', size=1.75) +
  geom_abline(size=3, linetype='dashed', color='darkgrey') +
  xlab('k600 via observed \nhydraulics [m/dy]') +
  ylab('BIKER k600 [m/dy]') +
  scale_y_log10(limits=c(10^-1,10^2),
                breaks=c(0.1, 1, 10, 100),
                labels=c('0.1', '1', '10', '100'))+
  scale_x_log10(limits=c(10^-1,10^2),
                breaks=c(0.1, 1, 10, 100),
                labels=c('0.1', '1', '10', '100'))+
  scale_fill_discrete_qualitative(palette = 'Dark3') +
  annotate("text", x = 10^-0.4, y = 10^1.5, label = '(Colors correspond \nto rivers)', size=6)+
  theme(legend.position = "none",
        axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'))
ggsave('cache/validation/validation.jpg', valPlot, width=8, height=8)

########################
##CALCULATE BY-RIVER ERROR METRICS---------------------------------------------
########################
stats_by_river <- group_by(results, river, errFlag) %>%
  summarise(r = sqrt(summary(lm(kest_mean~kobs))$r.squared),
            NRMSE = sqrt(mean((kobs - kest_mean)^2, na.rm=T)) / mean(kobs, na.rm=T),
            rBIAS =   mean(kest_mean- kobs, na.rm=T) / mean(kobs, na.rm=T),
            KGE = KGE(kest_mean, kobs),
            CV = sd(kobs)/mean(kobs),
            meanKobs = mean(kobs, na.rm=T),
            meanWobs = mean(Wobs, na.rm=T),
            meanSobs = mean(Sobs, na.rm=T),
            meanDobs = mean(Dobs, na.rm=T),
            n_data=n())

plot_stats <- gather(stats_by_river, key=key, value=value, c('NRMSE', 'rBIAS', 'KGE', 'r'))
plot_stats <- filter(plot_stats, key %in% c('NRMSE', 'KGE', 'rBIAS', 'r'))

########################
##PLOT BY-RIVER ERROR METRICS---------------------------------------
########################
lowKGE_0 <- nrow(plot_stats[plot_stats$key == 'KGE' & plot_stats$value < -2 & plot_stats$errFlag==0,])
lowKGE_1 <- nrow(plot_stats[plot_stats$key == 'KGE' & plot_stats$value < -2 & plot_stats$errFlag==1,])
plotSWOTreaches <- ggplot(plot_stats, aes(x=key, y=value, fill=factor(errFlag))) +
  geom_boxplot(size=1, alpha=0.75) +
  geom_hline(yintercept=1, linetype='dashed') +
  geom_hline(yintercept=0, linetype='dashed') +
  geom_hline(yintercept=0.50, linetype='dashed') +
  xlab('Metric') +
  ylab('Value') +
  coord_cartesian(ylim=c(-2,2))+
  annotate("text", x = 'KGE', y = -2, label = paste0('# rivers \n< -2: ', lowKGE_0), size=5, color='#7fc97f')+
  annotate("text", x = 'KGE', y = -1.6, label = paste0('# rivers \n< -2: ', lowKGE_1), size=5, color='#beaed4')+
  scale_fill_brewer(palette = 'Accent', name='', labels=c('No Error', 'SWOT Measurement Error')) +
  theme(legend.position = "bottom",
        axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'))

############
##r VS AMOUNT OF DATA-----------------------
############
r_vs_n <- ggplot(stats_by_river[stats_by_river$errFlag == 0,]) +
  geom_point(aes(y=r, x=n_data), color='#7fc97f',size=6) +
  geom_hline(aes(yintercept = median(stats_by_river[stats_by_river$errFlag==0,]$r)), linetype='dashed', color='darkblue', size=1.5)+
  ylab('r') +
  xlab('Num. SWOT Observations') +
  theme(legend.position = "bottom",
        axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'))

############
##r VS AMOUNT OF DATA-----------------------
############
r_vs_w <- ggplot(stats_by_river[stats_by_river$errFlag==0,]) +
  geom_point(aes(y=r, x=meanWobs),color='#7fc97f', size=6) +
  geom_hline(aes(yintercept = median(stats_by_river[stats_by_river$errFlag==0,]$r)), linetype='dashed', color='darkblue', size=1.5)+
  xlab('Mean observed width [m]') +
  ylab('r') +
  scale_x_log10()+
  theme(legend.position = "none",
        axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'))

temp <- plot_grid(r_vs_n, r_vs_w, ncol=1, labels=c('b', 'c'), label_size=18)
valPlot_byRiver <- plot_grid(plotSWOTreaches, temp, ncol=2, labels=c('a', NA), label_size=18)
ggsave('cache/validation/validation_by_river.jpg', valPlot_byRiver, width=12, height=8)

#################
##EXAMPLE TIMESERIES PLOTS FOR 3 RIVERS------------------------------------------------------------------------------------
#################
set.seed(113)

# #Randomly select 3 rivers representative of good, eh, and bad KGE scores
# kge_bins <- quantile(stats_by_river[stats_by_river$errFlag ==0,]$KGE, c(0.33, 0.66), na.rm = T)
# kge_bad <- filter(stats_by_river[stats_by_river$errFlag ==0,], KGE <= kge_bins[1]) %>% select(river)
# kge_good <- filter(stats_by_river[stats_by_river$errFlag ==0,], KGE >= kge_bins[2]) %>% select(river)
# kge_eh <- filter(stats_by_river[stats_by_river$errFlag ==0,], KGE >= kge_bins[1] & KGE <= kge_bins[2]) %>% select(river)
# 
# #bad river plot
# badRiver <- filter(full_output, river == sample(kge_bad$river, 1)) %>%
#   gather(key=key, value=value, c(kobs, kest_mean))
# badRiver$kest_low <- ifelse(badRiver$key == 'kobs', NA, badRiver$kest_low)
# badRiver$kest_high <- ifelse(badRiver$key == 'kobs', NA, badRiver$kest_high)
# riv <- as.character(badRiver[1,]$river)
# badRiverPlot <- ggplot(badRiver, aes(x=time, y=value, color=key)) + #model
#   geom_ribbon(aes(ymin = kest_low, ymax = kest_high), alpha=0.75, fill='grey')+
#   geom_line(size=1.5) +
#   ylab('k600 [m/dy]') +
#   xlab('Timestep') +
#   scale_color_brewer(palette='Set2') +
#   theme(legend.position = "none") +
#   ggtitle(paste0('Bad River: ', riv))
# 
# #eh river plot
# ehRiver <- filter(full_output, river == sample(kge_eh$river, 1)) %>%
#   gather(key=key, value=value, c(kobs, kest_mean))
# ehRiver$kest_low <- ifelse(ehRiver$key == 'kobs', NA, ehRiver$kest_low)
# ehRiver$kest_high <- ifelse(ehRiver$key == 'kobs', NA, ehRiver$kest_high)
# riv <- as.character(ehRiver[1,]$river)
# ehRiverPlot <- ggplot(ehRiver, aes(x=time, y=value, color=key)) + #model
#   geom_ribbon(aes(ymin = kest_low, ymax = kest_high), alpha=0.75, fill='grey')+
#   geom_line(size=1.5) +
#   ylab('k600 [m/dy]') +
#   xlab('Timestep') +
#   scale_color_brewer(palette='Set2') +
#   theme(legend.position = "none") +
#   ggtitle(paste0('OK River: ', riv))
# 
# #good river plot
# goodRiver <- filter(full_output, river == sample(kge_good$river, 1)) %>%
#   gather(key=key, value=value, c(kobs, kest_mean))
# goodRiver$kest_low <- ifelse(goodRiver$key == 'kobs', NA, goodRiver$kest_low)
# goodRiver$kest_high <- ifelse(goodRiver$key == 'kobs', NA, goodRiver$kest_high)
# riv <- as.character(goodRiver[1,]$river)
# goodRiverPlot <- ggplot(goodRiver, aes(x=time, y=value, color=key)) + #model
#   geom_ribbon(aes(ymin = kest_low, ymax = kest_high), alpha=0.75, fill='grey')+
#   geom_line(size=1.5) +
#   ylab('k600 [m/dy]') +
#   xlab('Timestep') +
#   scale_color_brewer(palette='Set2', name='k600 [m/dy]', labels=c('BIKER', 'Model using \nobserved hydraulics')) +
#   ggtitle(paste0('Good River: \n', riv))
# 
# # extract the legend from one of the plots
# legend <- get_legend(goodRiverPlot + theme(legend.box.margin = margin(0, 0, 0, 100)))
# 
# #bring it alllllll together via cowplot package
# timeseriesPlot <- plot_grid(goodRiverPlot + theme(legend.position = 'none'), ehRiverPlot, badRiverPlot, legend, ncol=2, labels=c('b','c','d',NA), label_size = 18)
#plot2 <- plot_grid(valPlot, timeseriesPlot, ncol=2, labels=c('a', NA), label_size = 18)
#ggsave('cache/validation/validation.jpg', valPlot, width=8, height=8)

####################
##SAVE ALL OTHER TIMESERIES FOR THE SUPPLEMENT
####################
files <- list.files('cache/validation/by_river', pattern="*.csv", full.names = FALSE)
file_paths <- list.files('cache/validation/by_river', pattern="*.csv", full.names = TRUE)

#Save all timeseries plots as a list of ggplot objects
pltList <- list()
for (i in 1:63){ #47 rivers with no measurement error, 16 with errors modeled via Frasson et al. 2021
  #get river name
  plot <- data.frame()
  file = files[i]
  river <- substr(file, 9, nchar(file)-4)

  #get results
  results <- read.csv(file_paths[i], header=TRUE, sep=",")
  results <- select(results, c('time', 'kobs', 'kest_mean', 'kest_low', 'kest_high')) %>%
    gather(key=key, value=value, c(kobs, kest_mean))
  results$kest_low <- ifelse(results$key == 'kobs', NA, results$kest_low/max(results$value))
  results$kest_high <- ifelse(results$key == 'kobs', NA, results$kest_high/max(results$value))
  riv <- as.character(results[1,]$river)

  pltList[[ river ]]  <- ggplot(results, aes(x=time, y=value/max(value), color=key)) + #model
    geom_ribbon(aes(ymin = kest_low, ymax = kest_high), alpha=0.75, fill='grey')+
    geom_line(size=1.5) +
    ylab('') +
    xlab('') +
    ylim(0,1.5)+
    scale_color_brewer(palette='Set2', name='', labels=c('BIKER', 'Model using \nobserved hydraulics')) +
    ggtitle(river)
  }

#gather and plot all ggplot objects (for no measurement error) as a single image via the cowplot package
plotgrid <- plot_grid(pltList$AshSlough + theme(legend.position='none'),
                      pltList$BerendaSlough + theme(legend.position='none'),
                      pltList$Brahmaputra + theme(legend.position='none'),
                      pltList$ChowchillaCanal + theme(legend.position='none'),
                      pltList$Connecticut + theme(legend.position='none'),
                      pltList$Cumberland + theme(legend.position='none'),
                      pltList$FresnoRiver + theme(legend.position='none'),
                      pltList$Ganges + theme(legend.position='none'),
                      pltList$GaronneDownstream + theme(legend.position='none'),
                      pltList$GaronneUpstream + theme(legend.position='none'),
                      pltList$GrantLineCanal + theme(legend.position='none'),
                      pltList$IowaRiver + theme(legend.position='none'),
                      pltList$Jamuna + theme(legend.position='none'),
                      pltList$Kanawha + theme(legend.position='none'),
                      pltList$Kushiyara + theme(legend.position='none'),
                      pltList$MariposaBypass + theme(legend.position='none'),
                      pltList$MercedRiver + theme(legend.position='none'),
                      pltList$MiddleRiver + theme(legend.position='none'),
                      pltList$MississippiDownstream + theme(legend.position='none'),
                      pltList$MississippiIntermediate + theme(legend.position='none'),
                      pltList$MississippiUpstream + theme(legend.position='none'),
                      pltList$MissouriDownstream + theme(legend.position='none'),
                      pltList$MissouriMidsection + theme(legend.position='none'),
                      pltList$MissouriUpstream + theme(legend.position='none'),
                      pltList$Ohio + theme(legend.position='none'),
                      pltList$OhioSection1 + theme(legend.position='none'),
                      pltList$OhioSection2 + theme(legend.position='none'),
                      pltList$OhioSection3 + theme(legend.position='none'),
                      pltList$OhioSection4 + theme(legend.position='none'),
                      pltList$OhioSection5 + theme(legend.position='none'),
                      pltList$OhioSection7 + theme(legend.position='none'),
                      pltList$OhioSection8 + theme(legend.position='none'),
                      pltList$Olentangy + theme(legend.position='none'),
                      pltList$Padma + theme(legend.position='none'),
                      pltList$Platte + theme(legend.position='none'),
                      pltList$Po + theme(legend.position='none'),
                      pltList$SacramentoDownstream + theme(legend.position='none'),
                      pltList$SacramentoUpstream + theme(legend.position='none'),
                      pltList$SanJoaquin + theme(legend.position='none'),
                      pltList$SanJoaquinRiver2 + theme(legend.position='none'),
                      pltList$Seine + theme(legend.position='none'),
                      pltList$SeineDownstream + theme(legend.position='none'),
                      pltList$SeineUpstream + theme(legend.position='none'),
                      pltList$Severn + theme(legend.position='none'),
                      pltList$StanislausRiver + theme(legend.position='none'),
                      pltList$TuolumneRiver + theme(legend.position='none'),
                      pltList$Wabash + theme(legend.position='none'),
                      ncol=4)

#grab legend from one of these plots
legend <- get_legend(
  # create some space to the left of the legend
  pltList$Wabash + theme(legend.text=element_text(size=18))
)

#draw legend in remaining empty space in figure
plotTimeseries <- plotgrid + draw_grob(legend, 0.8, -0.49, 1.1, 1.1)

#create x and y axis labels for entire figure
yTitleCombo <- textGrob(expression(k[600]/Max~k[600]~(m/dy)), gp=gpar(fontface="bold", col="black", fontsize=30), rot=90)
xTitleCombo <- textGrob(expression(Timestep), gp=gpar(fontface="bold", col="black", fontsize=30))
plotTimeseries <- gridExtra::grid.arrange(gridExtra::arrangeGrob(plotTimeseries, left = yTitleCombo, bottom = xTitleCombo))

#write to file
ggsave('cache/validation/timeseries_noerr.jpg', plotTimeseries, width = 18, height = 24)

####################
##SAVE RESULTS TO FILE---------------------------------------------------
###################
#results_3_2 <- data.frame('rmse'=rmse, 'r2'=r2)
#write.csv(results_3_2, 'cache/validation/results_all_riv.csv')
write.csv(stats_by_river, 'cache/validation/results_by_riv.csv')
