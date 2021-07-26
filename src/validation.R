########################
#Description: Produce figures for BIKER runs on SWOT-simulated rivers
#Creator: Craig Brinkerhoff
#Date: Summer 2021
######################

print('validating BIKER...')

`%notin%` <- Negate(`%in%`)

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
  xlab('k600 via observed \nshear velocity [m/dy]') +
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
  annotate("text", label = paste0('RMSE: ', rmse, ' m/dy'), x = 10^0, y = 10^1, size = 7, colour = "black")+
  annotate("text", label = paste0('r2: ', r2), x = 10^0, y = 10^0.8, size = 7, colour = "black")

###################
##CDFS OF BIKER VS OBSERVED---------------------
##################
t <- gather(full_output, key=key, value=value, c(kobs, kest_mean, kest_high, kest_low))
t$flag <- ifelse(t$key == 'kobs', 1, 0)
k_cdfs <- ggplot(t, aes(x=value, color=key, linetype=factor(flag))) +
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

#valPlot <- plot_grid(valPlot, k_cdfs, ncol=2, labels=c('a', 'b'), label_size=18)
#ggsave('cache/validation/validation.jpg', valPlot, width=12, height=7)

########################
##CALCULATE BY-RIVER ERROR METRICS---------------------------------------------
########################
stats_by_reach <- group_by(results, river, errFlag) %>%
  summarise(r2 = summary(lm(kest_mean~kobs))$r.squared,
            nrmse = sqrt(mean((kobs - kest_mean)^2, na.rm=T)) / mean(kobs, na.rm=T),
            rBIAS =   mean(kest_mean- kobs, na.rm=T) / mean(kobs, na.rm=T),
            rrmse =   sqrt(mean((kobs- kest_mean)^2 / kobs^2, na.rm=T)),
            meanKobs = mean(kobs, na.rm=T),
            meanWobs = mean(Wobs, na.rm=T),
            meanSobs = mean(Sobs, na.rm=T),
            meanDobs = mean(Dobs, na.rm=T),
            n_data=n())

plot_stats <- gather(stats_by_reach, key=key, value=value, c('nrmse', 'rBIAS', 'rrmse', 'r2'))
plot_stats <- filter(plot_stats, key %in% c('nrmse', 'rrmse', 'rBIAS', 'r2'))

########################
##PLOT BY-RIVER ERROR METRICS---------------------------------------
########################
plotSWOTreaches <- ggplot(plot_stats, aes(x=key, y=value, fill=factor(errFlag))) +
  geom_boxplot(size=1, alpha=0.75, outlier.shape=NA) +
  geom_point(size=3, position=position_jitterdodge(), aes(pch=factor(errFlag)),show.legend=FALSE)+
  geom_hline(yintercept=1, linetype='dashed') +
  geom_hline(yintercept=0, linetype='dashed') +
  geom_hline(yintercept=0.50, linetype='dashed') +
  xlab('Metric') +
  ylab('Value') +
#  coord_cartesian(ylim = c(-1.5,1))+
  scale_fill_brewer(palette = 'Accent', name='', labels=c('No Error', 'SWOT Measurement Error')) +
  theme(legend.position = "bottom",
        axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'))

############
##NRMSE VS AMOUNT OF DATA-----------------------
############
nrmse_vs_n <- ggplot(stats_by_reach, aes(y=nrmse, x=n_data, color=factor(errFlag))) +
  geom_point(size=5, aes(pch=factor(errFlag)),show.legend=FALSE) +
  geom_hline(yintercept=1, linetype='dashed') +
  geom_hline(yintercept=0, linetype='dashed') +
  geom_hline(yintercept=0.50, linetype='dashed') +
  ylab('NRMSE') +
  xlab('Num. SWOT Observations') +
  scale_color_brewer(palette = 'Accent', name='', labels=c('No Error', 'SWOT Measurement Error')) +
  coord_cartesian(ylim=c(-0.5, 2.3))+
  theme(legend.position = "bottom",
        axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'))

valPlot_byRiver <- plot_grid(plotSWOTreaches, nrmse_vs_n, ncol=2, labels=c('a', 'b'), label_size=18)
ggsave('cache/validation/validation_by_river.jpg', valPlot_byRiver, width=13, height=7)

#################
##PLOT EXAMPLE TIMESERIES PLOTS FOR 3 RIVERS------------------------------------------------------------------------------------
#################
set.seed(1141)

#Randomly select 3 rivers representative of good, eh, and bad KGE scores
nrmse_bins <- quantile(stats_by_reach[stats_by_reach$errFlag ==0,]$nrmse, c(0.33, 0.66), na.rm = T)
nrmse_bad <- filter(stats_by_reach[stats_by_reach$errFlag ==0,], nrmse >= nrmse_bins[2]) %>% select(river)
nrmse_good <- filter(stats_by_reach[stats_by_reach$errFlag ==0,], nrmse <= nrmse_bins[1]) %>% select(river)
nrmse_eh <- filter(stats_by_reach[stats_by_reach$errFlag ==0,], nrmse >= nrmse_bins[1] & nrmse <= nrmse_bins[2]) %>% select(river)

#bad river plot
badRiver <- filter(full_output, river == sample(nrmse_bad$river, 1)) %>%
  gather(key=key, value=value, c(kobs, kest_mean))
badRiver$kest_low <- ifelse(badRiver$key == 'kobs', NA, badRiver$kest_low)
badRiver$kest_high <- ifelse(badRiver$key == 'kobs', NA, badRiver$kest_high)
riv <- as.character(badRiver[1,]$river)
badRiverPlot <- ggplot(badRiver, aes(x=time, y=value, color=key)) + #model
  geom_ribbon(aes(ymin = kest_low, ymax = kest_high), alpha=0.75, fill='grey')+
  geom_line(size=1.5) +
  ylab('k600 [m/dy]') +
  xlab('Timestep') +
  scale_color_brewer(palette='Set2') +
  theme(legend.position = "none") +
  ggtitle(paste0('Bad River: ', riv))

#eh river plot
ehRiver <- filter(full_output, river == sample(nrmse_eh$river, 1)) %>%
  gather(key=key, value=value, c(kobs, kest_mean))
ehRiver$kest_low <- ifelse(ehRiver$key == 'kobs', NA, ehRiver$kest_low)
ehRiver$kest_high <- ifelse(ehRiver$key == 'kobs', NA, ehRiver$kest_high)
riv <- as.character(ehRiver[1,]$river)
ehRiverPlot <- ggplot(ehRiver, aes(x=time, y=value, color=key)) + #model
  geom_ribbon(aes(ymin = kest_low, ymax = kest_high), alpha=0.75, fill='grey')+
  geom_line(size=1.5) +
  ylab('k600 [m/dy]') +
  xlab('Timestep') +
  scale_color_brewer(palette='Set2') +
  theme(legend.position = "none") +
  ggtitle(paste0('OK River: ', riv))

#good river plot
goodRiver <- filter(full_output, river == sample(nrmse_good$river, 1)) %>%
  gather(key=key, value=value, c(kobs, kest_mean))
goodRiver$kest_low <- ifelse(goodRiver$key == 'kobs', NA, goodRiver$kest_low)
goodRiver$kest_high <- ifelse(goodRiver$key == 'kobs', NA, goodRiver$kest_high)
riv <- as.character(goodRiver[1,]$river)
goodRiverPlot <- ggplot(goodRiver, aes(x=time, y=value, color=key)) + #model
  geom_ribbon(aes(ymin = kest_low, ymax = kest_high), alpha=0.75, fill='grey')+
  geom_line(size=1.5) +
  ylab('k600 [m/dy]') +
  xlab('Timestep') +
  scale_color_brewer(palette='Set2', name='k600 [m/dy]', labels=c('BIKER', 'Model using \nobserved shear \nvelocity')) +
  ggtitle(paste0('Good River: \n', riv))

# extract the legend from one of the plots
legend <- get_legend(goodRiverPlot + theme(legend.box.margin = margin(0, 0, 0, 100)))

#bring it alllllll together via cowplot
timeseriesPlot <- plot_grid(goodRiverPlot + theme(legend.position = 'none'), ehRiverPlot, badRiverPlot, legend, ncol=2, labels=c('b','c','d',NA), label_size = 18)
plot2 <- plot_grid(valPlot, timeseriesPlot, ncol=2, labels=c('a', NA), label_size = 18)
ggsave('cache/validation/validation.jpg', plot2, width=14, height=8)

####################
##SAVE ALL OTHER TIMESERIES FOR THE SUPPLEMENT
####################
files <- list.files('cache/validation/by_river', pattern="*.csv", full.names = FALSE)
file_paths <- list.files('cache/validation/by_river', pattern="*.csv", full.names = TRUE)

pltList <- list()
for (i in 1:63){
  #get river name
  plot <- data.frame()
  file = files[i]
  river <- substr(file, 9, nchar(file)-4)

  #get results
  results <- read.csv(file_paths[i], header=TRUE, sep=",")
  results <- select(results, c('time', 'kobs', 'kest_mean', 'kest_low', 'kest_high')) %>%
    gather(key=key, value=value, c(kobs, kest_mean))
  results$kest_low <- ifelse(results$key == 'kobs', NA, results$kest_low)
  results$kest_high <- ifelse(results$key == 'kobs', NA, results$kest_high)
  riv <- as.character(results[1,]$river)

  pltList[[ river ]]  <- ggplot(results, aes(x=time, y=value, color=key)) + #model
    geom_ribbon(aes(ymin = kest_low, ymax = kest_high), alpha=0.75, fill='grey')+
    geom_line(size=1.5) +
    ylab('') +
    xlab('') +
    scale_color_brewer(palette='Set2', name='k600 [m/dy]', labels=c('BIKER', 'Model using \nobserved shear \nvelocity')) +
    ggtitle(river)
  }

#gather and plot all images for no measurement error
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

legend <- get_legend(
  # create some space to the left of the legend
  pltList$Wabash + theme(legend.text=element_text(size=18))
)

plotTimeseries <- plotgrid + draw_grob(legend, 0.8, -0.49, 1.1, 1.1)

yTitleCombo <- textGrob(expression(k600~(m/dy)), gp=gpar(fontface="bold", col="black", fontsize=22), rot=90)
xTitleCombo <- textGrob(expression(Timestep), gp=gpar(fontface="bold", col="black", fontsize=22))

plotTimeseries <- gridExtra::grid.arrange(gridExtra::arrangeGrob(plotTimeseries, left = yTitleCombo, bottom = xTitleCombo))
ggsave('cache/validation/timeseries_noerr.jpg', plotTimeseries, width = 18, height = 24)




####################
##SAVE RESULTS FOR MS---------------------------------------------------
###################
results_3_2 <- data.frame('rmse'=rmse, 'r2'=r2)
write.csv(results_3_2, 'cache/validation/results_all_riv.csv')
write.csv(stats_by_reach, 'cache/validation/results_by_riv.csv')

#########################
##NRMSE VS HYDRAULIC PROPERTIES FOR SUPPLEMENT-------------------------------------------------------
#########################
riverPropertiesPlot_k <- ggplot(stats_by_reach[stats_by_reach$errFlag==0,], aes(y=nrmse, x=meanKobs)) +
  geom_point(size=3) +
  geom_hline(yintercept = median(stats_by_reach[stats_by_reach$errFlag==0,]$nrmse), linetype='dashed')+
  xlab('Mean observed k [m/dy]') +
  ylab('NRMSE') +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'))
riverPropertiesPlot_W <- ggplot(stats_by_reach[stats_by_reach$errFlag==0,], aes(y=nrmse, x=meanWobs)) +
  geom_point(size=3) +
  geom_hline(yintercept = median(stats_by_reach[stats_by_reach$errFlag==0,]$nrmse), linetype='dashed')+
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
riverPropertiesPlot_S <- ggplot(stats_by_reach[stats_by_reach$errFlag==0,], aes(y=nrmse, x=meanSobs)) +
  geom_point(size=3) +
  geom_hline(yintercept = median(stats_by_reach[stats_by_reach$errFlag==0,]$nrmse), linetype='dashed')+
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

#bring it alllllll together via cowplot
riverPropertiesPlot <- plot_grid(riverPropertiesPlot_k + theme(legend.position = 'none'), riverPropertiesPlot_W, riverPropertiesPlot_S, ncol=2, label_size = 18)
ggsave('cache/validation/riverProperties_nrmse.jpg', riverPropertiesPlot, width=14, height=8)
