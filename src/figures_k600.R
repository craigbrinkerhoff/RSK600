########################
#Description: Produce figures for BIKER runs on SWOT-simulated rivers
#Creator: Craig Brinkerhoff
#Date: Winter 2022
######################

print('validating BIKER...')

`%notin%` <- Negate(`%in%`)

theme_set(theme_classic())

################
##READ IN RESULTS----------------------------------
###############
results <- read.csv('cache/validation/BIKER_validation_results.csv')
full_output <- filter(results, errFlag == 0) #remove results with SWOT measurement error

########################
##CALCULATE BY-RIVER ERROR METRICS---------------------------------------------
########################
stats_by_river <- group_by(results, river, errFlag) %>%
  summarise(r = sqrt(summary(lm(kest_mean~kobs))$r.squared),
            NRMSE = sqrt(mean((kobs - kest_mean)^2, na.rm=T)) / mean(kobs, na.rm=T),
            rBIAS =   mean(kest_mean- kobs, na.rm=T) / mean(kobs, na.rm=T),
            KGE = KGE(kest_mean, kobs),
            obsCV = sd(kobs)/mean(kobs),
            meanKobs = mean(kobs, na.rm=T),
            meanWobs = mean(Wobs, na.rm=T),
            meanSobs = mean(Sobs, na.rm=T),
            meanDobs = mean(Dobs, na.rm=T),
            n_data=n(),
            priorCV = sd(kprior)/mean(kprior),
            posteriorCV = sd(kest_mean)/mean(kest_mean))

plot_stats <- gather(stats_by_river, key=key, value=value, c('NRMSE', 'rBIAS', 'KGE', 'r'))
plot_stats <- filter(plot_stats, key %in% c('NRMSE', 'KGE', 'rBIAS', 'r'))

########################
##PLOT BY-RIVER ERROR METRICS---------------------------------------
########################
# lowKGE_0 <- nrow(plot_stats[plot_stats$key == 'KGE' & plot_stats$value < -2 & plot_stats$errFlag==0,])
# lowKGE_1 <- nrow(plot_stats[plot_stats$key == 'KGE' & plot_stats$value < -2 & plot_stats$errFlag==1,])
# plotSWOTreaches <- ggplot(plot_stats, aes(x=key, y=value, fill=factor(errFlag))) +
#   geom_boxplot(size=1, alpha=0.75) +
#   geom_hline(yintercept=1, linetype='dashed') +
#   geom_hline(yintercept=0, linetype='dashed') +
#   geom_hline(yintercept=0.50, linetype='dashed') +
#   xlab('Metric') +
#   ylab('Value') +
#   coord_cartesian(ylim=c(-2,2))+
#   annotate("text", x = 'KGE', y = -2, label = paste0('# rivers \n< -2: ', lowKGE_0), size=5, color='#7fc97f')+
#   annotate("text", x = 'KGE', y = -1.6, label = paste0('# rivers \n< -2: ', lowKGE_1), size=5, color='#beaed4')+
#   scale_fill_brewer(palette = 'Accent', name='', labels=c('No Error', 'SWOT Measurement Error')) +
#   theme(legend.position = "bottom",
#         axis.text=element_text(size=20),
#         axis.title=element_text(size=24,face="bold"),
#         legend.text = element_text(size=17),
#         legend.title = element_text(size=17, face='bold'))

#do it with CDFs...
kge_cdf <- ggplot(stats_by_river[stats_by_river$errFlag==0,], aes(x=KGE)) +
  stat_ecdf(size=3, color='#1c9099') +
  geom_hline(yintercept = 0.5, linetype='dashed', size=1) +
  xlab('')+
  ylab('Percentile')+
  ggtitle('KGE')+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        plot.title = element_text(size = 30, face = "bold"))

r_cdf <- ggplot(stats_by_river[stats_by_river$errFlag==0,], aes(x=r)) +
  stat_ecdf(size=3, color='#1c9099') +
  geom_hline(yintercept = 0.5, linetype='dashed', size=1) +
  xlab('Value')+
  ylab('Percentile')+
  ggtitle('r')+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        plot.title = element_text(size = 30, face = "bold"))

nrmse_cdf <- ggplot(stats_by_river[stats_by_river$errFlag==0,], aes(x=NRMSE)) +
  stat_ecdf(size=3, color='#1c9099') +
  geom_hline(yintercept = 0.5, linetype='dashed', size=1) +
  xlab('')+
  ylab('')+
  ggtitle('NRMSE')+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        plot.title = element_text(size = 30, face = "bold"))

rbias_cdf <- ggplot(stats_by_river[stats_by_river$errFlag==0,], aes(x=rBIAS)) +
  stat_ecdf(size=3, color='#1c9099') +
  geom_hline(yintercept = 0.5, linetype='dashed', size=1) +
  xlab('Value')+
  ylab('')+
  ggtitle('rBIAS')+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        plot.title = element_text(size = 30, face = "bold"))

statsPlot <- plot_grid(kge_cdf, nrmse_cdf, r_cdf, rbias_cdf, ncol=2, labels = NA)
ggsave('cache/validation/validation_by_river.jpg', statsPlot, width=10, height=10)

########################
##ERROR/NO ERROR COMPARISON------------------------------------
#########################
forError <- filter(stats_by_river, errFlag==1)
forError <- stats_by_river[stats_by_river$river %in% forError$river,]
forError <- select(forError, c('river', 'errFlag', 'NRMSE', 'rBIAS', 'KGE', 'r'))

#KGE
forError_KGE <- pivot_wider(forError, names_from=errFlag, values_from = KGE) %>%
  group_by(river) %>%
  summarise(KGE_no_err = sum(`0`, na.rm=T),
            KGE_err = sum(`1`, na.rm=T))

kge_scatter <- ggplot(forError_KGE, aes(x=KGE_no_err, y=KGE_err)) +
  geom_point(size=8, color='#2ca25f') +
  geom_smooth(method='lm', se=F, size=3, color='black')+
  geom_abline(linetype='dashed', color='darkgrey', size=2) +
  xlim(-5,1)+
  ylim(-5,1) +
  xlab('')+
  ylab('Error')+
  ggtitle('KGE')+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        plot.title = element_text(size = 30, face = "bold"))

#r
forError_r <- pivot_wider(forError, names_from=errFlag, values_from = r) %>%
  group_by(river) %>%
  summarise(r_no_err = sum(`0`, na.rm=T),
            r_err = sum(`1`, na.rm=T))

r_scatter <- ggplot(forError_r, aes(x=r_no_err, y=r_err)) +
  geom_point(size=8, color='#2ca25f') +
  geom_smooth(method='lm', se=F, size=3, color='black')+
  geom_abline(linetype='dashed', color='darkgrey', size=2) +
  xlim(0,1)+
  ylim(0,1) +
  xlab('No Error')+
  ylab('Error')+
  ggtitle('r')+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        plot.title = element_text(size = 30, face = "bold"))

#NRMSE
forError_NRMSE <- pivot_wider(forError, names_from=errFlag, values_from = NRMSE) %>%
  group_by(river) %>%
  summarise(NRMSE_no_err = sum(`0`, na.rm=T),
            NRMSE_err = sum(`1`, na.rm=T))

NRMSE_scatter <- ggplot(forError_NRMSE, aes(x=NRMSE_no_err, y=NRMSE_err)) +
  geom_point(size=8, color='#2ca25f') +
  geom_smooth(method='lm', se=F, size=3, color='black')+
  geom_abline(linetype='dashed', color='darkgrey', size=2) +
  xlim(0,2)+
  ylim(0,2) +
  xlab('')+
  ylab('')+
  ggtitle('NRMSE')+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'),
        plot.title = element_text(size = 30, face = "bold"))

#rBIAS
forError_rbias <- pivot_wider(forError, names_from=errFlag, values_from = rBIAS) %>%
  group_by(river) %>%
  summarise(rbias_no_err = sum(`0`, na.rm=T),
            rbias_err = sum(`1`, na.rm=T))

rbias_scatter <- ggplot(forError_rbias, aes(x=rbias_no_err, y=rbias_err)) +
  geom_point(size=8, color='#2ca25f') +
  geom_smooth(method='lm', se=F, size=3, color='black')+
  geom_abline(linetype='dashed', color='darkgrey', size=2) +
  xlim(-1,2)+
  ylim(-1,2) +
  xlab('No Error')+
  ylab('')+
  ggtitle('rBIAS')+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'),
        plot.title = element_text(size = 30, face = "bold"))

errorPlot <- plot_grid(kge_scatter, NRMSE_scatter, r_scatter, rbias_scatter, ncol=2, labels=NA)
ggsave('cache/validation/errorPlot.jpg', errorPlot, width=10, height=10)

###########################
## PLOT SCALING DYNAMICS-------------------------------
###########################
forPlot <- group_by(full_output, river) %>%
  do(model=lm((kest_mean)~(kobs), data=.)) %>%
  mutate(Slope = model$coefficients[2]) %>%
  select('river', 'Slope')
write.csv(forPlot, 'cache/validation/results_dynamics.csv')

bw <- 0.3
dynamicsPlot_stats <- ggplot(forPlot) +
  geom_histogram(aes(x=Slope),color='black', size=1, binwidth=bw) +
  geom_vline(xintercept=1, size=3, color='darkblue', linetype='dotted') +
  ylab("Count")+
  xlab("Regression Slope") +
  annotate("text", x = 4, y = 10, label = 'Slope of 1: correctly\ninferred temporal dynamics', size=6, color='darkblue')+
  theme(legend.position = "none",
        axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'))

#################
## PLOT TEMPORAL DYNAMICS---------------------
#################
forPlot <- gather(stats_by_river[stats_by_river$errFlag == 0,], key=key, value=value, c('obsCV', 'priorCV', 'posteriorCV'))
forPlot$id <- ifelse(forPlot$key == 'obsCV', 'obs', 'notobs')
CVPlot <- ggplot(forPlot, aes(x=value, color=key)) +
  stat_ecdf(size=2) +
  scale_color_manual(values=c('#756bb1', '#006d2c', '#74c476'), name='', labels=c('Observed k', 'Posterior k', 'Prior k')) +
  coord_cartesian(xlim = c(0,1)) +
  ylab('Density') +
  xlab('Coefficient of variation')+
  theme(legend.position=c(.75, .3),
        axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'))

temporalPlot <- plot_grid(dynamicsPlot_stats, CVPlot, labels='auto', label_size=18)
ggsave('cache/validation/validation_temporal.jpg', temporalPlot, width=12, height=7)
break
####################
##SAVE ALL TIMESERIES FOR THE SUPPLEMENT
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
    ylim(0,1.7)+
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
plotTimeseries <- plotgrid + draw_grob(legend, 0.65, -0.2, 0.5, 0.5)

#create x and y axis labels for entire figure
yTitleCombo <- textGrob(expression(k[600]/Max~k[600]~(m/dy)), gp=gpar(fontface="bold", col="black", fontsize=30), rot=90)
xTitleCombo <- textGrob(expression(Timestep), gp=gpar(fontface="bold", col="black", fontsize=30))
plotTimeseries <- gridExtra::grid.arrange(gridExtra::arrangeGrob(plotTimeseries, left = yTitleCombo, bottom = xTitleCombo))

#write to file
ggsave('cache/validation/timeseries_noerr.jpg', plotTimeseries, width = 18, height = 24)

################
##TIMERSERIES PLOT WITH SWOT ERRORS-----------------------------------
#################
#gather and plot all ggplot objects (for no measurement error) as a single image via the cowplot package
plotgrid_err <- plot_grid(pltList$MissouriDownstream_err + theme(legend.position='none'),
                      pltList$MissouriUpstream_err + theme(legend.position='none'),
                      pltList$OhioSection8_err+ theme(legend.position='none'),
                      pltList$OhioSection7_err+ theme(legend.position='none'),
                      pltList$OhioSection5_err+ theme(legend.position='none'),
                      pltList$OhioSection4_err+ theme(legend.position='none'),
                      pltList$OhioSection3_err+ theme(legend.position='none'),
                      pltList$OhioSection2_err+ theme(legend.position='none'),
                      pltList$OhioSection1_err+ theme(legend.position='none'),
                      pltList$SeineDownstream_err+ theme(legend.position='none'),
                      pltList$Brahmaputra_err+ theme(legend.position='none'),
                      pltList$IowaRiver_err+ theme(legend.position='none'),
                      pltList$SeineUpstream_err+ theme(legend.position='none'),
                      pltList$Kushiyara_err+ theme(legend.position='none'),
                      pltList$Jamuna_err+ theme(legend.position='none'),
                      pltList$Padma_err+ theme(legend.position='none'),
                      ncol=3)
                      
#draw legend in remaining empty space in figure
plotTimeseries <- plotgrid_err + draw_grob(legend, 0.6, -0.15, 0.5, 0.5)
plotTimeseries <- gridExtra::grid.arrange(gridExtra::arrangeGrob(plotTimeseries, left = yTitleCombo, bottom = xTitleCombo))

#write to file
ggsave('cache/validation/timeseries_err.jpg', plotTimeseries, width = 13, height = 13)

#################
##GRAB REPRESENTIVE TIMESERIES FOR THE MAIN TEXT-----------------------
##################
set.seed(165)
quantiles <- quantile(stats_by_river[stats_by_river$errFlag == 0,]$KGE, c(0.33, 0.66))

stats_by_river$river <- as.character(stats_by_river$river)

temp <- stats_by_river[stats_by_river$errFlag==0 & stats_by_river$KGE < quantiles[1],]
first_tertile <- temp[sample(nrow(temp), 2), ]$river

temp <- stats_by_river[stats_by_river$errFlag==0 & stats_by_river$KGE >= quantiles[1] & stats_by_river$KGE < quantiles[2],]
second_tertile <- temp[sample(nrow(temp), 2), ]$river

temp <- stats_by_river[stats_by_river$errFlag==0 & stats_by_river$KGE >= quantiles[2],]
third_tertile <- temp[sample(nrow(temp), 2), ]$river

plotGrid_short <- plot_grid(
  pltList[[first_tertile[1]]]+ ggtitle('Sacramento Upstream')+ theme(legend.position='none',
                                     axis.text=element_text(size=20),
                                     axis.title=element_text(size=24,face="bold"),
                                     legend.text = element_text(size=17),
                                     legend.title = element_text(size=17, face='bold')),
  pltList[[first_tertile[2]]]+ ggtitle('Ohio')+ theme(legend.position='none',
                                      axis.text=element_text(size=20),
                                      axis.title=element_text(size=24,face="bold"),
                                      legend.text = element_text(size=17),
                                      legend.title = element_text(size=17, face='bold')),
  textGrob('Low KGE', y=0.6, gp=gpar(fontface="bold", col="black", fontsize=30)),
  pltList[[second_tertile[1]]]+ ggtitle('Po') + theme(legend.position='none',
                                     axis.text=element_text(size=20),
                                     axis.title=element_text(size=24,face="bold"),
                                     legend.text = element_text(size=17),
                                     legend.title = element_text(size=17, face='bold')),
  pltList[[second_tertile[2]]]+ggtitle('Chowchilla Canal')+ theme(legend.position='none',
                                      axis.text=element_text(size=20),
                                      axis.title=element_text(size=24,face="bold"),
                                      legend.text = element_text(size=17),
                                      legend.title = element_text(size=17, face='bold')),
  textGrob('Middle KGE', y=0.6, gp=gpar(fontface="bold", col="black", fontsize=30)),
  pltList[[third_tertile[1]]]+ ggtitle('Tuolumne River')+theme(legend.position='none',
                                        axis.text=element_text(size=20),
                                        axis.title=element_text(size=24,face="bold"),
                                        legend.text = element_text(size=17),
                                        legend.title = element_text(size=17, face='bold')),
  pltList[[third_tertile[2]]]+ggtitle('Sacramento Downstream')+theme(legend.position='none',
                                        axis.text=element_text(size=20),
                                        axis.title=element_text(size=24,face="bold"),
                                        legend.text = element_text(size=17),
                                        legend.title = element_text(size=17, face='bold')),
  textGrob('High KGE', y=0.6, gp=gpar(fontface="bold", col="black", fontsize=30)),
  ncol=3,
  label_size = 18,
  labels=c('a', 'b', NA, 'c', 'd', NA, 'e', 'f', NA))

#draw legend in remaining empty space in figure
plotTimeseries <- plotGrid_short + draw_grob(legend, 0.6, -0.2, 0.5, 0.5)
plotTimeseries <- gridExtra::grid.arrange(gridExtra::arrangeGrob(plotTimeseries, left = yTitleCombo, bottom = xTitleCombo))

#write to file
ggsave('cache/validation/timeseries_noerr_short.jpg', plotTimeseries, width = 13, height = 10)

####################
##SAVE RESULTS TO FILE---------------------------------------------------
###################
write.csv(stats_by_river, 'cache/validation/results_by_riv.csv')
