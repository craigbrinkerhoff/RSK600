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

plotSWOTreaches_cdfs_kge <- ggplot(plot_stats[plot_stats$key=='KGE',], aes(x=value, color=factor(errFlag))) +
  stat_ecdf(size=2) +
  xlab('KGE') +
  ylab('Density') +
  scale_color_brewer(palette = 'Accent', name='', labels=c('No Error', 'SWOT Measurement Error')) +
  theme(legend.position='none',
        axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'))

plotSWOTreaches_cdfs_r <- ggplot(plot_stats[plot_stats$key=='r',], aes(x=value, color=factor(errFlag))) +
  stat_ecdf(size=2) +
  xlab('r') +
  ylab('Density') +
  scale_color_brewer(palette = 'Accent', name='', labels=c('No Error', 'SWOT Measurement Error')) +
  theme(legend.position='none',
        axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'))

temp <- plot_grid(plotSWOTreaches_cdfs_kge, plotSWOTreaches_cdfs_r, ncol=1, labels=c('b', 'c'), label_size = 18)
plotOverall <- plot_grid(plotSWOTreaches, temp, labels=c('a', NA), label_size = 18, ncol=2)
ggsave('cache/validation/validation_by_river.jpg', plotOverall, width=12, height=7)

###########################
## PLOT SCALING DYNAMICS-------------------------------
###########################
forPlot <- group_by(full_output, river) %>%
  do(model=lm(log(kest_mean)~log(kobs), data=.)) %>%
  mutate(Slope = model$coefficients[2]) %>%
  select('river', 'Slope')
write.csv(forPlot, 'cache/validation/results_dynamics.csv')

bw <- 0.5
dynamicsPlot_stats <- ggplot(forPlot) +
  geom_histogram(aes(x=Slope),color='black', size=1, binwidth=bw) +
  geom_vline(xintercept=1, size=3, color='darkblue', linetype='dotted') +
  ylab("Count")+
  xlab("Regression Slope") +
  annotate("text", x = 3, y = 10, label = 'Slope of 1: correctly\ninferred temporal dynamics', size=6, color='darkblue')+
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
  xlab('Coefficient of variation (by river)')+
  theme(legend.position=c(.75, .3),
        axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'))

temporalPlot <- plot_grid(dynamicsPlot_stats, CVPlot, labels='auto', label_size=18)
ggsave('cache/validation/validation_temporal.jpg', temporalPlot, width=12, height=7)
break
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
write.csv(stats_by_river, 'cache/validation/results_by_riv.csv')
