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
            CV = sd(kobs)/mean(kobs),
            meanKobs = mean(kobs, na.rm=T),
            meanWobs = mean(Wobs, na.rm=T),
            meanSobs = mean(Sobs, na.rm=T),
            meanDobs = mean(Dobs, na.rm=T),
            n_data=n(),
            priorKbias = (mean(kprior, na.rm=T)-mean(kobs, na.rm=T))/mean(kobs, na.rm=T),
            BIKERbias = (sum(kest_mean - kobs, na.rm=T)/n())/(mean(kobs, na.rm=T)))

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

###########################
## PLOT SCALING DYNAMICS-------------------------------
###########################
forPlot <- group_by(full_output, river) %>%
  do(model=lm(log(kest_mean)~log(kobs), data=.)) %>%
  mutate(Slope = model$coefficients[2]) %>%
  select('river', 'Slope')
write.csv(forPlot, 'cache/validation/results_dynamics.csv')

dynamicsPlot_stats <- ggplot(forPlot, aes(x=Slope)) +
  geom_histogram(color='black', size=1) +
  geom_vline(xintercept=1, size=3, color='darkblue', linetype='dotted') +
  ylab("Count")+
  xlab("Regression Slope") +
  annotate("text", x = 3, y = 12, label = 'Slope of 1 indicates\nperfectly captured dynamics', size=6, color='darkblue')+
  theme(legend.position = "none",
        axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'))

byRiverPlot <- plot_grid(plotSWOTreaches, dynamicsPlot_stats, ncol=2, labels=c('a', 'b'), label_size = 18)
ggsave('cache/validation/validation_by_river.jpg', byRiverPlot, width=13.5, height=7)

############
##r VS AMOUNT OF DATA-----------------------
############
r_vs_n <- ggplot(stats_by_river[stats_by_river$errFlag == 0,]) +
  geom_point(aes(y=r, x=n_data), color='#7fc97f',size=6) +
  geom_hline(aes(yintercept = median(stats_by_river[stats_by_river$errFlag==0,]$r)), linetype='dashed', color='darkblue', size=1.5)+
  ylab('r') +
  xlab('') +
  theme(legend.position = "bottom",
        axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'))

nBiasPlot <- ggplot(stats_by_river[stats_by_river$errFlag==0,], aes(y=BIKERbias, x=priorKbias)) +
  geom_point(,color='#7fc97f', size=6) +
  geom_abline(size=2, linetype='dashed', color='darkgrey')+
  geom_smooth(method='lm', se=F, color='#7fc97f')+
  xlab('Prior bias') +
  ylab('BIKER bias') +
  theme(legend.position = "none",
        axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'))

hydraulicsPlot <- plot_grid(r_vs_n, nBiasPlot, ncol=1, labels=c('auto'), label_size=18)
ggsave('cache/validation/validation_hydraulics.jpg', hydraulicsPlot, width=12, height=12)

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
