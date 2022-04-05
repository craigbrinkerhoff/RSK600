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
full_output <- filter(results, errFlag == 0) #remove results with SWOT measurement error (for now)

########################
##CALCULATE BY-RIVER ERROR METRICS---------------------------------------------
########################
stats_by_river <- group_by(results, river, errFlag) %>%
  summarise(r = cor(kobs, kest_mean, method='pearson'),
            NRMSE = sqrt(mean((kest_mean - kobs)^2, na.rm=T)) / mean(kobs, na.rm=T),
            KGE = KGE(kest_mean, kobs),
            NMAE = mean(abs(kest_mean-kobs)/mean(kobs, na.rm=T), na.rm=T),
            NMAE_prior = mean(abs(kprior - kobs), na.rm=T)/mean(kobs, na.rm=T),
            meanKobs = mean(kobs, na.rm=T),
            meanWobs = mean(Wobs, na.rm=T),
            meanSobs = mean(Sobs, na.rm=T),
            meanDobs = mean(Dobs, na.rm=T),
            n_data=n(),
            cvObs = sd(kobs, na.rm=T)/mean(kobs))

plot_stats <- gather(stats_by_river, key=key, value=value, c('NRMSE', 'NMAE', 'KGE', 'r'))
plot_stats <- filter(plot_stats, key %in% c('NRMSE', 'KGE', 'NMAE', 'r'))

########################
##PLOT BY-RIVER ERROR METRICS w/ eCDFs---------------------------------------
########################
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

nmae_cdf <- ggplot(stats_by_river[stats_by_river$errFlag==0,], aes(x=NMAE)) +
  stat_ecdf(size=3, color='#1c9099') +
  scale_x_log10()+
  geom_hline(yintercept = 0.5, linetype='dashed', size=1) +
  xlab('Value')+
  ylab('')+
  ggtitle('NMAE')+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        plot.title = element_text(size = 30, face = "bold"))

statsPlot <- plot_grid(kge_cdf, nrmse_cdf, r_cdf, nmae_cdf, ncol=2, labels = NA)
ggsave('cache/validation/fig4.jpg', statsPlot, width=10, height=10)

########################
##ERROR/NO ERROR COMPARISON------------------------------------
#########################
forError <- filter(stats_by_river, errFlag==1)
forError <- stats_by_river[stats_by_river$river %in% forError$river,]
forError <- select(forError, c('river', 'errFlag', 'NRMSE', 'NMAE', 'KGE', 'r'))

#KGE
forError_KGE <- pivot_wider(forError, names_from=errFlag, values_from = KGE) %>%
  group_by(river) %>%
  summarise(KGE_no_err = sum(`0`, na.rm=T),
            KGE_err = sum(`1`, na.rm=T))

kge_scatter <- ggplot(forError_KGE, aes(x=KGE_no_err, y=KGE_err)) +
  geom_polygon(data = data.frame(x = c(-5, 1, -5), y = c(-5,1,1)), aes(x = x, y = y), fill = "#8dd3c7" )+  # better triangle
  geom_polygon(data = data.frame(x = c(-5, 1, 1), y = c(-5,1,-5)), aes(x = x, y = y), fill = "#bebada" )+  # worse triangle
  geom_point(size=8, color='#386cb0') +
  geom_smooth(method='lm', se=F, color='black', size=3)+
  #geom_abline(linetype='dashed', color='darkgrey', size=2) +  xlim(-5,1)+
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
  geom_polygon(data = data.frame(x = c(-1, 1, -1), y = c(-1,1,1)), aes(x = x, y = y), fill = "#8dd3c7" )+  # better triangle
  geom_polygon(data = data.frame(x = c(-1, 1, 1), y = c(-1,1,-1)), aes(x = x, y = y), fill = "#bebada" )+  # worse triangle
  geom_point(size=8, color='#386cb0') +
  geom_smooth(method='lm', se=F, color='black', size=3)+
  #geom_abline(linetype='dashed', color='darkgrey', size=2) +
  xlim(-1,1)+
  ylim(-1,1) +
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
  geom_polygon(data = data.frame(x = c(2, 0, 2), y = c(2,0,0)), aes(x = x, y = y), fill = "#8dd3c7" )+  # better triangle
  geom_polygon(data = data.frame(x = c(2, 0, 0), y = c(2,0,2)), aes(x = x, y = y), fill = "#bebada" )+  # worse triangle
  geom_point(size=8, color='#386cb0') +
  geom_smooth(method='lm', se=F, color='black', size=3)+
  #geom_abline(linetype='dashed', color='darkgrey', size=2) +  xlim(2,0)+
  ylim(2,0) +
  xlim(2,0)+
  xlab('')+
  ylab('')+
  ggtitle('NRMSE')+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'),
        plot.title = element_text(size = 30, face = "bold"))

#NMAE
forError_nmae <- pivot_wider(forError, names_from=errFlag, values_from = NMAE) %>%
  group_by(river) %>%
  summarise(nmae_no_err = sum(`0`, na.rm=T),
            nmae_err = sum(`1`, na.rm=T))

nmae_scatter <- ggplot(forError_nmae, aes(x=abs(nmae_no_err), y=abs(nmae_err))) +
  geom_polygon(data = data.frame(x = c(1, 0, 1), y = c(1,0,0)), aes(x = x, y = y), fill = "#8dd3c7" )+  # better triangle
  geom_polygon(data = data.frame(x = c(1, 0, 0), y = c(1,0,1)), aes(x = x, y = y), fill = "#bebada" )+  # worse triangle
  geom_point(size=8, color='#386cb0') +
  geom_smooth(method='lm', se=F, color='black', size=3)+
  xlim(1,0)+
  ylim(1,0) +
  xlab('No Error')+
  ylab('')+
  ggtitle('NMAE')+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'),
        plot.title = element_text(size = 30, face = "bold"))

fig6 <- plot_grid(kge_scatter, NRMSE_scatter, r_scatter, nmae_scatter, ncol=2, labels=NA)
ggsave('cache/validation/fig5.jpg', fig6, width=10, height=10)

###########################
## PLOT PRIOR/POSTEROR COMPARISON-------------------------------
###########################
#Prior/Posterior bias
forPlot <- stats_by_river[stats_by_river$errFlag == 0,]
forPlot2 <- gather(forPlot, key=key, value=value, c('NMAE', 'NMAE_prior'))

priorPosterPlot_all <- ggplot(forPlot2, aes(x=value, color=key)) +
  stat_ecdf(size=3) +
  scale_color_brewer(palette = 'Accent', name='')+
  xlab('')+
  ylab('Percentile')+
  xlim(0,2.2)+
  annotate("text", x = 1.25, y = 0.5, label = 'All Rivers\n(n = 47 rivers)', size=8, color='darkblue')+
  theme(axis.text=element_text(size=22),
        axis.title=element_text(size=26,face="bold"),
        plot.title = element_text(size = 30, face = "bold"),
        legend.position = 'none')

priorPosterPlot_2 <- ggplot(forPlot2[forPlot2$cvObs >= 0.10,], aes(x=value, color=key)) +
  stat_ecdf(size=3) +
  scale_color_brewer(palette = 'Accent', name='')+
  xlab('')+
  ylab('')+
  xlim(0,2.2)+
  annotate("text", x = 1.25, y = 0.5, label = 'Observed CV > 0.10\n(n = 38 rivers)', size=8, color='darkblue')+
  theme(axis.text=element_text(size=22),
        axis.title=element_text(size=26,face="bold"),
        plot.title = element_text(size = 30, face = "bold"),
        legend.position = 'none')

priorPosterPlot_3 <- ggplot(forPlot2[forPlot2$cvObs >= 0.20,], aes(x=value, color=key)) +
  stat_ecdf(size=3) +
  scale_color_brewer(palette = 'Accent', name='')+
  xlab('NMAE')+
  ylab('Percentile')+
  xlim(0,2.2)+
  annotate("text", x = 1.25, y = 0.5, label = 'Observed CV > 0.20\n(n = 18 rivers)', size=8, color='darkblue')+
  theme(axis.text=element_text(size=22),
        axis.title=element_text(size=26,face="bold"),
        plot.title = element_text(size = 30, face = "bold"),
        legend.position = 'none')

priorPosterPlot_4 <- ggplot(forPlot2[forPlot2$cvObs >= 0.30,], aes(x=value, color=key)) +
  stat_ecdf(size=3) +
  scale_color_brewer(palette = 'Accent', name='', labels=c('Posterior', 'Prior'))+
  xlab('NMAE')+
  ylab('')+
  xlim(0,2.2)+
  annotate("text", x = 1.25, y = 0.5, label = 'Observed CV > 0.30\n(n = 10 rivers)', size=8, color='darkblue')+
  theme(axis.text=element_text(size=22),
        axis.title=element_text(size=26,face="bold"),
        legend.text = element_text(size=28),
        legend.title = element_text(size=28, face='bold'),
        plot.title = element_text(size = 30, face = "bold"),
        legend.position = c(0.7, 0.8))


fig7 <- plot_grid(priorPosterPlot_all, priorPosterPlot_2, priorPosterPlot_3, priorPosterPlot_4, labels='auto', label_size = 18, ncol = 2)
ggsave('cache/validation/fig6.jpg', fig7, width=12, height=12)

####################
##SAVE ALL TIMESERIES IN A SINGLE FINGURE FOR THE SUPPLEMENT
####################
files <- list.files('cache/validation/by_river', pattern=".csv", full.names = FALSE)
file_paths <- list.files('cache/validation/by_river', pattern=".csv", full.names = TRUE)

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
  
  #make ggplot object
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
plotgrid <- plot_grid(pltList$AshSlough + theme(legend.position='none',
                                                plot.title = element_text(size = 25, face = "bold"),
                                                axis.text=element_text(size=17),
                                                axis.title=element_text(size=17,face="bold")),
                      pltList$BerendaSlough + theme(legend.position='none',
                                                    plot.title = element_text(size = 25, face = "bold"),
                                                    axis.text=element_text(size=17),
                                                    axis.title=element_text(size=17,face="bold")),
                      pltList$Brahmaputra + theme(legend.position='none',
                                                  plot.title = element_text(size = 25, face = "bold"),
                                                  axis.text=element_text(size=17),
                                                  axis.title=element_text(size=17,face="bold")),
                      pltList$ChowchillaCanal + theme(legend.position='none',
                                                      plot.title = element_text(size = 25, face = "bold"),
                                                      axis.text=element_text(size=17),
                                                      axis.title=element_text(size=17,face="bold")),
                      pltList$Connecticut + theme(legend.position='none',
                                                  plot.title = element_text(size = 25, face = "bold"),
                                                  axis.text=element_text(size=17),
                                                  axis.title=element_text(size=17,face="bold")),
                      pltList$Cumberland + theme(legend.position='none',
                                                 plot.title = element_text(size = 25, face = "bold"),
                                                 axis.text=element_text(size=17),
                                                 axis.title=element_text(size=17,face="bold")),
                      pltList$FresnoRiver + theme(legend.position='none',
                                                  plot.title = element_text(size = 25, face = "bold"),
                                                  axis.text=element_text(size=17),
                                                  axis.title=element_text(size=17,face="bold")),
                      pltList$Ganges + theme(legend.position='none',
                                             plot.title = element_text(size = 25, face = "bold"),
                                             axis.text=element_text(size=17),
                                             axis.title=element_text(size=17,face="bold")),
                      pltList$GaronneDownstream + theme(legend.position='none',
                                                        plot.title = element_text(size = 25, face = "bold"),
                                                        axis.text=element_text(size=17),
                                                        axis.title=element_text(size=17,face="bold")),
                      pltList$GaronneUpstream + theme(legend.position='none',
                                                      plot.title = element_text(size = 25, face = "bold"),
                                                      axis.text=element_text(size=17),
                                                      axis.title=element_text(size=17,face="bold")),
                      pltList$GrantLineCanal + theme(legend.position='none',
                                                     plot.title = element_text(size = 25, face = "bold"),
                                                     axis.text=element_text(size=17),
                                                     axis.title=element_text(size=17,face="bold")),
                      pltList$IowaRiver + theme(legend.position='none',
                                                plot.title = element_text(size = 25, face = "bold"),
                                                axis.text=element_text(size=17),
                                                axis.title=element_text(size=17,face="bold")),
                      pltList$Jamuna + theme(legend.position='none',
                                             plot.title = element_text(size = 25, face = "bold"),
                                             axis.text=element_text(size=17),
                                             axis.title=element_text(size=17,face="bold")),
                      pltList$Kanawha + theme(legend.position='none',
                                              plot.title = element_text(size = 25, face = "bold"),
                                              axis.text=element_text(size=17),
                                              axis.title=element_text(size=17,face="bold")),
                      pltList$Kushiyara + theme(legend.position='none',
                                                plot.title = element_text(size = 25, face = "bold"),
                                                axis.text=element_text(size=17),
                                                axis.title=element_text(size=17,face="bold")),
                      pltList$MariposaBypass + theme(legend.position='none',
                                                     plot.title = element_text(size = 25, face = "bold"),
                                                     axis.text=element_text(size=17),
                                                     axis.title=element_text(size=17,face="bold")),
                      pltList$MercedRiver + theme(legend.position='none',
                                                  plot.title = element_text(size = 25, face = "bold"),
                                                  axis.text=element_text(size=17),
                                                  axis.title=element_text(size=17,face="bold")),
                      pltList$MiddleRiver + theme(legend.position='none',
                                                  plot.title = element_text(size = 25, face = "bold"),
                                                  axis.text=element_text(size=17),
                                                  axis.title=element_text(size=17,face="bold")),
                      pltList$MississippiDownstream + theme(legend.position='none',
                                                            plot.title = element_text(size = 25, face = "bold"),
                                                            axis.text=element_text(size=17),
                                                            axis.title=element_text(size=17,face="bold")),
                      pltList$MississippiIntermediate + theme(legend.position='none',
                                                              plot.title = element_text(size = 25, face = "bold"),
                                                              axis.text=element_text(size=17),
                                                              axis.title=element_text(size=17,face="bold")),
                      pltList$MississippiUpstream + theme(legend.position='none',
                                                          plot.title = element_text(size = 25, face = "bold"),
                                                          axis.text=element_text(size=17),
                                                          axis.title=element_text(size=17,face="bold")),
                      pltList$MissouriDownstream + theme(legend.position='none',
                                                         plot.title = element_text(size = 25, face = "bold"),
                                                         axis.text=element_text(size=17),
                                                         axis.title=element_text(size=17,face="bold")),
                      pltList$MissouriMidsection + theme(legend.position='none',
                                                         plot.title = element_text(size = 25, face = "bold"),
                                                         axis.text=element_text(size=17),
                                                         axis.title=element_text(size=17,face="bold")),
                      pltList$MissouriUpstream + theme(legend.position='none',
                                                       plot.title = element_text(size = 25, face = "bold"),
                                                       axis.text=element_text(size=17),
                                                       axis.title=element_text(size=17,face="bold")),
                      pltList$Ohio + theme(legend.position='none',
                                           plot.title = element_text(size = 25, face = "bold"),
                                           axis.text=element_text(size=17),
                                           axis.title=element_text(size=17,face="bold")),
                      pltList$OhioSection1 + theme(legend.position='none',
                                                   plot.title = element_text(size = 25, face = "bold"),
                                                   axis.text=element_text(size=17),
                                                   axis.title=element_text(size=17,face="bold")),
                      pltList$OhioSection2 + theme(legend.position='none',
                                                   plot.title = element_text(size = 25, face = "bold"),
                                                   axis.text=element_text(size=17),
                                                   axis.title=element_text(size=17,face="bold")),
                      pltList$OhioSection3 + theme(legend.position='none',
                                                   plot.title = element_text(size = 25, face = "bold"),
                                                   axis.text=element_text(size=17),
                                                   axis.title=element_text(size=17,face="bold")),
                      pltList$OhioSection4 + theme(legend.position='none',
                                                   plot.title = element_text(size = 25, face = "bold"),
                                                   axis.text=element_text(size=17),
                                                   axis.title=element_text(size=17,face="bold")),
                      pltList$OhioSection5 + theme(legend.position='none',
                                                   plot.title = element_text(size = 25, face = "bold"),
                                                   axis.text=element_text(size=17),
                                                   axis.title=element_text(size=17,face="bold")),
                      pltList$OhioSection7 + theme(legend.position='none',
                                                   plot.title = element_text(size = 25, face = "bold"),
                                                   axis.text=element_text(size=17),
                                                   axis.title=element_text(size=17,face="bold")),
                      pltList$OhioSection8 + theme(legend.position='none',
                                                   plot.title = element_text(size = 25, face = "bold"),
                                                   axis.text=element_text(size=17),
                                                   axis.title=element_text(size=17,face="bold")),
                      pltList$Olentangy + theme(legend.position='none',
                                                plot.title = element_text(size = 25, face = "bold"),
                                                axis.text=element_text(size=17),
                                                axis.title=element_text(size=17,face="bold")),
                      pltList$Padma + theme(legend.position='none',
                                            plot.title = element_text(size = 25, face = "bold"),
                                            axis.text=element_text(size=17),
                                            axis.title=element_text(size=17,face="bold")),
                      pltList$Platte + theme(legend.position='none',
                                             plot.title = element_text(size = 25, face = "bold"),
                                             axis.text=element_text(size=17),
                                             axis.title=element_text(size=17,face="bold")),
                      pltList$Po + theme(legend.position='none',
                                         plot.title = element_text(size = 25, face = "bold"),
                                         axis.text=element_text(size=17),
                                         axis.title=element_text(size=17,face="bold")),
                      pltList$SacramentoDownstream + theme(legend.position='none',
                                                           plot.title = element_text(size = 25, face = "bold"),
                                                           axis.text=element_text(size=17),
                                                           axis.title=element_text(size=17,face="bold")),
                      pltList$SacramentoUpstream + theme(legend.position='none',
                                                         plot.title = element_text(size = 25, face = "bold"),
                                                         axis.text=element_text(size=17),
                                                         axis.title=element_text(size=17,face="bold")),
                      pltList$SanJoaquin + theme(legend.position='none',
                                                 plot.title = element_text(size = 25, face = "bold"),
                                                 axis.text=element_text(size=17),
                                                 axis.title=element_text(size=17,face="bold")),
                      pltList$SanJoaquinRiver2 + theme(legend.position='none',
                                                       plot.title = element_text(size = 25, face = "bold"),
                                                       axis.text=element_text(size=17),
                                                       axis.title=element_text(size=17,face="bold")),
                      pltList$Seine + theme(legend.position='none',
                                            plot.title = element_text(size = 25, face = "bold"),
                                            axis.text=element_text(size=17),
                                            axis.title=element_text(size=17,face="bold")),
                      pltList$SeineDownstream + theme(legend.position='none',
                                                      plot.title = element_text(size = 25, face = "bold"),
                                                      axis.text=element_text(size=17),
                                                      axis.title=element_text(size=17,face="bold")),
                      pltList$SeineUpstream + theme(legend.position='none',
                                                    plot.title = element_text(size = 25, face = "bold"),
                                                    axis.text=element_text(size=17),
                                                    axis.title=element_text(size=17,face="bold")),
                      pltList$Severn + theme(legend.position='none',
                                             plot.title = element_text(size = 25, face = "bold"),
                                             axis.text=element_text(size=17),
                                             axis.title=element_text(size=17,face="bold")),
                      pltList$StanislausRiver + theme(legend.position='none',
                                                      plot.title = element_text(size = 25, face = "bold"),
                                                      axis.text=element_text(size=17),
                                                      axis.title=element_text(size=17,face="bold")),
                      pltList$TuolumneRiver + theme(legend.position='none',
                                                    plot.title = element_text(size = 25, face = "bold"),
                                                    axis.text=element_text(size=17),
                                                    axis.title=element_text(size=17,face="bold")),
                      pltList$Wabash + theme(legend.position='none',
                                             plot.title = element_text(size = 25, face = "bold"),
                                             axis.text=element_text(size=17),
                                             axis.title=element_text(size=17,face="bold")),
                      ncol=4)

#grab legend from one of these plots
legend <- get_legend(pltList$Wabash + theme(legend.text=element_text(size=25)))

#draw legend in remaining empty space in figure
plotTimeseries <- plotgrid + draw_grob(legend, 0.65, -0.2, 0.5, 0.5)

#create x and y axis labels for entire figure
yTitleCombo <- textGrob(expression(k[600]/Max~k[600]~(m/dy)), gp=gpar(fontface="bold", col="black", fontsize=30), rot=90)
xTitleCombo <- textGrob(expression(Timestep), gp=gpar(fontface="bold", col="black", fontsize=30))
plotTimeseries <- gridExtra::grid.arrange(gridExtra::arrangeGrob(plotTimeseries, left = yTitleCombo, bottom = xTitleCombo))

#write to file
ggsave('cache/validation/figS4.jpg', plotTimeseries, width = 20, height = 26)

################
##TIMERSERIES PLOT WITH SWOT ERRORS-----------------------------------
#################
#gather and plot all ggplot objects (for no measurement error) as a single image via the cowplot package
plotgrid_err <- plot_grid(pltList$MissouriDownstream_err + theme(legend.position='none',
                                                                 plot.title = element_text(size = 25, face = "bold"),
                                                                 axis.text=element_text(size=17),
                                                                 axis.title=element_text(size=17,face="bold")),
                      pltList$MissouriUpstream_err + theme(legend.position='none',
                                                           plot.title = element_text(size = 25, face = "bold"),
                                                           axis.text=element_text(size=17),
                                                           axis.title=element_text(size=17,face="bold")),
                      pltList$OhioSection8_err+ theme(legend.position='none',
                                                      plot.title = element_text(size = 25, face = "bold"),
                                                      axis.text=element_text(size=17),
                                                      axis.title=element_text(size=17,face="bold")),
                      pltList$OhioSection7_err+ theme(legend.position='none',
                                                      plot.title = element_text(size = 25, face = "bold"),
                                                      axis.text=element_text(size=17),
                                                      axis.title=element_text(size=17,face="bold")),
                      pltList$OhioSection5_err+ theme(legend.position='none',
                                                      plot.title = element_text(size = 25, face = "bold"),
                                                      axis.text=element_text(size=17),
                                                      axis.title=element_text(size=17,face="bold")),
                      pltList$OhioSection4_err+ theme(legend.position='none',
                                                      plot.title = element_text(size = 25, face = "bold"),
                                                      axis.text=element_text(size=17),
                                                      axis.title=element_text(size=17,face="bold")),
                      pltList$OhioSection3_err+ theme(legend.position='none',
                                                      plot.title = element_text(size = 25, face = "bold"),
                                                      axis.text=element_text(size=17),
                                                      axis.title=element_text(size=17,face="bold")),
                      pltList$OhioSection2_err+ theme(legend.position='none',
                                                      plot.title = element_text(size = 25, face = "bold"),
                                                      axis.text=element_text(size=17),
                                                      axis.title=element_text(size=17,face="bold")),
                      pltList$OhioSection1_err+ theme(legend.position='none',
                                                      plot.title = element_text(size = 25, face = "bold"),
                                                      axis.text=element_text(size=17),
                                                      axis.title=element_text(size=17,face="bold")),
                      pltList$SeineDownstream_err+ theme(legend.position='none',
                                                         plot.title = element_text(size = 25, face = "bold"),
                                                         axis.text=element_text(size=17),
                                                         axis.title=element_text(size=17,face="bold")),
                      pltList$Brahmaputra_err+ theme(legend.position='none',
                                                     plot.title = element_text(size = 25, face = "bold"),
                                                     axis.text=element_text(size=17),
                                                     axis.title=element_text(size=17,face="bold")),
                      pltList$IowaRiver_err+ theme(legend.position='none',
                                                   plot.title = element_text(size = 25, face = "bold"),
                                                   axis.text=element_text(size=17),
                                                   axis.title=element_text(size=17,face="bold")),
                      pltList$SeineUpstream_err+ theme(legend.position='none',
                                                       plot.title = element_text(size = 25, face = "bold"),
                                                       axis.text=element_text(size=17),
                                                       axis.title=element_text(size=17,face="bold")),
                      pltList$Kushiyara_err+ theme(legend.position='none',
                                                   plot.title = element_text(size = 25, face = "bold"),
                                                   axis.text=element_text(size=17),
                                                   axis.title=element_text(size=17,face="bold")),
                      pltList$Jamuna_err+ theme(legend.position='none',
                                                plot.title = element_text(size = 25, face = "bold"),
                                                axis.text=element_text(size=17),
                                                axis.title=element_text(size=17,face="bold")),
                      pltList$Padma_err+ theme(legend.position='none',
                                               plot.title = element_text(size = 25, face = "bold"),
                                               axis.text=element_text(size=17),
                                               axis.title=element_text(size=17,face="bold")),
                      ncol=3)

#draw legend in remaining empty space in figure
plotTimeseries <- plotgrid_err + draw_grob(legend, 0.6, -0.15, 0.5, 0.5)
plotTimeseries <- gridExtra::grid.arrange(gridExtra::arrangeGrob(plotTimeseries, left = yTitleCombo, bottom = xTitleCombo))

#write to file
ggsave('cache/validation/figS5.jpg', plotTimeseries, width = 13, height = 13)

#################
##PLOT REPRESENTIVE TIMESERIES FOR THE MAIN TEXT FIGURE-----------------------
  #Use KGE quantiles to inform timerseries selection
##################
set.seed(165)
quantiles <- quantile(stats_by_river[stats_by_river$errFlag == 0,]$KGE, c(0.33, 0.66))

stats_by_river$river <- as.character(stats_by_river$river)

temp <- stats_by_river[stats_by_river$errFlag==0 & stats_by_river$KGE < quantiles[1],] #low skill
temp2 <- stats_by_river[stats_by_river$errFlag==0 & stats_by_river$KGE >= quantiles[1] & stats_by_river$KGE < quantiles[2],] #med skill
temp3 <- stats_by_river[stats_by_river$errFlag==0 & stats_by_river$KGE >= quantiles[2],] #high skill

plotGrid_short <- plot_grid(
  pltList[[temp[13,]$river]]+ ggtitle('Sacramento River \n(upstream)')+ theme(legend.position='none',
                                     axis.text=element_text(size=20),
                                     axis.title=element_text(size=24,face="bold"),
                                     legend.text = element_text(size=17),
                                     legend.title = element_text(size=17, face='bold'),
                                     plot.title = element_text(size = 30, face = "bold")),
  pltList[[temp[9,]$river]]+ ggtitle('Ohio')+ theme(legend.position='none',
                                      axis.text=element_text(size=20),
                                      axis.title=element_text(size=24,face="bold"),
                                      legend.text = element_text(size=17),
                                      legend.title = element_text(size=17, face='bold'),
                                      plot.title = element_text(size = 30, face = "bold")),
  textGrob('Low Skill', y=0.6, gp=gpar(fontface="bold", col="darkblue", fontsize=34)),
  pltList[[temp2[7,]$river]]+ ggtitle('Iowa River') + theme(legend.position='none',
                                     axis.text=element_text(size=20),
                                     axis.title=element_text(size=24,face="bold"),
                                     legend.text = element_text(size=17),
                                     legend.title = element_text(size=17, face='bold'),
                                     plot.title = element_text(size = 30, face = "bold")),
  pltList[[temp2[4,]$river]]+ggtitle('Connecticut River')+ theme(legend.position='none',
                                      axis.text=element_text(size=20),
                                      axis.title=element_text(size=24,face="bold"),
                                      legend.text = element_text(size=17),
                                      legend.title = element_text(size=17, face='bold'),
                                      plot.title = element_text(size = 30, face = "bold")),
  textGrob('Middle Skill', y=0.6, gp=gpar(fontface="bold", col="darkblue", fontsize=34)),
  pltList[[temp3[5,]$river]]+ ggtitle('Missouri River\n(midsection)')+theme(legend.position='none',
                                        axis.text=element_text(size=20),
                                        axis.title=element_text(size=24,face="bold"),
                                        legend.text = element_text(size=17),
                                        legend.title = element_text(size=17, face='bold'),
                                        plot.title = element_text(size = 30, face = "bold")),
  pltList[[temp3[12,]$river]]+ggtitle('Sacramento River\n(downstream)')+theme(legend.position='none',
                                        axis.text=element_text(size=20),
                                        axis.title=element_text(size=24,face="bold"),
                                        legend.text = element_text(size=17),
                                        legend.title = element_text(size=17, face='bold'),
                                        plot.title = element_text(size = 30, face = "bold")),
  textGrob('High Skill', y=0.6, gp=gpar(fontface="bold", col="darkblue", fontsize=34)),
  ncol=3,
  label_size = 18,
  labels=c('a', 'b', NA, 'c', 'd', NA, 'e', 'f', NA))

#draw legend in remaining empty space in figure
plotTimeseries <- plotGrid_short + draw_grob(legend, 0.6, -0.2, 0.5, 0.5)
plotTimeseries <- gridExtra::grid.arrange(gridExtra::arrangeGrob(plotTimeseries, left = yTitleCombo, bottom = xTitleCombo))

#write to file
ggsave('cache/validation/fig3.jpg', plotTimeseries, width = 15, height = 12)

####################
##SAVE RESULTS TO FILE---------------------------------------------------
###################
write.csv(stats_by_river, 'cache/validation/results_by_riv.csv')
