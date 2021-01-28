#scrap-------------------------------------------------

# #Is dataset representative for SWOT observable rivers?
# ggplot(data, aes(x=width)) +
#   stat_ecdf(size=2)+
#   geom_vline(xintercept = 100, linetype='dashed', size=1.25) +
#   geom_vline(xintercept = 50, linetype='dashed', size=1.25) +
#   annotate("text", x = 800, y = 0.8, label = 'SWOT Observable')+
#   annotate("text", x = 1, y = 0.7, label = 'Not SWOT Observable')+
#   scale_x_log10(
#     breaks = scales::trans_breaks("log10", function(x) 10^x),
#     labels = scales::trans_format("log10", scales::math_format(10^.x))
#   ) +
#   xlab('Width [m]') +
#   ylab('Density')


# #plot roughness vs k600
# df_temp1 <- filter(data, slope >= 0.001)
# plot1 <- ggplot(df_temp1, aes(y=(k600), x=(keulegan_yr), color=data)) +
#   geom_point(size=3) +
#   geom_smooth(se=F, method='lm', color='black', size=2) +
#   annotate("text", x = 0.3, y = 1000, label = "Slope >= 0.001")+
#   annotate("text", x = 0.3, y = 2000, label = paste0('r2: ', round(summary(lm(log10(k600) ~ log10(keulegan_yr), data=df_temp1))$r.squared, 2)))+
#   scale_color_brewer(palette = 'Dark2', name='') +
#   scale_y_log10(
#     breaks = scales::trans_breaks("log10", function(x) 10^x),
#     labels = scales::trans_format("log10", scales::math_format(10^.x))
#   ) +
#   scale_x_log10(
#     breaks = scales::trans_breaks("log10", function(x) 10^x),
#     labels = scales::trans_format("log10", scales::math_format(10^.x))
#   )+
#   ylab('log K600 [m/dy]') +
#   xlab('log Effective Roughness Height [m]')
#
# df_temp2 <- filter(data, slope >= 0.01)
# plot2 <- ggplot(df_temp2, aes(y=(k600), x=(keulegan_yr), color=data)) +
#   geom_point(size=3) +
#   geom_smooth(se=F, method='lm', color='black', size=2) +
#   annotate("text", x = 0.3, y = 1000, label = "Slope >= 0.01")+
#   annotate("text", x = 0.3, y = 2000, label = paste0('r2: ', round(summary(lm(log10(k600) ~ log10(keulegan_yr), data=df_temp2))$r.squared, 2)))+
#   scale_color_brewer(palette = 'Dark2', name='') +
#   scale_y_log10(
#     breaks = scales::trans_breaks("log10", function(x) 10^x),
#     labels = scales::trans_format("log10", scales::math_format(10^.x))
#   ) +
#   scale_x_log10(
#     breaks = scales::trans_breaks("log10", function(x) 10^x),
#     labels = scales::trans_format("log10", scales::math_format(10^.x))
#   )+
#   ylab('log K600 [m/dy]') +
#   xlab('log Effective Roughness Height [m]')
#
# df_temp3 <- filter(data, slope >= 0.1)
# plot3 <- ggplot(df_temp3, aes(y=(k600), x=(keulegan_yr), color=data)) +
#   geom_point(size=3) +
#   geom_smooth(se=F, method='lm', color='black', size=2) +
#   annotate("text", x = 0.3, y = 1000, label = "Slope >= 0.1")+
#   annotate("text", x = 0.3, y = 2000, label = paste0('r2: ', round(summary(lm(log10(k600) ~ log10(keulegan_yr), data=df_temp3))$r.squared, 2)))+
#   scale_color_brewer(palette = 'Dark2', name='') +
#   scale_y_log10(
#     breaks = scales::trans_breaks("log10", function(x) 10^x),
#     labels = scales::trans_format("log10", scales::math_format(10^.x))
#   ) +
#   scale_x_log10(
#     breaks = scales::trans_breaks("log10", function(x) 10^x),
#     labels = scales::trans_format("log10", scales::math_format(10^.x))
#   )+
#   ylab('log K600 [m/dy]') +
#   xlab('log Effective Roughness Height [m]')
#
# df_temp4 <- filter(data, slope >= 0.15)
# plot4 <- ggplot(df_temp4, aes(y=(k600), x=(keulegan_yr), color=data)) +
#   geom_point(size=3) +
#   geom_smooth(se=F, method='lm', color='black', size=2) +
#   annotate("text", x = 0.3, y = 1000, label = "Slope >= 0.15")+
#   annotate("text", x = 0.3, y = 2000, label = paste0('r2: ', round(summary(lm(log10(k600) ~ log10(keulegan_yr), data=df_temp4))$r.squared, 2)))+
#   scale_color_brewer(palette = 'Dark2', name='') +
#   scale_y_log10(
#     breaks = scales::trans_breaks("log10", function(x) 10^x),
#     labels = scales::trans_format("log10", scales::math_format(10^.x))
#   ) +
#   scale_x_log10(
#     breaks = scales::trans_breaks("log10", function(x) 10^x),
#     labels = scales::trans_format("log10", scales::math_format(10^.x))
#   )+
#   ylab('log K600 [m/dy]') +
#   xlab('log Effective Roughness Height [m]')
#
# total_plot <- plot_grid(plot1, plot2, plot3, plot4, ncol=1)


#raymond et al 2012
# data <- read.csv('raymond_2012.csv')
# data <- select(data, c('Mean.width..ft.', 'Mean.depth..ft.', 'Water.surface.slope', 'Mean.stream.flow.velocity..ft.s.','Mean.discharge..cfs.', 'k600', 'Subreach.Name', 'k2.t1'))
# colnames(data) <- c('width_ft', 'depth_ft', 'slope', 'mean_velocity_ft_s','discharge_cfs', 'k600', 'name', 'K2_t1' )
# #data <- filter(data, k600 > 0, K2_t1 >= 0.3) #per Melching & Flores 1999

# data$discharge_cms <- data$discharge_cfs * 0.2083 #cfs to cms
# data$mean_velocity_m_s <- data$mean_velocity_ft_s * 0.3048 #ft to m
# data$mean_width_m <- data$width_ft * 0.3048 #ft to m
# data$mean_depth_m <- data$depth_ft * 0.3048 #ft to m
# data$area_m2 <- data$discharge_cms / (data$mean_velocity_m_s)
# data$hydraulic_radius <- data$area_m2 /(data$mean_width_m + 2*(data$mean_depth_m))
# data$mannings_n <- (data$mean_depth_m^(2/3) * data$slope^(1/2))/data$mean_velocity_m_s
# data$chezys_c <- data$mean_velocity_m_s * data$slope^(-1/2) * data$mean_depth_m^(-1/2) * 0.552
# data$friction_f <- (8*9.8*data$mean_depth_m*data$slope) / data$mean_velocity_m_s^2







# t <- goodRiver
# Qdf <- data.frame(NA, NA, c(1:15), NA, NA, NA, NA, NA, NA, NA, 'Qobs', colMeans(Q_obs)/1000)
# colnames(Qdf) <- colnames(t)
# t <- rbind(t, Qdf)
# temp <- ggplot(t, aes(x=time, y=value, color=key)) + #model
#   geom_pointrange(aes(ymin = kest_low, ymax = kest_high), fatten=10) +
#   geom_line(size=1) +
#   ylab('k600 [m/dy]') +
#   xlab('Timestep') +
#   scale_color_brewer(palette='Set2', name='', labels=c('Modeled k600', 'Observed k600', 'Discharge\n/1000')) +
#   theme(legend.position = 'bottom') +
#   ggtitle(riv)
# temp
# #Compare BIKER uncertainties against MC uncertanties-------------------------------------------------------
# prior_k600_uncertainity <- 1.029
# mean_modelSD_by_riv <- group_by(full_output, river) %>%
#   summarise(meanSD = mean(kest_sd, na.rm=T))
# uncertainity_comparison_alltimesteps <- ggplot(full_output, aes(x=kest_sd)) +
#   geom_histogram(size=1, color='black', fill='darkgreen', bins=30) +
#   geom_vline(xintercept = prior_k600_uncertainity, linetype='dashed', size=1.3, color='darkblue') +
#   geom_vline(xintercept = mean(full_output$kest_sd), size=1.3, color='darkred') +
#   xlab('BIKER k600 ln(sd)') +
#   ylab('Count')
# uncertainity_comparison_byRiv <- ggplot(mean_modelSD_by_riv, aes(x=meanSD)) +
#   geom_histogram(size=1, color='black', fill='darkgreen', bins=30) +
#   geom_vline(xintercept = prior_k600_uncertainity, linetype='dashed', size=1.3, color='darkblue') +
#   geom_vline(xintercept = mean(mean_modelSD_by_riv$meanSD), size=1.3, color='darkred') +
#   xlab('BIKER k600 ln(sd)') +
#   ylab('Count')
# 
# RS_uncertainity <- plot_grid(uncertainity_comparison_alltimesteps, uncertainity_comparison_byRiv, ncol=2, labels = 'auto')
# ggsave('outputs//validation/validation_uncertainity.jpg', RS_uncertainity)