#Monte Carlo analysis for quantifying uncertainty of k600 model and BIKER's assumptions
#Creator: Craig Brinkerhoff
#Summer 2020

print('staring MC analysis...')

set.seed(455)

#Pull in hydraulic measurements from Brinkerhoff etal 2019 to test parameter uncertainty------------------------------------------
if(munge == 1){
  hydraulics_init <- read.csv('cache/MonteCarlo/field_measurements.csv')
  nhd_med <- read.csv('cache/MonteCarlo/NHD_join_table.csv')

  #add reach slopes from NHD
  hydraulics <- left_join(hydraulics_init, select(nhd_med, SLOPE, SOURCE_FEA, LatSite, LonSite), by = c("site_no" = "SOURCE_FEA"))

  #intial cleaning of field measurements
  hydraulics <- filter(hydraulics, is.finite(chan_width) ==1) %>%
    filter(is.finite(chan_velocity)==1) %>%
    filter(is.finite(chan_discharge)==1) %>%
    filter(is.finite(chan_area)==1) %>%
    filter(chan_width > 0) %>%
    filter(chan_velocity > 0) %>%
    filter(chan_discharge > 0) %>%
    filter(chan_area > 0) %>%
    filter(measured_rating_diff != 'Poor') %>%
    filter(measured_rating_diff != 'POOR') %>%
    filter(is.finite(SLOPE)==1)

    #imperial to metric
    hydraulics$area <- hydraulics$chan_area * 0.092903 #ft2 to m2
    hydraulics$width <- hydraulics$chan_width*0.305 #m
    hydraulics$velocity <- hydraulics$chan_velocity*0.305 #m/s
    hydraulics$slope <- hydraulics$SLOPE
    hydraulics$depth <- hydraulics$area / hydraulics$width #m
    hydraulics$n <- (hydraulics$depth^(2/3)*hydraulics$SLOPE^(1/2))/hydraulics$velocity

    #many slopes were saved with floor values of 1e-5, so I just remove these because they influence the k600 model to provide artificial results
    hydraulics <- filter(hydraulics, slope > 1e-5)

    #randomly grab M hydraulic sets from the Brinkerhoff etal 2019 dataset (which is 500,000+...)
    ids <- runif(M, 1, nrow(hydraulics))
    hydraulics <- hydraulics[ids,]
    write.csv(hydraulics, 'inputs/MonteCarlo/mc_hydraulics.csv')
}

#Gather measurement ids
ids <- c(1:M)

#Function that runs M Monte Carlo simulations (of 10,000 runs each) in parallel--------------------------------------------------------------
func_MC <- function(id, mannings_flag, slopeFlag) {
  hydraulics <- read.csv('data/MonteCarlo/mc_hydraulics.csv')

  s <- hydraulics$slope[id]
  v <- hydraulics$velocity[id]
  w <- hydraulics$width[id]

  #different tests to parse out influence of different assumptions
  vel <- ifelse(mannings_flag == 1, exp(rnorm(10000, log(v), mannings_uncertainity)), v) #Include mannings error?
  s <- ifelse(slopeFlag == 1, rnorm(10000, s, 1.7e-5), s) #include swot slope error?

  if(w < 10 & s < 0.05){ #< 10m wide and slope < 0.05
    exp <- rnorm(10000, 0.6488131, 0.173) #Sample from normal distribution. sigma taken from my rule-based regression model
    int <- rnorm(10000, 111.58121, 1.038) #Sample from normal distribution. sigma taken from my rule-based regression model
  }
  else if(w < 10 & s >= 0.05){ #< 10m wide and slope >= 0.05
    exp <- rnorm(10000, 1.3160065, 0.177) #Sample from normal distribution. sigma taken from my rule-based regression model
    int <- rnorm(10000, 792.63149, 1.087) #Sample from normal distribution. sigma taken from my rule-based regression model
  }
  else if(w < 50 & w >= 10){ #10-50m wide
    exp <- rnorm(10000, 0.6624354, 0.241) #Sample from normal distribution. sigma taken from my rule-based regression model
    int <- rnorm(10000, 109.04977, 1.046) #Sample from normal distribution. sigma taken from my rule-based regression model
  }
  else if(w < 100 & w >= 50){ #50-100m wide
    exp <- rnorm(10000, 0.4366114, 0.359) #Sample from normal distribution. sigma taken from my rule-based regression model
    int <- rnorm(10000, 31.84344, 1.059) #Sample from normal distribution. sigma taken from my rule-based regression model
  }
  else{ #> 100m wide
    exp <- rnorm(10000, 0.2724834, 0.201) #Sample from normal distribution. sigma taken from my rule-based regression model
    int <- rnorm(10000, 14.16939, 1.031) #Sample from normal distribution. sigma taken from my rule-based regression model
  }

  #k600 model
  log_eD <- log(s * 9.8 * vel)
  logk <- exp*log_eD + log(int)

  #save a few random k600 pdfs for plotting----------------
  if (id == 100) {
    logk_1 <- logk
    if(mannings_flag == 0 & slopeFlag == 0){
      write.csv(logk_1, 'cache/MonteCarlo/pdf_a_0_0.csv')
    }
    if(mannings_flag == 1 & slopeFlag == 0){
      write.csv(logk_1, 'cache/MonteCarlo/pdf_a_1_0.csv')
    }
    if(mannings_flag == 1 & slopeFlag == 1){
      write.csv(logk_1, 'cache/MonteCarlo/pdf_a_1_1.csv')
    }
  }
  if (id == 1000) {
    logk_2 <- logk
    if(mannings_flag == 0 & slopeFlag == 0){
      write.csv(logk_2, 'cache/MonteCarlo/pdf_b_0_0.csv')
    }
    if(mannings_flag == 1 & slopeFlag == 0){
      write.csv(logk_2, 'cache/MonteCarlo/pdf_b_1_0.csv')
    }
    if(mannings_flag == 1 & slopeFlag == 1){
      write.csv(logk_2, 'cache/MonteCarlo/pdf_b_1_1.csv')
    }
  }
  if (id == 3400) {
    logk_3 <- logk
    if(mannings_flag == 0 & slopeFlag == 0){
      write.csv(logk_3, 'cache/MonteCarlo/pdf_c_0_0.csv')
    }
    if(mannings_flag == 1 & slopeFlag == 0){
      write.csv(logk_3, 'cache/MonteCarlo/pdf_c_1_0.csv')
    }
    if(mannings_flag == 1 & slopeFlag == 1){
      write.csv(logk_3, 'cache/MonteCarlo/pdf_c_1_1.csv')
    }
  }

  #get M different model estimate uncertainties and return that------------------------------
  temp <- data.frame('sigma_logk600' = sd(logk, na.rm=T), 'mean_logk600'= mean(logk, na.rm = T), 'id'=id)
}

#Run MC simulations in parallel--------------------------
if(munge2 == 1){
  #No manning's or slope error
  print('MC simulation upscaling error')
  system.time(
    output <- mclapply(ids, func_MC, 0, 0, mc.cores=cores)
  )
  output <- data.frame(t(sapply(output, function(x) x[1:max(lengths(output))])))
  output$sigma_logk600 <- as.numeric(output$sigma_logk600)
  output$mean_logk600 <- as.numeric(output$mean_logk600)
  output$id <- as.numeric(output$id)
  write.csv(output, 'cache/MonteCarlo/mc_0_0.csv')

  #No manning's or slope error
  print('MC simulation upscaling + mannings errors')
  system.time(
    output <- mclapply(ids, func_MC, 1, 0, mc.cores=cores)
  )
  output <- data.frame(t(sapply(output, function(x) x[1:max(lengths(output))])))
  output$sigma_logk600 <- as.numeric(output$sigma_logk600)
  output$mean_logk600 <- as.numeric(output$mean_logk600)
  output$id <- as.numeric(output$id)
  write.csv(output, 'cache/MonteCarlo/mc_1_0.csv')

  #No manning's or slope error
  print('MC simulation upscaling + mannings and swot errors')
  system.time(
    output <- mclapply(ids, func_MC, 1, 1, mc.cores=cores)
  )
  output <- data.frame(t(sapply(output, function(x) x[1:max(lengths(output))])))
  output$sigma_logk600 <- as.numeric(output$sigma_logk600)
  output$mean_logk600 <- as.numeric(output$mean_logk600)
  output$id <- as.numeric(output$id)
  write.csv(output, 'cache/MonteCarlo/mc_1_1.csv')
}

#plot results-------------------------------
#ERROR SECNARIO 1
df_0_0 <- read.csv('cache/MonteCarlo/mc_0_0.csv')
medianSD <- round(median(df_0_0$sigma_logk600, na.rm=T),2)
print(medianSD)
plot_0_0 <- ggplot(df_0_0, aes(x=sigma_logk600)) +
  geom_histogram(color='black', fill='darkgreen', size=1, bins=50) +
  xlab('')+
  ylab('Count') +
  geom_vline(xintercept = medianSD, linetype='dashed', size=1.2, color='blue') +
  geom_text(x=3, y=300, label=paste0('median \u03c3: ', medianSD)) +
  ggtitle('Upscaling Uncertainty')

df_1 <- read.csv('cache/MonteCarlo/pdf_a_0_0.csv')
df_2 <- read.csv('cache/MonteCarlo/pdf_b_0_0.csv')
df_3 <- read.csv('cache/MonteCarlo/pdf_c_0_0.csv')
df <- data.frame('a'= as.numeric(df_1$x), 'b'= as.numeric(df_2$x), 'c' = as.numeric(df_3$x))
df <- gather(df, key=key, value=value)
plot2_0_0 <- ggplot(df, aes(x=value, color=key)) +
  geom_density(size=2) +
  coord_cartesian(ylim=c(0,0.5))+
  ylab('Density') +
  scale_color_brewer(palette = 'Dark2') +
  xlab('') +
  theme(legend.position = 'none')

#ERROR SCENARIO 2
df_1_0<- read.csv('cache/MonteCarlo/mc_1_0.csv')
medianSD <- round(median(df_1_0$sigma_logk600, na.rm=T),2)
print(medianSD)
plot_1_0 <- ggplot(df_1_0, aes(x=sigma_logk600)) +
  geom_histogram(color='black', fill='darkgreen', size=1, bins=50) +
  xlab('') +
  ylab('Count') +
  geom_vline(xintercept = medianSD, linetype='dashed', size=1.2, color='blue') +
  geom_text(x=3, y=300, label=paste0('median \u03c3: ', medianSD)) +
  ggtitle('Upscaling + \nMannings Uncertainty')

df_1 <- read.csv('cache/MonteCarlo/pdf_a_1_0.csv')
df_2 <- read.csv('cache/MonteCarlo/pdf_b_1_0.csv')
df_3 <- read.csv('cache/MonteCarlo/pdf_c_1_0.csv')
df <- data.frame('a'= as.numeric(df_1$x), 'b'= as.numeric(df_2$x), 'c' = as.numeric(df_3$x))
df <- gather(df, key=key, value=value)
plot2_1_0 <- ggplot(df, aes(x=value, color=key)) +
    geom_density(size=2) +
    coord_cartesian(ylim=c(0,0.5))+
    ylab('Density') +
    scale_color_brewer(palette = 'Dark2') +
    xlab('') +
    theme(legend.position = 'none')

#ERROR SCENARIO 3
df_1_1 <- read.csv('cache/MonteCarlo/mc_1_1.csv')
medianSD <- round(median(df_1_1$sigma_logk600, na.rm=T),2)
print(medianSD)
plot_1_1 <- ggplot(df_1_1, aes(x=sigma_logk600)) +
  geom_histogram(color='black', fill='darkgreen', size=1, bins=50) +
  xlab('lnSD of MC Simulated Estimates [m/dy]') +
  ylab('Count') +
  geom_vline(xintercept = medianSD, linetype='dashed', size=1.2, color='blue')  +
  geom_text(x=3, y=300, label=paste0('median \u03c3: ', medianSD)) +
  ggtitle('Upscaling + \nMannings + SWOT \nUncertainty')

df_1 <- read.csv('cache/MonteCarlo/pdf_a_1_1.csv')
df_2 <- read.csv('cache/MonteCarlo/pdf_b_1_1.csv')
df_3 <- read.csv('cache/MonteCarlo/pdf_c_1_1.csv')
df <- data.frame('a'= as.numeric(df_1$x), 'b'= as.numeric(df_2$x), 'c' = as.numeric(df_3$x))
df <- gather(df, key=key, value=value)
plot2_1_1 <- ggplot(df, aes(x=value, color=key)) +
    geom_density(size=2) +
    coord_cartesian(ylim=c(0,0.5))+
    ylab('Density') +
    scale_color_brewer(palette = 'Dark2') +
    xlab('ln k600 [m/dy]') +
    theme(legend.position = 'none')

#combine subplots-------------------------------------------------
combinedPlot <- plot_grid(plot_0_0, plot2_0_0, plot_1_0, plot2_1_0, plot_1_1, plot2_1_1, ncol=2, labels = 'auto', label_size = 18)
ggsave('cache/MonteCarlo/MC_results.jpg', combinedPlot, width=8, height=9)

#Save uncertainty to file for ms---------------------------------------
results <- data.frame('medianSigma'=medianSD)
write.csv(results, 'cache/MonteCarlo/results.csv')
