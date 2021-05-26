#Monte Carlo analysis for quantifying uncertainty of k02 model and BIKER's assumptions
#Creator: Craig Brinkerhoff
#Summer 2021

print('staring MC analysis...')

set.seed(455)

#Pull in hydraulic measurements from Brinkerhoff etal 2019 to test parameter uncertainty------------------------------------------
if(munge == 1){
  hydraulics_init <- read.csv('data/MonteCarlo/field_measurements.csv')
  nhd_med <- read.csv('data/MonteCarlo/NHD_join_table.csv')

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

    #Get uncertainty from Rh~D only (for validation runs in the paper)----------------------------
    hydraulics$Rh <- hydraulics$area / (hydraulics$width + 2*hydraulics$depth)
    temp <- hydraulics #filter(hydraulics, width/depth > 10)
    lm <- lm((temp$Rh)~(temp$depth))
    print(summary(lm))

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
func_MC <- function(id) {
  hydraulics <- read.csv('data/MonteCarlo/mc_hydraulics.csv')

  s <- hydraulics$slope[id]
  d <- hydraulics$depth[id]

  #different tests to parse out influence of different assumptions
  Rh <- rnorm(10000, d, 0.20) #0.20 comes from the Rh ~ A/W comparison performed in munge above

  ######## 0.25 IS A DUMMY VALUE FOR NOW. THIS IS NOT REAL
  alpha <- rnorm(10000, 48, 0.25) #Sample from normal distribution. sigma taken from my rule-based regression model

  #k600 model
  Ustar <- sqrt(9.8*s*Rh)
  ko2 <- alpha*Ustar

  #save a few random k600 pdfs for plotting----------------
  if (id == 100) {
    ko2_1 <- ko2
    write.csv(ko2_1, 'cache/MonteCarlo/pdf_a.csv')
  }
  if (id == 1000) {
    ko2_2 <- ko2
    write.csv(ko2_2, 'cache/MonteCarlo/pdf_b.csv')
  }
  if (id == 3400) {
    ko2_3 <- ko2
    write.csv(ko2_3, 'cache/MonteCarlo/pdf_c.csv')
  }

  #get M different model estimate uncertainties and return that------------------------------
  temp <- data.frame('sigma_ko2' = sd(ko2, na.rm=T), 'mean_ko2'= mean(ko2, na.rm = T), 'id'=id)
}

#Run MC simulations in parallel--------------------------
if(munge2 == 1){
  print('MC simulation alpha + Rh errors')
  system.time(
    output <- mclapply(ids, func_MC, mc.cores=cores)
  )
  output <- data.frame(t(sapply(output, function(x) x[1:max(lengths(output))])))
  output$sigma_ko2 <- as.numeric(output$sigma_ko2)
  output$mean_ko2 <- as.numeric(output$mean_ko2)
  output$id <- as.numeric(output$id)
  write.csv(output, 'cache/MonteCarlo/mc_results.csv')
}

#plot results-------------------------------
df <- read.csv('cache/MonteCarlo/mc_results.csv')
medianSD <- round(median(df$sigma_ko2, na.rm=T),2)
print(medianSD)
plot <- ggplot(df, aes(x=sigma_ko2)) +
  geom_histogram(color='black', fill='darkgreen', size=1, bins=50) +
  xlab('SD of MC Simulated Estimates [m/dy]') +
  ylab('Count') +
  geom_vline(xintercept = medianSD, linetype='dashed', size=1.2, color='blue')  +
  geom_text(x=3, y=300, label=paste0('median \u03c3: ', medianSD))

df_1 <- read.csv('cache/MonteCarlo/pdf_a.csv')
df_2 <- read.csv('cache/MonteCarlo/pdf_b.csv')
df_3 <- read.csv('cache/MonteCarlo/pdf_c.csv')
df <- data.frame('a'= as.numeric(df_1$x), 'b'= as.numeric(df_2$x), 'c' = as.numeric(df_3$x))
df <- gather(df, key=key, value=value)
plot2 <- ggplot(df, aes(x=value, color=key)) +
    geom_density(size=2) +
    ylab('Density') +
    scale_color_brewer(palette = 'Dark2') +
    xlab('ko2 [m/dy]') +
    theme(legend.position = 'none')

#combine subplots-------------------------------------------------
combinedPlot <- plot_grid(plot, plot2, ncol=1, labels = 'auto', label_size = 18)
ggsave('cache/MonteCarlo/MC_results.jpg', combinedPlot, width=6, height=11)

#Save uncertainty to file for ms---------------------------------------
results <- data.frame('medianSigma'=medianSD)
write.csv(results, 'cache/MonteCarlo/results.csv')
