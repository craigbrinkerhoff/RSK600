###########
##SCRAPS-----------------------------------
######
#test plot------------------------
#prepare ulseth data for this plot...
#ulseth_data <- read.csv('data/k600/Ulseth_etal_2019.csv', fileEncoding="UTF-8-BOM") #contains data from 5 studies
#ulseth_data <- ulseth_data[,-8]
#ulseth_data$data <- ifelse(ulseth_data$data == 'This Study', 'Ulseth et al. 2019', as.character(ulseth_data$data)) #fix some labeling
#ulseth_data$data <- ifelse(ulseth_data$data == 'Raymond et al. 2012', 'Raymond et al. 2012 [5 studies]', as.character(ulseth_data$data))

#data <- ulseth_data
#data <- filter(data, is.na(width)==0) #filter out no width measurements

#calculate eD and other hydraulics---------------------------------------------------------
#data$eD <- g * data$slope * data$Vms #dissipation rate of surface turbulence originating from depth-scale form drag [W/kg]
#data$ustar <- sqrt(g*data$slope*data$depth)
#data$eS <- (sqrt(data$depth*g*data$slope))^3/data$depth #dissipation rate of surface turbulence originating from bed friction [W/kg]
#data$Rh <- (data$depth*data$width)/(data$width + 2*data$depth) #hydraulic radius [m]
#data$area <- data$width*data$depth #channel area [m2]
#data$keulegan <- 11*data$depth*(1/exp(data$Vms/(2.5*(9.8*data$depth*data$slope)^(1/2)))) #[m]

#log transform some variables
#data$log_eD <- log(data$eD)
#data$log_k600 <- log(data$k600)
#data$log_slope <- log(data$slope)
#data$log_area <- log(data$area)

#dummy values with regression models
#lseq <- function(from, to, length.out) {
  # logarithmic spaced sequence
  # blatantly stolen from library("emdbook"), because need only this
#  exp(seq(log(from), log(to), length.out = length.out))
#}

#test <- data.frame('ustar'=seq(0,0.8, 0.1)) #data.frame('eD'=lseq(10^-5, 10^1, 100))
#test$k600_1d <- 48*test$ustar #exp(4.44383)*test$eD^0.59957
#test$k600_2d <- ifelse(test$eD < exp(-6.494), exp(0.6277)*test$eD^(0.030075),
#                          ifelse(test$eD< exp(-2.220), exp(4.7048)*test$eD^(0.65790), exp(7.1266)*test$eD^(1.7489)))

#weirdPlot <- ggplot() +
#  geom_point(data=data, aes(y=k600, x=ustar), alpha=0.4, size=3, color='blue')+
#  geom_pointrange(data=full_output, mapping=aes(y=kest_mean, x=(kobs/48), ymin = kest_low, ymax = kest_high), fatten=10, fill='#1b9e77', pch=21, color='black') +
#  geom_line(data=test, aes(x=ustar, y=k600_1d), color='darkred', size=3) +
#  #geom_line(data=test, aes(x=eD, y=k600_2d), color='darkblue', size=3) +
#  xlab('Ustar [m/s]') +
#  ylab('BIKER ko2 [m/dy]') +
#  scale_y_log10(
#      breaks = scales::trans_breaks("log10", function(x) 10^x),
#      labels = scales::trans_format("log10", scales::math_format(10^.x)),
#      limits=c(10^-1,10^4))+
#  scale_x_log10(
#      breaks = scales::trans_breaks("log10", function(x) 10^x),
#      labels = scales::trans_format("log10", scales::math_format(10^.x)),
#      limits=c(10^-2.5, 10^0.5))+
#    scale_color_discrete_qualitative(palette = 'Harmonic') +
#  theme(legend.position = "none",
#        axis.text=element_text(size=20),
#        axis.title=element_text(size=24,face="bold"),
#        legend.text = element_text(size=17),
#        legend.title = element_text(size=17, face='bold'))
#ggsave('cache/validation/test.jpg', weirdPlot, width=9, height=7)
