###########################
#Description: Assessing the k600~eD relationship using data from Ulseth etal 2019
#Creator: Craig Brinkerhoff
#Summer 2021
###########################

print('starting eD~k600 relationship...')

options(scipen = 999)

#Ulseth etal 2019------------------------------------------------------------
ulseth_data <- read.csv('data/Ulseth_etal_2019.csv', fileEncoding="UTF-8-BOM") #contains data from 5 studies
ulseth_data <- ulseth_data[,-8]
ulseth_data$data <- ifelse(ulseth_data$data == 'This Study', 'Ulseth et al. 2019', as.character(ulseth_data$data)) #fix some labeling
ulseth_data$data <- ifelse(ulseth_data$data == 'Raymond et al. 2012', 'Raymond et al. 2012 [5 studies]', as.character(ulseth_data$data))

data <- ulseth_data
data <- filter(data, is.na(width)==0) #filter out no width measurements

#calculate eD and other hydraulics---------------------------------------------------------
data$eD <- g * data$slope * data$Vms #dissipation rate of surface turbulence originating from depth-scale form drag [W/kg]

#log transform some variables-----------------------------------------------------------
data$log_eD <- log(data$eD)
data$log_k600 <- log(data$k600)

#Generate regression line from Ulseth etal 2019---------------
data$ulseth_k600 <- ifelse(data$eD < 0.02, exp(3.10+0.35*log(data$eD)), exp(6.43+1.18*log(data$eD)))

data$SWOTflag <- ifelse(data$width > 100, 'SWOT-observable', 'Small')

t <- data[data$SWOTflag == 'SWOT-observable',]
lm_swot <- lm(t$log_k600 ~ t$log_eD)
r2_swot <- round(summary(lm_swot)$r.squared,2)

t <- data[data$SWOTflag != 'SWOT-observable',]
lm_small <- lm(t$log_k600 ~ t$log_eD)
r2_small <- round(summary(lm_small)$r.squared,2)

################
##K600 VS ED PLOT-------------------
#################
eD_k600 <- ggplot(data, aes(x=eD, y=k600, color=SWOTflag)) +
  geom_point(size=4, alpha=0.50) +
  geom_smooth(method='lm', se=F, size=3)+
  #geom_line(aes(x=eD, y=ulseth_k600), color='darkgreen', size=3) +
  xlab('eD [W/kg]') +
  ylab('k600 [m/dy]') +
  scale_color_brewer(palette='Dark2', name='')+
  annotate("text", label = paste0('SWOT r2: ', r2_swot), x = 10^-4, y = 10^3, size = 5, colour = "#1b9e77")+
  annotate("text", label = paste0('Small r2: ', r2_small), x = 10^-4, y = 10^2.5, size = 5, colour = "#d95f02")+
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'),
        legend.position = 'bottom')
ggsave('cache/k_eD/eD_k600.jpg', eD_k600, width=9, height=7)
