###########################
#Description: Modeling k600 from hydraulics. Comparing upscaling models for k600 from eD
#Creator: Craig Brinkerhoff
#Winter 2021
###########################

#load model settings and packages-------------------------------
source(here::here('scripts' , 'inputs.R'))

#Ulseth etal 2019------------------------------------------------------------
ulseth_data <- read.csv('inputs//k600//Ulseth_etal_2019.csv', fileEncoding="UTF-8-BOM") #contains data from 5 studies
ulseth_data <- ulseth_data[,-8]
ulseth_data$data <- ifelse(ulseth_data$data == 'This Study', 'Ulseth et al. 2019', as.character(ulseth_data$data)) #fix some labeling
ulseth_data$data <- ifelse(ulseth_data$data == 'Raymond et al. 2012', 'Raymond et al. 2012 [5 studies]', as.character(ulseth_data$data))

data <- ulseth_data
data <- filter(data, is.na(width)==0) #filter out no width measurements

data$eD <- g * data$slope * data$Vms #water column turbulent dissipation rate m2/s3

#width regime models for k600-----------------------------------------------------------------------
data$widthRegime <- ifelse(data$width < 10 & data$slope < 0.05, '1', #meters
                            ifelse(data$width < 10 & data$slope > 0.05, '1.5',
                                ifelse(data$width < 50, '2',
                                   ifelse(data$width < 100, '3','4'))))

#log transform some variables-----------------------------------------------------------
data$log_eD <- log(data$eD)
data$log_k600 <- log(data$k600)
data$log_slope <- log(data$slope)
data$raymond <- data$slope*data$Vms #for raymonds 2012 model

#Check linearity assumptions within each width class. If you cycle through the different classes, everything holds
#lm <- lm(log_k600~log_eD, data=filter(data, widthRegime==1))
#autoplot(lm)

#80/20 training and test split----------------------------------------
  #This is just for validation of method. In implementation I train model on entire dataset (see below)
set.seed(234)
require(caTools)
sample = sample.split(data$k600, SplitRatio = .80)
data_train = subset(data, sample == TRUE)
data_test  = subset(data, sample == FALSE)

#Other models: Ulseth and raymond models using only V and S
lm_raymond <- lm(k600~raymond, data=data_train) #raymond model- linear, eq 5
lm_raymond_2 <- lm(log_k600~log(raymond), data=data_train)#raymond model 2- power law, eq 4
lm_raymond_3 <- lm(log_k600~log(slope)+log(Vms), data=data_train)#raymond model 3- multi power law, eq 3

#Ulseth model
linearModel <- lm(log_k600~log_eD, data=data_train)
lm_ulseth <- segmented(linearModel, npsi = 1)

#Our model
#model k600 using rule-based linear regression
k_model_train <- group_by(data_train, widthRegime) %>%
  do(model = lm(log(k600) ~ log(eD), data=.)) %>%
  mutate(a = exp(model$coefficients[1])) %>%
  mutate(b = model$coefficients[2]) %>%
  mutate(a_se = round(exp(summary(model)$coefficients[2,2]), 3)) %>% #confidence intervals for later uncertainity analysis
  mutate(b_se = round(summary(model)$coefficients[1,2], 3)) %>% #confidence intervals for later uncertainity analysis
  mutate(k600_eD_r2 = round(summary(model)$r.squared, 2))
k_model_groups_train <- group_by(data_train, widthRegime) %>%
  summarise(n = n())

k_model_train <- left_join(k_model_train, k_model_groups_train, by='widthRegime')
data_train <- left_join(data_train,k_model_train, by="widthRegime")
data_test <- left_join(data_test,k_model_train, by="widthRegime")

#predict k600
data_test$k600_pred_me <- data_test$a * data_test$eD^data_test$b
data_test$k600_pred_raymond <- predict(lm_raymond, data_test)
data_test$k600_pred_raymond_2 <- exp(predict(lm_raymond_2, data_test))
data_test$k600_pred_raymond_3 <- exp(predict(lm_raymond_3, data_test))
data_test$k600_pred_ulseth <- exp(predict(lm_ulseth, data_test))

#Plot validation - me
r2_me <- round(summary(lm(log10(data_test$k600) ~ log10(data_test$k600_pred_me)))$r.squared,2)
rmse_me <- round(10^(Metrics::rmse(log10(data_test$k600), log10(data_test$k600_pred_me))),2)
test_vali_me <- ggplot(data_test, aes(x=(k600), y=(k600_pred_me), color=data)) +
  geom_abline(size=2, color='grey', linetype='dashed') +
  geom_point(size=5, alpha=0.75) +
  geom_smooth(color='black', method='lm', size=2, se=F) +
  xlab('Observed k600 [m/dy]') +
  ylab('Predicted k600 [m/dy]') +
  geom_richtext(aes(x=4, y=300), label=paste0('r<sup>2</sup>: ', r2_me, '0'), color='black')+
  geom_richtext(aes(x=4, y=170), label=paste0('RMSE: ', rmse_me, ' m/dy'), color='black')+
  scale_color_brewer(palette = 'Dark2', name='Model Validation (80/20 Split)')+
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  )

#Plot validation - raymond
# No error metrics b/c model produces negatives and then NAs in log space
test_vali_raymond <- ggplot(data_test, aes(x=(k600), y=(k600_pred_raymond), color=data)) +
  geom_abline(size=2, color='grey', linetype='dashed') +
  geom_point(size=5, alpha=0.75) +
  geom_smooth(color='black', method='lm', size=2, se=F) +
  xlab('Observed k600 [m/dy]') +
  ylab('Predicted k600 [m/dy]') +
  geom_richtext(aes(x=4, y=800), label='r<sup>2</sup>: NA', color='black')+
  geom_richtext(aes(x=4, y=500), label='RMSE: NA', color='black')+
  scale_color_brewer(palette = 'Dark2')+
  theme(legend.position='none')+
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  )

#Plot validation - raymond 2
r2_raymond2 <- round(summary(lm(log10(data_test$k600) ~ log10(data_test$k600_pred_raymond_2)))$r.squared,2)
rmse_raymond2 <- round(10^(Metrics::rmse(log10(data_test$k600), log10(data_test$k600_pred_raymond_2))),2)
test_vali_raymond_2 <- ggplot(data_test, aes(x=(k600), y=(k600_pred_raymond_2), color=data)) +
  geom_abline(size=2, color='grey', linetype='dashed') +
  geom_point(size=5, alpha=0.75) +
  geom_smooth(color='black', method='lm', size=2, se=F) +
  xlab('Observed k600 [m/dy]') +
  ylab('Predicted k600 [m/dy]') +
  geom_richtext(aes(x=4, y=200), label=paste0('r<sup>2</sup>: ', r2_raymond2), color='black')+
  geom_richtext(aes(x=4, y=100), label=paste0('RMSE: ', rmse_raymond2, ' m/dy'), color='black')+
  scale_color_brewer(palette = 'Dark2')+
  theme(legend.position='none')+
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  )

#Plot validation - raymond 3
r2_raymond3 <- round(summary(lm(log10(data_test$k600) ~ log10(data_test$k600_pred_raymond_3)))$r.squared,2)
rmse_raymond3 <- round(10^(Metrics::rmse(log10(data_test$k600), log10(data_test$k600_pred_raymond_3))),2)
test_vali_raymond_3 <- ggplot(data_test, aes(x=(k600), y=(k600_pred_raymond_3), color=data)) +
  geom_abline(size=2, color='grey', linetype='dashed') +
  geom_point(size=5, alpha=0.75) +
  geom_smooth(color='black', method='lm', size=2, se=F) +
  xlab('Observed k600 [m/dy]') +
  ylab('Predicted k600 [m/dy]') +
  geom_richtext(aes(x=4, y=200), label=paste0('r<sup>2</sup>: ', r2_raymond3), color='black')+
  geom_richtext(aes(x=4, y=100), label=paste0('RMSE: ', rmse_raymond3, ' m/dy'), color='black')+
  scale_color_brewer(palette = 'Dark2')+
  theme(legend.position='none')+
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  )

#Plot validation - Ulseth
r2_ulseth <- round(summary(lm(log10(data_test$k600) ~ log10(data_test$k600_pred_ulseth)))$r.squared,2)
rmse_ulseth <- round(10^(Metrics::rmse(log10(data_test$k600), log10(data_test$k600_pred_ulseth))),2)
test_vali_ulseth <- ggplot(data_test, aes(x=(k600), y=(k600_pred_ulseth), color=data)) +
  geom_abline(size=2, color='grey', linetype='dashed') +
  geom_point(size=5, alpha=0.75) +
  geom_smooth(color='black', method='lm', size=2, se=F) +
  xlab('Observed k600 [m/dy]') +
  ylab('Predicted k600 [m/dy]') +
  geom_richtext(aes(x=4, y=300), label=paste0('r<sup>2</sup>: ', r2_ulseth), color='black')+
  geom_richtext(aes(x=4, y=170), label=paste0('RMSE: ', rmse_ulseth, ' m/dy'), color='black')+
  scale_color_brewer(palette = 'Dark2', name='')+
  theme(legend.position='none')+
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  )

legend <- get_legend(test_vali_me + theme(legend.box.margin = margin(0, 0, 0, 90)))
validation <- plot_grid(test_vali_me + theme(legend.position='none') + xlab(''),
                        test_vali_ulseth + xlab('') + ylab(''),
                        test_vali_raymond_2 + xlab(''),
                        test_vali_raymond_3 + ylab(''),
                        test_vali_raymond,
                        legend, ncol=2,
                        labels=c('    RS-able Ulseth etal 2019', '   Ulseth etal 2019', 'Raymond etal 2012- Eq 4','Raymond etal 2012- Eq. 3', 'Raymond etal 2012- Eq 5', NA))
ggsave('outputs/k600/validation.jpg', validation, width=10, height=12)

#retrain model on all data for actual implementation in BIKER algorithm (OOBtrain/test split)----------------------------------------
theory_plot <- ggplot(data, aes(x=(eD), y=(k600), color=factor(widthRegime))) +
  geom_point(size=4, alpha=0.50) +
  geom_smooth(size=3, method='lm', se=F) +
  xlab('eD [m2/s3]') +
  ylab('k600 [m/dy]') +
  scale_color_discrete_qualitative(palette = 'Harmonic', name='River Width [m]', label=c('0-10 (slope < 0.05)', '0-10 (slope > 0.05)','10-50', '50-100', '100+'))+
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'))
ggsave('outputs//k600//theory_plot.jpg', theory_plot, width=9, height=7)

#save model coefficients to be implemented within BIKER algorithm
k_model <- group_by(data, widthRegime) %>%
  do(model = lm(log(k600) ~ log(eD), data=.)) %>%
  mutate(a = exp(model$coefficients[1])) %>%
  mutate(b = model$coefficients[2]) %>%
  mutate(a_se = round(exp(summary(model)$coefficients[2,2]), 3)) %>% #confidence intervals for later uncertainity analysis
  mutate(b_se = round(summary(model)$coefficients[1,2], 3)) %>% #confidence intervals for later uncertainity analysis
  mutate(k600_eD_r2 = round(summary(model)$r.squared, 2))
k_model_groups <- group_by(data, widthRegime) %>%
  summarise(n = n())

k_model <- left_join(k_model, k_model_groups, by='widthRegime')
data <- left_join(data,k_model, by="widthRegime")

#Other models trained on full dataset-------------------------------------
#Other models: Ulseth and raymond models using only V and S
lm_raymond <- lm(k600~raymond, data=data) #raymond model- linear, eq 5
lm_raymond_2 <- lm(log_k600~log(raymond), data=data)#raymond model 2- power law, eq 4
lm_raymond_3 <- lm(log_k600~log(slope)+log(Vms), data=data)#raymond model 3- multi power law, eq 3

#Ulseth model
linearModel <- lm(log_k600~log_eD, data=data)
lm_ulseth <- segmented(linearModel, npsi = 1)

#k600 prior pdfs used within BIKER algorithm--------------------------------------
#Hagemann et al. (2017)-style prior using slope. IMPLEMENTED.
lmSlope <- lm(log_k600~log_slope, data=data)
segLmSlope <- segmented(lmSlope, npsi=2) #r2 0.55
slope(segLmSlope)
intercept(segLmSlope)
summary(segLmSlope)

#geoBAM-style prior. NOT IMPLEMENTED.
data$slopeRegime <- ifelse(data$slope <= 0.0005, 1,
                           ifelse(data$slope <= 0.005, 2,
                                  ifelse(data$slope <= 0.05, 3,4)))

biker_k_prior <- group_by(data, slopeRegime) %>%
  summarise(k600_hat = median(k600, na.rm = T),
            k600_hat_mean = mean(k600, na.rm=T),
            k600_sd = sd(k600, na.rm=T),
            lowerbound_k600 = min(k600, na.rm=T),
            upperbound_k600 = max(k600, na.rm=T),
            n=n())
write.csv(biker_k_prior, 'outputs//k600//k_priors.csv')

#Figure S4: bed roughness vs. k600-------------------------------------------
data$keulegan <- 11*data$depth*(1/exp(data$Vms/(2.5*(9.8*data$depth*data$slope)^(1/2))))

lm_roughness <- lm(log10(data$k600)~log10(data$keulegan))
r2_roughness <- round(summary(lm_roughness)$r.squared, 2)
roughness_0 <- ggplot(data, aes(x=keulegan, y=k600)) +
  geom_point() +
  geom_smooth(method='lm', se=F)+
  geom_richtext(aes(x=10^-1.5, y=10^3), label=paste0('r<sup>2</sup>: ', r2_roughness), color='black') +
  coord_cartesian(x=c(10^-2.5, 10^1.5), c(10^-2, 10^5)) +
  ylab('K600 [m/dy]') +
  xlab('') +
  ggtitle('All Data')+
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  )

lm_roughness <- lm(log10(k600)~log10(keulegan), data = filter(data, slope > 0.001))
r2_roughness <- round(summary(lm_roughness)$r.squared, 2)
roughness_1 <- ggplot(filter(data, slope > 0.001), aes(x=keulegan, y=k600)) +
  geom_point() +
  geom_smooth(method='lm', se=F)+
  geom_richtext(aes(x=10^-1.5, y=10^3), label=paste0('r<sup>2</sup>: ', r2_roughness), color='black') +
  coord_cartesian(x=c(10^-2.5, 10^1.5), c(10^-2, 10^5)) +
  ylab('') +
  xlab('') +
  ggtitle('Slope > 0.001')+
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  )

lm_roughness <- lm(log10(k600)~log10(keulegan), data = filter(data, slope > 0.1))
r2_roughness <- round(summary(lm_roughness)$r.squared, 2)
roughness_2 <- ggplot(filter(data, slope > 0.1), aes(x=keulegan, y=k600)) +
  geom_point() +
  geom_smooth(method='lm', se=F)+
  geom_richtext(aes(x=10^-1.5, y=10^3), label=paste0('r<sup>2</sup>: ', r2_roughness), color='black') +
  coord_cartesian(x=c(10^-2.5, 10^1.5), c(10^-2, 10^5)) +
  ylab('K600 [m/dy]') +
  xlab('Effective Bed Roughness height [m]') +
  ggtitle('Slope > 0.1')+
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  )

lm_roughness <- lm(log10(k600)~log10(keulegan), data = filter(data, slope > 0.2))
r2_roughness <- round(summary(lm_roughness)$r.squared, 2)
roughness_3 <- ggplot(filter(data, slope > 0.2), aes(x=keulegan, y=k600)) +
  geom_point() +
  geom_smooth(method='lm', se=F)+
  geom_richtext(aes(x=10^-1.5, y=10^3), label=paste0('r<sup>2</sup>: ', r2_roughness), color='black') +
  coord_cartesian(x=c(10^-2.5, 10^1.5), c(10^-2, 10^5)) +
  ylab('') +
  xlab('Effective Bed Roughness height [m]') +
  ggtitle('Slope > 0.2')+
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  )

roughness_plot <- plot_grid(roughness_0, roughness_1, roughness_2, roughness_3, ncol=2)
ggsave('outputs/k600/roughness.jpg', roughness_plot, width=9, height=8)

#Save results to file for manuscript-----------------------------------
results <- data.frame('r2'=c(r2_me, NA, r2_raymond2, r2_raymond3, r2_ulseth),
                      'rmse'=c(rmse_me, NA, rmse_raymond2, rmse_raymond3,rmse_ulseth),
                      'Model'=c('RS-able', 'Raymond_1', 'Raymond_2', 'Raymond_3','Ulseth'))
write.csv(results, 'outputs//k600//results.csv')


#non linear function tests
k_model_poly <- lm(log(k600) ~ poly(log(eD), 2), data=data)
data$k600_pred_poly <- predict(k_model_poly, data)
theory_plot2 <- ggplot(data, aes(x=(eD), y=(k600))) +
  geom_point(size=4, alpha=0.50) +
  #geom_smooth(size=3, se=F, method='lm') +
  geom_line(aes(x=eD, y=exp(k600_pred_poly)), color='blue', size=3) +
  xlab('eD [m2/s3]') +
  ylab('k600 [m/dy]') +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'))
ggsave('outputs//k600//theory_plot2.jpg', theory_plot2, width=9, height=7)

print(summary(k_model_poly))
