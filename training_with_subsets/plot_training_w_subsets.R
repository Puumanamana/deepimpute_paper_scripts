library("ggplot2")
library("scales")

data <- read.csv('../results/training_w_subsets/metrics_with_increasing_fraction.csv',
                 check.names=F)

m_MSE <- min(data$MSE)
M_MSE <- max(data$MSE)
m_pears <- min(data$Pearson)
M_pears <- max(data$Pearson)

data$Pearson <- m_MSE+(M_MSE-m_MSE)*(data$Pearson-m_pears)/(M_pears-m_pears)

df <- do.call(data.frame,aggregate(. ~ Fraction, data, function(x) c(mean=mean(x), sd=sd(x))))

ggplot(data=df, aes(x=Fraction)) +
  geom_line(aes(y=MSE.mean,colour='MSE')) +
  geom_errorbar(aes(ymin=MSE.mean-MSE.sd, ymax=MSE.mean+MSE.sd,colour='MSE')) +
  
  geom_line(aes(y=Pearson.mean,colour='Pearson')) +
  geom_errorbar(aes(ymin=Pearson.mean-Pearson.sd, ymax=Pearson.mean+Pearson.sd,colour='Pearson')) +
  scale_y_continuous(sec.axis = sec_axis(~m_pears+(.-m_MSE)*(M_pears-m_pears)/(M_MSE-m_MSE), name = "Pearson")) +
  
  guides(linetype = "none", shape = "none", line = "none", point = "none") +
  theme_classic(base_size = 16) +
  theme(legend.title=element_blank()) +
  labs(x="Fraction of the data for model training", y="Mean Squared Error", color='')

