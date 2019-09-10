library("ggplot2")
library("scales")

data <- read.csv('../results/speed_memory/speed_memory.csv',row.names=1)
colors <- c("#3E63B2","#FF9966","#9933FF","#FF0000","#009900","#FF66FF","#996633")
names(colors) <- c("deepImpute","DCA","scImpute","SAVER","MAGIC","DrImpute","VIPER")

## Speed chart
speed_data <- data[data$variable=='runTime',]
ggplot(data=speed_data, aes(x=cellCount,y=value,colour=speed_data$packageName)) +
  geom_point(size=3) +
  geom_line(size=1) +
  guides(linetype = "none", shape = "none", line = "none", point = "none") +
  ## scale_y_log10(limits=c(0.06,1000),labels=rescale_none) +
  scale_y_log10(limits=c(0.06,1000),labels=rescale_none) +
  scale_x_continuous(labels=comma) +
  theme_classic(base_size = 16) +
  theme(legend.title=element_blank()) +
  labs(x="Cell count", y="Duration in minutes", color='') +
  scale_color_manual(breaks=names(colors),values=unname(colors)) 

  ## Memory chart
  memory_data <- data[data$variable=='peakMem',]
  ggplot(data=memory_data, aes(x=cellCount,y=value,colour=memory_data$packageName)) +
    geom_point(size=3) +
    geom_line(size=1) +
    guides(linetype = "none", shape = "none", line = "none", point = "none") +
    theme_classic(base_size = 16) +
    theme(legend.title=element_blank()) +
    labs(x="Cell count", y="Peak memory", color='') +
    scale_x_continuous(labels=comma) +
    scale_y_continuous(labels=comma) +
    scale_color_manual(breaks=names(colors),values=unname(colors)) +
    coord_cartesian(ylim=c(0, 30)) +
    geom_hline(yintercept = 30, linetype = "dotted") +
    # theme(legend.title=element_blank()) +
    annotate("text", 24000, 30, vjust = -.5, label = "Test Machine Memory Limit")
  
