
# plot ECG averages

setwd("D:/Psycho/cikkek/Saját/ANT_HEP")

ecg_data=read.csv("ecg_data.csv",sep=";")



library(ggplot2)


# Raincloud plot example
library(ggdist)
pt=ggplot(ecg_data, aes(x = cond_pt, y = pt, fill = cond_pt)) +
  geom_boxplot(
    width = 0.15,
    alpha = 0.5
  ) +
  scale_fill_manual(values = c("blue", "red")) +  # Set fill colors for the boxplot
  geom_point(aes(color = cond_pt), size = 3) +  # Set point colors
  scale_color_manual(values = c("blue", "red")) +  # Define colors for points
  geom_line(aes(group = ID), 
            color = "black", 
            alpha = 0.5) +  # Connect points by ID with lines
  theme_bw() +
  ylab("ECG Amplitude (424-600 ms)") +
  xlab("") +
  theme(axis.title = element_text(size = 28)) +
  theme(legend.title = element_blank()) + 
  theme(legend.text = element_text(size = 25)) +
  ylim(-25, 25) +
  theme(axis.text.x = element_blank())+theme(axis.text.y=element_text(size=15))


pw=ggplot(ecg_data, aes(x = cond_pw, y = pw, fill = cond_pw)) +
  geom_boxplot(
    width = 0.15,
    alpha = 0.5
  ) +
  scale_fill_manual(values = c("blue", "yellow")) +  # Set fill colors for the boxplot
  geom_point(aes(color = cond_pw), size = 3) +  # Set point colors
  scale_color_manual(values = c("blue", "yellow")) +  # Define colors for points
  geom_line(aes(group = ID), 
            color = "black", 
            alpha = 0.5) +  # Connect points by ID with lines
  theme_bw() +
  ylab("ECG Amplitude (373-555 ms)") +
  xlab("") +
  theme(axis.title = element_text(size = 28)) +
  theme(legend.title = element_blank()) + 
  theme(legend.text = element_text(size = 25)) +
  ylim(-25, 25) +
  theme(axis.text.x = element_blank())+theme(axis.text.y=element_text(size=15))
 
require(gridExtra)

ecg_plots=grid.arrange(pt,pw,ncol=2)

ggsave("ecg_avg.pdf",ecg_plots,dpi=600, width = 20, height = 10)

######################## ant-hep plots

ant_data=read.csv("ant_data.csv",sep=";")

pt_ant=ggplot(ant_data, aes(x = cond_pt, y = pt, fill = cond_pt)) +
  geom_boxplot(
    width = 0.15,
    alpha = 0.5
  ) +
  scale_fill_manual(values = c("blue", "red")) +  # Set fill colors for the boxplot
  geom_point(aes(color = cond_pt), size = 3) +  # Set point colors
  scale_color_manual(values = c("blue", "red")) +  # Define colors for points
  geom_line(aes(group = ID), 
            color = "black", 
            alpha = 0.5) +  # Connect points by ID with lines
  theme_bw() +
  ylab("HEP Amplitude (424-600 ms)") +
  xlab("") +
  theme(axis.title = element_text(size = 28)) +
  theme(legend.title = element_blank()) + 
  theme(legend.text = element_text(size = 25)) +
  ylim(-1, 5) +
  theme(axis.text.x = element_blank())+theme(axis.text.y=element_text(size=15))


pw_ant=ggplot(ant_data, aes(x = cond_pw, y = pw, fill = cond_pw)) +
  geom_boxplot(
    width = 0.15,
    alpha = 0.5
  ) +
  scale_fill_manual(values = c("blue", "yellow")) +  # Set fill colors for the boxplot
  geom_point(aes(color = cond_pw), size = 3) +  # Set point colors
  scale_color_manual(values = c("blue", "yellow")) +  # Define colors for points
  geom_line(aes(group = ID), 
            color = "black", 
            alpha = 0.5) +  # Connect points by ID with lines
  theme_bw() +
  ylab("HEP Amplitude (373-555 ms)") +
  xlab("") +
  theme(axis.title = element_text(size = 28)) +
  theme(legend.title = element_blank()) + 
  theme(legend.text = element_text(size = 25)) +
  ylim(-1, 5) +
  theme(axis.text.x = element_blank())+theme(axis.text.y=element_text(size=15))

require(gridExtra)

ant_plots=grid.arrange(pt_ant,pw_ant,ncol=2)

ggsave("hep_avg.pdf",ant_plots,dpi=600, width = 20, height = 10)

