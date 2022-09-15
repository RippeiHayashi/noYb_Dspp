### this script is to plot the in trans pingpong linkage values

### R libraries
library(ggplot2)
library(reshape2)
pdf.options(useDingbats = FALSE)

### simulation
table <- read.table("simulated.numbers.1Mio.stats", header=F)
colnames(table) <- c("linkage","library")

table$library <- factor(table$library,levels = c("0percent", "1percent", "3percent", "5percent", "10percent"))

pdf(file="in-trans-pingpong_simulation_2022-09-08.pdf", height=6, width=6)
ggplot(table, aes(x=library, y=linkage)) +
  ylim(c(0,20))+
  geom_jitter(size=2,shape=1, position=position_jitter(0.2)) +
  labs(title="in-trans ping-pong, simulation", y="linkage")
dev.off()

### small RNA data
table <- read.table("in-trans-pingpong.stats", header=F)
colnames(table) <- c("type","linkage","library")

pdf(file="in-trans-pingpong_mouse_Deug_Dpse_w1118_2022-09-08.pdf", height=7, width=10)
ggplot(table, aes(x=library, y=linkage, fill=type)) +
  geom_bar(stat="identity", position="dodge")+
  ylim(c(0,25))+
  labs(title="in-trans ping-pong, all-mappers and genome-unique mappers", y="linkage")+
  theme(axis.text.x=element_text(angle = -90, hjust = 0))
dev.off()
