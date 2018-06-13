#Plotting Python Sim results
#Load Sim Results
library(readr)
library(ggplot2)
#simresults <- read_csv("~/git_repo/PyBASC/SimResults/SimResults_TR400.csv")


#IBS_GBS NIFTI ANALYSIS
IBS_GBS_Sim<- read_csv('/Users/aki.nikolaidis/git_repo/PyBASC/Updated_IBS_GBS_Jovo4_20_2018.csv')

IBS_GBS_Sim$IBS<- as.factor(IBS_GBS_Sim$IBS)
ggplot(IBS_GBS_Sim, aes(x=GBS)) + 
  #scale_x_reverse() +
  geom_line(aes(y=ism_gsm_corrstd, col=IBS), size=1) + 
  scale_y_continuous(limits = c(0,1)) +
  #  ggtitle("Clustering with Bootstraps vs None") +
  #  ggsubtitle("Bootstrap stability approaches true adjacency matrix")+
  #  ylab("Accuracy") +
  labs(title="IBS & GBS Improve Reduce Variance in Individual-Group Coherence", 
       #subtitle="200 vols- Bootstraps most important for lower SNR ranges" , 
       x="GBS",
       y="ism_gsm_corrstd")


#IBS_Sim Plotting
IBS_Sim <- read_csv('/Users/aki.nikolaidis/git_repo/PyBASC/SimResultsIBS/PaperSimResults_IBS_Noise.csv')
IBS_Sim$SNR<- as.factor(IBS_Sim$SNR)

ggplot(subset(IBS_Sim, TRs==200), aes(x=bootstraps)) + 
  #scale_x_reverse() +
  geom_line(aes(y=AvgAcc, col=SNR), size=1) + 
  scale_y_continuous(limits = c(0,1)) +
  #  ggtitle("Clustering with Bootstraps vs None") +
  #  ggsubtitle("Bootstrap stability approaches true adjacency matrix")+
  #  ylab("Accuracy") +
  labs(title="Individual Bootstraps improves individual cluster accuracy", 
       subtitle="200 vols- Bootstraps most important for lower SNR ranges" , 
       x="bootstraps",
       y="AvgAcc")

ggplot(subset(IBS_Sim, TRs==400), aes(x=bootstraps)) + 
  #scale_x_reverse() +
  geom_line(aes(y=AvgAcc, col=SNR), size=1) + 
  scale_y_continuous(limits = c(0,1)) +
  #  ggtitle("Clustering with Bootstraps vs None") +
  #  ggsubtitle("Bootstrap stability approaches true adjacency matrix")+
  #  ylab("Accuracy") +
  labs(title="Individual Bootstraps improves individual cluster accuracy", 
       subtitle="400 Vols- Bootstraps most important for lower SNR ranges" , 
       x="bootstraps",
       y="AvgAcc")

ggplot(subset(IBS_Sim, TRs==800), aes(x=bootstraps)) + 
  #scale_x_reverse() +
  geom_line(aes(y=AvgAcc, col=SNR), size=1) + 
  scale_y_continuous(limits = c(0,1)) +
  #  ggtitle("Clustering with Bootstraps vs None") +
  #  ggsubtitle("Bootstrap stability approaches true adjacency matrix")+
  #  ylab("Accuracy") +
  labs(title="Individual Bootstraps improves individual cluster accuracy", 
       subtitle="800 Vols- Bootstraps most important for lower SNR ranges" , 
       x="bootstraps",
       y="AvgAcc")

multiplot(tr200, tr400, tr800, cols=3)
# plot

ggplot(simresults, aes(x=GBS)) + 
  #scale_x_reverse() +
  geom_line(aes(y=group_label_acc, col=IBS), size=1) + 
  scale_y_continuous(limits = c(0,1)) +
#  ggtitle("Clustering with Bootstraps vs None") +
#  ggsubtitle("Bootstrap stability approaches true adjacency matrix")+
#  ylab("Accuracy") +
  labs(title="Without bootstrapping- less accurate cluster labels", 
     #subtitle="Across all SNR levels- Individual Bootstrapping helps" , 
     x="GBS",
     y="group_label_acc")
    
  

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}





  labs(title="Clustering with Bootstraps vs None", 
     subtitle="Showing stability approaches the true adjacency matrix", 
     caption="Source: PyBASC", 
     color=NULL) +  # title and caption
  #scale_x_date(labels = lbls, breaks = brks) +  # change to monthly ticks and labels
  scale_color_manual(labels = c("psavert", "uempmed"), 
                     values = c("psavert"="#00ba38", "uempmed"="#f8766d")) +  # line color
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, size = 8),  # rotate x axis text
        panel.grid.minor = element_blank())  # turn off minor grid