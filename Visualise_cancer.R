#run with Rscript Visualise.R basenames NSAMS

rm(list=ls())
graphics.off()
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(zoo)
args = commandArgs(trailingOnly=TRUE)

input = args[1]
file_list <- unlist(read.table(input, header=FALSE, as.is=T))
basename_files <- c()
genolike_files <- c()
ploids_files <- c()
out_plot <- c()


#mode function
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

# Extract file names and build output names
for(i in 1:length(file_list)){ ###
  genolike_files[i] <- paste(file_list[i],".genolikes.gz",sep="")
  splittedName <- unlist(strsplit(file_list[i],split="/"))
  basename_files[i] <- splittedName[length(splittedName)]
  directory <- paste(splittedName[-length(splittedName)],collapse = '/')
  out_plot[i] <- paste(file_list[i],"_plot.pdf",sep="")
  ploids_files[i] <- paste(file_list[i],".ploids",sep="")
} 
###### Note NSAMS fixed for project work only
NSAMS=strtoi(args[2])
sam_depths <- vector("list",NSAMS)
ploidy <- vector("list",NSAMS)
length_of_samples <- c()

q=0
Exp_Ploidies <- vector("list",NSAMS*length(file_list))
for(i in 1:length(file_list)){  ###
  q=q+1
  
  ploidies = ploids_files[i]
  output = out_plot[i]
  

  Ploidies=read.csv(ploidies,sep="\t",header=FALSE) # load the ploidies
  sams<-head(Ploidies,1)
  NSAMS=length(Ploidies)-sum(is.na(sams))
  Expected_Ploidies<-tail(Ploidies,NSAMS)
  Expected_Ploidies <- Expected_Ploidies[,colSums(is.na(Expected_Ploidies))<nrow(Expected_Ploidies)]
  number_of_windows <- length(Expected_Ploidies[3,])
  Expected_Ploidies <- melt(as.matrix(Expected_Ploidies))
  Expected_Ploidies$Var2<-as.numeric(Expected_Ploidies$Var2)
  Ploidies<-head(Ploidies,1)
  Ploidies <- Ploidies[,colSums(is.na(Ploidies))<nrow(Ploidies)]
  
  

  
  plot <- ggplot(data=Expected_Ploidies, aes(x = Var2, y = Var1)) + ggtitle("Window analysis of ploidies") +
    geom_tile(aes(fill=factor(value))) + 
    scale_fill_manual(values=c("blue","red","green","orange","magenta","cyan","purple","yellow","brown","black","red")) +
    labs(x="Base position", y="Sample") + scale_y_continuous('Sample',limits=c(2.5,(NSAMS+2.5)),breaks = c(3:(NSAMS+2)),labels = c(1:NSAMS))+
    guides(fill=guide_legend("Window Ploidy"))+
    scale_x_continuous('Window',limits = c(0.5,(number_of_windows+0.5)),breaks=seq(0.5,number_of_windows+0.5,length.out=11),labels=round(seq(0,number_of_windows,length.out=11)))+
    theme(legend.position = "right",plot.title = element_text(hjust = 0.5,size = 20,face="bold"),axis.text=element_text(size=14),axis.title=element_text(size=16,face="bold"))
  

  pdf(output, 11.7, 8.3)
  plot(plot)
  graphics.off()
}
