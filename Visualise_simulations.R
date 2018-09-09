#run with Rscript Visualise.R basenames

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
win<-strtoi(args[2])

# Extract file names and build output names
for(i in 1:length(file_list)){
  genolike_files[i] <- paste(file_list[i],".genolikes.gz",sep="")
  splittedName <- unlist(strsplit(file_list[i],split="/"))
  basename_files[i] <- splittedName[length(splittedName)]
  out_plot[i] <- paste(file_list[i],"_plot.pdf",sep="")
  ploids_files[i] <- paste(file_list[i],".ploids",sep="")
  
} 
directory<-paste0(paste(head(splittedName,-1),collapse = "/"),"/")
#sample_list<-levels(unlist(read.table(paste0(directory,"bamlist.txt"))))
#for(s in 1:length(sample_list)){
#  sample_list[s]<- unlist(strsplit(sample_list[s],split=".indel"))[1]
#}


for(i in 1:length(file_list)){                 
  genolikes = genolike_files[i]
  ploidies = ploids_files[i]
  output = out_plot[i]
  
  Genolikes=read.csv(gzfile(genolikes),sep="\t",header=FALSE) # Load the data
  Ploidies=read.csv(ploidies,sep="\t",header=FALSE) # load the ploidies
  #extract the number of smaples and extract expected ploidies removing NA values
  sams<-head(Ploidies,1)
  NSAMS=length(Ploidies)-sum(is.na(sams))
  Expected_Ploidies<-tail(Ploidies,NSAMS)
  Expected_Ploidies <- Expected_Ploidies[,colSums(is.na(Expected_Ploidies))<nrow(Expected_Ploidies)]
  Ploidies<-head(Ploidies,2)
  Ploidies <- Ploidies[,colSums(is.na(Ploidies))<nrow(Ploidies)]
  
  Depths=Genolikes$V5
  length_of_sample=length(Depths)/NSAMS
  Ploidies_to_sum=Ploidies[1,]
  Ploidies_to_sum[is.na(Ploidies_to_sum)]<-0
  Meandepth=sum(Depths)/(length_of_sample*sum(Ploidies_to_sum))
  # create a vector of lists of depths for each sample
  positions<-c()
  sample_depths <- vector("list",NSAMS)
  for(j in 0:(NSAMS-1)){
    p=1
    for(i in 1:length(Depths)){
      positions[p]<-Genolikes$V2[i]
      if(i%%NSAMS==j){
        if(j==0){
          sample_depths[[NSAMS]][p]<-Depths[i]
          
          p<-p+1
          
        }else{
          sample_depths[[j]][p]<-Depths[i]
          p<-p+1}
      }
    }
  }
  sample_mean_depth<-c()
  for(sample in 1:NSAMS){
    sample_depths[[sample]]<-sample_depths[[sample]][sample_depths[[sample]]!=0]
    sample_depths[[sample]]<-sample_depths[[sample]][sample_depths[[sample]]<quantile(sample_depths[[sample]],1)]
    #sample_depths[[sample]]<-rollmean(sample_depths[[sample]],length_of_sample/100,fill=list(0,0,0)) # take rolling average 
    sample_depths[[sample]]<-sample_depths[[sample]][seq(from=1,to = length(sample_depths[[sample]]) ,length.out = 1000)] # take every len/1000  value
    sample_mean_depth[sample]<-mean(sample_depths[[sample]])  
  }
  positions<-signif((seq(from=0,to=max(positions),length.out = 11)/100),2)
  depths <- as.data.frame(sample_depths[[1]])
  
  for(i in c(2:NSAMS)){
    depths<-cbind(depths,sample_depths[[i]])
  }
  col_names=c(1:NSAMS)
  colnames(depths)<-col_names
  
  
  depths <-melt(depths)
  
  ploidy<-c()
  
  for(j in Ploidies){
    ploidy<- c(ploidy,rep(j[1],1000))
    
    
  }
  depths <- cbind(depths,ploidy)
  
  # Add the window by window analysis
  number_of_windows<-length(Expected_Ploidies[1,])
  normalised_window_length<-floor(1000/number_of_windows)
  ex<-1000-(number_of_windows*normalised_window_length)
  expected_ploidy<-c()
  ovr_mean_hap=sum(depths$value)/sum(depths$ploidy)
  ovr_mean=mean(depths$value)
  plot_mean_depths<-c()
  for(n in 1:NSAMS){
    
    for(j in Expected_Ploidies[n,]){
      expected_ploidy<- c(expected_ploidy,rep(j,normalised_window_length))
    }
    expected_ploidy<-c(expected_ploidy,rep(Expected_Ploidies[n,length(Expected_Ploidies)],ex))
    plot_mean_depths<-c(plot_mean_depths,c(ploidy[(((n-1)*1000)+1):(n*1000)]*ovr_mean_hap))
  }
  
  depths <- cbind(depths,expected_ploidy)
  depths <- cbind(depths,plot_mean_depths)
  sample_mean_depth_long<-rep(sample_mean_depth,each=1000)
  norm_value<-depths$value*(sample_mean_depth_long/ovr_mean)
  depths<-cbind(depths,norm_value)
  
  
  myColors <- mycols <- colors()[c(12,414,576,573,524,436,106,74,75,86,99,137,627,656,367,419,81,410,512,402,468,592,535,429,404,477,50,79,102,20,101,52,51,24,134,616)]
  names(myColors) <- levels(depths$variable)
  
  new_colours=c("blue","green","orange","magenta","cyan","purple","yellow","brown","black","red")
  
  

  pdf(paste0(output), 8, 5)
  temp=length(unique(depths$expected_ploidy))
  colScale<-scale_colour_manual(values = new_colours[c(c(1:temp),9,10)],
                                  guide = guide_legend(override.aes = list(
                                    linetype = c(rep("blank", temp), "solid","dashed"),
                                    shape = c(rep(1, temp), NA,NA))))
    
    
  plot<-ggplot(data = depths[c(1:(NSAMS*1000)),]) + xlim(0,1000*NSAMS) +ylim(0,quantile(Depths,0.8)) # plot axis
  plot <- plot + theme(legend.position = "top",plot.title = element_text(hjust = 0.5,size = 20,face="bold"),axis.text=element_text(size=14),axis.title=element_text(size=16,face="bold"),legend.text = element_text(size=10),legend.title = element_text(size=12)) # sprt out axis
    
  #plot <- plot + geom_point(data=depths[c((((sample-1)*1000)+1):(sample*1000)),],aes(x=c(1:(1000)),y=value,colour=factor(expected_ploidy)),size=2,alpha=1/3)  # plot depths
  plot <- plot + geom_point(data=depths[c(1:(NSAMS*1000)),],aes(x=c(1:(1000*NSAMS)),y=value,colour=factor(expected_ploidy)),size=2,alpha=1/6)  # plot depths
  plot <- plot + geom_line(data= depths[c(1:(NSAMS*1000)),],aes(x=c(1:(1000*NSAMS)),y=depths$plot_mean_depths,colour="Chromosome inferred ploidy"),size=1)
  plot <- plot + geom_line(data= depths[c(1:(NSAMS*1000)),],aes(x=c(1:(1000*NSAMS)),y=rep(sample_mean_depth,each=1000),colour="Mean Depth"),size=0.4,linetype="dashed") +
    guides(color=guide_legend("Localised Ploidy",override.aes = list(linetype = c(rep("blank", temp), "solid","dashed"),size=c(rep(2,temp),1,0.5),shape = c(rep(16, temp), NA,NA),alpha=2/3)))+colScale
  plot <- plot + facet_grid(.~depths$variable,scales = "free")#,shrink=FALSE,drop=FALSE)
  plot <- plot + ggtitle(paste0("Inferred ploidies vs depth ")) + ylab("Depth") # add titles
    
  #plot <- plot + geom_line(data = depths[c((((sample-1)*100)+1):(sample*100)),],aes(x=c(1:100),y=expected_ploidy*sample_mean_depth[sample]), size=0.5,colour='red') # plot inferred ploidy
  plot <- plot +scale_y_continuous('Depth',limits=c(0,ovr_mean_hap*(max(depths$ploidy)+4)),sec.axis = sec_axis(~./ovr_mean_hap, name = 'Inferred Ploidy',breaks = c(0:(max(ploidy)+1)))) # add 2nd y axis
  plot <- plot + scale_x_continuous(name = "Sample",breaks = (c(1:NSAMS)*1000)-500,labels = c(1:NSAMS)) # change scale on x axis to samples
    
  plot(plot)
  graphics.off()

}

