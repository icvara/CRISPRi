###############################################################################
###############################################################################
##    this script is made to analyse time course data from plate reader      ##
## 											           ##
## 				made by Içvara v 2.0			  	           ##
###############################################################################
###############################################################################


#######################
#####Variable PART#####
#######################
path<-"C:\\Users\\Administrator\\Desktop\\Modeling\\CRISPRi\\data"

rawDataFileName = c("OD_A","GFP_A","mko2_A")
rawDataFileName2 = c("OD_B","GFP_B","mko2_B")
rawDataFileName3 = c("OD_C","GFP_C","mko2_C")



timeDirection="Y"



#######################
#####FUNCTION PART#####
#######################
library("plyr")


GetNiceData<-function(rawData,rawLayout,timeDirection){
	niceLayout<-GetLayout(rawLayout)
	niceData<-GetData(rawData)
	dataReady<-AssembleData(rawData,niceData,niceLayout,timeDirection)
	return(dataReady)
}

#put the layout in one columns to fit the dataframe structure

GetLayout<-function(rawLayout){
	mLayout<-as.matrix(rawLayout)
	splitLayout<-strsplit(mLayout,",")
	niceLayout<-list()
	name<-c()
	l<-c()
	c=FALSE
	for(y in 1:length(mLayout)){
		if(splitLayout[[y]][1]=="LAYOUT"){
			if(c==TRUE){
				niceLayout[[name]]<-l
				name=splitLayout[[y]][2]
				l<-c()
			}
			if(c==FALSE){
				name=splitLayout[[y]][2]
				c=TRUE
				l<-c()
			}	
		}		
		if(splitLayout[[y]][1]!="LAYOUT"){
			l<-c(l,splitLayout[[y]])
			if( y==length(mLayout)){
				niceLayout[[name]]<-l
				l<-c()
			}
		}	
	}
	return(niceLayout)
}

GetData<-function(rawData){
	niceData<-list()
	for(di in 1:length(rawData)){
		lineData<-c()
		for(i in 1:ncol(rawData[[di]])){
		 lineData<-c(lineData,rawData[[di]][,i])
		}
		niceData[[names(rawData)[di]]]<-lineData
	}
	return(niceData)
}


AssembleData<-function(rawData,niceData,niceLayout,timeDirection){
	adjustedLayout<-list()
	if(timeDirection=="X"){
		nTime=ncol(rawData[[1]])
		for(i in 1:length(niceLayout)){
			adjustedLayout[[names(niceLayout)[i]]]<-rep(niceLayout[[i]],time=nTime)
		}
		time=rep(1:nTime,each=nrow(rawData[[1]]))
	}
	if(timeDirection=="Y"){
		nTime=nrow(rawData[[1]])
		for(i in 1:length(niceLayout)){
			adjustedLayout[[names(niceLayout)[i]]]<-rep(niceLayout[[i]],each=nTime)
		}
		time=rep(1:nTime,time=ncol(rawData[[1]]))

	}
	return(finalData<-data.frame(time,adjustedLayout,niceData))	
}

LoadLayout<-function(path,file){
	rawLayout<-read.table(paste(path, file,sep="\\"),h=F)
	return(rawLayout)
}

LoadData<-function(path,rawDataFileName ){
	rawData<-list()
	for( name in rawDataFileName ){
		rd<-read.table(paste(paste(path, name ,sep="\\"),".txt",sep=""),h=F)
		rawData[[name]]<-data.frame(rd)
	}
	return(rawData)
}

GetLine<-function(data,index){
	l<-ncol(data)
	f<-data
	for( i in 1:l){
		if(index[i] != "NA"){
			f<-subset(f, f[,i]==index[i] )
		}	
	}
	return(f)
}

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

#######################
#####Call Function#####
#######################

rawLayout<-LoadLayout(path,"Layout.csv")
##rawLayout<-LoadLayout(path,"Layout_corrected.csv")

rawData<-LoadData(path,rawDataFileName )
data1<-GetNiceData(rawData,rawLayout,timeDirection)

rawData<-LoadData(path,rawDataFileName2 )
data2<-GetNiceData(rawData,rawLayout,timeDirection)

rawData<-LoadData(path,rawDataFileName3 )
data3<-GetNiceData(rawData,rawLayout,timeDirection)


names(data1)<-c("time","sample","arabinose","OD","GFP","mko2")
names(data2)<-c("time","sample","arabinose","OD","GFP","mko2")
names(data3)<-c("time","sample","arabinose","OD","GFP","mko2")

data1$rep="A"
data2$rep="B"
data3$rep="C"


fulldata<-rbind(data1,data2,data3)

##################3
#quick visualisation

ggplot(fulldata, aes(y=OD,x= time, linetype=rep, col=arabinose)) + geom_line()+
facet_grid(sample~.)

###########################


##cas 3 neg sub
D<-fulldata

blank<-subset(D,sample=="blank")

neg<-subset(D,sample=="empty"  )
neg$GFP<-(neg$GFP-mean(blank$GFP))/(neg$OD-mean(blank$OD)) 
neg$mko2<-(neg$mko2-mean(blank$mko2))/(neg$OD-mean(blank$OD)) 



meanNEG<-ddply(neg, c("sample", "arabinose","time"), summarise,
      mean_GFP = mean(GFP), sd_GFP = sd(GFP),
      mean_mko2 = mean(mko2), sd_mko2 = sd(mko2),
      mean_OD = mean(OD), sd_OD = sd(OD) )


D<-subset(D,time==100)
neg<-subset(meanNEG,time==100)

D$GFP<-(D$GFP-mean(blank$GFP))/(D$OD-mean(blank$OD)) - neg$mean_GFP # (neg$GFP-mean(blank$GFP))/(neg$OD-mean(blank$OD))
D$mko2<-(D$mko2-mean(blank$mko2))/(D$OD-mean(blank$OD))

yl="Fluo-blk / OD-blk "
#yl="Fluo-blk / OD-blk - empty-blk / OD -blk "
 "blank" & sample != "empty" & sample != "blank2")

#mean,SD
meanD<-ddply(D, c("sample", "arabinose","time"), summarise,
      mean_GFP = mean(GFP), sd_GFP = sd(GFP),
      mean_mko2 = mean(mko2), sd_mko2 = sd(mko2),
     # mean_mkate2 = mean(mkate2), sd_mkate2 = sd(mkate2),
      mean_OD = mean(OD), sd_OD = sd(OD) )

names(D)
   


###################

#######################
#####Graph Part   #####
#######################

library("ggplot2")




d2<-subset(D, sample != "blank" & sample!= "empty" & sample !="NA")
#plot
ggplot(d2,aes(y=GFP,x=time, linetype=rep, col=arabinose))+
geom_line(size=1)+
ylab(yl)+
xlab("time")+
theme(text = element_text(size=20), axis.text.x = element_text(angle=90, hjust=1))+
facet_grid(sample ~ .)


d2<-subset(meanD, sample != "blank" & sample!= "empty" & sample !="NA")
#plot
ggplot(d2,aes(y=mean_GFP,x=time, col=arabinose))+
geom_line(size=1)+
ylab(yl)+
xlab("time")+
theme(text = element_text(size=20), axis.text.x = element_text(angle=90, hjust=1))+
facet_grid(sample ~ .)



+
#theme_classic()+
facet_grid(scRNA ~ bsA)+
theme(strip.text.x = element_text(size=8))



d2<-subset(meanD, sample != "blank" & sample!= "empty" & sample !="NA" & time==100)

ggplot(d2,aes(y=mean_GFP,x=sample, fill=arabinose))+ 
geom_bar(stat="identity", position=position_dodge(),size=0.2)+
geom_errorbar( aes(x=sample, ymin=mean_GFP-sd_GFP, ymax=mean_GFP+sd_GFP), width=0.4, alpha=0.9, size=1, position=position_dodge(0.9))

ggplot(d2,aes(y=mean_mko2,x=sample, fill=arabinose))+ 
geom_bar(stat="identity", position=position_dodge(),size=0.2)+
geom_errorbar( aes(x=sample, ymin=mean_mko2-sd_mko2, ymax=mean_mko2+sd_mko2), width=0.4, alpha=0.9, size=1, position=position_dodge(0.9))

ggplot(d2,aes(y=mean_OD,x=sample, fill=arabinose))+ 
geom_bar(stat="identity", position=position_dodge(),size=0.2)+
geom_errorbar( aes(x=sample, ymin=mean_OD-sd_OD, ymax=mean_OD+sd_OD), width=0.4, alpha=0.9, size=1, position=position_dodge(0.9))

d2$arabinose=as.numeric(d2$arabinose)

d2$arabinose[d2$arabinose==0] = 0.0000001
d2$arabinose[d2$arabinose==-1] =0.00000001


ggplot(d2,aes(y=mean_GFP,x=arabinose, col=as.factor(sample)))+  geom_line()+geom_point() +
geom_errorbar( aes(x=arabinose, ymin=mean_GFP-sd_GFP, ymax=mean_GFP+sd_GFP), width=0.1, alpha=0.9)+
scale_x_log10() + facet_grid(sample~.)





-------------------------------

d2<-subset(meanD, sample != "blank" & sample!= "empty" & sample !="NA" & time==100)
d2$arabinose[d2$arabinose==-1] ="off-target"

d3<-data.frame(d2$sample,d2$arabinose,d2$mean_GFP)
names(d3)<-c("sample","arabinose","GFP")


d3$GFP<-d3$GFP/max(d3$GFP)
ggplot(d3,aes(y=GFP,x=as.factor(arabinose), col=as.factor(sample)))+geom_point() +
facet_grid(sample~.)

write.table(d3,paste(path,"data.txt",sep="\\"),quote=FALSE,col.names=TRUE,sep="\t")

