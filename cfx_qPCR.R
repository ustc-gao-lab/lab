if(!require(magick)) install.packages("magick")
if(!require(data.table)) install.packages("data.table")
if(!require(clipr)) install.packages("clipr")
if(!require(psych)) install.packages("psych")
if(!require(ggpubr)) install.packages("ggpubr")
if(!require(ggsci)) install.packages("ggsci")
if(!require(scales)) install.packages("scales")
if(!require(lemon)) install.packages("lemon")
if(!require(xlsx)) install.packages("xlsx")
if(!require(cowplot)) install.packages("cowplot")
if(!require(ggplot2)) install.packages("ggplot2")
library(ggplot2)
library(data.table)
library(dplyr)
library(clipr)
library(psych)
library(reshape2)
library(ggpubr)
library(ggsci)
library(scales)
library(lemon)
library(xlsx)
library(cowplot)
#write_clip(as.data.frame(rep(paste("Sample",1:12),each=2)),col.names=F)
#write_clip(as.data.frame(rep(paste("Sample",13:24),each=2)),col.names=F)
#write_clip(as.data.frame(rep(paste("Sample",13:24),each=1)),col.names=F)
#write_clip(as.data.frame(rep(paste("Sample",1:12),each=1)),col.names=F)
#write_clip(as.data.frame(rep(paste("Sample",25:36),each=2)),col.names=F)
#write_clip(as.data.frame(rep(paste("Sample",seq(i,24+i,by=8)),each=4)),col.names=F)
#c(outer(8*(0:2),1:8,FUN="+"))
#write_clip(as.data.frame(rep(paste("Sample",c(outer(8*(0:2),1:8,FUN="+"))),each=1)))
#[1]  1  9 17  2 10 18  3 11 19  4 12 20  5 13 21  6 14 22  7 15 23  8 16 24

options(warn=-1)
setwd("C:/Users/gaoda/OneDrive - USTC/Gao lab/Scripts")##working folder 初次使用需要更改
input=readxl::read_excel("./qPCR analysis/qPCR.xlsx","Sheet1")
title=as.character(input[1,2])
file_number=as.numeric(input[2,2])
log_scale=as.character(input[3,2])
gene_number=as.numeric(input[4,2])
ylab=input[6,2:(1+gene_number)]
dt=as.character(input[1,4])
colnames(ylab)=as.character(input[5,2:(1+gene_number)])
miss_=1
remove_undetermine=as.character(input[7,2])
undetermined=as.numeric(input[8,2])
normalization=as.character(input[9,2])
##dir.create(("./qPCR analysis/qPCR results"))
dir.create(paste0("./qPCR analysis/qPCR results/",dt))
dir.create(paste0("./qPCR analysis/qPCR results/",dt,"/",format(Sys.time(), "%Y%m%d"),title))##save folder
folder_path=paste0("./qPCR analysis/qPCR results/",dt,"/",format(Sys.time(), "%Y%m%d"),title)
file.copy("./qPCR analysis/qPCR.xlsx",paste(folder_path,"/input",format(Sys.time(), "%Y%m%d"),".xlsx",sep=""))
Cq_file=paste0("Z:/高大兴/qPCR/",as.character(input[1,4]))### 初次使用需要更改
#for(i in 1:file_number){##read all qPCR data into ts
  t=read.xlsx(paste0(Cq_file,"/",list.files(Cq_file)[grep("Cq",list.files(Cq_file))]),"0")
  t=as_tibble(t)
  t[[1]]=NULL
  file_name=list.files(Cq_file)[grep("Cq",list.files(Cq_file))]##tail(strsplit(as.character(input[3,i+1]),"/")[[1]],n=1)
  file.copy(paste0(Cq_file,"/",file_name),paste(folder_path,"/",file_name,sep=""))
  ##colnames(t)=t[1,]
  t=t[!is.na(t$Sample),]
  t=t[1:nrow(t),c(3,5,7)]
  colnames(t)=c("gene","Sample Name","CT")
  ts=t
  #t=t[t[[1]]!=""&t[[2]]!="",]
  #if(i==1)ts=t
  #else ts=rbind(ts,t)
#}
###sample ID change
inputn=input[11:nrow(input),1:4]
colnames(inputn)=c("Sample_Name","SampleName","group","order")
if(any(is.na(inputn[4]))){
unique_sample=unique(inputn[2:3],sorting=F)
unique_sample$order=1:nrow(unique_sample)
inputn=merge(inputn[1:3],unique_sample,by=c("SampleName","group"))}
write.xlsx(merge(ts,inputn,by.x="Sample Name",by.y="Sample_Name",all.y=T),paste0(folder_path,"/Input_label",format(Sys.time(), "%Y%m%d"),".xlsx"))

if(miss_==1){ ##label missing values
  ts$miss_show=as.numeric(is.na(ts$CT))
}else ts$miss_show=0
ts=ts[ts$CT>10||is.na(ts$CT),]## 15 CT
if(remove_undetermine=="T"){
  ts=ts[!is.na(ts$CT),]
}

for(i in 1:nrow(ts)){##convert undetermined into value
  if(is.na(ts[i,3]))ts[i,3]=undetermined
}

ts[[3]]=as.numeric(ts[[3]])
ts=cbind(ts[,c(2,1,3)],ts$miss_show)#####Sample Name, Target Name, CT
colnames(ts)=c("id","gene","ct","miss")

if(input[7,4]=="F"){ts=ts%>%group_by(id,gene)%>%summarise_all(funs(mean))}##technical replicates
ts2=ts%>%group_by(id,gene)%>%mutate(count=sequence(n()))

tss=split(ts2,as.factor(ts2[[2]]))##reshape the tabl
if(unique(tss[[1]]$id)==""){tss[[1]]=NULL}

for(i in 1:length(tss)){#add names
  colnames(tss[[i]])[3]=as.character(unlist(tss[[i]][1,2]))
  colnames(tss[[i]])[4]=paste0("miss_",as.character(unlist(tss[[i]][1,2])))
  tss[[i]]=tss[[i]][,-2]
}
tr=tss[[1]]##reshape
if(length(tss)>1){
for(i in 2:length(tss)){
 tr=merge(tr,tss[[i]],by=c("id","count"),all=T,sort=F)
}
tr=tr[,-2]
for(i in 1:length(tss)){##internal control be first one
  
  if(colnames(tr)[2*i]==as.character(input[5,2])){
      a=tr[2:3]
      b=colnames(tr)[2:3]
      tr[2:3]=tr[(2*i):(2*i+1)]
      colnames(tr)[2:3]=colnames(tr)[(2*i):(2*i+1)]
      tr[(2*i):(2*i+1)]=a
      colnames(tr)[(2*i):(2*i+1)]=b
    break
  }
}
}

for(i in 1:(length(tss)-1)){
  tr[[2*i+2]]=as.numeric(tr[[2]])-as.numeric(tr[[2*i+2]])
  tr[[2*i+2]]=2^tr[[2*i+2]]
}


xlab=inputn
xlab$order=as.numeric(xlab$order)
result=merge(tr,xlab,by.x="id",by.y="Sample_Name",sort = F,all.y=T)
result=result[,-1]
##result=result%>%group_by(SampleName,genotype)%>%summarise_all(funs(mean))
result=result[order(result$order),]

if(normalization=="T"){
  for(i in 1:(length(tss)-1)){
    if(!is.na(result[1,2*i+1]))
      result[[2*i+1]]=as.numeric(result[[2*i+1]])/as.numeric(result[[2*i+1]][1])
  }
}
##result=result[order(result[[1]]),]
##result=result[order(sapply(as.character(result[[1]]),nchar)),]
if(nrow(unique(result[ncol(result)-1]))==1){pal="#80818099"}else pal=c("blue","red","purple","dodgerblue","darkorange","coral4","navy","green1",pal_aaas("default")(10),pal_npg("nrc")(10))
##http://applied-r.com/category/r-colors/ show_col("blue")

result$group=factor(result$group,levels=unique(result$group))#keep the original order
result$SampleName=factor(result$SampleName,levels=unique(result$SampleName))
write.xlsx(result,paste0(folder_path,"/result",format(Sys.time(), "%Y%m%d"),".xlsx"))
result2=result
result=result[c((1:length(tss))*2-1,(ncol(result)-2):(ncol(result)-1))]
result_miss=result2[(1:length(tss))*2]
result_miss=as.data.frame(ifelse(result_miss==0,"","#"))

if(log_scale=="T"){result[2:(ncol(result)-2)]=log10(result[2:(ncol(result)-2)])}

plot_genes=intersect(colnames(result),colnames(ylab))

plots1=ggbarplot(result,"SampleName",plot_genes,color=0,fill="group",palette =pal,##plot each graph
                 add = c("mean_sd"),error.plot = "upper_errorbar", 
                 remove=NA,
                 add.params = list(group = "group",position=position_dodge(preserve = "single",width = 0.7),color="black",fill="group"),#
                position=position_dodge(preserve = "single",width = 0.7),
                 x.text.angle = 45,xlab ="",title=title,ggtheme = theme_pubr(9,legend = "right"))###base font size


plots2=ggbarplot(result,"SampleName",plot_genes,color="group",fill="group",size=0.1,palette =pal, ##plot layout
                 add = c("mean_sd"##, "dotplot"
                         ),error.plot = "upper_errorbar", 
                #width = 0.5,
                remove=NA,
                add.params = list(group = "group",position=position_dodge(preserve = "single",width = 0.7),color="black",fill="group"),
                  position=position_dodge(preserve = "single",width = 0.7),
                  x.text.angle = 45,xlab ="",title=title,ggtheme=theme_pubr(9,legend = "none"))


lgd=get_legend(plots1[[2]]+ theme(legend.position="bottom",legend.title = element_text(size=9,face="bold"),
                                  legend.text = element_text(colour="black", size=9,face="bold"),
                                  legend.key.size=unit(0.2, "cm")))
ylab=ylab[,match(plot_genes,colnames(ylab))]#change o=rder of ylab alphabetical

#ylab=cbind(as.character(input[6,2]),ylab)
  for(i in 1:length(plots1)){
    if(all(is.na(select(result,plot_genes)[[i]])))next###position=position_dodge(width=0.9, preserve = "single"),
  plot=plots1[[i]]+geom_point(aes(group=result$group),size=.7,alpha=0.3,color="black",fill=0.01,position=position_dodge(width = 0.7))+scale_y_continuous(name=as.character(ylab[1,i]),expand=c(0,0),breaks= pretty_breaks())+coord_capped_cart(left='both',ylim=c(floor(min(0,min(select(result,plot_genes)[[i]],na.rm=T))),max(select(result,plot_genes)[[i]]*1.4,na.rm=T)))+theme(plot.margin = margin(30, 30, 40, 50))+geom_text(aes(x=result$SampleName,y=result[[i]],label=result_miss[[i]],group=result$group),position=position_dodge(width = 0.7),check_overlap =T,vjust=-2)
       #trans="log10",ylim=c(log10(min(select(result,plot_genes)[[i]]/1.4,na.rm=T))
  
  ggsave(paste0(format(Sys.time(), "%Y%m%d"),plot_genes[i],".png"),plot=plot,width = length(unique(xlab[[1]]))*0.5+length(unique(xlab[[2]]))*0.25+3,height = 4,path = folder_path,device="png")
  
  plots2[[i]]=plots2[[i]]+geom_point(aes(group=result$group),size=.7,alpha=0.3,color="black",fill=0.01,position=position_dodge(width=0.7))+scale_y_continuous(name=as.character(ylab[1,i]),expand=c(0,0),breaks= pretty_breaks())+coord_capped_cart(left='both',ylim=c(floor(min(0,min(select(result,plot_genes)[[i]],na.rm=T))),max(select(result,plot_genes)[[i]]*1.4,na.rm=T)))+theme(plot.margin = margin(30, 60, 50, 50))
  ##margin up right bottom left
  }
#i
for(i in 1:ceiling((length(plots1)+1)/4)){
  plot_layout=do.call(plot_grid,c(plots2[(4*i-3):min(4*i,length(plots2))],scale=1.4-nrow(result)/200,nrow=2,ncol=2))##scale enlarge the plot in graph
  plot_layout=plot_grid(plot_layout,lgd,rel_heights=c(1,.03),nrow=2)
  plot_layout=plot_layout+theme(plot.margin = margin(10, 10, 10, 5))
  save_plot(paste0(folder_path,"/",format(Sys.time(), "%Y%m%d"),"_layout",i,".png"),plot=plot_layout,base_width = length(unique(xlab[[1]]))*0.6+length(unique(xlab[[2]]))*0.3+6,base_height=nrow(result)*0.15+5,device="png")
}

result3=result %>% group_by(SampleName,group) %>% mutate(flag=sequence(n()))#对于重复的数据，标记flag
result4= reshape(as.data.frame(result3),idvar=c("SampleName","flag"),timevar = c("group"),direction="wide")#根据不同组reshape
result5= reshape(result4,idvar=c("SampleName"),timevar = c("flag"),direction="wide")#将重复的数据reshape
result6=result5[-1]
result7=cbind(result5[1],result6[order(colnames(result6))])#排序
write.xlsx(result7,paste0(folder_path,"/Graphpad",format(Sys.time(), "%Y%m%d"),".xlsx"))

options(warn=0)

