if(!require(magick)) install.packages("magick")
if(!require(readxl)) install.packages("readxl")
library(magick)
library(readxl)##setwd("/Users/daxinggao/scripts/")##mac
setwd("C:/Users/gaoda/OneDrive - USTC/Gao lab/Scripts")##windows
input=readxl::read_excel("./figure layout/figure layout.xlsx","Sheet1")
title=as.character(colnames(input)[1])
rn=nrow(input)/2
cn=ncol(input)
ex=4
figure_width=300*ex
figure_folder=as.character(colnames(input)[2])
dir.create(paste("./figure layout/layouts/",format(Sys.time(), "%Y%m%d"),title))
folder_path=paste("./figure layout/layouts/",format(Sys.time(), "%Y%m%d"),title)
file.copy("./figure layout/figure layout.xlsx",paste(folder_path,"/input.xlsx",sep=""))
for(j in 1:rn){##row
  for(i in 1:cn){##column
    path=paste0(figure_folder,"/",as.character(input[2*j-1,i]),".jpg")
    file_name=tail(strsplit(path,"/")[[1]],n=1)
    file.copy(path,paste(folder_path,"/",file_name,sep=""))
    if(is.na(as.character(input[2*j-1,i]))){
      blot=image_blank(figure_width+1*ex,figure_width+1*ex,"white")
    }
    else{blot=image_read(path)
    blot=image_scale(blot,geometry_area(width=figure_width))
    blot=image_border(blot,"black","4X4")}#change 4 fold
    blot=image_border(blot,"white","60x100")#change 4 fold
    if(!is.na(as.character(input[2*j-1,i])))blot=image_annotate(blot,as.character(input[2*j,i]),size = 20*ex,gravity="south")##title of western blot
      if(i==1){##not first graph but the first column
        blot_row=blot
      }
      else{

        blots=c(blot_row,blot)
        blot_row=image_append(blots)##append graph on each row
        if(i==cn){## end of the row
          if(j==1){blot_all=blot_row}
          else{

          blotsall=c(blot_all,blot_row)
          blot_all=image_append(blotsall,stack = T)

        }
      }
    }
  }
}
blot_all=image_border(blot_all,"white","40x120")
blot_all=image_annotate(blot_all,title,size = 30*ex,gravity="north")
image_write(path = paste0(folder_path,"/",format(Sys.time(), "%Y%m%d"),title,".png"),blot_all)
