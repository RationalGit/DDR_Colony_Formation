library("plyr")
library("stringr")
library("data.table")
library("doMC")
library("gtools")
library("plyr")
library("ggplot2")
library("ggrepel")
library("scales")
library("RColorBrewer")
library("gridExtra")
library("cowplot")
library("ggpubr")

#Clear out existing environment
rm(list = ls())

#Set the working directory to the local input folder
getwd()
dev.off()
setwd("/Users/webbert/Documents/Work/Experiments/LateSubDriv Large files:folders (unsynced)/Gene workup/DDR/Drug Dose-response/Colony formation/20181219 expt 2/R analysis")

#Get output file list
fileList = list.files()
fileList
df <- read.csv("DataTable_20181219_2_post.csv",na.strings = c("NA","NaN","","#DIV/0!","#VALUE!"))
df <- df[-grep(paste(c("Treatment","Cells.Plated"), collapse = "|"), colnames(df))]

df.key <- read.csv("DataKey_20181219_2_post.csv",na.strings = c("NA","NaN","","#DIV/0!","#VALUE!"))
df.colnames <- read.csv("Colnames.csv",na.strings = c("NA","NaN","","#DIV/0!","#VALUE!"))
df.colnames <- as.character(df.colnames$Colnames)

colnames(df)[7:23] <- df.colnames
data.df <- merge(df.key, df)

colnames(data.df)

df.median <- ddply(Filter(is.numeric, data.df), .(Sample.number), colwise(median))

df.final <- merge(df.key, df.median)
colnames(df.final)
cols <- colnames(df.final)[1:12]
df.final[,cols] <- lapply(df.final[,cols],as.factor)

df.final <- df.final[!is.na(df.final$siRNA), ]
df.final <- df.final[-grep("UT", df.final$siRNA), ]

Plate.treatment.order <- levels(df.final$IR.Dose)
Doseage <- levels(df.final$Dose.number)[c(1:3)]
Drugs <- levels(df.final$Drug)[c(8,14,19,5,7,9,10,15,16,17)]
Drugs_APH <- levels(df.final$Drug)[c(18,4,6,1,2,3)]
Drugs_HU <- levels(df.final$Drug)[c(18,11,13,12)]

Measurements <- colnames(df.final)[c(15:(length(colnames(df.final))-1))]
siRNA.order <- levels(df.final$siRNA)[c(3,1)]
siRNA.order.APH <- levels(df.final$siRNA)[c(3,1,2,4)]

df.Drugs <- df.final[df.final$Drug %in% Drugs,]
df.Drugs_APH <- df.final[df.final$Sample.number %in% 145:192,]
df.Drugs_HU <- df.final[df.final$Drug %in% Drugs_HU,]

length(unique(df.Drugs_APH$Drug))

Make.plots = function(save.image = FALSE, df, measurement){
  drug.list <- c()
  # Choice of label for plot title: 
  # Colour of graph based on IR or UV
  col.pal = brewer.pal(n=9, "PuBu")[6:3]
  
  # dataframe aspect selection
  if(length(unique(df$Drug)) > 6){
    drug.list = Drugs
    siRNA.list = siRNA.order
    assay = "Drug list"
  } else if(length(unique(df$Drug)) == 6){
    drug.list = Drugs_APH
    siRNA.list = siRNA.order.APH
    assay = "Rep Stress Test"
  } else if(length(unique(df$Drug)) < 6){
    drug.list = Drugs_HU
    siRNA.list = siRNA.order
    assay = "HU combinations"
  }
  
  # Graph title 
  graph.title = paste(as.character(measurement), assay, "Colony formation (2) treated post-plating", sep = " - ")
  
  # Plot function, takes dataframe, plots siRNA (plate order) on x axis by treatment time, against measurement on y axis
  working.plot = ggplot(df,
                        aes(x = factor(siRNA, level = siRNA.list), 
                            fill = factor(Dose.number, level = Doseage), y = df[[measurement]])) + 
    
    # labels
    labs(x = "siRNA", y = measurement, fill = "Dose Concentration", title = graph.title) + #S9.6 #pRPA2 S33    
    # graph type, stat+ position plot treatments alongside one another
    geom_bar(stat="identity", position = "dodge") + 
    facet_grid(cols = vars(factor(Drug, levels = drug.list)),
               rows = vars(Seeding.number)) +
    # SEM error bars
    # geom_errorbar(aes(x = reorder(siRNA, Order), ymin = (df[[feature]]-df[[gsub("MEAN", "SE", feature)]]), ymax = (df[[feature]]+df[[gsub("MEAN", "SE", feature)]])), width = .2, position = position_dodge(.9)) +
    
    # colour of bars
    scale_fill_manual(values = col.pal) +
    # text orientation
    # scale_y_continuous(breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1))))) +
    
    theme_light(base_size = 12) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank()) +
    theme(strip.background = element_rect(fill="gray90")) +
    theme(strip.text = element_text(colour = 'gray20')) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  
  ##Saving in function
  if(save.image){
    ggsave(
      paste(make.names(graph.title), ".pdf", sep=""),
      working.plot,
      width = 12,
      height = 12,
      scale = 1,
      dpi = 180
    )
    dev.off()
    print("saved")
  } else{
    print(working.plot)
  }
  
  #file_name = paste(gsub("/", "x", graph.title), ".png", sep="")
  #png(file_name, width = 960, height = 960, bg = "transparent", pointsize = 12)
  #print(working.plot)
  
  
}

Make.plots(FALSE, df.Drugs_APH, Measurements[1])
setwd("./Output/")
getwd()
lapply(Measurements, function(x) {Make.plots(TRUE, df.Drugs,x)})
lapply(Measurements, function(x) {Make.plots(TRUE, df.Drugs_APH,x)})
lapply(Measurements, function(x) {Make.plots(TRUE, df.Drugs_HU,x)})
