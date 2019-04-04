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
library("reshape2")

# Clear out existing environment
rm(list = ls())

# Set the working directory to the local input folder
getwd()
dev.off()
setwd("/Users/webbert/Documents/Work/Experiments/LateSubDriv Large files:folders (unsynced)/Gene workup/DDR/Drug Dose-response/Colony formation/201901 cgRNA/R analysis")

# Get output file list
fileList = list.files()
fileList

# Import data, data key and readable column names
df.cg5 <- read.csv("DataTable_201901_cg5.csv",na.strings = c("NA","NaN","","#DIV/0!","#VALUE!"))
df.cg5 <- df.cg5[-grep(paste(c("Treatment","Cells.Plated"), collapse = "|"), colnames(df.cg5))]

df.cg5.key <- read.csv("DataKey_201901_cg5.csv",na.strings = c("NA","NaN","","#DIV/0!","#VALUE!"))
df.cg5.colnames <- read.csv("Colnames.csv",na.strings = c("NA","NaN","","#DIV/0!","#VALUE!"))
df.cg5.colnames <- as.character(df.cg5.colnames$Colnames)

colnames(df.cg5)[7:23]
colnames(df.cg5)[7:23] <- df.cg5.colnames
data.df.cg5 <- merge(df.cg5.key, df.cg5)

#remove cgNT_5 (data not useful)
data.df.cg5 <- data.df.cg5[-grep("cgNT_5", data.df.cg5$siRNA),]

colnames(data.df.cg5)

# Import second set of data, data key and readable column names
df.cg6 <- read.csv("DataTable_201901_cg6.csv",na.strings = c("NA","NaN","","#DIV/0!","#VALUE!"))
df.cg6 <- df.cg6[-grep(paste(c("Treatment","Cells.Plated"), collapse = "|"), colnames(df.cg6))]

df.cg6.key <- read.csv("DataKey_201901_cg6.csv",na.strings = c("NA","NaN","","#DIV/0!","#VALUE!"))
df.cg6.colnames <- read.csv("Colnames.csv",na.strings = c("NA","NaN","","#DIV/0!","#VALUE!"))
df.cg6.colnames <- as.character(df.cg6.colnames$Colnames)

colnames(df.cg6)[7:23]
colnames(df.cg6)[7:23] <- df.cg6.colnames
data.df.cg6 <- merge(df.cg6.key, df.cg6)

colnames(data.df.cg6)

df.key <- rbind(df.cg5.key,df.cg6.key)
data.df <- rbind(data.df.cg5,data.df.cg6)

# remove unnecessary columns
filter.cols.out <- c(1:13, 16:17, 36:42)
filter.cols.in <- c(1:13, 16:17)

# Find mean and stdev of triplicate plate values
df.mean <- ddply(Filter(is.numeric, data.df), .(Sample.number,cg), colwise(mean))
df.sd <- ddply(Filter(is.numeric, data.df[, -filter.cols.out]), .(Sample.number,cg), colwise(sd))

colnames(df.sd) <- paste(colnames(df.sd), ".sd",sep = "")
colnames(df.mean)
View(df.sd)
colnames(df.sd)
df.SD <- cbind(df.mean, df.sd[,-c(1:2)])

df.final <- merge(df.key, df.SD)
df.final <- df.final[,-c(16, 35)]
colnames(df.final)
cols <- colnames(df.final)[1:15]

df.final[,cols] <- lapply(df.final[,cols],as.factor)

df.final <- df.final[!is.na(df.final$siRNA), ]

#String wrap function to generate line breaks in names for display on graphs 
swr = function(string, nwrap=8) {
  paste(strwrap(string, width=nwrap), collapse="\n")
}
swr = Vectorize(swr)

# Create line breaks in Drugs
levels(df.final$Drug) = swr(levels(df.final$Drug))

Plate.treatment.order <- levels(df.final$IR.Dose)
Doseage <- levels(df.final$Dose.number)[c(1:3)]
Drugs <- c("DMSO","KU55933\n(ATMi)", "VE821\n(ATRi)", "Camptothecin", "Cisplatin", "MMC", "DRB", "MK1775\n(Wee1i)", "Olaparib", "Olaparib+Wee1i", "Pladeinolide\nB")
Drugs.RS <- c("DMSO", "APH", "APH+MK1775\n(Wee1i)", "APH+CCT2355737\n(Chk1i)", "HU", "HU+MK1775\n(Wee1i)", "HU+CCT2355737\n(Chk1i)", "HU+KU55933+CCT2355737\n(+ATMi+Chk1i)")
Drugs.IR <- c("DMSO", "KU55933\n(ATMi)", "VE821\n(ATRi)", "Camptothecin", "DRB", "Olaparib", "Pladeinolide\nB", "MK1775\n(Wee1i)", "CCT2355737\n(Chk1i)", "APH", "HU")
Measurements <- colnames(df.final)[-grep(".sd",colnames(df.final))][c(16:33)]
siRNA.order <- levels(df.final$siRNA)[c(2,4,1,3)]

df.Drugs <- df.final[df.final$Sample.number %in% c(1:48, 145:192),]
df.Drugs.RS <- df.final[df.final$Drug %in% Drugs.RS & df.final$IR.Dose  %in% 0,]
df.Drugs.IR <- df.final[df.final$IR.Dose  %in% c(3,6),]
df.Drugs.IR.short.pre <- df.Drugs.IR[df.Drugs.IR$Pretreat.time]
df.Drugs.IR$Pretreat.time

length(unique(df.Drugs_RS$Drug))
"CCT2355737 (Chk1i)" %in% unique(df.Drugs_IR$Drug)

Make.plots = function(save.image = FALSE, df, measurement){
  drug.list <- c()
  # Choice of label for plot title: 
  # Colour of graph based on IR or UV
  col.pal = brewer.pal(n=9, "PuBu")[6:3]
  error.measure = paste(as.character(measurement),".sd",sep = "")
  print(error.measure)
  
  # dataframe aspect selection
  if("Cisplatin" %in% unique(df$Drug)){
    drug.list = Drugs
    assay = "Genomic Stress Inducers"
    drug.fill = factor(df$Dose.number, level = Doseage)
    fill.label = "Dose Concentration"
    col.pal = brewer.pal(n=9, "PuBu")[6:3]
    chart.height = 9
  } else if("HU+KU55933+CCT2355737\n(+ATMi+Chk1i)" %in% unique(df$Drug)){
    drug.list = Drugs.RS
    assay = "Replication Stress combinations"
    drug.fill = factor(df$Dose.number, level = Doseage)    
    fill.label = "Dose Concentration"
    col.pal = brewer.pal(n=9, "PuBu")[6:3]
    chart.height = 9
  } else if("CCT2355737\n(Chk1i)" %in% unique(df$Drug)){
    drug.list = Drugs.IR
    assay = "IR dose (Gy)"
    drug.fill = factor(df$Pretreat.time)
    fill.label = "Pretreat Time"
    col.pal = brewer.pal(n=9, "Oranges")[4:5]
    chart.height = 12
  }
  
  # Graph title 
  graph.title = paste(gsub(measurement, pattern = "\\.", replacement = " "), assay, "Colony formation (2) treated post-plating", sep = " - ")
  
  # Plot function, takes dataframe, plots siRNA (plate order) on x axis by treatment time, against measurement on y axis
  working.plot = ggplot(df,
                        aes(x = factor(siRNA, level = siRNA.order), 
                            fill = drug.fill, y = df[[measurement]])) + 
    
    # labels
    labs(x = "cgRNA", y =gsub(measurement, pattern = "\\.", replacement = " "), fill = fill.label, title = graph.title) + #S9.6 #pRPA2 S33    
    # graph type, stat+ position plot treatments alongside one another
    geom_bar(stat="identity", position = "dodge") + 
    
    facet_grid(cols = vars(factor(Drug, levels = drug.list)),
               rows = vars(IR.Dose, Treatment)) +
    geom_errorbar(data = df, 
                  mapping = aes(x = factor(siRNA, level = siRNA.order), 
                                ymin = (df[[measurement]] - df[[error.measure]]), 
                                ymax = (df[[measurement]] + df[[error.measure]])),
                  colour = "grey57",
                  width = .2, 
                  position = position_dodge(.9)) +
    # SEM error bars
    # geom_errorbar(aes(x = reorder(siRNA, Order), ymin = (df[[feature]]-df[[gsub("MEAN", "SE", feature)]]), ymax = (df[[feature]]+df[[gsub("MEAN", "SE", feature)]])), width = .2, position = position_dodge(.9)) +
    
    # colour of bars
    scale_fill_manual(values = col.pal) +
    # text orientation
    # scale_y_continuous(breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1))))) +
    
    theme_light(base_size = 11) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank()) +
    theme(strip.background = element_rect(fill="gray90")) +
    theme(strip.text = element_text(size = 7, colour = 'gray20')) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  
  
  ##Saving in function
  if(save.image){
    ggsave(
      paste(make.names(graph.title), ".pdf", sep=""),
      working.plot,
      width = 12,
      height = chart.height,
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

Make.plots(FALSE, df.Drugs.IR, Measurements[1])
setwd("./Output/")
getwd()
lapply(Measurements, function(x) {Make.plots(TRUE, df.Drugs,x)})
lapply(Measurements, function(x) {Make.plots(TRUE, df.Drugs.RS,x)})
lapply(Measurements, function(x) {Make.plots(TRUE, df.Drugs.IR,x)})
