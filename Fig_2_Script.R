
setwd("E:/Dropbox/Git_Files/FAME02/")

## Libraries
require(reshape2)
require(ggplot2)
library(RColorBrewer)
library(Rmisc)


## Nudge parameter
nudge <- 0.15
a <- 1

## Read the Lectin Data
lec.fname <- "FAME02_LW_data.txt"
raw.data <- read.table(lec.fname, sep = "\t", header = T)

## Read the Microbiome Data
mb.v2 <- read.delim("FAME02_vag_isolates_visit2.txt")     # Visit 1
mb.v3 <- read.delim("FAME02_vag_isolates_visit3.txt")     # Visit 2
mb.v4 <- read.delim("FAME02_vag_isolates_visit4.txt")     # Visit 3

## Melt the raw data into the long format
lec.data <- melt(raw.data, id.vars = c("Sample", "Visit", "Group", "Formulation", "Drug"))
names(lec.data)[6] <- "Lectin"

## Remove CVN from the data as that lectin was reported to be acting weird by Lara and agreed by Brian. 
lec.data <- subset(lec.data, !(Lectin %in% "CVN"))


## Subset the data from gel related samples
lec.data.gel <- subset(lec.data, Formulation %in% "Gel")
lec.data.GP <- subset(lec.data, Group %in% "Gel_Placebo")
lec.data.GD <- subset(lec.data, Group %in% "Gel_Dapivirine")

## Subset the data from film related samples
lec.data.film <- subset(lec.data, Formulation %in% "Film")
lec.data.FP <- subset(lec.data, Group %in% "Film_Placebo")
lec.data.FD <- subset(lec.data, Group %in% "Film_Dapivirine")
# 
# ################################### Gel (Placebo + Dapivirine) #############################################################
# 
# ## Gel samples :: data processing
# lec.data.gel2 <- dcast(lec.data.gel, Lectin~Visit, value.var = "value", fun.aggregate = mean)
# 
# lec.data.gel3 <- dcast(lec.data.gel, Lectin~Visit+Sample, value.var = "value", fun.aggregate = mean)
# n.samples.gel <- length(unique(lec.data.gel$Sample))
# 
# 
# lec.data.gel2$V2V1.pval <- apply(lec.data.gel3[,-1], 1, 
#                                  function(x){ t.test(as.numeric(x)[1:n.samples.gel], as.numeric(x)[(n.samples.gel+1):(2*n.samples.gel)],
#                                                      alternative = "two.sided", paired = T, var.equal = F)$p.val})
#                                  
# lec.data.gel2$V3V1.pval <- apply(lec.data.gel3[,-1], 1, 
#                                  function(x){ t.test(as.numeric(x)[1:n.samples.gel], as.numeric(x)[((2*n.samples.gel)+1):(3*n.samples.gel)],
#                                                      alternative = "two.sided", paired = T, var.equal = F)$p.val})
# 
# lec.data.gel2$log2FC.V2V1 <- lec.data.gel2$Visit.2 - lec.data.gel2$Visit.1
# lec.data.gel2$log2FC.V3V1 <- lec.data.gel2$Visit.3 - lec.data.gel2$Visit.1
# 
# lec.data.gel2$log10pval.V2V1 <- -log10(lec.data.gel2$V2V1.pval)
# lec.data.gel2$log10pval.V3V1 <- -log10(lec.data.gel2$V3V1.pval)
# 
# lec.data.gel2$Threshold <- as.factor(lec.data.gel2$log10pval.V2V1 >= -log10(0.05))
# 
# 
# Gel.volcanoplot <- ggplot(lec.data.gel2, aes(log2FC.V2V1, log10pval.V2V1, label = Lectin)) + 
#   geom_point() + 
#   xlim(c(-1,1)) +
#   geom_hline(yintercept = -log10(0.05), col = "red", lty = 2) + 
#   geom_text(data = subset(lec.data.gel2, log10pval.V2V1 >= -log10(0.05))) + 
#   theme_bw()
# 

########################################### Gel Placebo #############################################################

## Gel Placebo :: Data Processing
lec.data.GP2 <- dcast(lec.data.GP, Lectin~Visit, value.var = "value", fun.aggregate = mean)

lec.data.GP3 <- dcast(lec.data.GP, Lectin~Visit+Sample, value.var = "value", fun.aggregate = mean)
n.samples.GP <- length(unique(lec.data.GP$Sample))


lec.data.GP2$V2V1.pval <- apply(lec.data.GP3[,-1], 1, 
                                 function(x){ t.test(as.numeric(x)[1:n.samples.GP], as.numeric(x)[(n.samples.GP+1):(2*n.samples.GP)],
                                                     alternative = "two.sided", paired = T, var.equal = F)$p.val})

lec.data.GP2$V3V1.pval <- apply(lec.data.GP3[,-1], 1, 
                                 function(x){ t.test(as.numeric(x)[1:n.samples.GP], as.numeric(x)[((2*n.samples.GP)+1):(3*n.samples.GP)],
                                                     alternative = "two.sided", paired = T, var.equal = F)$p.val})

lec.data.GP2$log2FC.V2V1 <- lec.data.GP2$Visit.2 - lec.data.GP2$Visit.1
lec.data.GP2$log2FC.V3V1 <- lec.data.GP2$Visit.3 - lec.data.GP2$Visit.1

lec.data.GP2$log10pval.V2V1 <- -log10(lec.data.GP2$V2V1.pval)
lec.data.GP2$log10pval.V3V1 <- -log10(lec.data.GP2$V3V1.pval)

lec.data.GP2$Threshold <- as.factor(lec.data.GP2$log10pval.V2V1 >= -log10(0.05))

## Manualy adding color to the significant lectin labels.
lec.data.GP2$Lectin.Color <- "lg"
lec.data.GP2[lec.data.GP2$log10pval.V2V1 >= -log10(0.05), "Lectin.Color"] <- "dg"
lec.data.GP2[lec.data.GP2$Lectin %in% c("SVN","SVN.2","GRFT"), "Lectin.Color"] <- "g"
lec.data.GP2[lec.data.GP2$Lectin %in% "a.Sialyl.Le.X", "Lectin.Color"] <- "r"


## Gel Placebo:: Volcano Plot
GP.volcanoplot <- ggplot(lec.data.GP2, aes(log2FC.V2V1, log10pval.V2V1, label = Lectin, colour = as.factor(Lectin.Color))) + 
                    geom_point(show.legend = F, size = 2, shape = 21) + 
                    labs(x = "Gel Placebo [log2(V3/V1)]", y = "-log10(p-value)") +
                    xlim(c(-2,2)) + ylim(c(0,3)) + 
                    scale_color_manual(values = c("lg" = "lightgray", "dg" = "darkgray", "g" = "chartreuse3", "r" = "red")) +
                    geom_hline(yintercept = -log10(0.05), col = "red", lty = 2) + 
                    geom_text(data = subset(lec.data.GP2, log10pval.V2V1 >= -log10(0.05)),
                              nudge_y = nudge, show.legend = F, fontface = "bold", alpha = a) + 
                    theme_bw()




# 
# ########################################### Gel Dapivirine #############################################################
# 
# ## Gel Dapivirine :: Data Processing
# lec.data.GD2 <- dcast(lec.data.GD, Lectin~Visit, value.var = "value", fun.aggregate = mean)
# 
# lec.data.GD3 <- dcast(lec.data.GD, Lectin~Visit+Sample, value.var = "value", fun.aggregate = mean)
# n.samples.GD <- length(unique(lec.data.GD$Sample))
# 
# 
# lec.data.GD2$V2V1.pval <- apply(lec.data.GD3[,-1], 1, 
#                                 function(x){ t.test(as.numeric(x)[1:n.samples.GD], as.numeric(x)[(n.samples.GD+1):(2*n.samples.GD)],
#                                                     alternative = "two.sided", paired = T, var.equal = F)$p.val})
# 
# lec.data.GD2$V3V1.pval <- apply(lec.data.GD3[,-1], 1, 
#                                 function(x){ t.test(as.numeric(x)[1:n.samples.GD], as.numeric(x)[((2*n.samples.GD)+1):(3*n.samples.GD)],
#                                                     alternative = "two.sided", paired = T, var.equal = F)$p.val})
# 
# lec.data.GD2$log2FC.V2V1 <- lec.data.GD2$Visit.2 - lec.data.GD2$Visit.1
# lec.data.GD2$log2FC.V3V1 <- lec.data.GD2$Visit.3 - lec.data.GD2$Visit.1
# 
# lec.data.GD2$log10pval.V2V1 <- -log10(lec.data.GD2$V2V1.pval)
# lec.data.GD2$log10pval.V3V1 <- -log10(lec.data.GD2$V3V1.pval)
# 
# lec.data.GD2$Threshold <- as.factor(lec.data.GD2$log10pval.V2V1 >= -log10(0.05))
# 
# ## Manualy adding color to the significant lectin labels.
# lec.data.GD2$Lectin.Color <- "lg"
# lec.data.GD2[lec.data.GD2$log10pval.V2V1 >= -log10(0.05), "Lectin.Color"] <- "dg"
# lec.data.GD2[lec.data.GD2$Lectin %in% c("SVN","SVN.2","GRFT"), "Lectin.Color"] <- "g"
# lec.data.GD2[lec.data.GD2$Lectin %in% "a.Le.X", "Lectin.Color"] <- "r"
# 
# 
# ## Gel Dapivirine:: Volcano Plot
# GD.volcanoplot <- ggplot(lec.data.GD2, aes(log2FC.V2V1, log10pval.V2V1, label = Lectin, colour = as.factor(Lectin.Color))) + 
#   geom_point(show.legend = F) + 
#   xlim(c(-2,2)) + ylim(c(0,3)) + 
#   scale_color_manual(values = c("lg" = "lightgray", "dg" = "darkgray", "g" = "chartreuse3", "r" = "red")) +
#   geom_hline(yintercept = -log10(0.05), col = "red", lty = 2) + 
#   geom_text(data = subset(lec.data.GD2, log10pval.V2V1 >= -log10(0.05)),
#             nudge_y = 0.25, show.legend = F, fontface = "bold") + 
#   theme_bw()
# 
            


########################################### Film Placebo #############################################################

## Film Placebo :: Data Processing
lec.data.FP2 <- dcast(lec.data.FP, Lectin~Visit, value.var = "value", fun.aggregate = mean)

lec.data.FP3 <- dcast(lec.data.FP, Lectin~Visit+Sample, value.var = "value", fun.aggregate = mean)
n.samples.FP <- length(unique(lec.data.FP$Sample))


lec.data.FP2$V2V1.pval <- apply(lec.data.FP3[,-1], 1, 
                                function(x){ t.test(as.numeric(x)[1:n.samples.FP], as.numeric(x)[(n.samples.FP+1):(2*n.samples.FP)],
                                                    alternative = "two.sided", paired = T, var.equal = F)$p.val})

lec.data.FP2$V3V1.pval <- apply(lec.data.FP3[,-1], 1, 
                                function(x){ t.test(as.numeric(x)[1:n.samples.FP], as.numeric(x)[((2*n.samples.FP)+1):(3*n.samples.FP)],
                                                    alternative = "two.sided", paired = T, var.equal = F)$p.val})

lec.data.FP2$log2FC.V2V1 <- lec.data.FP2$Visit.2 - lec.data.FP2$Visit.1
lec.data.FP2$log2FC.V3V1 <- lec.data.FP2$Visit.3 - lec.data.FP2$Visit.1

lec.data.FP2$log10pval.V2V1 <- -log10(lec.data.FP2$V2V1.pval)
lec.data.FP2$log10pval.V3V1 <- -log10(lec.data.FP2$V3V1.pval)

lec.data.FP2$Threshold <- as.factor(lec.data.FP2$log10pval.V2V1 >= -log10(0.05))

## Manualy adding color to the significant lectin labels.
lec.data.FP2$Lectin.Color <- "lg"
lec.data.FP2[lec.data.FP2$log10pval.V2V1 >= -log10(0.05), "Lectin.Color"] <- "dg"
lec.data.FP2[lec.data.FP2$Lectin %in% c("CVN","UDA","GNA_EY"), "Lectin.Color"] <- "b"


## Film Placebo:: Volcano Plot
FP.volcanoplot <- ggplot(lec.data.FP2, aes(log2FC.V2V1, log10pval.V2V1, label = Lectin, colour = as.factor(Lectin.Color))) + 
  geom_point(show.legend = F, size = 2, shape = 21) +  
  labs(x = "Film Placebo [log2(V3/V1)]", y = "-log10(p-value)") + 
  xlim(c(-2,2)) + ylim(c(0,3)) + 
  scale_color_manual(values = c("lg" = "lightgray", "dg" = "darkgray", "b" = "blue")) +
  geom_hline(yintercept = -log10(0.05), col = "red", lty = 2) + 
  geom_text(data = subset(lec.data.FP2, log10pval.V2V1 >= -log10(0.05)),
            nudge_y = nudge, show.legend = F, fontface = "bold", alpha = a) + 
  theme_bw()



# 
# ########################################### Film Dapivirine #############################################################
# 
# ## Film Dapivirine :: Data Processing
# lec.data.FD2 <- dcast(lec.data.FD, Lectin~Visit, value.var = "value", fun.aggregate = mean)
# 
# lec.data.FD3 <- dcast(lec.data.FD, Lectin~Visit+Sample, value.var = "value", fun.aggregate = mean)
# n.samples.FD <- length(unique(lec.data.FD$Sample))
# 
# 
# lec.data.FD2$V2V1.pval <- apply(lec.data.FD3[,-1], 1, 
#                                 function(x){ t.test(as.numeric(x)[1:n.samples.FD], as.numeric(x)[(n.samples.FD+1):(2*n.samples.FD)],
#                                                     alternative = "two.sided", paired = T, var.equal = F)$p.val})
# 
# lec.data.FD2$V3V1.pval <- apply(lec.data.FD3[,-1], 1, 
#                                 function(x){ t.test(as.numeric(x)[1:n.samples.FD], as.numeric(x)[((2*n.samples.FD)+1):(3*n.samples.FD)],
#                                                     alternative = "two.sided", paired = T, var.equal = F)$p.val})
# 
# lec.data.FD2$log2FC.V2V1 <- lec.data.FD2$Visit.2 - lec.data.FD2$Visit.1
# lec.data.FD2$log2FC.V3V1 <- lec.data.FD2$Visit.3 - lec.data.FD2$Visit.1
# 
# lec.data.FD2$log10pval.V2V1 <- -log10(lec.data.FD2$V2V1.pval)
# lec.data.FD2$log10pval.V3V1 <- -log10(lec.data.FD2$V3V1.pval)
# 
# lec.data.FD2$Threshold <- as.factor(lec.data.FD2$log10pval.V2V1 >= -log10(0.05))
# 
# ## Manualy adding color to the significant lectin labels.
# lec.data.FD2$Lectin.Color <- "lg"
# lec.data.FD2[lec.data.FD2$log10pval.V2V1 >= -log10(0.05), "Lectin.Color"] <- "dg"
# lec.data.FD2[lec.data.FD2$Lectin %in% c("CVN","UDA","GNA_EY"), "Lectin.Color"] <- "b"
# 
# 
# ## Film Dapivirine:: Volcano Plot
# FD.volcanoplot <- ggplot(lec.data.FD2, aes(log2FC.V2V1, log10pval.V2V1, label = Lectin, colour = as.factor(Lectin.Color))) + 
#   geom_point(show.legend = F, size = 2, shape = 21) +  
#   xlab ("Film Dapivirine [log2(V3/V1)]") + 
#   ylab("-log10(p-value)") + 
#   xlim(c(-2,2)) + ylim(c(0,3)) + 
#   scale_color_manual(values = c("lg" = "lightgray", "dg" = "darkgray", "b" = "blue")) +
#   geom_hline(yintercept = -log10(0.05), col = "red", lty = 2) + 
#   geom_text(data = subset(lec.data.FD2, log10pval.V2V1 >= -log10(0.05)),
#             nudge_y = 0.25, show.legend = F, fontface = "bold") + 
#   theme_bw()


#################################################################################################

## Gel Placebo vs Film Placebo Plot


FP.data <- dcast(lec.data.FP, Lectin+Sample~Visit, value.var = "value", fun.aggregate = mean)
FP.data$V2V1 <- FP.data$Visit.2-FP.data$Visit.1
FP.data$V3V1 <- FP.data$Visit.3-FP.data$Visit.1

GP.data <- dcast(lec.data.GP, Lectin+Sample~Visit, value.var = "value", fun.aggregate = mean)
GP.data$V2V1 <- GP.data$Visit.2-GP.data$Visit.1
GP.data$V3V1 <- GP.data$Visit.3-GP.data$Visit.1

FP.data.2 <- melt(FP.data, id.vars = c("Lectin", "Sample"), measure.vars = c("V2V1","V3V1"))
GP.data.2 <- melt(GP.data, id.vars = c("Lectin", "Sample"), measure.vars = c("V2V1","V3V1"))
FP.data.2$Formulation <- "Film"
GP.data.2$Formulation <- "Gel"
data.2 <- rbind(FP.data.2, GP.data.2)

data.3 <- dcast(data.2, Lectin~variable+Formulation+Sample, value.var = "value", fun.aggregate = mean)
data.4 <- dcast(data.2, Lectin~variable+Formulation, value.var = "value", fun.aggregate = mean)

data.4$V2V1.pval <- apply(data.3[,-1], 1, 
                                function(x){ t.test(as.numeric(x)[1:n.samples.FP], as.numeric(x)[(n.samples.FP+1):(n.samples.FP+n.samples.GP)],
                                                    alternative = "two.sided", paired = F, var.equal = F)$p.val})

data.4$V3V1.pval <- apply(data.3[,-1], 1, 
                                function(x){ t.test(as.numeric(x)[(n.samples.FP+n.samples.GP+1):((2*n.samples.FP)+n.samples.GP)], as.numeric(x)[(((2*n.samples.FP)+n.samples.GP)+1):(2*(n.samples.FP+n.samples.GP))],
                                                    alternative = "two.sided", paired = F, var.equal = F)$p.val})

data.4$Color.Label = "lg"
data.4[data.4$Lectin %in% c("SVN","SVN.2","GRFT", "UDA"), "Color.Label"] <- "g"
data.4[data.4$Lectin %in% c("UEA.II"), "Color.Label"] <- "dg"

data.4$Threshold <- as.factor(data.4$V2V1.pval <= 0.05)


## Gel Placebo vs Film Placebo Scatter Plot
GPvsFP.plot <- ggplot(data.4, aes(V2V1_Film, V2V1_Gel, shape = Threshold, label = Lectin, col = Color.Label)) + 
                geom_point(show.legend = F, size = 2) +  
                labs(x = "Film Placebo [log2(V3/V1)]", y = "Gel Placebo [log2(V3/V1)]") + 
                xlim(c(-2,1)) + ylim(c(-2,1)) +
                scale_color_manual(values = c("lg" = "lightgray", "dg" = "darkgray", "g" = "chartreuse3")) +
                geom_text(data = subset(data.4, V2V1.pval <= 0.05 | Lectin %in% c("SVN","SVN..old.","GRFT")),
                          nudge_y = nudge, show.legend = F, fontface = "bold", alpha = a) + 
                geom_abline(slope = 1, intercept = 0, col = "red", lty = 2) +
                # geom_abline(slope = 1.25, intercept = 0, col = "darkgray", lty = 2) + 
                # geom_abline(slope = 0.75, intercept = 0, col = "darkgray", lty = 2) + 
                theme_bw()



###################################################################################################
############################ Microbiome Analysis #################################################

lec.data$Sample.2 <- sapply(lec.data$Sample, function(x){return(as.numeric(strsplit(as.character(x), "S")[[1]][2]))})
lec.data$Uniq.ID <- with(lec.data, paste(Visit, Sample.2, sep = "_"))

## Replace the NA values with 0's
mb.v2[is.na(mb.v2)] <- 0
mb.v3[is.na(mb.v3)] <- 0
mb.v4[is.na(mb.v4)] <- 0

# Replace the lengthy column name for the data sets with Sample ID
colnames(mb.v2)[1] <- "Sample.ID"
colnames(mb.v3)[1] <- "Sample.ID"
colnames(mb.v4)[1] <- "Sample.ID"

## Extract the unique sample IDs from the three data sets from the three visits.
uniq.samples <- intersect(mb.v2$Sample.ID, intersect(mb.v3$Sample.ID, mb.v4$Sample.ID))

## Isolate the data only for the common samples.
mbv2.1 <- mb.v2[mb.v2$Sample.ID %in% uniq.samples, ]
mbv3.1 <- mb.v3[mb.v3$Sample.ID %in% uniq.samples, ]
mbv4.1 <- mb.v4[mb.v4$Sample.ID %in% uniq.samples, ]

## Melt the data sets into long format to add the appropriate visit annotations. 
mbv2.melt <- melt(mbv2.1, id.vars = "Sample.ID")
mbv3.melt <- melt(mbv3.1, id.vars = "Sample.ID")
mbv4.melt <- melt(mbv4.1, id.vars = "Sample.ID")

mbv2.melt$Visit <- "Visit.1"
mbv3.melt$Visit <- "Visit.2"
mbv4.melt$Visit <- "Visit.3"

## combine into a master data frame and generate a unique ID to merge with lec.data annotations
mb.DF <- rbind(mbv2.melt, mbv3.melt, mbv4.melt)
mb.DF$Uniq.ID <- with(mb.DF, paste(Visit, Sample.ID, sep = "_"))

## Merge the microbiome with other annotations (takes about 20 seconds for me)
mb.DF2 <- merge(mb.DF, unique(lec.data[, !(names(lec.data) %in% c("Lectin","value"))]), 
                  by.x = "Uniq.ID", by.y = "Uniq.ID", all.x = F, all.y = F)

## Convert the presence of microbial species to binary units to count them later to generate diversity profiles
mb.DF3 <- mb.DF2
mb.DF3[mb.DF3$value != 0, "value"] <- 1

## Aggregate the microbiome data to add the number of microbial species per sample.
mb.DF4 <- aggregate(value ~ Sample+Group+Formulation+Drug+Sample.ID+Visit.x, data = mb.DF3, FUN = sum)

## Subset the samples with Placebo formulation only
mb.DF5 <- subset(mb.DF4, Drug %in% "Placebo")

## Changing the order of the formulation such that Gel is plotted before Film.
mb.DF5$Formulation.Order <- with(mb.DF5, relevel(Formulation, "Gel"))

## Box plot showing change in # of Microbial Species with each visit in Film and Gel samples.
mb.boxplots <- ggplot(mb.DF5, aes(Formulation.Order, value, col = Visit.x)) +
                geom_jitter(position = position_jitterdodge(jitter.width = 0.1)) +
                geom_boxplot(notch = T, alpha = 0) +
                labs(x = "Placebo", y = "# of Microbial Species") +
                theme_bw()

## Combined Plot
multiplot(GP.volcanoplot, GPvsFP.plot, FP.volcanoplot, mb.boxplots, cols = 2)
