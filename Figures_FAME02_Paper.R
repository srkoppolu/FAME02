library(ggplot2)
library(reshape2)

# WD <- "E:/Dropbox/Mahal_lab/Bioinformatics_Project/CVL_project/FAME02/Microbiome_Data/"
# setwd(WD)

lec.fname <- "LectinData_meltedformat.txt"
lec.data <- read.table(lec.fname, sep = "\t", header = T)

# fname <- "Microbiome_Data_mbvDF2.txt"
# wide.data <- read.table(fname, sep = "\t", header = T)


setwd("E:/Dropbox/Git_Files/FAME02/")

mb.v2 <- read.delim("FAME02_vag_isolates_visit2.txt")
mb.v3 <- read.delim("FAME02_vag_isolates_visit3.txt")
mb.v4 <- read.delim("FAME02_vag_isolates_visit4.txt")

sample.treat <- read.delim("sample_treatment_list.txt")


mb.v2[is.na(mb.v2)] <- 0
mb.v3[is.na(mb.v3)] <- 0
mb.v4[is.na(mb.v4)] <- 0


colnames(mb.v2)[1] <- "Sample.ID"
colnames(mb.v3)[1] <- "Sample.ID"
colnames(mb.v4)[1] <- "Sample.ID"






uniq.samples <- intersect(unique(lec.data$Pt.ID.2), 
                            intersect(mb.v2$Sample.ID, 
                                      intersect(mb.v3$Sample.ID, mb.v4$Sample.ID)))



mbv2.1 <- mb.v2[mb.v2$Sample.ID %in% uniq.samples, -dim(mb.v2)[2]]
mbv3.1 <- mb.v3[mb.v3$Sample.ID %in% uniq.samples, -dim(mb.v3)[2]]
mbv4.1 <- mb.v4[mb.v4$Sample.ID %in% uniq.samples, -dim(mb.v4)[2]]

mbv2.melt <- melt(mbv2.1, id.vars = "Sample.ID")
mbv3.melt <- melt(mbv3.1, id.vars = "Sample.ID")
mbv4.melt <- melt(mbv4.1, id.vars = "Sample.ID")

mbv2.melt$Visit <- "Visit 1"
mbv3.melt$Visit <- "Visit 2"
mbv4.melt$Visit <- "Visit 3"

mb.DF <- rbind(mbv2.melt, mbv3.melt, mbv4.melt)
mb.DF2 <- unique(merge(mb.DF, lec.data[,5:6], 
                by.x = "Sample.ID", by.y = "Pt.ID.2",
                all.x = F, all.y = F))

mb.DF3 <- unique(merge(wide.data[,1:4], mb.DF2, 
                by.x = "Pt.ID", by.y = "Pt.ID",
                all.x = F, all.y = F))


mb.DF4 <- mb.DF3
mb.DF4[mb.DF4$value != 0, "value"] <- 1

mb.DF5 <- aggregate(value ~ Pt.ID+Treatment+Platform+Drug+Sample.ID+Visit, data = mb.DF4, FUN = sum)

# 
# 
# ggplot(mb.DF5, aes(Platform, value, fill = Visit)) +
#   geom_jitter(aes(col = Visit), position = position_jitterdodge(jitter.width = 0.1)) +
#   geom_boxplot(aes(col = Visit), notch = T, alpha = 0) + 
#   theme_bw() + 
#   facet_wrap(~Drug)
# 
# 
ggplot(subset(mb.DF5, Drug %in% "Placebo"), aes(Platform, value, col = Visit)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1)) +
  geom_boxplot(notch = T, alpha = 0) +
  theme_bw()


mb.DF5.pl <- subset(mb.DF5, Drug %in% "Placebo")


x.film.v1 <- subset(mb.DF5.pl, Platform %in% "Film" & Visit %in% "Visit 1", "value")
x.film.v2 <- subset(mb.DF5.pl, Platform %in% "Film" & Visit %in% "Visit 2", "value")
x.film.v3 <- subset(mb.DF5.pl, Platform %in% "Film" & Visit %in% "Visit 3", "value")

x.gel.v1 <- subset(mb.DF5.pl, Platform %in% "Gel" & Visit %in% "Visit 1", "value")
x.gel.v2 <- subset(mb.DF5.pl, Platform %in% "Gel" & Visit %in% "Visit 2", "value")
x.gel.v3 <- subset(mb.DF5.pl, Platform %in% "Gel" & Visit %in% "Visit 3", "value")


p.film.v1v2 <- t.test(x.film.v1, x.film.v2,
                      alternative = "two.sided", paired = F,
                      var.equal = F, conf.level = 0.95)$p.value
p.film.v1v3 <- t.test(x.film.v1, x.film.v3,
                      alternative = "two.sided", paired = F,
                      var.equal = F, conf.level = 0.95)$p.value
p.film.v2v3 <- t.test(x.film.v2, x.film.v3,
                      alternative = "two.sided", paired = F,
                      var.equal = F, conf.level = 0.95)$p.value


p.gel.v1v2 <- t.test(x.gel.v1, x.gel.v2,
                      alternative = "two.sided", paired = F,
                      var.equal = F, conf.level = 0.95)$p.value
p.gel.v1v3 <- t.test(x.gel.v1, x.gel.v3,
                      alternative = "two.sided", paired = F,
                      var.equal = F, conf.level = 0.95)$p.value
p.gel.v2v3 <- t.test(x.gel.v2, x.gel.v3,
                      alternative = "two.sided", paired = F,
                      var.equal = F, conf.level = 0.95)$p.value





###############################################


mb.DF6 <- subset(mb.DF5, !(Drug %in% "Dapivirine" & Visit %in% c("Visit 2", "Visit 3")))

# mean(mb.DF6$value)
# sd(mb.DF6$value)

sd.low <- mean(mb.DF6$value) - sd(mb.DF6$value)
sd.high <- mean(mb.DF6$value) + sd(mb.DF6$value)
hsd.low <- mean(mb.DF6$value) - sd(mb.DF6$value)*0.5
hsd.high <- mean(mb.DF6$value) + sd(mb.DF6$value)*0.5


lec.data[lec.data$Visit %in% "Visit 2", "Visit.2"] <- "Visit 1"
lec.data[lec.data$Visit %in% "Visit 3", "Visit.2"] <- "Visit 2"
lec.data[lec.data$Visit %in% "Visit 4", "Visit.2"] <- "Visit 3"

lec.data$Uniq.ID.2 <- paste(lec.data$Pt.ID, lec.data$Visit.2, sep = "_")

mb.DF5$Uniq.ID <- paste(mb.DF5$Pt.ID, mb.DF5$Visit, sep = "_")



lec.data.2 <- dcast(lec.data, Uniq.ID.2+variable~Name, value.var = "value")

mblec.DF1 <- merge(mb.DF5, lec.data.2[,c("Uniq.ID.2","GRFT","SVN","SVN (old)")],
                   by.x = "Uniq.ID", by.y = "Uniq.ID.2",
                   all.x = F, all.y = F)
names(mblec.DF1)[9:11] <- c("GRFT","rSVN","SVN")


mblec.DF2 <- mblec.DF1
mblec.DF2$Div.label <- "b.Diverse"
mblec.DF2$Div.label.2 <- "b.Diverse"
mblec.DF2[mblec.DF2$value <= sd.low, "Div.label"] <- "a.Less Diverse"
mblec.DF2[mblec.DF2$value >= sd.high, "Div.label"] <- "c.More Diverse"
mblec.DF2[mblec.DF2$value <= hsd.low, "Div.label.2"] <- "a.Less Diverse"
mblec.DF2[mblec.DF2$value >= hsd.high, "Div.label.2"] <- "c.More Diverse"

names(mblec.DF2)[8] <- "Species.Count"

mblec.melt <- melt(mblec.DF2, id.vars = c("Uniq.ID","Pt.ID","Treatment","Platform","Drug",
                                          "Sample.ID","Visit","Species.Count","Div.label","Div.label.2"))

ggplot(mblec.melt, aes(Div.label.2, value, col = Div.label.2)) + 
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1)) +
  geom_boxplot(notch = T, alpha = 0) +
  theme_bw() + 
  facet_wrap(~variable)
    
    
ggplot(mblec.melt, aes(Div.label, value, col = Div.label)) + 
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1)) +
  geom_boxplot(notch = T, alpha = 0) +
  theme_bw() + 
  facet_wrap(~variable)    
    

mblec.melt2 <- mblec.melt[!(mblec.melt$Div.label %in% "b.Diverse"),]
mblec.melt3 <- mblec.melt[!(mblec.melt$Div.label.2 %in% "b.Diverse"),]


ggplot(mblec.melt2, aes(Div.label, value, col = Div.label)) + 
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1)) +
  geom_boxplot(notch = T, alpha = 0) +
  theme_bw() + 
  facet_wrap(~variable)    

ggplot(mblec.melt3, aes(Div.label.2, value, col = Div.label.2)) + 
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1)) +
  geom_boxplot(notch = T, alpha = 0) +
  theme_bw() + 
  facet_wrap(~variable)    




x.1sd.grft.ld <- subset(mblec.melt2, Div.label %in% "a.Less Diverse" & 
                          variable %in% "GRFT", "value")

x.1sd.grft.md <- subset(mblec.melt2, Div.label %in% "c.More Diverse" & 
                          variable %in% "GRFT", "value")


x.1sd.rsvn.ld <- subset(mblec.melt2, Div.label %in% "a.Less Diverse" & 
                          variable %in% "rSVN", "value")

x.1sd.rsvn.md <- subset(mblec.melt2, Div.label %in% "c.More Diverse" & 
                          variable %in% "rSVN", "value")


x.1sd.svn.ld <- subset(mblec.melt2, Div.label %in% "a.Less Diverse" & 
                          variable %in% "SVN", "value")

x.1sd.svn.md <- subset(mblec.melt2, Div.label %in% "c.More Diverse" & 
                          variable %in% "SVN", "value")


p.1sd.grft <- t.test(x = x.1sd.grft.ld, y = x.1sd.grft.md,
                     alternative = "two.sided", paired = F,
                     var.equal = F, conf.level = 0.95)$p.value

p.1sd.rsvn <- t.test(x = x.1sd.rsvn.ld, y = x.1sd.rsvn.md,
                     alternative = "two.sided", paired = F,
                     var.equal = F, conf.level = 0.95)$p.value

p.1sd.svn <- t.test(x = x.1sd.svn.ld, y = x.1sd.svn.md,
                     alternative = "two.sided", paired = F,
                     var.equal = F, conf.level = 0.95)$p.value
 
# p.aov.1sd.grft <- summary(aov(value~Div.label, subset(mblec.melt2, variable %in% "GRFT")))[[1]]["Pr(>F)"][[1]][1]
# p.aov.1sd.rsvn <- summary(aov(value~Div.label, subset(mblec.melt2, variable %in% "rSVN")))[[1]]["Pr(>F)"][[1]][1]
# p.aov.1sd.svn <- summary(aov(value~Div.label, subset(mblec.melt2, variable %in% "SVN")))[[1]]["Pr(>F)"][[1]][1]



x.hsd.grft.ld <- subset(mblec.melt3, Div.label.2 %in% "a.Less Diverse" & 
                          variable %in% "GRFT", "value")

x.hsd.grft.md <- subset(mblec.melt3, Div.label.2 %in% "c.More Diverse" & 
                          variable %in% "GRFT", "value")


x.hsd.rsvn.ld <- subset(mblec.melt3, Div.label.2 %in% "a.Less Diverse" & 
                          variable %in% "rSVN", "value")

x.hsd.rsvn.md <- subset(mblec.melt3, Div.label.2 %in% "c.More Diverse" & 
                          variable %in% "rSVN", "value")


x.hsd.svn.ld <- subset(mblec.melt3, Div.label.2 %in% "a.Less Diverse" & 
                         variable %in% "SVN", "value")

x.hsd.svn.md <- subset(mblec.melt3, Div.label.2 %in% "c.More Diverse" & 
                         variable %in% "SVN", "value")


p.hsd.grft <- t.test(x = x.hsd.grft.ld, y = x.hsd.grft.md,
                     alternative = "two.sided", paired = F,
                     var.equal = F, conf.level = 0.95)$p.value

p.hsd.rsvn <- t.test(x = x.hsd.rsvn.ld, y = x.hsd.rsvn.md,
                     alternative = "two.sided", paired = F,
                     var.equal = F, conf.level = 0.95)$p.value

p.hsd.svn <- t.test(x = x.hsd.svn.ld, y = x.hsd.svn.md,
                    alternative = "two.sided", paired = F,
                    var.equal = F, conf.level = 0.95)$p.value








