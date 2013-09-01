########
# All data analysis for Lloyd et al
########


require(ggplot2)
require(reshape2)
require(plyr)
require(grid) 
require(reshape2)
require(gridExtra)
require(xlsx)
require(boot)
require(RColorBrewer)
require(segmented)

# remove this line before publication
setwd("~/Dropbox/Metadata Analysis/")

#====================
# global figure settings
#====================
clm <- 3.42
clm2 <- 7.08
#ps <- 2 point size should be diffeerent in different figures

theme_set(theme_bw() + 
            theme(panel.grid.major=element_blank(),
                  panel.grid.minor=element_blank(),
                  panel.background=element_rect(colour="black"),
                  text=element_text(size=8)))



na.strings <- c("")

#======================
# Read & munge sediment data
#=====================

raw_seds <- read.xlsx("SI_database_sediments.xlsx", sheetName="Compilation")

#=================
# Fix column names
#=================
names(raw_seds)[names(raw_seds)=="Depth..m."] <- "depth"
attr(raw_seds$depth, "units") <- "meters"

names(raw_seds)[names(raw_seds)=="Cells.per.cc"] <- "totalcells"
attr(raw_seds$totalcells, "units") <- "cells per cc"

names(raw_seds)[names(raw_seds)=="qPCR.Total.per.cc"] <- "qPCRtotal"
attr(raw_seds$qPCRtotal, "units") <- "qPCR total copies per cc"

names(raw_seds)[names(raw_seds)=="Fraction.Arc.qPCR"] <- "percentqPCR"
attr(raw_seds$percentqPCR, "units") <- "fraction of 1"

#names(raw_seds)[names(raw_seds)=="X.2.mM.sulfate."] <- "greater_2mM_sulfate"
#raw_seds$greater_2mM_sulfate <- as.logical(raw_seds$greater_2mM_sulfate)

names(raw_seds)[names(raw_seds)=="Sulfate..mM."] <- "sulfate"
attr(raw_seds$sulfate, "units") <- "mM"

names(raw_seds)[names(raw_seds)=="qPCR.with.universal.primers.per.cc"] <- "qPCRuniversal"
attr(raw_seds$qPCRuniversal, "units") <- "qPCR with universal primers, per cc"

names(raw_seds)[names(raw_seds)=="CARDFISH.Arc.below.detection.limit.or.not.measured."] <- "CARDFISH.Arc.nd"


#============
# Convert columns into the correct data type
#============
#raw_seds$greater_2mM_sulfate <- as.logical(raw_seds$greater_2mM_sulfate)
#raw_seds$sulfate <- as.numeric(as.character(raw_seds$sulfate))

raw_seds$Mud.volcano.or.seep.[raw_seds$Mud.volcano.or.seep.==FALSE | is.na(raw_seds$Mud.volcano.or.seep.)] <- FALSE
raw_seds$Environment.Type <- as.character(raw_seds$Environment.Type)
raw_seds$Environment.Type[is.na(raw_seds$Environment.Type)] <- ""



corrected_seds <- raw_seds

#corrected_seds$paper <- as.character(corrected_seds$paper)
names(corrected_seds)[names(corrected_seds)=="Fraction.Arc.qPCR.out.of.bac.plus.arc"] <- "percentqPCR"
#corrected_seds$location <- as.character(corrected_seds$location)


#==========
# Create some useful columns
#==========

# Convert all zero depths to 1 cm & log-transform the depth
corrected_seds$depth[(corrected_seds$depth<=0.01)&(is.na(corrected_seds$depth)==FALSE)] <- 0.01
attr(corrected_seds$depth, "units") <- "m"
corrected_seds$depth_log10 <- log10(corrected_seds$depth)


# Append a number corresponding to each paper, to indicate the paper the core came from
corrected_seds$paperNumber <- as.numeric(corrected_seds$paper)
corrected_seds$core <- paste(as.character(corrected_seds$core), " [", corrected_seds$paperNumber, "]", sep="")
corrected_seds$Uses.516.for.Arc[is.na(corrected_seds$Uses.516.for.Arc)] <- FALSE

paperTable <- ddply(corrected_seds, .(paper), function(x) x[1, c("paper", "paperNumber")])

savePaperTable <- FALSE

if(savePaperTable) {
  pdf("sed_refs_table.pdf", height=20, width=10)
  grid.table(paperTable)
  dev.off()
}


#=============
# Read & munge seawater data
#=============

# div/0! and blanks to be read as NA  - what is source of 'div/0'?
na.strings <- c("#DIV/0!", "")
water_raw_data <- read.xlsx("SI_database_seawater.xlsx", sheetName="Compilation")
corrected_sw <- water_raw_data

names(corrected_sw)[names(corrected_sw)=="Depth..m."] <- "depth"
# What is 'Sample'?
names(corrected_sw)[names(corrected_sw)=="qPCR.Bacteria..copies.mL.water."] <- "qPCRbac"
attr(corrected_sw$qPCRbac, "units") <- "copiesPerMl"

names(corrected_sw)[names(corrected_sw)=="qPCR.Archaea..copies.mL.water."] <- "qPCRarc"
attr(corrected_sw$qPCRarc, "units") <- "copiesPerMl"

names(corrected_sw)[names(corrected_sw)=="Total.qPCR.copies.mL.water."] <- "qPCRtotal"
attr(corrected_sw$qPCRarc, "units") <- "copiesPerMl"

names(corrected_sw)[names(corrected_sw)=="Cells.per.cc"] <- "totalcells"
attr(corrected_sw$qPCRarc, "units") <- "cellsPerMl"

#names(corrected_sw)[names(corrected_sw)=="Fraction.Arc.CARDFISH"] <- "percentArcCARD"
#attr(corrected_sw$qPCRarc, "units") <- "percent"

# Convert all seawater "0" depths less than 1 m (mostly zeros) to 1 m to prevent log problems
corrected_sw$depth[(corrected_sw$depth<=0)&(is.na(corrected_sw$depth)==FALSE)] <- 1

# Convert all seawater Archaea "0" counts to 1 to prevent log problems
corrected_sw$CARDFISH.Arc.per.cc[corrected_sw$CARDFISH.Arc.per.cc<=0 & !is.na(corrected_sw$CARDFISH.Arc.per.cc)] <- 1

# Log-transform the depths
corrected_sw$depth_log10 <- log10(corrected_sw$depth)

# Convert all of the Schattenhoffer permeabilization data to "proteinase K"
corrected_sw$Arc.permeabilization[corrected_sw$paper=="Schattenhofer 2009"] <- "proteinase K"

#===========
# Prep sed and sw data frames for merge
#===========
seds_for_merge <- corrected_seds[ , c("paper", "core", "depth", "totalcells", "CARDFISH.Bac.per.cc", "CARDFISH.Arc.per.cc",
                                      "CARDFISH.Total.per.cc", "FISH.yield", "qPCR.Bac.per.cc", "qPCR.Arc.per.cc", "qPCRtotal", 
                                      "qPCRuniversal", "Fraction.Arc.CARDFISH", "percentqPCR", "Fish.or.cardFish",
                                      "Cell.stain", "Fixative", "Bac.permeabilization", "Arc.permeabilization",
                                      "SYBR.vs.Taqman", "DNA.extraction.procedure")]
#names(seds_for_merge)[names(seds_for_merge)=="depths"] <- "depth"
#names(seds_for_merge)[names(seds_for_merge)=="CARDFISH.Bac.per.cc"] <- "CARDbac"
#names(seds_for_merge)[names(seds_for_merge)=="CARDFISH.Arc.per.cc"] <- "CARDarc"
#names(seds_for_merge)[names(seds_for_merge)=="qPCR.Bac.per.cc"] <- "qPCRbac"
#names(seds_for_merge)[names(seds_for_merge)=="qPCR.Arc.per.cc"] <- "qPCRarc"
#names(seds_for_merge)[names(seds_for_merge)=="X.Arc.CARFDISH"] <- "percentCARDFISH"

seds_for_merge$environment <- "sediments"


sw_for_merge <- corrected_sw[ , c("paper", "depth", "totalcells", "qPCRbac", "qPCRarc", "qPCRtotal", "Cell.stain", 
                                  "CARDFISH.Bac.per.cc", "CARDFISH.Arc.per.cc", "CARDFISH.Total.per.cc", 
                                  "FISH.yield", "Fraction.Arc.CARDFISH", "Fish.or.cardFish", 
                                  "Cell.stain", "Fixative", 
                                  "Bac.permeabilization", "Arc.permeabilization", 
                                  "Bac.probe", "Arc.probe")]
sw_for_merge$environment <- "seawater"

allData <- merge(sw_for_merge, seds_for_merge, all=TRUE)

# Reorder environment vector so that sw always comes out in front of seds
allData$environment <- factor(allData$environment, levels=c("seawater", "sediments"), ordered=TRUE)

#========
# Fix up allData a little bit
#========
#allData$paper <- factor(as.character(allData$paper))
#allData <- allData[-which(is.na(allData$paper)), ]

#===========
# Clean up workspace
#===========
rm(na.strings)
rm(raw_seds)
rm(water_raw_data)

##############
# Figure 1: total prokaryotes, by CARDFISH, vs. total cells, by direct counts, in sw and sediments
##############

# Notes: change labels (sediments & seawater), change order of sed vs sw
# Fig 1 excludes polyribonucleotide-FISH
fig1data <- allData[!is.na(allData$totalcells) & 
                      !is.na(allData$CARDFISH.Total.per.cc) &
                      !is.na(allData$Fish.or.cardFish) & 
                      !is.na(allData$Arc.permeabilization), ]

# Create vector of permeabilization method and/or FISH 
fig1data$permForFig1 <- NA
fig1data$permForFig1[fig1data$Fish.or.cardFish == "FISH"] <- "FISH"
fig1data$permForFig1[fig1data$Fish.or.cardFish == "polyribonucFISH"] <- "polyribonucFISH"
fig1data$permForFig1[fig1data$Fish.or.cardFish == "CARDFISH" & fig1data$Arc.permeabilization == "proteinase K"] <- "CARD_w_protK"
fig1data$permForFig1[fig1data$Fish.or.cardFish == "CARDFISH" & fig1data$Arc.permeabilization != "proteinase K"] <- "CARD_other_perm"
fig1data <- fig1data[fig1data$totalcells!=0 &
                       fig1data$CARDFISH.Total.per.cc != 0, ]

fig1data$permForFig1 <- factor(fig1data$permForFig1, 
                               levels=c("FISH", "CARD_w_protK", "CARD_other_perm"), ordered=TRUE)

# Create a 1:1 line
oneOne <- data.frame(x=c(1e4, 1e11), y=c(1e4, 1e11))

# How many data points from each method and environment
nPoints <- ddply(fig1data, c("environment", "permForFig1"), nrow)
names(nPoints)[names(nPoints)=="V1"] <- "nPoints"
nPoints$nPoints <- paste("n=", nPoints$nPoints, sep="")
nPoints$x <- 1e6
nPoints$y <- rep(c(1e10, 1e9, 1e8), 2)

fig1col <- brewer.pal(n=3, name="Paired")
fig1col[3] <- brewer.pal(n=3, name="OrRd")[3]
#fig1col <- fig1col[c(1, 3, 2)]


fig1 <- ggplot() + 
  geom_point(data=fig1data, aes(x=totalcells, y=CARDFISH.Total.per.cc, colour=permForFig1), size=0.75, alpha=1) +
  geom_line(data=oneOne, aes(x=x, y=y), colour="black") +
  geom_text(data=nPoints, aes(x=x, y=y, label=nPoints, colour=permForFig1), size=2, hjust=0) +
  scale_x_log10(expression(paste("Total Cells (cells ", ml^{-1}, ")")), breaks=10^(4:11)) + 
  scale_y_log10(expression(paste("FISH or CARD-FISH\nSum of Bac + Arc (cells ", ml^{-1}, ")")), breaks=10^(4:11)) +
  scale_colour_manual(values=fig1col, 
                      name="Method",
                      breaks=c("FISH", "CARD_w_protK", "CARD_other_perm"),
                      labels=c("FISH", "CARD-FISH,\nproteinase K", "CARD-FISH,\nother permeabilization")
  ) +
  coord_fixed() +
  facet_wrap(~environment) +
  guides(colour = guide_legend(override.aes = list(size=3, alpha=1))) +
  theme(legend.position="top",
        axis.title.y=element_text(hjust=0.5),
        axis.text.x=element_text(angle=-45, hjust=0))
print(fig1)
#ggsave("fig1_tosubmit.tiff", height=3, width=clm, dpi=900, compression="lzw")
tiff("fig1_tosubmit.tif", height=3, width=clm, units="in", res=900, compression="lzw")
print(fig1)
dev.off()


################
# Fig 2: FISH / CARDFISH yield in sediments, by permeabilization method and core
###############

yield_label <- expression(frac("Bacteria + Archaea (FISH or CARD-FISH)", "Total Cell Count (DAPI, AO or SYBR Green)"))
oldTheme <- theme_get()
theme_set(theme_get()) #+ theme(axis.title=element_text(size=4)))
yaxsize <- 7

fig2a_data <- corrected_seds[!is.na(corrected_seds$Arc.permeabilization) &
                               !is.na(corrected_seds$FISH.yield) &
                               !is.na(corrected_seds$Fish.or.cardFish) &
                               !corrected_seds$Arc.permeabilization == "CARD-FISH, archaea not measured" &
                               !corrected_seds$Arc.permeabilization == "FISH, archaea not measured" &
                               (corrected_seds$Fish.or.cardFish == "CARDFISH" | corrected_seds$Fish.or.cardFish=="FISH") &
                               corrected_seds$Environment.Type != "Intertidal" &
                               corrected_seds$Environment.Type != "Salt marsh", ]


# Make a permeabilization column, to reflect that FISH doesn't need a permeabilization method
fig2a_data$perm <- as.character(fig2a_data$Arc.permeabilization)
fig2a_data$perm[fig2a_data$Fish.or.cardFish=="FISH"] <- "FISH"

# Re-level the permeabilization methods
fig2a_data$perm <- factor(fig2a_data$perm, 
                          levels=c("proteinase K", "FISH", 
                                   "detergent", "unknown",  
                                   "lysozyme/achromopeptidase", "lysozyme"), 
                          labels=c("proteinase K", "none (FISH)", "detergent", 
                                   "unknown", "lysozyme/\nachromopeptidase", "lysozyme"),
                          ordered=TRUE)

### FIg 2a reduced: cut out all the cores with only one data point
# Count the # points in the core, remove all < 1
fig2a_data2 <- ddply(fig2a_data, .(core), transform, count=length(FISH.yield))
fig2a_data2 <- fig2a_data2[fig2a_data2$count > 1, ]
fig2a_data2$core <- as.factor(as.character(fig2a_data2$core))

# Order the papers in order of decreasing median yield
fig2a_dataMedianYields <- ddply(fig2a_data2, .(core), summarise, medianYield=median(FISH.yield, na.rm=TRUE))
fig2a_dataMedianYields <- fig2a_dataMedianYields[order(fig2a_dataMedianYields$medianYield, decreasing=TRUE), ]
fig2a_data2$core <- factor(fig2a_data2$core, levels=fig2a_dataMedianYields$core, ordered=TRUE)


# Determine the quartiles of yield from SW data; overplot that on sediments plot
#swQuant <- quantile(allData$FISH.yield[allData$environment=="seawater"], probs=c(0.5-0.341, 0.5, 0.5+0.341), na.rm=TRUE)
swQuant <- quantile(corrected_sw$FISH.yield, probs=c(0.25, 0.5, 0.75), na.rm=TRUE)
swYieldDF <- data.frame(ymin=swQuant[1], ymed=swQuant[2], ymax=swQuant[3])
rownames(swYieldDF) <- NULL
swYieldDF <- data.frame(x=c(0, length(unique(fig2a_data2$core))+0.5), rbind(swYieldDF, swYieldDF))

## Count number of points in each core
nPoints <- ddply(fig2a_data2, .(core), function(x) data.frame(nrow=nrow(x), perm=x$perm[1]))
nPoints$core <- factor(nPoints$core, levels=levels(fig2a_data2$core), ordered=TRUE)
nPoints$perm <- factor(nPoints$perm, levels=levels(fig2a_data2$perm), ordered=TRUE)

envDF <- ddply(fig2a_data2, .(core), function(x) as.character(x$Environment.Type[1]))

# Set the text position
colors <- brewer.pal(n=6, name="RdBu")
colors[4] <- "#999999"
textPos <- 1.5

fig2a <- ggplot() +
  geom_blank(data=fig2a_data2, aes(x=core, y=FISH.yield, fill=perm)) +
  geom_ribbon(data=swYieldDF, aes(x=x, ymin=ymin, ymax=ymax), fill="skyblue", alpha=1) +
  geom_line(data=swYieldDF, aes(x=x, y=ymed), colour="black") +
  geom_boxplot(data=fig2a_data2, aes(x=core, y=FISH.yield, fill=perm), colour="black", outlier.size=1) +
  geom_rect(data=fig2a_data2, aes(xmin=as.numeric(core)-0.5, xmax=as.numeric(core)+0.5, ymin=textPos-0.05, ymax=textPos+0.05, fill=perm)) +
  geom_text(data=nPoints, aes(x=core, y=textPos, label=nrow), size=2) +
  geom_text(data=envDF, aes(x=core, y=textPos-0.075, label=V1), size=2, angle=-90, hjust=0) +
  ylab(yield_label) +
  coord_cartesian(ylim=c(0,1.6)) + 
  scale_fill_manual(values=colors, name="Archaeal\nPermeabilization Method") +
  theme(axis.text.x=element_text(angle=-45, hjust=0, size=6),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=yaxsize),
        legend.position="top",
        legend.title.align=1,
        plot.margin=unit(c(0, 0.9, 0, 0), "in")) #default margin is c(1, 1, 0.5, 0.5) 
#print(fig2a)

##################
# Fig 2b : Boxplot of yields by permeabilization method
#################

# Spacing around the edges of panels 2b and 2c
fig2bc_spacing = c(0, 0.1, 0.1, 0.1)

# Y-axis label
yl <- expression(frac("Archaea (CARD-FISH)", "Bacteria + Archaea (CARD-FISH)"))


# Want to replace with boxplot of all cardfish % Archaea + Fish, colored by permeabilization method, for sediments only
fig2b_data <- corrected_seds[!is.na(corrected_seds$Arc.permeabilization) &
                               !is.na(corrected_seds$Fish.or.cardFish) &
                               !is.na(corrected_seds$Fraction.Arc.CARDFISH) &
                               !corrected_seds$Arc.permeabilization == "unknown" &
                               !corrected_seds$Arc.permeabilization == "not measured"  &
                               !corrected_seds$Fish.or.cardFish == "FISH" &
                               corrected_seds$Environment.Type != "Intertidal" &
                               corrected_seds$Environment.Type != "Salt marsh", ]

fig2b_data$perm <- as.character(fig2b_data$Arc.permeabilization)
#fig2b_data$perm[fig2b_data$Fish.or.cardFish=="FISH"] <- "FISH"
fig2b_data$perm <- factor(fig2b_data$perm, 
                          levels=c("proteinase K", "detergent", 
                                   "lysozyme/achromopeptidase", "lysozyme"),
                          labels=c("proteinase K", "detergent", 
                                   "lysozyme/\nachromopeptidase", "lysozyme"),
                          ordered=TRUE)

fig2b_medians <- ddply(fig2b_data, .(perm), summarise, med=median(Fraction.Arc.CARDFISH))
fig2b_data$perm <- factor(fig2b_data$perm, levels=fig2b_medians$perm[order(fig2b_medians$med, decreasing=TRUE)])

pointsAndStudies <- ddply(fig2b_data, .(perm),
                          function(x) data.frame(nPoints=nrow(x), 
                                                 nStudies=length(unique(x$paper))))
pointsAndStudies$label <- paste(pointsAndStudies$nPoints, " Points\n", pointsAndStudies$nStudies, " Studies", sep="")

# Do ANOVA on rank-transformed data (b/c is obviously vy non-normal)
fig2b_data$rank.Fraction.Arc.CARDFISH <-rank(fig2b_data$Fraction.Arc.CARDFISH)
fig2b_rm <- aov(rank.Fraction.Arc.CARDFISH ~ perm, data=fig2b_data)
summary(fig2b_rm)
TukeyHSD(fig2b_rm)
# proteinase K: A
# Lysozyme/achromopeptidase: B
# Detergent: B
# Lysozyme: C

# Add the letters to a data frame in order to display them
pointsAndStudies$letter <- NA
pointsAndStudies[pointsAndStudies$perm=="proteinase K", "letter"] <- "a"
pointsAndStudies[pointsAndStudies$perm=="lysozyme/\nachromopeptidase", "letter"] <- "b"
pointsAndStudies[pointsAndStudies$perm=="detergent", "letter"] <- "b"
pointsAndStudies[pointsAndStudies$perm=="lysozyme", "letter"] <- "c"



fig2b <- ggplot() + 
  geom_boxplot(data=fig2b_data, aes(x=perm, y=Fraction.Arc.CARDFISH, fill=perm), outlier.size=1) +
  geom_text(data=pointsAndStudies, aes(x=perm, y=1.2, label=label), vjust=1, size=2) +
  geom_text(data=pointsAndStudies, aes(x=perm, y=-0.03, label=letter), vjust=1, size=2.5) +
  scale_y_continuous(yl, limits=c(-0.1, 1.25), breaks=seq(from=0, to=1, by=0.25)) +
  scale_x_discrete("All Data") +
  scale_fill_manual(values=colors[c(1, 5, 3, 6)], name="permeabilization\nmethod") +
  theme(legend.position="none",
        text=element_text(size=8), 
        axis.text.x=element_text(angle=0, hjust=0.5),
        axis.title.y=element_text(size=yaxsize),
        plot.margin=unit(fig2bc_spacing, "in")) 
print(fig2b)
#ggsave("fig2b_boxplot.png", height=6, width=8, units="in", dpi=300)

#=========
# Fig 2c: boxplots of percent ARC for only seds < 1 m deep
#=========

fig2c_data <- fig2b_data[fig2b_data$depth <= 1, ]

# Rank-order ANOVA tests
fig2c_data$rank.Fraction.Arc.CARDFISH <-rank(fig2c_data$Fraction.Arc.CARDFISH)
fig2c_rm <- aov(rank.Fraction.Arc.CARDFISH ~ perm, data=fig2c_data)
summary(fig2c_rm)
TukeyHSD(fig2c_rm)
# Proteinase K: A
# Lysozyme: B
# Lysozyme/achromopeptidase: B
# Detergent: B


# Calculate the number of data points and studies
pointsAndStudies1m <- ddply(fig2c_data, .(perm),
                            function(x) data.frame(nPoints=nrow(x), 
                                                   nStudies=length(unique(x$paper))))
pointsAndStudies1m$label <- paste(pointsAndStudies1m$nPoints, " Points\n", pointsAndStudies1m$nStudies, " Studies", sep="")

# Add letter labels for signifciantly different permeabilization methods
pointsAndStudies1m$letter <- NA
pointsAndStudies1m[pointsAndStudies1m$perm=="proteinase K", "letter"] <- "a"
pointsAndStudies1m[pointsAndStudies1m$perm=="lysozyme/\nachromopeptidase", "letter"] <- "b"
pointsAndStudies1m[pointsAndStudies1m$perm=="detergent", "letter"] <- "b"
pointsAndStudies1m[pointsAndStudies1m$perm=="lysozyme", "letter"] <- "b"




fig2c <- ggplot() + 
  geom_boxplot(data=fig2c_data, aes(x=perm, y=Fraction.Arc.CARDFISH, fill=perm), outlier.size=1) +
  geom_text(data=pointsAndStudies1m, aes(x=perm, y=1.2, label=label), vjust=1, size=2) +
  geom_text(data=pointsAndStudies1m, aes(x=perm, y=-0.03, label=letter), vjust=1, size=2.5) +
  scale_y_continuous(yl,limits=c(-0.1, 1.25), breaks=seq(from=0, to=1, by=0.25)) +
  scale_x_discrete("Shallower than 1 m") +
  scale_fill_manual(values=colors[c(1, 5, 3, 6)]) +
  theme(legend.position="none",
        text=element_text(size=8), 
        axis.text.x=element_text(angle=0, hjust=0.5),
        axis.text.x=element_blank(),
        axis.title.y=element_text(size=yaxsize),
        plot.margin=unit(fig2bc_spacing, "in")) 
print(fig2c)
#ggsave("fig2c_boxplot.png", height=6, width=8, units="in", dpi=300)



#png("fig2_integrated.png", width=clm2, height=8, res=300, units="in")
#print(fig2a, vp=vp1)
#print(fig2b, vp=vp2b)
#print(fig2c, vp=vp2c)
#dev.off()




#====================
# Fig 2d and 2e: same as Fig 2b and c, but with yield
#====================

fig2d_data <- fig2b_data <- corrected_seds[!is.na(corrected_seds$Arc.permeabilization) &
                                             !is.na(corrected_seds$Fish.or.cardFish) &
                                             !is.na(corrected_seds$FISH.yield) &
                                             !corrected_seds$Arc.permeabilization=="unknown" &
                                             !corrected_seds$Arc.permeabilization=="not measured"  &
                                             !corrected_seds$Fish.or.cardFish=="FISH" &
                                             corrected_seds$Environment.Type != "Intertidal" &
                                             corrected_seds$Environment.Type != "Salt marsh", ]

#fig2d_data <- fig2b_data
fig2d_data$perm <- as.character(fig2d_data$Arc.permeabilization)
#fig2b_data$perm[fig2b_data$Fish.or.cardFish=="FISH"] <- "FISH"
fig2d_data$perm <- factor(fig2d_data$perm, 
                          levels=c("proteinase K", "lysozyme/achromopeptidase", "detergent", "lysozyme"), 
                          labels=c("proteinase K", "lysozyme/\nachromopeptidase", "detergent", "lysozyme"),
                          ordered=TRUE)

fig2e_data <- fig2d_data[!is.na(fig2d_data$depth) &
                           fig2d_data$depth <= 1, ]

# ANOVA on rank-transformed data
fig2d_data$rank.FISH.yield <-rank(fig2b_data$FISH.yield)
fig2d_rm <- aov(rank.FISH.yield ~ perm, data=fig2d_data)
summary(fig2d_rm)
TukeyHSD(fig2d_rm)

# proteinase K gets A
# lysozyme-achromopeptidase gets B
# detergent gets B
# lysozyme gets C
pointsAndStudiesYield <- ddply(fig2d_data, .(perm),
                               function(x) data.frame(nPoints=nrow(x), 
                                                      nStudies=length(unique(x$paper))))

pointsAndStudiesYield$label <- paste(pointsAndStudies$nPoints, " Points\n", pointsAndStudies$nStudies, " Studies", sep="")
pointsAndStudiesYield$letter <- NA
pointsAndStudiesYield[pointsAndStudiesYield$perm=="proteinase K", "letter"] <- "a"
pointsAndStudiesYield[pointsAndStudiesYield$perm=="lysozyme/\nachromopeptidase", "letter"] <- "b"
pointsAndStudiesYield[pointsAndStudiesYield$perm=="detergent", "letter"] <- "b"
pointsAndStudiesYield[pointsAndStudiesYield$perm=="lysozyme", "letter"] <- "c"

swDF <- data.frame(perm=c(0, 1, 2, 5), #perm = levels(fig2d_data$perm), #
                   ymin=rep(swQuant[1], nlevels(fig2d_data$perm)), 
                   ymed=rep(swQuant[2], nlevels(fig2d_data$perm)),
                   ymax=rep(swQuant[3], nlevels(fig2d_data$perm)))

fig2d <- ggplot() + 
  geom_blank(data=fig2d_data, aes(x=perm, y=FISH.yield, fill=perm)) +
  geom_ribbon(data=swDF, aes(x=perm, ymin=ymin, ymax=ymax), fill="skyblue") +
  geom_line(data=swDF, aes(x=perm, y=ymed)) +
  geom_boxplot(data=fig2d_data, aes(x=perm, y=FISH.yield, fill=perm), outlier.size=1) +
  geom_text(data=pointsAndStudiesYield, aes(x=perm, y=1.4, label=label), vjust=1, size=2) +
  geom_text(data=pointsAndStudiesYield, aes(x=perm, y=0, label=letter), vjust=1, size=2.5) +
  #scale_y_continuous("yield", limits=c(-0.1, 1.25)) +
  coord_cartesian(xlim=c(0.25, 4.75), ylim=c(-0.1, 1.5)) +
  xlab("All data") +
  ylab(yield_label) +
  scale_fill_manual(values=colors[c(1, 5, 3, 6)], name="permeabilization\nmethod") +
  theme(legend.position="none",
        text=element_text(size=8), 
        axis.text.x=element_text(angle=0, hjust=0.5),
        axis.text.x=element_blank(),
        axis.title.y=element_text(size=yaxsize),
        plot.margin=unit(fig2bc_spacing, "in"),
        axis.title.x=element_blank())  
#print(fig2d)

# ANOVA on rank-transformed data (shallower than 1 m)
fig2e_data$rank.FISH.yield <-rank(fig2e_data$FISH.yield)
fig2e_rm <- aov(rank.FISH.yield ~ perm, data=fig2e_data)
summary(fig2e_rm)
TukeyHSD(fig2e_rm)

# proteinase K gets A
# lysozyme-achromopeptidase gets B
# detergent gets B
# lysozyme gets C
pointsAndStudiesYield1m <- ddply(fig2e_data, .(perm),
                                 function(x) data.frame(nPoints=nrow(x), 
                                                        nStudies=length(unique(x$paper))))

pointsAndStudiesYield1m$label <- paste(pointsAndStudies$nPoints, " Points\n", pointsAndStudies$nStudies, " Studies", sep="")
pointsAndStudiesYield1m$letter <- NA
pointsAndStudiesYield1m[pointsAndStudiesYield1m$perm=="proteinase K", "letter"] <- "a"
pointsAndStudiesYield1m[pointsAndStudiesYield1m$perm=="lysozyme/\nachromopeptidase", "letter"] <- "b"
pointsAndStudiesYield1m[pointsAndStudiesYield1m$perm=="detergent", "letter"] <- "b"
pointsAndStudiesYield1m[pointsAndStudiesYield1m$perm=="lysozyme", "letter"] <- "c"


fig2e <- ggplot() + 
  geom_blank(data=fig2e_data, aes(x=perm, y=FISH.yield, fill=perm)) +
  geom_ribbon(data=swDF, aes(x=perm, ymin=ymin, ymax=ymax), fill="skyblue") +
  geom_line(data=swDF, aes(x=perm, y=ymed)) +
  geom_boxplot(data=fig2e_data, aes(x=perm, y=FISH.yield, fill=perm), outlier.size=1) +
  geom_text(data=pointsAndStudies1m, aes(x=perm, y=1.4, label=label), vjust=1, size=2) +
  geom_text(data=pointsAndStudies1m, aes(x=perm, y=0, label=letter), vjust=1, size=2.5) +
  #scale_y_continuous("yield", limits=c(-0.1, 1.25)) +
  coord_cartesian(xlim=c(0.25, 4.75), ylim=c(-0.1, 1.5)) +
  #xlab("Shallower than 1 m") +
  ylab(yield_label) +
  scale_fill_manual(values=colors[c(1, 5, 3, 6)], name="permeabilization\nmethod") +
  theme(legend.position="none",
        text=element_text(size=8), 
        axis.text.x=element_text(angle=0, hjust=0.5),
        axis.text.x=element_blank(),
        axis.title.y=element_text(size=yaxsize),
        plot.margin=unit(fig2bc_spacing, "in"),
        plot.margin=unit(c(0, 0, 0, 0), "in"),
        axis.title.x=element_blank())  
#print(fig2e)


# Print 2d and 2e as supplemental figure
#vp2d <- viewport(x=0.25, y=0.5, width=0.5, height=1)
#vp2e <- viewport(x=0.75, y=0.5, width=0.5, height=1)

# png("fig2d_and_e_yield.png", height=4, width=clm2, units="in", res=300)
# print(fig2d, vp=vp2d)
# print(fig2e, vp=vp2e)
# dev.off()

# Print all the figures as 1 giant fig
vp1 <- viewport(width=1, height=0.5, x=0.5, y=0.75)
vp2b <- viewport(width=0.5, height=0.25, x=0.25, y=1/8)
vp2c <- viewport(width=0.5, height=0.25, x=0.75, y=1/8)
vp2d <- viewport(width=0.5, height=0.25, x=0.25, y=0.375)
vp2e <- viewport(width=0.5, height=0.25, x=0.75, y=3/8)

tiff("fig2_tosubmit.tif", width=clm2, height=8, 
     res=300, units="in", compression="lzw")
print(fig2a, vp=vp1)
print(fig2b, vp=vp2b)
print(fig2c, vp=vp2c)
print(fig2d, vp=vp2d)
print(fig2e, vp=vp2e)
dev.off()

theme_set(oldTheme)

######
# Fig 3: qPCR quantifications
######
ps <- 0.75

fig3a_df <- corrected_seds[!is.na(corrected_seds$qPCR.Bac.per.cc) &
                             !is.na(corrected_seds$qPCR.Arc.per.cc) &
                             !is.na(corrected_seds$qPCRuniversal) ,]

oneOne <- data.frame(x=10^(3:9), y=10^(3:9))

# Not colored
fig3a <- ggplot(fig3a_df, aes(x=qPCRuniversal, y=(qPCR.Bac.per.cc + qPCR.Arc.per.cc))) + 
  geom_point(size=ps) +
  geom_smooth(method="lm", se=TRUE, colour="black", linetype=2) +
  geom_abline(slope=1, intercept=0) +
  scale_x_log10(breaks=10^c(3:9)) +
  scale_y_log10(breaks=10^c(3:9)) +
  ylab(expression(atop(paste("Gene Copies ", ml^{-1}), "Bacterial + Archaeal primers"))) +
  #xlab(expression(atop(paste("Gene Copies ", ml^{-1}), "Universal Primers"))) +
  xlab(expression(paste("Gene Copies ", ml^{-1}, ", Universal Primers"))) +
  coord_fixed() +
  theme(text=element_text(size=8),
        axis.text.x=element_text(angle=-45, hjust=0))
print(fig3a)
#ggsave("fig3a_black_and_white.png", height=3, width=clm2/3, units="in", dpi=300)
# expression(over(paste("Gene Copies, ",  ml^{-1}), ",/nSum of Bacterial + Archaeal primers"))



fig3a_lm <- lm(log10(qPCR.Bac.per.cc + qPCR.Arc.per.cc) ~ log10(qPCRuniversal), data=fig3a_df)
summary(fig3a_lm)

fig3a_lm2 <- lm(log10(qPCR.Bac.per.cc + qPCR.Arc.per.cc) ~ log10(qPCRuniversal) + core, data=fig3a_df)
summary(fig3a_lm2)

#############
# 3b
############

oneOneLim <- 10^seq(from=3, to=11, by=0.1)
reasonableCopyNumber <- 3.04 # average of bacterial copy numbers in sequenced genomes
oneOne <- data.frame(x=oneOneLim, y=oneOneLim, ymax=24*oneOneLim, yReasonable=reasonableCopyNumber*oneOneLim)
fig3b_data <- allData[!is.na(allData$qPCRtotal) &
                        !is.na(allData$totalcells) &
                        !is.na(allData$core) &
                        allData$environment=="sediments" &
                        allData$qPCRtotal!=0 & allData$totalcells != 0, ]

#For some reason the papers are out of alphabetical order.
fig3b_data$core <- factor(fig3b_data$core, levels=sort(unique(fig3b_data$core)), ordered=TRUE)

fig3bm <- lm(log10(qPCRtotal) ~ log10(totalcells), data=fig3b_data)
rng <- seq(from=min(log10(fig3b_data$totalcells)), to=max(log10(fig3b_data$totalcells)), length.out=100)
rngDF <- data.frame(totalcells=rng)
fig3b_pred <- predict(fig3bm,
                      #newdata=rngDF,
                      interval="prediction")
fig3b_predDF <- as.data.frame(fig3b_pred)
fig3b_predDF$x <- fig3b_data$totalcells
fig3b_predDF$span <- fig3b_predDF$upr - fig3b_predDF$lwr
#ggplot(fig3b_predDF, aes(x=x, y=span)) + geom_line() + scale_y_continuous(limits=c(0, 5))
mean(fig3b_predDF$span, na.omit=TRUE)


fig3b_predDF$lwr <- 10^fig3b_predDF$lwr
fig3b_predDF$fit <- 10^fig3b_predDF$fit
fig3b_predDF$upr <- 10^fig3b_predDF$upr
#fig3b_predDF$lwr[fig3b_predDF$lwr < 1000] <- 1000

fig3b_predDF
fig3b_predDF$span <- fig3b_predDF$upr / fig3b_predDF$lwr
#ggplot(fig3b_predDF, aes(x=x, y=span)) + geom_line()

# Make the various data frame have max values at 1e11
oneOne[oneOne$ymax > 1e11, "ymax"] <- 1e11
oneOne[oneOne$yReasonable > 1e11, "yReasonable"] <- 1e11

fig3b_predDF[fig3b_predDF$lwr < 1e4, "lwr"] <- 1e4
fig3b_predDF[fig3b_predDF$upr > 1e11, "upr"] <- 1e11

# lm of qPCRtotal vs totalcells
summary(lm(log10(qPCRtotal) ~ log10(totalcells), data=fig3b_data))
ggplot(fig3b_data, aes(x=totalcells, y=qPCRtotal)) + geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  geom_smooth(method="lm")

# Big hairy figure showing qPCR total vs totalcells
fig3b <- ggplot() + 
  geom_blank(data=fig3b_data, aes(x=totalcells, y=qPCRtotal, colour=core, shape=core)) +
  geom_line(data=oneOne, aes(x=x, y=y)) +
  geom_ribbon(data=fig3b_predDF, aes(x=x, ymin=lwr, ymax=upr), fill="black", alpha=0.1) +
  geom_line(data=fig3b_predDF, aes(x=x, y=fit), linetype=2) +
  geom_ribbon(data=oneOne, aes(x=x, ymin=y, ymax=ymax), fill="green", alpha=0.2) +
  geom_ribbon(data=oneOne, aes(x=x, ymin=y, ymax=yReasonable), fill="darkgreen", alpha=0.4) +
  geom_point(data=fig3b_data, aes(x=totalcells, y=qPCRtotal, colour=core, shape=core), size=1) +
  geom_smooth(data=fig3b_data, aes(x=totalcells, y=qPCRtotal, colour=core), method="lm", se=FALSE, size=0.2) +
  scale_x_log10(limits = c(1e5, 1e10), breaks=10^(5:10)) +
  scale_y_log10(limits=c(1e4, 1e11), breaks=10^(4:11)) +
  scale_shape_manual(values=rep(15:18, nlevels(fig3b_data$core))) +
  coord_fixed() +
  guides(colour = guide_legend(ncol = 2), shape=guide_legend(ncol=2)) +
  theme(axis.text.x=element_text(angle=-45, hjust=0),
        legend.position="right")
print(fig3b)
#ggsave("fig3b_test.png", height=9, width=9, units="in", dpi=300)

####
# Fig 3c: variance in qPCR and CARDFISH vs variance in total cells
#####

# Way 1 to pull out desired data. Problem is, the qPCR data set includes CARDFISH data with bad methods
varVvar_data_qPCR <- corrected_seds[!is.na(corrected_seds$totalcells) &
                                      !is.na(corrected_seds$qPCRtotal), ]
varVvar_data_CARDFISH <- corrected_seds[!is.na(corrected_seds$totalcells) &
                                          !is.na(corrected_seds$Fish.or.cardFish) &
                                          !is.na(corrected_seds$CARDFISH.Total.per.cc) &
                                          corrected_seds$Fish.or.cardFish == "CARDFISH" &
                                          corrected_seds$Arc.permeabilization == "proteinase K", ]
varVvar_data <- rbind(varVvar_data_qPCR, varVvar_data_CARDFISH)

# Way 2 to pull out desired data, which is any row with EITHER valid qPCR and totalcells data, 
#    OR any row with valid CARDFISH data AND the right methods (Arc.permeabilization=="proteinase K")
varVvar_data <- corrected_seds[ #get rid of bad qPCR data
  (!is.na(corrected_seds$totalcells) & !is.na(corrected_seds$qPCRtotal) |
     # get rid of bad CARDFISH data
     !is.na(corrected_seds$totalcells) &
     !is.na(corrected_seds$CARDFISH.Total.per.cc) &
     !is.na(corrected_seds$Arc.permeabilization) &
     corrected_seds$Arc.permeabilization != "proteinase K"), ]

varVvar_data_CARDFISH <- corrected_seds[!is.na(corrected_seds$totalcells) &
                                          !is.na(corrected_seds$CARDFISH.Total.per.cc) &
                                          !is.na(corrected_seds$Arc.permeabilization) &
                                          corrected_seds$Arc.permeabilization == "proteinase K", ]
vVv_CF <- ddply(varVvar_data_CARDFISH, .(core), summarise, # there are 4 NAs in here b/c 4 cores only have 1 data point
                cellsSD=sd(log10(totalcells), na.rm=TRUE), 
                CARDsd=sd(log10(CARDFISH.Total.per.cc), na.rm=TRUE),
                nSamples = length(CARDFISH.Total.per.cc))

nrow(vVv_CF[!is.na(vVv_CF$cellsSD) & !is.na(vVv_CF$CARDsd), ])
cfm <- lm(CARDsd ~ cellsSD, data=vVv_CF)
summary(cfm)

# Variance in total cells vs variance in qPCR
varVvar_data_qPCR <- corrected_seds[!is.na(corrected_seds$totalcells) &
                                      !is.na(corrected_seds$qPCRtotal), ]
vVv_q <- ddply(varVvar_data_qPCR, .(core), summarise,
               cellsSD=sd(log10(totalcells), na.rm=TRUE), 
               qPCRsd=sd(log10(qPCRtotal), na.rm=TRUE),
               nSamples = length(qPCRtotal))

nrow(vVv_q[!is.na(vVv_q$cellsSD) & !is.na(vVv_q$qPCRsd), ])
qm <- lm(qPCRsd ~ cellsSD, data=vVv_q)
summary(qm)

vVv_CF$qPCRsd <- NA
vVv_q$CARDsd <- NA

vVv_df <- rbind(vVv_CF, vVv_q)
vVv_m <- melt(vVv_df, measure.vars=c("CARDsd", "qPCRsd"), variable.name="enumeration_method", value.name="sd")

fig3c <- ggplot(vVv_m, aes(x=cellsSD, y=sd, colour=enumeration_method)) + geom_point(pointsize=ps) +
  geom_smooth(method="lm", se=FALSE, linetype=2) +
  geom_abline(slope=1, intercept=0) +
  scale_colour_manual(values=c("red", "blue"),
                      labels=c("CARD-FISH\nWith Proteinase K", "qPCR"),
                      name="Enumeration Method") +
  xlab("Std. dev. of log of total cells") +
  ylab("Std. dev. of log of qPCR or CARD-FISH") +
  coord_cartesian(xlim=seq(0, 1.25, by=0.25), ylim=seq(0, 1.25, by=0.25)) +
  theme(legend.position="top")
print(fig3c)

#vp3a <- viewport(x=1/6, y=0.5, height=1, width=1/3)
#vp3b <- viewport(x=3/6, y=0.5, height=1, width=1/3)
vpA <-viewport(x=0.75, y=0.25, height=0.5, width=0.5)
vpB <- viewport(x=0.5, y=0.75, height=0.5, width=1)
vpC <- viewport(x=0.25, y=0.25, height=0.5, width=0.5)

tiff("fig3_tosubmit.tif", height=8, width=clm2, 
     units="in", res=900, compression="lzw")
print(fig3a, vp=vpA)
print(fig3b, vp=vpB)
print(fig3c, vp=vpC)
dev.off()

##########
# Fig 4: Trends in archaea & bacteria (relative & absolute) by depth 
##########

# Size of axis titles
ax.tit.sz <- 10
ps4 <- 0.75
aMargin <- unit(c(0, 0, 0, 0), "in")
bcMargin <- unit(c(0, 0.1, 0, 0.2), "in")

# Set special theme for this specific plot
oldTheme <- theme_get()
theme_set <- oldTheme + theme(text=element_text(size=9))

qPCR_depth_df <- corrected_seds[!is.na(corrected_seds$percentqPCR) &
                                  !is.na(corrected_seds$depth) &
                                  corrected_seds$Environment.Type != "Intertidal", ]
qPCR_depth_df$Uses.516.for.Arc[is.na(qPCR_depth_df$Uses.516.for.Arc)] <- FALSE



qPCR_depth_df$qual_depth <- cut(qPCR_depth_df$depth_log10, breaks=c(-3:3), right=TRUE, 
                                labels=c("0-1cm", "1-10cm", "10cm-1m", "1m-10m", "10m-100m", ">100m")) 
qPCR_depth_df$qual_depth <- factor(qPCR_depth_df$qual_depth, levels=rev(levels(qPCR_depth_df$qual_depth)), ordered=TRUE)

# What is the purpose of this? It is overwritten later
nPoints <- ddply(qPCR_depth_df, c("qual_depth"), summarise, 
                 nPointsText=paste(length(qual_depth[Uses.516.for.Arc==TRUE]), 
                                   "\n", length(qual_depth[Uses.516.for.Arc==FALSE])))



#######
# I seem to use this data frame later, but not (directly) for fig b
######
fig4b_data <- corrected_seds[!is.na(corrected_seds$depth) &
                               !is.na(corrected_seds$Fraction.Arc.CARDFISH) &
                               !is.na(corrected_seds$Arc.permeabilization) &
                               corrected_seds$Arc.permeabilization == "proteinase K", ]




########## 
# Experimental fig 4a: Include CARDFISH
##########
fig4b_data$qual_depth <- cut(fig4b_data$depth_log10, breaks=c(-3:3), right=TRUE, 
                             labels=c("0-1cm", "1-10cm", "10cm-1m", "1m-10m", "10m-100m", ">100m")) 
fig4b_data$qual_depth <- factor(fig4b_data$qual_depth, levels=rev(levels(fig4b_data$qual_depth)), ordered=TRUE)

names(fig4b_data)[names(fig4b_data)=="Fraction.Arc.CARDFISH"] <- "percentArc"
fig4b_data$method <- "CARD-FISH (prot-K)"
names(qPCR_depth_df)[names(qPCR_depth_df)=="percentqPCR"] <- "percentArc"
qPCR_depth_df$method <- NA
qPCR_depth_df$method[qPCR_depth_df$Uses.516.for.Arc == TRUE] <- "qPCR (516)"
qPCR_depth_df$method[qPCR_depth_df$Uses.516.for.Arc == FALSE] <- "qPCR (not 516)"

fig4aEx_df <- rbind(qPCR_depth_df[ , c("depth_log10", "qual_depth", "percentArc", "method")], 
                    fig4b_data[ , c("depth_log10", "qual_depth", "percentArc", "method")])

# #This approach to make all the bars the same width is supposed to work, ubt apparently does not
#  xtabs(formula = ~ as.factor(qual_depth) + as.factor(method), data=fig4aEx_df)
#  newRows <- data.frame(qual_depth=c(">100m", "1m-10m", "10m-100m"), percentArc=rep(NA, 3), method=rep("CARD-FISH (prot-K)"))
#  fig4aEx_df <- rbind(fig4aEx_df, newRows)
# xtabs(formula = ~ as.factor(qual_depth) + as.factor(method), data=fig4aEx_df)

# Set the color palette for this figure
fig1col <- brewer.pal(n=3, name="Paired")
fig1col[3] <- brewer.pal(n=3, name="OrRd")[3]
fig1col <- fig1col[c(3, 1, 2)]

nPoints <- ddply(fig4aEx_df, c("qual_depth"), summarise, 
                 nPointsText=paste(length(qual_depth[method=="qPCR (not 516)"]), 
                                   "\n", length(qual_depth[method=="qPCR (516)"]),
                                   if(length(qual_depth[method=="CARD-FISH (prot-K)"]) != 0) {
                                     paste("\n", length((qual_depth[method=="CARD-FISH (prot-K)"])))
                                   } else {
                                     paste("\n", "")
                                   },
                                   sep=""))


fig4a <- ggplot() + 
  geom_boxplot(data=fig4aEx_df, aes(x=qual_depth, fill=method, y=percentArc), outlier.size=ps4) +
  geom_text(data=nPoints, aes(x=qual_depth, y=-0.05, label=nPointsText), size=1.75, hjust=1) +
  xlab("depth range") +
  scale_y_continuous(expression(paste(over("archaea", "archaea + bacteria"))), lim=c(-0.1, 1)) +
  scale_fill_manual(values=fig1col, name="Method",
                    label=c("CARD-FISH\n(prot-K)", "qPCR\n(516)", "qPCR\n(not 516)")
  ) +
  coord_flip() +
  theme(legend.position="top", 
        text=element_text(size=ax.tit.sz),
        plot.margin=aMargin,
        legend.justification=c(0, 0))
print(fig4a)


#ggsave("fig4a with CARDFISH.png", height=5, width=6)

# Calculate medians & interquartile ranges of each bar in fig 4a

fig4a_quant <- ddply(fig4aEx_df, c("qual_depth", "method"), summarise, 
                     median=round(quantile(percentArc, .5, na.rm=TRUE), 3),
                     twentyFifth=round(quantile(percentArc, .25, na.rm=TRUE), 3),
                     seventyFifth=round(quantile(percentArc, .75, na.rm=TRUE), 3),
                     IQR = seventyFifth - twentyFifth)
#pdf("fig4a_medians_and_IQRs.pdf", height=8.5, width=11)
#grid.table(fig4a_quant[nrow(fig4a_quant):1, ])
#dev.off()

# Make figs b and c, depth profiles of absolute quantities of Archaea and Bacteria
# by direct card counts and estimation from qPCR

# parameters for plotting
protK_yield <- median(corrected_seds$FISH.yield[corrected_seds$Arc.permeabilization=="proteinase K"], na.rm=TRUE)
xl <- c(1e4, 1e11)


fig4bc_df_Arc_directCounts <- corrected_seds[!is.na(corrected_seds$Environment.Type) &
                                               !is.na(corrected_seds$depth_log10) &
                                               !is.na(corrected_seds$totalcells) &
                                               !is.na(corrected_seds$CARDFISH.Arc.per.cc) &
                                               !is.na(corrected_seds$Arc.permeabilization) &
                                               #corrected_seds$CARDFISH.Arc.per.cc != 0 &
                                               corrected_seds$Environment.Type != "Intertidal" &
                                               corrected_seds$Environment.Type != "Salt marsh" &
                                               corrected_seds$Arc.permeabilization=="proteinase K", 
                                             
                                             c("paper", "depth_log10", "CARDFISH.Arc.per.cc")]
fig4bc_df_Arc_directCounts$method <- "measured"

fig4bc_df_Arc_directCounts[fig4bc_df_Arc_directCounts==0] <- 1 #So as not to cause problems with log10

mArcMeas <- lm(log10(CARDFISH.Arc.per.cc) ~ depth_log10, data=fig4bc_df_Arc_directCounts)
summary(mArcMeas)

fig4bc_df_Arc_est <- corrected_seds[!is.na(corrected_seds$Environment.Type) &
                                      !is.na(corrected_seds$depth_log10) &
                                      !is.na(corrected_seds$totalcells) &
                                      !is.na(corrected_seds$percentqPCR) &
                                      corrected_seds$Uses.516.for.Arc==FALSE &
                                      corrected_seds$Environment.Type != "Intertidal" &
                                      corrected_seds$Environment.Type != "Salt marsh" ,
                                    c("paper", "depth_log10", "totalcells", "percentqPCR")]
fig4bc_df_Arc_est$method <- "estimated"
fig4bc_df_Arc_est$estArc <- fig4bc_df_Arc_est$totalcells * fig4bc_df_Arc_est$percentqPCR * protK_yield

mArcEst <- lm(log10(estArc) ~ depth_log10, data=fig4bc_df_Arc_est)
summary(mArcEst)

fig4bc_df_Arc_directCounts <- rename(fig4bc_df_Arc_directCounts, c("CARDFISH.Arc.per.cc" = "Arc.per.cc"))
fig4bc_df_Arc_est <- rename(fig4bc_df_Arc_est, c("estArc" = "Arc.per.cc"))

# Make the data frame from which the plot will be made
fig4b_df <- rbind(fig4bc_df_Arc_directCounts[ , c("paper", "depth_log10", "Arc.per.cc", "method")], 
                  fig4bc_df_Arc_est[ , c("paper", "depth_log10", "Arc.per.cc", "method")])
fig4b_df <- fig4b_df[sample(1:nrow(fig4b_df)), ] #prevents systematic pattern in overplotting

# Absolute number of Archaea (measured and estimated)
fig4b <- ggplot(fig4b_df, aes(x=depth_log10, y=Arc.per.cc, colour=method)) + geom_point(size=ps4) +
  #geom_smooth(method="lm", se=TRUE, linetype=2) +
  geom_smooth(data=subset(fig4b_df, method=="estimated"), 
              aes(x=depth_log10, y=Arc.per.cc, colour=method), 
              method="lm", 
              linetype=2) +
  scale_x_reverse("log of depth (m)") +
  scale_y_log10(expression(paste("Archaea (cells ", ml^{-1}, ")")),
                limits=xl,
                breaks=10^(log10(xl[1]):log10(xl[2]))) +
  scale_colour_manual(values=fig1col[c(3, 1)]) +
  coord_flip() +
  theme(legend.position="top",
        axis.text.x=element_text(angle=-45, hjust=0),
        plot.margin=bcMargin)
print(fig4b)
#ggsave("fig4b.png", height=4, width=clm, units="in", dpi=300)

# Do the same thing, but with Bacteria
fig4bc_df_Bac_directCounts <- corrected_seds[!is.na(corrected_seds$Environment.Type) &
                                               !is.na(corrected_seds$depth_log10) &
                                               !is.na(corrected_seds$totalcells) &
                                               !is.na(corrected_seds$CARDFISH.Bac.per.cc) &
                                               corrected_seds$CARDFISH.Bac.per.cc != 0 &
                                               corrected_seds$Environment.Type != "Intertidal" &
                                               corrected_seds$Environment.Type != "Salt marsh",
                                             #corrected_seds$Arc.permeabilization=="proteinase K"
                                             c("paper", "depth_log10", "CARDFISH.Bac.per.cc")]
fig4bc_df_Bac_directCounts$method <- "measured"

mBacMeas <- lm(log10(CARDFISH.Bac.per.cc) ~ depth_log10, data=fig4bc_df_Bac_directCounts)
summary(mBacMeas)

fig4bc_df_Bac_est <- corrected_seds[!is.na(corrected_seds$Environment.Type) &
                                      !is.na(corrected_seds$depth_log10) &
                                      !is.na(corrected_seds$totalcells) &
                                      !is.na(corrected_seds$percentqPCR) &
                                      corrected_seds$Uses.516.for.Arc==FALSE &
                                      corrected_seds$Environment.Type != "Intertidal" &
                                      corrected_seds$Environment.Type != "Salt marsh" ,
                                    c("paper", "depth_log10", "totalcells", "percentqPCR")]
fig4bc_df_Bac_est$method <- "estimated"
fig4bc_df_Bac_est$estBac <- fig4bc_df_Arc_est$totalcells * (1-fig4bc_df_Arc_est$percentqPCR) * protK_yield

mBacEst <- lm(log10(estBac) ~ depth_log10, data=fig4bc_df_Bac_est)
summary(mBacEst)

fig4bc_df_Bac_directCounts <- rename(fig4bc_df_Bac_directCounts, c("CARDFISH.Bac.per.cc" = "Bac.per.cc"))
fig4bc_df_Bac_est <- rename(fig4bc_df_Bac_est, c("estBac" = "Bac.per.cc"))

# Make the data frame from which the plot will be made
fig4c_df <- rbind(fig4bc_df_Bac_directCounts[ , c("paper", "depth_log10", "Bac.per.cc", "method")], 
                  fig4bc_df_Bac_est[ , c("paper", "depth_log10", "Bac.per.cc", "method")])

fig4c_df$Bac.per.cc[fig4c_df$Bac.per.cc==0] <- 1 #Set the 0 
fig4c_df <- fig4c_df[sample(1:nrow(fig4c_df)), ]
fig4c <- ggplot(fig4c_df, aes(x=depth_log10, y=Bac.per.cc, colour=method)) + geom_point(size=ps4) +
  geom_smooth(method="lm", se=TRUE, linetype=2) +
  scale_x_reverse("log of depth (m)"
  ) +
  scale_y_log10(expression(paste("Bacteria (cells ", ml^{-1}, ")")), 
                limits=xl,
                breaks=10^(log10(xl[1]):log10(xl[2]))) +
  scale_colour_manual(values=fig1col[c(3, 1)]) +
  coord_flip() +
  theme(legend.position="top",
        axis.text.x=element_text(angle=-45, hjust=0),
        plot.margin=bcMargin)
print(fig4c)
#ggsave("fig4c.png", height=4, width=clm, units="in", dpi=300)


vp1 <- viewport(height=1, width=0.4, x=0.2, y=0.5)
vp2 <- viewport(height=1, width=0.3, x=0.55, y=0.5)
vp3 <- viewport(height=1, width=0.3, x=0.85, y=0.5)

tiff("fig4_tosubmit.tif", height=4, width=clm2+0.5, 
     units="in", res=900, compression="lzw")
print(fig4b, vp=vp2)
print(fig4c, vp=vp3)
print(fig4a, vp=vp1) # Have to print this one last because the legend impinges on the other plots' space
dev.off()

theme_set(oldTheme)

# ########
# # Miscellaneous stats about depth trends; I think these didn't make the final manuscript
# ########
# 
# 
# qPCR_npoints <- ddply(qPCR_depth_df, .(core), summarise, length(unique(depth_log10)))
# names(qPCR_npoints)[2] <- "nDepths"
# 
# qPCR_depth_df <- merge(qPCR_depth_df, qPCR_npoints, by="core")
# qPCR_depth_df_forTrends <- qPCR_depth_df[qPCR_depth_df$nDepths > 3, ]
# 
# #ddply(corrected_seds, .(core), transform, nPoints = function(x) nrow(subset(x, !is.na(percentqPCR)))
# 
# 
# # WHat fraction have significant decreasing, significant increasing, or no trend?
# depthTrends <- dlply(qPCR_depth_df_forTrends, .(core, Uses.516.for.Arc), function(x) lm(percentqPCR ~ depth_log10, data=x))
# 
# # Function to get the slope, slope error, and significance fo lms
# cor_result <- function(m) {
#   pVal <- round(summary(m)$coefficients[2, 4], digits=3)
#   slope <- round(summary(m)$coefficients[2, 1], digits=2)
#   slopeErr <- round(summary(m)$coefficients[2, 2], digits=2)
#   data.frame(pVal = pVal, slope=slope, slopeErr=slopeErr)
# }
# 
# 
# percentArcqPCRvDepth <- ldply(depthTrends, cor_result)
# percentArcqPCRvDepth$isSig <- FALSE
# percentArcqPCRvDepth$isSig[percentArcqPCRvDepth$pVal < 0.05] <- TRUE
# percentArcqPCRvDepth$trend <- "no sig. trend"
# percentArcqPCRvDepth$trend[percentArcqPCRvDepth$isSig == TRUE & percentArcqPCRvDepth$slope >0] <- "increasing"
# percentArcqPCRvDepth$trend[percentArcqPCRvDepth$isSig == TRUE & percentArcqPCRvDepth$slope <0] <- "decreasing"
# percentArcqPCRvDepth$trend <- as.factor(percentArcqPCRvDepth$trend)
# percentArcqPCRvDepth
# 
# # pdf("table_depth_trends_of_percent_archaea_by_qPCR.pdf", height=12, width=8,
# #     pointsize=10)
# # grid.table(percentArcqPCRvDepth)
# # dev.off()
# 
# # What fraction are increasing (whether or not significantly)
# nrow(percentArcqPCRvDepth[percentArcqPCRvDepth$trend == "increasing", ]) / nrow(percentArcqPCRvDepth)
# 
# # What fraction are decreasing
# nrow(percentArcqPCRvDepth[percentArcqPCRvDepth$trend == "decreasing" &
#                             !is.na(percentArcqPCRvDepth$trend), ]) / nrow(percentArcqPCRvDepth)
# nrow(percentArcqPCRvDepth[is.na(percentArcqPCRvDepth$trend), ]) / nrow(percentArcqPCRvDepth)

##############
# Miscellaneous analysis that goes into the text 
#    Note: this should be run after all the figure scripts, as it relies on some of the 
#          data frames generated by those scripts
##############

# Abstract
# Number of papers
nlevels(allData$paper)

#########################
# Results and discussion
########################

# Want the total number of studies with any seawater data
summary(corrected_sw$paper)
nlevels(corrected_sw$paper)
nrow(corrected_sw)

# minimum yield
corrected_sw_FISH.yield <- corrected_sw[!is.na(corrected_sw$FISH.yield), ]

# Total # studies with FISH yield data
nrow(corrected_sw_FISH.yield)
min(corrected_sw_FISH.yield$FISH.yield, na.rm=TRUE) #0.02
nrow(corrected_sw_FISH.yield[corrected_sw_FISH.yield$FISH.yield >= 0.1, ]) / nrow(corrected_sw_FISH.yield) #99%
quantile(corrected_sw_FISH.yield$FISH.yield, c(.25, .5, .75))

# Total number of studies and samples that have FISH yield
swC <- corrected_sw[!is.na(corrected_sw$FISH.yield), ]
nrow(swC)
swC$paper <- as.factor(as.character(swC$paper))
summary(swC$paper)
nlevels(swC$paper)
quantile(swC$FISH.yield, c(.25, .50, .75))

# Total number of studies and samples that HAVE FISH yield
sedsCount <- corrected_seds[!is.na(corrected_seds$FISH.yield), ]
nrow(sedsCount)
sedsCount$paper <- as.factor(as.character(sedsCount$paper))
summary(sedsCount$paper)
nlevels(sedsCount$paper)
median(sedsCount$FISH.yield)
quantile(sedsCount$FISH.yield, c(.25, .50, .75))
sd(sedsCount$FISH.yield)

# Total number of studies that have qPCR data
sedsQ <- corrected_seds[!is.na(corrected_seds$percentqPCR), ]
nrow(sedsQ)
sedsQ$paper <- as.factor(as.character(sedsQ$paper))
summary(sedsQ$paper)
nlevels(sedsQ$paper)

# Line 249: MEdian yields of 15 lowest-yielding cores
medYds <- ddply(fig2a_data, .(core), summarise, medYd = median(FISH.yield))
medYds <- medYds[order(medYds$medYd, decreasing=FALSE), ]
medYds

# Line 284: low yields of seawater studies with archaea
swLysX.CARDFISH <- corrected_sw[!is.na(corrected_sw$Arc.permeabilization) &
                                  !is.na(corrected_sw$Fraction.Arc.CARDFISH) &
                                  corrected_sw$Arc.permeabilization=="lysozyme", ]
ddply(swLysX.CARDFISH, .(paper), summarise, medX.CARDFISH = median(Fraction.Arc.CARDFISH))
# There's only 1 paper, the yield is 4.2%

# Line 288: Yields from intertidal
#    (This is copied from file fig2aS_intertidal.R)
fig2aS_data <- corrected_seds[!is.na(corrected_seds$Arc.permeabilization) &
                                !is.na(corrected_seds$FISH.yield) &
                                !is.na(corrected_seds$Fish.or.cardFish) &
                                #!corrected_seds$Arc.permeabilization == "CARD-FISH, archaea not measured" &
                                #!corrected_seds$Arc.permeabilization == "FISH, archaea not measured" &
                                (corrected_seds$Fish.or.cardFish == "CARDFISH" | corrected_seds$Fish.or.cardFish=="FISH") &
                                corrected_seds$Environment.Type == "Intertidal", ]
ddply(fig2aS_data, .(core), summarise, medYield = median(FISH.yield))
ddply(fig2aS_data, .(core), summarise, medX.CARDFISH= median(Fraction.Arc.CARDFISH, na.rm=TRUE))

# Number of papers & samples with qPCRtotal values and totalcells
nrow(fig3b_data)
length(unique(fig3b_data$paper))
summary(lm(log10(qPCRtotal) ~ log10(totalcells), data=fig3b_data))

# Number of samples with sdVar > totalcellsVar
sum((vVv_q$qPCRsd > vVv_q$cellsSD), na.rm=TRUE)
sum((vVv_q$qPCRsd < vVv_q$cellsSD), na.rm=TRUE)
summary(lm(qPCRsd ~ cellsSD, data=vVv_q))

# Slope, intercept and p value of relationship between sd of total cells and sd of CARDFISH
summary(lm(CARDsd ~ cellsSD, data=vVv_CF))


ddply(fig4aEx_df, c("qual_depth", "method"), summarise, med=round(median(percentArc), 2), 
      err = round((quantile(percentArc, .75) - quantile(percentArc, .25))/2, 2))

t.test(subset(fig4aEx_df, qual_depth==">100m" & method=="qPCR (not 516)")$percentArc, 
       subset(fig4aEx_df, qual_depth=="10m-100m" & method=="qPCR (not 516)")$percentArc)

# SW comparison t-tests
sz <- 0.75
bts <- TRUE # run the bootstrap analysis?
ylsw <- c(0, 1.3)
xlsw <- c(0.25, 2.75)
al <- 0.1 #alpha (for rug plots)
lgps <- "none" # legend position

# Function to mimic bootstrapped analog of t test
tFun <- function(x, i=NA, splitVar, dataCol) {
  
  # Can be used in resampling functions, or not
  if(!is.numeric(i)){
    xResamp <- x
  } else {
    xResamp <- x[i, ]
  }
  
  xResamp <- xResamp[!is.na(xResamp[ , splitVar]), ]
  
  # Check that the splitting variable has exactly 2 levels 
  if(nlevels(as.factor(x[ , splitVar])) != 2) {
    warning("the split variable splitVar must have exactly 2 levels")
    NA
  } else {
    #print("its cool, you've got 2 or fewer levels")
    #browser()
    medians <- ddply(xResamp, splitVar, function(x) median(x[ , dataCol], na.rm=TRUE))
    medDiff <- medians[2, "V1"] - medians[1, "V1"]
    medDiff
  }
  
}

##################
# Fish vs cardfish BOXPLOTS
##################
swYield_1 <- corrected_sw[!is.na(corrected_sw$Fish.or.cardFish) &
                            !is.na(corrected_sw$FISH.yield), ]
swYield_1$Fish.or.cardFish <- factor(swYield_1$Fish.or.cardFish, levels=c("FISH", "CARDFISH"), ordered=TRUE)

p_fish_card_box_sw <- ggplot(swYield_1, aes(y=FISH.yield, x=Fish.or.cardFish)) +
  geom_boxplot(outlier.size=sz) +
  geom_rug(data=subset(swYield_1, as.numeric(Fish.or.cardFish)==1), aes(x=FISH.yield), sides="l", alpha=0.2) +
  geom_rug(data=subset(swYield_1, as.numeric(Fish.or.cardFish)==2), aes(x=FISH.yield), sides="r", alpha=0.2) +
  #scale_colour_manual(values=c("black", "black"), name="method") + 
  scale_y_continuous("yield") +
  scale_x_discrete("hybridization method") +
  coord_cartesian(xlim=xlsw, ylim=ylsw) +
  theme(legend.position=lgps)
print(p_fish_card_box_sw)


##################
# Cell stain BOXPLOTS
##################
swYield_2 <- corrected_sw[!is.na(corrected_sw$Cell.stain) &
                            !is.na(corrected_sw$FISH.yield), ]
swYield_2$stain_reduced <- NA
swYield_2$stain_reduced[swYield_2$Cell.stain=="DAPI"] <- "DAPI"
swYield_2$stain_reduced[swYield_2$Cell.stain=="AO" | swYield_2$Cell.stain=="AO/isopropanol"] <- "AO"
swYield_2$stain_reduced <- factor(swYield_2$stain_reduced)

p_stain_box_sw <- ggplot(swYield_2, aes(x=stain_reduced, y=FISH.yield)) + 
  geom_boxplot(outlier.size=sz) +
  geom_rug(data=subset(swYield_2, as.numeric(stain_reduced)==1), aes(x=FISH.yield), sides="l", alpha=0.2) +
  geom_rug(data=subset(swYield_2, as.numeric(stain_reduced)==2), aes(x=FISH.yield), sides="r", alpha=0.2) +
  scale_y_continuous("yield") +
  scale_x_discrete("cell stain") +
  coord_cartesian(xlim=xlsw, ylim=ylsw) +
  theme(legend.position=lgps)
print(p_stain_box_sw)

##################
# Counting method BOXPLOTS
##################
swYield_3 <- corrected_sw[!is.na(corrected_sw$Counting.method) & !is.na(corrected_sw$FISH.yield), ]


p_automated_box_sw <- ggplot(swYield_3, aes(x=Counting.method, y=FISH.yield)) + #geom_density() +
  geom_boxplot(outlier.size=sz) +
  geom_rug(data=subset(swYield_3, as.numeric(Counting.method)==1), aes(x=FISH.yield), sides="l", alpha=0.2) +
  geom_rug(data=subset(swYield_3, as.numeric(Counting.method)==2), aes(x=FISH.yield), sides="r", alpha=0.2) +
  scale_y_continuous("yield") +
  scale_x_discrete("counting method") +
  coord_cartesian(xlim=xlsw, ylim=ylsw) +
  theme(legend.position=lgps)
print(p_automated_box_sw)

##################
# Counting method BOXPLOTS
##################
swYield_4 <- corrected_sw[!is.na(corrected_sw$Arc.probe) & 
                            !is.na(corrected_sw$FISH.yield) &
                            #(corrected_sw$Arc.probe == "ARCH915" |
                            #   corrected_sw$Arc.probe == "CREN554, EURY806"), ]
                            corrected_sw$Arc.probe != "not tried", ]


# Test whether it is reasonable to collapse the 3 'CREN' probe methods into one
#     First make the probes into a new factor & reorder the levels
probeMed <- ddply(swYield_4, .(Arc.probe), summarise, medYield = median(FISH.yield))
swYield_4$corr_probe <- factor(swYield_4$Arc.probe, 
                               levels=probeMed$Arc.probe[order(probeMed$medYield, decreasing=TRUE)], 
                               ordered=TRUE)

# Make a boxplot of all the methods, just to see if there are big differences between the various CREN probes
ggplot(swYield_4, aes(x=corr_probe, y=FISH.yield)) + geom_boxplot(outlier.size=0) + 
  geom_point(position=position_jitter(width=0.2), alpha=0.2) +
  ggtitle("no sig difference between the cren probes, which are both different from ARCH915")
ggsave("fig comparison of all sw Arc probes_tosubmit.png", height=5, width=6, units="in")


# Make the linear model
probeAOV <- aov(FISH.yield ~ corr_probe, data=swYield_4)
summary(probeAOV)

# Perform post-hoc analysis
TukeyHSD(probeAOV)
# So there is no difference among the CREN probes, but there is a difference between the CREN and ARCH probes.

# Therefore we'll lump the two CREN probes together
swYield_4$corr_probe_condense <- NA
swYield_4$corr_probe_condense[swYield_4$corr_probe=="ARCH915"] <- "ARCH915"
swYield_4$corr_probe_condense[is.na(swYield_4$corr_probe_condense)] <- "var. CREN probes"


p_probe_box_sw <- ggplot(swYield_4, aes(x=corr_probe_condense, y=FISH.yield)) + 
  geom_boxplot(outlier.size=sz) +
  geom_rug(data=subset(swYield_4, as.numeric(corr_probe)==1), aes(x=FISH.yield), sides="l", alpha=0.2) +
  geom_rug(data=subset(swYield_4, as.numeric(corr_probe)==2), aes(x=FISH.yield), sides="r", alpha=0.2) +
  scale_y_continuous("yield") +
  scale_x_discrete("probe") +
  coord_cartesian(xlim=xlsw, ylim=ylsw) +
  theme(legend.position=lgps)
print(p_probe_box_sw)

# Make seawater plots
vp1 <- viewport(width = 0.5, height = 0.5, x = 0.25, y = 0.75)
vp2 <- viewport(width= 0.5, height = 0.5, x=0.75, y= 0.75)
vp3 <- viewport(width= 0.5, height = 0.5, x=0.25, y= 0.25)
vp4 <- viewport(width= 0.5, height = 0.5, x=0.75, y= 0.25)


tiff("figS1_sw_boxplots_tosubmit.tif", width=clm2, height=clm2/1.4, units="in",res=600, compression="lzw")
print(p_fish_card_box_sw, vp=vp1)
print(p_stain_box_sw, vp=vp2)
print(p_automated_box_sw, vp=vp3)
print(p_probe_box_sw, vp=vp4) # The probe analysis is currently up in the air
dev.off()



######## 
# Bootsrapping effect sizes
########

if(bts == TRUE) {
  # Create a data frame to hold the bootstrap results
  
  
  
  # bootstrap test of whether median yields are different for FISH or CARDFISH
  Fish.or.cardFish.b <- boot(swYield_1, tFun, R=1000, splitVar="Fish.or.cardFish", dataCol="FISH.yield")
  print(Fish.or.cardFish.b)
  print(median(Fish.or.cardFish.b$t))
  plot(Fish.or.cardFish.b)
  Fish.or.cardFish.ci <- boot.ci(Fish.or.cardFish.b)
  df1 <- data.frame(comp = "FISH or CARDFISH", med = median(Fish.or.cardFish.b$t), lowConf = Fish.or.cardFish.ci$bca[4], hiConf = Fish.or.cardFish.ci$bca[5])
  
  # bootstrap whether median yields are different for stain (reduced) 
  stain.b <- boot(swYield_2, tFun, R=1000, splitVar="stain_reduced", dataCol="FISH.yield")
  print(stain.b)
  print(median(stain.b$t))
  plot(stain.b)
  stain.ci <- boot.ci(stain.b)
  df2 <- data.frame(comp = "stain", med = median(stain.b$t), lowConf = stain.ci$bca[4], hiConf = stain.ci$bca[5])
  
  # bootstrap whether median yields are different by counting method (automatic or manual)
  count_method.b <- boot(swYield_3, tFun, R=1000, splitVar="Counting.method", dataCol="FISH.yield")
  count_method.b
  print(median(count_method.b$t))
  plot(count_method.b)
  count_method.ci <- boot.ci(count_method.b)
  df3 <- data.frame(comp = "count method", med = median(count_method.b$t), lowConf = count_method.ci$bca[4], hiConf = count_method.ci$bca[5])
  
  # Bootstrap whether median yields are different by archaeal probe
  # Note that all the CREN probes have been lumped together
  arc.probe.b <- boot(swYield_4, tFun, R=1000, splitVar="corr_probe_condense", dataCol="FISH.yield")
  arc.probe.b
  print(median(arc.probe.b$t))
  plot(arc.probe.b)
  arc.probe.ci <- boot.ci(arc.probe.b)
  df4 <- data.frame(comp = "Arc probe", med = median(arc.probe.b$t), lowConf = arc.probe.ci$bca[4], hiConf = arc.probe.ci$bca[5])
  
  boot_results_df <- rbind(df1, df2, df3, df4)
  pdf("sw_effect_sizes.pdf")
  grid.table(boot_results_df)
  dev.off()
}

### Fig 2a, but for intertidal
###   (Will actually go in supplemental somewhere)
### Fig 2
fig2aS_data <- corrected_seds[!is.na(corrected_seds$Arc.permeabilization) &
                                !is.na(corrected_seds$FISH.yield) &
                                !is.na(corrected_seds$Fish.or.cardFish) &
                                #!corrected_seds$Arc.permeabilization == "CARD-FISH, archaea not measured" &
                                #!corrected_seds$Arc.permeabilization == "FISH, archaea not measured" &
                                (corrected_seds$Fish.or.cardFish == "CARDFISH" | corrected_seds$Fish.or.cardFish=="FISH") &
                                corrected_seds$Environment.Type == "Intertidal", ]


# Make a permeabilization column, to reflect that FISH doesn't need a permeabilization method
fig2aS_data$perm <- as.character(fig2aS_data$Arc.permeabilization)
fig2aS_data$perm[fig2aS_data$perm=="none"] <- "FISH"

# Re-level the permeabilization methods
fig2aS_data$perm <- factor(fig2aS_data$perm, 
                           levels=c("proteinase K", 
                                    "FISH", 
                                    "FISH, archaea not measured",
                                    "CARD-FISH, archaea not measured",
                                    #"detergent", 
                                    #"unknown",  
                                    #"lysozyme",
                                    "lysozyme/achromopeptidase"), 
                           #labels=c("proteinase K", "none (FISH)", "detergent", 
                           #         "unknown", "lysozyme/achromopeptidase", "lysozyme"),
                           ordered=TRUE)

### FIg 2a reduced: cut out all the cores with only one data point
# Count the # points in the core, remove all < 1
fig2aS_data2 <- ddply(fig2aS_data, .(core), transform, count=length(FISH.yield))
fig2aS_data2 <- fig2aS_data2[fig2aS_data2$count > 1, ]
fig2aS_data2$core <- as.factor(as.character(fig2aS_data2$core))

# Order the papers in order of decreasing median yield
fig2aS_dataMedianYields <- ddply(fig2aS_data2, .(core), summarise, medianYield=median(FISH.yield, na.rm=TRUE))
fig2aS_dataMedianYields <- fig2aS_dataMedianYields[order(fig2aS_dataMedianYields$medianYield, decreasing=TRUE), ]
fig2aS_data2$core <- factor(fig2aS_data2$core, levels=fig2aS_dataMedianYields$core, ordered=TRUE)


# Determine the quartiles of yield from SW data; overplot that on sediments plot
#swQuant <- quantile(allData$FISH.yield[allData$environment=="seawater"], probs=c(0.5-0.341, 0.5, 0.5+0.341), na.rm=TRUE)
swQuant <- quantile(corrected_sw$FISH.yield, probs=c(0.25, 0.5, 0.75), na.rm=TRUE)
swYieldDF <- data.frame(ymin=swQuant[1], ymed=swQuant[2], ymax=swQuant[3])
rownames(swYieldDF) <- NULL
swYieldDF <- data.frame(x=c(0, length(unique(fig2aS_data2$core))+0.5), rbind(swYieldDF, swYieldDF))

## Count number of points in each core
nPoints <- ddply(fig2aS_data2, .(core), function(x) data.frame(nrow=nrow(x), perm=x$perm[1]))
nPoints$core <- factor(nPoints$core, levels=levels(fig2aS_data2$core), ordered=TRUE)
nPoints$perm <- factor(nPoints$perm, levels=levels(fig2aS_data2$perm), ordered=TRUE)

envDF <- ddply(fig2aS_data2, .(core), function(x) as.character(x$Environment.Type[1]))

# Set the text position
colors <- brewer.pal(n=6, name="RdBu")
colors[4] <- "#999999"
textPos <- 1.4

fig2aS_2 <- ggplot() +
  geom_blank(data=fig2aS_data2, aes(x=core, y=FISH.yield, fill=perm)) +
  geom_ribbon(data=swYieldDF, aes(x=x, ymin=ymin, ymax=ymax), fill="skyblue", alpha=1) +
  geom_line(data=swYieldDF, aes(x=x, y=ymed), colour="black") +
  geom_boxplot(data=fig2aS_data2, aes(x=core, y=FISH.yield, fill=perm), colour="black", outlier.size=1) +
  geom_rect(data=fig2aS_data2, aes(xmin=as.numeric(core)-0.5, xmax=as.numeric(core)+0.5, ymin=textPos-0.05, ymax=textPos+0.05, fill=perm)) +
  geom_text(data=nPoints, aes(x=core, y=textPos, label=nrow), size=2) +
  geom_text(data=envDF, aes(x=core, y=textPos-0.075, label=V1), size=2, angle=-90, hjust=0) +
  ylab(expression(frac("Bacteria + Archaea (FISH or CARD-FISH)", "Total Cell Count (DAPI, AO or SYBR Green)"))) +
  coord_cartesian(ylim=c(0,1.5)) + 
  scale_fill_manual(values=colors, name="Archaeal\nPermeabilization\nMethod") +
  theme(axis.text.x=element_text(angle=-45, hjust=0, size=6),
        axis.title.x=element_blank(),
        legend.position="top",
        legend.title.align=1,
        plot.margin=unit(c(0, 0.9, 0, 0), "in")) #default margin is c(1, 1, 0.5, 0.5) 
print(fig2aS_2)
ggsave("fig2aS_intertidal_only_yield_tosubmit.png", height=5, width=clm2, dpi=300)


#############
# Fig S3: Absolute quantities of bacteria & Archaea, in sw and seds, quantified by best practices
############

# Point size for all figures
ps <- 0.75

# Set up seawater data frames:
# Bacteria in seawater, by CARDFISH
s3bac_df <- corrected_sw[!is.na(corrected_sw$depth_log10) &
                           !is.na(corrected_sw$CARDFISH.Bac.per.cc), 
                         c("paper", "depth_log10", "CARDFISH.Bac.per.cc")]
s3bac_df$CARDFISH.Bac.per.cc[s3bac_df$CARDFISH.Bac.per.cc == 0] <- 1 # set zeros to one, to avoid problems with log fits
s3bac_df$logBac <- log10(s3bac_df$CARDFISH.Bac.per.cc)

# Archaea in seawater, by CARDFISH
s3arc_df <- corrected_sw[!is.na(corrected_sw$depth_log10) &
                           !is.na(corrected_sw$CARDFISH.Arc.per.cc), 
                         c("paper", "depth_log10", "CARDFISH.Arc.per.cc")]
s3arc_df$CARDFISH.Bac.per.cc[s3arc_df$CARDFISH.Arc.per.cc == 0] <- 1 # set zeros to one, to avoid problems with log fits
s3arc_df$logArc <- log10(s3arc_df$CARDFISH.Arc.per.cc)

# % Archaea (not exactly the same data points as above, since some studies report % Archaea but not actual numbers)
s3percent_df <- corrected_sw[!is.na(corrected_sw$depth_log10) &
                               !is.na(corrected_sw$Fraction.Arc.CARDFISH), 
                             c("paper", "depth_log10", "Fraction.Arc.CARDFISH")]

#source("breakpoint analysis.R")
# Breakpoint analysis: Bacteria
bac_depthM <- lm(logBac ~ depth_log10, data=s3bac_df)
seg_bac_depthM <- segmented(bac_depthM, seg.Z = ~ depth_log10, psi=2)
summary(seg_bac_depthM)
bacBreak <- summary(seg_bac_depthM)$psi[2] #1.81

# Fakey way to give linear models (really I should get them from summary(seg_bac_depthM), but I don't totally understand the output)
shallowBac <- lm(logBac ~ depth_log10, data=subset(s3bac_df, depth_log10 <= bacBreak))
summary(shallowBac)
deepBac <- lm(logBac ~ depth_log10, data=subset(s3bac_df, depth_log10 > bacBreak))
summary(deepBac)


# Breakpoint analysis: Archaea
arc_depthM <- lm(logArc ~ depth_log10, data=s3arc_df)
seg_arc_depthM <- segmented(arc_depthM, seg.Z = ~ depth_log10, psi=2)
summary(seg_arc_depthM)
arcBreak <- summary(seg_arc_depthM)$psi[2] #2.57

# Fakey way to give linear models (really I should get them from summary(seg_bac_depthM), but I don't totally understand the output)
shallowArc <- lm(logArc ~ depth_log10, data=subset(s3arc_df, depth_log10 <= arcBreak))
summary(shallowArc)
deepArc <- lm(logArc ~ depth_log10, data=subset(s3arc_df, depth_log10 > arcBreak))
summary(deepArc)

# Breakpoint analysis: Percent Archaea
percent_depthM <- lm(Fraction.Arc.CARDFISH ~ depth_log10, s3percent_df)
seg_percent_depthM <- segmented(percent_depthM, seg.Z= ~ depth_log10, psi=2)
summary(seg_percent_depthM)
percentBreak <- summary(seg_percent_depthM)$psi[2] #1.71

# Fakey way to give linear models (really I should get them from summary(seg_bac_depthM), but I don't totally understand the output)
shallowPer <- lm(Fraction.Arc.CARDFISH ~ depth_log10, data=subset(s3percent_df, depth_log10 <= percentBreak))
summary(shallowPer)
deepPer <- lm(Fraction.Arc.CARDFISH ~ depth_log10, data=subset(s3percent_df, depth_log10 > percentBreak))
summary(deepPer)

# Seawater figure limits and breaks
sw_cell_limits <- 10^c(2.8, 6.3)
sw_cell_breaks <- 10^(3:6)
sw_depth_limits <- c(4, 0)

sw_depth_breaks <- sw_depth_limits[1]:sw_depth_limits[2]
fig5sw_bac <- ggplot() + 
  geom_point(data=fig5sw_bac_df, aes(x=depth_log10, y=CARDFISH.Bac.per.cc), size=ps) +
  geom_smooth(data=subset(fig5sw_bac_df, depth_log10 <=bacBreak), 
              aes(x=depth_log10, y=CARDFISH.Bac.per.cc), method="lm", colour="black", linetype=2) +
  geom_smooth(data=subset(fig5sw_bac_df, depth_log10 > bacBreak), 
              aes(x=depth_log10, y=CARDFISH.Bac.per.cc), method="lm",  colour="black", linetype=2) +
  scale_y_log10(expression(paste("Bacteria direct counts, cells ", ml^{-1})), limits=sw_cell_limits, breaks=sw_cell_breaks) +
  scale_x_reverse("log depth, m", limits=sw_depth_limits) +
  coord_flip()
print(fig5sw_bac)

fig5sw_arc <- ggplot(fig5sw_arc_df, aes(x=depth_log10, y=CARDFISH.Arc.per.cc)) + geom_point(size=ps) +
  geom_smooth(data=subset(fig5sw_arc_df, depth_log10 <=arcBreak), 
              aes(x=depth_log10, y=CARDFISH.Arc.per.cc), method="lm", colour="black", linetype=2) +
  geom_smooth(data=subset(fig5sw_arc_df, depth_log10 > arcBreak), 
              aes(x=depth_log10, y=CARDFISH.Arc.per.cc), method="lm",  colour="black", linetype=2) +
  scale_y_log10(expression(paste("Archaea direct counts, cells ", ml^{-1})), limits=sw_cell_limits, breaks=sw_cell_breaks) +
  scale_x_reverse("log depth, m", limits=sw_depth_limits) +
  coord_flip()
print(fig5sw_arc)

fig5sw_percent <- ggplot(fig5sw_percent_df, aes(x=depth_log10, y=Fraction.Arc.CARDFISH)) + geom_point(size=ps) +
  geom_smooth(data=subset(fig5sw_percent_df, depth_log10 <=percentBreak), 
              aes(x=depth_log10, y=Fraction.Arc.CARDFISH), method="lm", colour="black", linetype=2) +
  geom_smooth(data=subset(fig5sw_percent_df, depth_log10 > percentBreak), 
              aes(x=depth_log10, y=Fraction.Arc.CARDFISH), method="lm",  colour="black", linetype=2) +
  scale_y_continuous("fraction Archaea", limits=c(-0.02, 1)) +
  scale_x_reverse("log depth, m", limits=sw_depth_limits) +
  coord_flip()
print(fig5sw_percent)

# Make integrated figure
vp1 <- viewport(width=1/3, height=1, x=1/6, y=0.5)
vp2 <- viewport(width=1/3, height=1, x=3/6, y=0.5)
vp3 <- viewport(width=1/3, height=1, x=5/6, y=0.5)

tiff("figS3sw_absolute_quantifications_v_depth_tosubmit.tif", 
     height=3, width=clm2, units="in", res=300, compression="lzw")
print(fig5sw_bac, vp=vp1)
print(fig5sw_arc, vp=vp2)
print(fig5sw_percent, vp=vp3)
dev.off()
