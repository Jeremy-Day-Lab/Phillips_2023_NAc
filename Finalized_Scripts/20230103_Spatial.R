#NDF MCN paper 2022
#20221229
library(scales)
library(dplyr)
#Read in the ROIs #LIST on stardist
#Remember to set emasurements with centroid on for this
G1F_ROIs       <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G1F_NAcfull_Dapi_Results.csv",header = TRUE)
G1F_ROIs$Index <- 1:nrow(G1F_ROIs)
G2F_ROIs       <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G2F_NAcfull_Dapi_Results.csv",header = TRUE)
G2F_ROIs$Index <- 1:nrow(G2F_ROIs)
G3F_ROIs       <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G3F_NAcfull_Dapi_Results.csv",header = TRUE)
G3F_ROIs$Index <- 1:nrow(G3F_ROIs)
G4M_ROIs       <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G4M_NAcfull_Dapi_Results.csv",header = TRUE)
G4M_ROIs$Index <- 1:nrow(G4M_ROIs)
G5M_ROIs       <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G5M_NAcfull_Dapi_Results.csv",header = TRUE)
G5M_ROIs$Index <- 1:nrow(G5M_ROIs)
G6M_ROIs       <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G6M2_NAcfull_Dapi_Results.csv",header = TRUE)
G6M_ROIs$Index <- 1:nrow(G6M_ROIs)


#G1F
#Drd1
Drd1_G1F <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G1F_NAcfull_Drd1_Results.csv",header = TRUE) #MEASURE on stardist
Drd1_G1F$Category <- ifelse(Drd1_G1F$Max >= 83, "Pos", "Neg")
Drd1_G1F$Midline_x  <- 14.50
Drd1_G1F$Midline_y  <- 4.56
Drd1_G1F$Position_x <- G1F_ROIs$X
Drd1_G1F$Position_y <- G1F_ROIs$Y
Drd1_G1F$Difference_x <- (Drd1_G1F$Midline_x - Drd1_G1F$Position_x)*403.39
Drd1_G1F$Difference_y <- (Drd1_G1F$Midline_y - Drd1_G1F$Position_y)*403.39
#Ebf1
Ebf1_G1F <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G1F_NAcfull_Ebf1_Results.csv",header = TRUE)
Ebf1_G1F$Category <- ifelse(Ebf1_G1F$Max >= 80, "Pos", "Neg")
Ebf1_G1F$Midline_x  <- 14.50
Ebf1_G1F$Midline_y  <- 4.56
Ebf1_G1F$Position_x <- G1F_ROIs$X
Ebf1_G1F$Position_y <- G1F_ROIs$Y
Ebf1_G1F$Difference_x <- (Ebf1_G1F$Midline_x - Ebf1_G1F$Position_x)*403.39
Ebf1_G1F$Difference_y <- (Ebf1_G1F$Midline_y - Ebf1_G1F$Position_y)*403.39
#Htr4
Htr4_G1F <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G1F_NAcfull_Htr4_Results.csv",header = TRUE)
Htr4_G1F$Category <- ifelse(Htr4_G1F$Max >= 30, "Pos", "Neg")
Htr4_G1F$Midline_x  <- 14.50
Htr4_G1F$Midline_y  <- 4.56
Htr4_G1F$Position_x <- G1F_ROIs$X
Htr4_G1F$Position_y <- G1F_ROIs$Y
Htr4_G1F$Difference_x <- (Htr4_G1F$Midline_x - Htr4_G1F$Position_x)*403.39
Htr4_G1F$Difference_y <- (Htr4_G1F$Midline_y - Htr4_G1F$Position_y)*403.39


#G2F
#Drd1
Drd1_G2F <- read.csv(file ="/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G2F_NAcfull_Drd1_Results.csv",header = TRUE)
Drd1_G2F$Category <- ifelse(Drd1_G2F$Max >= 100, "Pos", "Neg")
Drd1_G2F$Midline_x  <- 12.12
Drd1_G2F$Midline_y  <- 2.40
Drd1_G2F$Position_x <- G2F_ROIs$X
Drd1_G2F$Position_y <- G2F_ROIs$Y
Drd1_G2F$Difference_x <- (Drd1_G2F$Midline_x - Drd1_G2F$Position_x)*403.39
Drd1_G2F$Difference_y <- (Drd1_G2F$Midline_y - Drd1_G2F$Position_y)*403.39
#Ebf1
Ebf1_G2F <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G2F_NAcfull_Ebf1_Results.csv",header = TRUE)
Ebf1_G2F$Category <- ifelse(Ebf1_G2F$Max >= 80, "Pos", "Neg")
Ebf1_G2F$Midline_x  <- 12.12
Ebf1_G2F$Midline_y  <- 2.40
Ebf1_G2F$Position_x <- G2F_ROIs$X
Ebf1_G2F$Position_y <- G2F_ROIs$Y
Ebf1_G2F$Difference_x <- (Ebf1_G2F$Midline_x - Ebf1_G2F$Position_x)*403.39
Ebf1_G2F$Difference_y <- (Ebf1_G2F$Midline_y - Ebf1_G2F$Position_y)*403.39
#Htr4
Htr4_G2F <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G2F_NAcfull_Htr4_Results.csv",header = TRUE)
Htr4_G2F$Category <- ifelse(Htr4_G2F$Max >= 30, "Pos", "Neg")
Htr4_G2F$Midline_x  <- 12.12
Htr4_G2F$Midline_y  <- 2.40
Htr4_G2F$Position_x <- G2F_ROIs$X
Htr4_G2F$Position_y <- G2F_ROIs$Y
Htr4_G2F$Difference_x <- (Htr4_G2F$Midline_x - Htr4_G2F$Position_x)*403.39
Htr4_G2F$Difference_y <- (Htr4_G2F$Midline_y - Htr4_G2F$Position_y)*403.39


#G3F
#Drd1
Drd1_G3F <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G3F_NAcfull_Drd1_Results.csv", header = TRUE)
Drd1_G3F$Category <- ifelse(Drd1_G3F$Max >= 100,"Pos","Neg")
Drd1_G3F$Midline_x  <- 1122
Drd1_G3F$Midline_y  <- 368
Drd1_G3F$Position_x <- G3F_ROIs$X
Drd1_G3F$Position_y <- G3F_ROIs$Y
Drd1_G3F$Difference_x <- (Drd1_G3F$Midline_x - Drd1_G3F$Position_x)*4.202
Drd1_G3F$Difference_y <- (Drd1_G3F$Midline_y - Drd1_G3F$Position_y)*4.202
#Ebf1
Ebf1_G3F <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G3F_NAcfull_Ebf1_Results.csv",header = TRUE)
Ebf1_G3F$Category <- ifelse(Ebf1_G3F$Max >= 80,"Pos","Neg")
Ebf1_G3F$Midline_x  <- 1122
Ebf1_G3F$Midline_y  <- 368
Ebf1_G3F$Position_x <- G3F_ROIs$X
Ebf1_G3F$Position_y <- G3F_ROIs$Y
Ebf1_G3F$Difference_x <- (Ebf1_G3F$Midline_x - Ebf1_G3F$Position_x)*4.202
Ebf1_G3F$Difference_y <- (Ebf1_G3F$Midline_y - Ebf1_G3F$Position_y)*4.202
#Htr4
Htr4_G3F <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G3F_NAcfull_Htr4_Results.csv",header = TRUE)
Htr4_G3F$Category <- ifelse(Htr4_G3F$Max >= 30, "Pos", "Neg")
Htr4_G3F$Midline_x  <- 1122
Htr4_G3F$Midline_y  <- 368
Htr4_G3F$Position_x <- G3F_ROIs$X
Htr4_G3F$Position_y <- G3F_ROIs$Y
Htr4_G3F$Difference_x <- (Htr4_G3F$Midline_x - Htr4_G3F$Position_x)*4.202
Htr4_G3F$Difference_y <- (Htr4_G3F$Midline_y - Htr4_G3F$Position_y)*4.202


#G4M 
#Drd1
Drd1_G4M <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G4M_NAcfull_Drd1_Results.csv",header = TRUE)
Drd1_G4M$Category <- ifelse(Drd1_G4M$Max >= 100, "Pos", "Neg")
Drd1_G4M$Midline_x  <- 1208
Drd1_G4M$Midline_y  <- 396
Drd1_G4M$Position_x <- G4M_ROIs$X
Drd1_G4M$Position_y <- G4M_ROIs$Y
Drd1_G4M$Difference_x <- (Drd1_G4M$Midline_x - Drd1_G4M$Position_x)*4.202
Drd1_G4M$Difference_y <- (Drd1_G4M$Midline_y - Drd1_G4M$Position_y)*4.202
#Ebf1
Ebf1_G4M <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G4M_NAcfull_Ebf1_Results.csv",header = TRUE)
Ebf1_G4M$Category <- ifelse(Ebf1_G4M$Max >= 80, "Pos", "Neg")
Ebf1_G4M$Midline_x  <- 1208
Ebf1_G4M$Midline_y  <- 396
Ebf1_G4M$Position_x <- G4M_ROIs$X
Ebf1_G4M$Position_y <- G4M_ROIs$Y
Ebf1_G4M$Difference_x <- (Drd1_G4M$Midline_x - Drd1_G4M$Position_x)*4.202
Ebf1_G4M$Difference_y <- (Drd1_G4M$Midline_y - Drd1_G4M$Position_y)*4.202
#Htr4
Htr4_G4M <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G4M_NAcfull_Htr4_Results.csv",header = TRUE)
Htr4_G4M$Category <- ifelse(Htr4_G4M$Max >= 30, "Pos","Neg")
Htr4_G4M$Midline_x  <- 1208
Htr4_G4M$Midline_y  <- 396
Htr4_G4M$Position_x <- G4M_ROIs$X
Htr4_G4M$Position_y <- G4M_ROIs$Y
Htr4_G4M$Difference_x <- (Htr4_G4M$Midline_x - Htr4_G4M$Position_x)*4.202
Htr4_G4M$Difference_y <- (Htr4_G4M$Midline_y - Htr4_G4M$Position_y)*4.202


#G5M 
#Drd1
Drd1_G5M <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G5M_NAcfull_Drd1_Results.csv",header = TRUE)
Drd1_G5M$Category <- ifelse(Drd1_G5M$Max >= 100, "Pos", "Neg")
Drd1_G5M$Midline_x  <- 1142
Drd1_G5M$Midline_y  <- 560
Drd1_G5M$Position_x <- G5M_ROIs$X
Drd1_G5M$Position_y <- G5M_ROIs$Y
Drd1_G5M$Difference_x <- (Drd1_G5M$Midline_x - Drd1_G5M$Position_x)*4.202
Drd1_G5M$Difference_y <- (Drd1_G5M$Midline_y - Drd1_G5M$Position_y)*4.202
#Ebf1
Ebf1_G5M <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G5M_NAcfull_Ebf1_Results.csv",header = TRUE)
Ebf1_G5M$Category <- ifelse(Ebf1_G5M$Max >= 80, "Pos", "Neg")
Ebf1_G5M$Midline_x  <- 1142
Ebf1_G5M$Midline_y  <- 560
Ebf1_G5M$Position_x <- G5M_ROIs$X
Ebf1_G5M$Position_y <- G5M_ROIs$Y
Ebf1_G5M$Difference_x <- (Drd1_G5M$Midline_x - Drd1_G5M$Position_x)*403.39
Ebf1_G5M$Difference_y <- (Drd1_G5M$Midline_y - Drd1_G5M$Position_y)*403.39
#Htr4
Htr4_G5M <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G5M_NAcfull_Htr4_Results.csv",header = TRUE)
Htr4_G5M$Category <- ifelse(Htr4_G5M$Max >= 30, "Pos","Neg")
Htr4_G5M$Midline_x  <- 1142
Htr4_G5M$Midline_y  <- 560
Htr4_G5M$Position_x <- G5M_ROIs$X
Htr4_G5M$Position_y <- G5M_ROIs$Y
Htr4_G5M$Difference_x <- (Htr4_G5M$Midline_x - Htr4_G5M$Position_x)*403.39
Htr4_G5M$Difference_y <- (Htr4_G5M$Midline_y - Htr4_G5M$Position_y)*403.39


#G6M 
#Drd1
Drd1_G6M <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G6M_NAcfull_Drd1_Results.csv",header = TRUE)
Drd1_G6M$Category <- ifelse(Drd1_G6M$Max >= 100, "Pos", "Neg")
Drd1_G6M$Midline_x  <- 13.77
Drd1_G6M$Midline_y  <- 3.48
Drd1_G6M$Position_x <- G6M_ROIs$X
Drd1_G6M$Position_y <- G6M_ROIs$Y
Drd1_G6M$Difference_x <- (Drd1_G6M$Midline_x - Drd1_G6M$Position_x)*403.39 #G6M wasn't using pixels?
Drd1_G6M$Difference_y <- (Drd1_G6M$Midline_y - Drd1_G6M$Position_y)*403.39
#Ebf1
Ebf1_G6M <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G6M_NAcfull_Ebf1_Results.csv",header = TRUE)
Ebf1_G6M$Category <- ifelse(Ebf1_G6M$Max >= 80, "Pos", "Neg")
Ebf1_G6M$Midline_x  <- 13.77
Ebf1_G6M$Midline_y  <- 3.48
Ebf1_G6M$Position_x <- G6M_ROIs$X
Ebf1_G6M$Position_y <- G6M_ROIs$Y
Ebf1_G6M$Difference_x <- (Ebf1_G6M$Midline_x - Ebf1_G6M$Position_x)*403.39
Ebf1_G6M$Difference_y <- (Ebf1_G6M$Midline_y - Ebf1_G6M$Position_y)*403.39
#Htr4
Htr4_G6M <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G6M_NAcfull_Htr4_Results.csv",header = TRUE)
Htr4_G6M$Category <- ifelse(Htr4_G6M$Max >=30, "Pos","Neg")
Htr4_G6M$Midline_x  <- 13.77
Htr4_G6M$Midline_y  <- 3.48
Htr4_G6M$Position_x <- G6M_ROIs$X
Htr4_G6M$Position_y <- G6M_ROIs$Y
Htr4_G6M$Difference_x <- (Htr4_G6M$Midline_x - Htr4_G6M$Position_x)*403.39
Htr4_G6M$Difference_y <- (Htr4_G6M$Midline_y - Htr4_G6M$Position_y)*403.39


#combo all
Drd1_all <- rbind(Drd1_G4M,Drd1_G2F,Drd1_G1F,Drd1_G3F,Drd1_G5M,Drd1_G6M)
Ebf1_all <- rbind(Ebf1_G4M,Ebf1_G2F,Ebf1_G1F,Ebf1_G3F,Ebf1_G5M,Ebf1_G6M)
Htr4_all <- rbind(Htr4_G4M,Htr4_G2F,Htr4_G1F,Htr4_G3F,Htr4_G5M,Htr4_G6M)
All <- rbind(Htr4_G4M,Htr4_G2F,Htr4_G1F,Htr4_G3F,Htr4_G5M,Htr4_G6M,Drd1_G4M,Drd1_G2F,Drd1_G1F,Drd1_G3F,Drd1_G5M,Drd1_G6M,Ebf1_G4M,Ebf1_G2F,Ebf1_G1F,Ebf1_G3F,Ebf1_G5M,Ebf1_G6M)

All <- data.frame(Drd1_Max = Drd1_all$Max,
                  Drd1_Raw = Drd1_all$RawIntDen,
                  Drd1_Cat = Drd1_all$Category,
                  Ebf1_Max = Ebf1_all$Max,
                  Ebf1_Raw = Ebf1_all$RawIntDen,
                  Ebf1_Cat = Ebf1_all$Category,
                  Htr4_Max = Htr4_all$Max,
                  Htr4_Raw = Htr4_all$RawIntDen,
                  Htr4_Cat = Htr4_all$Category,
                  Midline_x = Drd1_all$Midline_x,
                  Midline_y = Drd1_all$Midline_y,
                  Position_x = Drd1_all$Position_x,
                  Position_y = Drd1_all$Position_y,
                  Difference_x = Drd1_all$Difference_x,
                  Difference_y = Drd1_all$Difference_y,
                  Category = NA)

All[which(All$Drd1_Cat == "Pos" & All$Ebf1_Cat == "Neg" & All$Htr4_Cat == "Neg"), "Category"] <- "Drd1+ only"
All[which(All$Drd1_Cat == "Pos" & All$Ebf1_Cat == "Pos" & All$Htr4_Cat == "Neg"), "Category"] <- "Drd1+/Ebf1+"
All[which(All$Drd1_Cat == "Pos" & All$Ebf1_Cat == "Neg" & All$Htr4_Cat == "Pos"), "Category"] <- "Drd1+/Htr4+"
All[which(All$Drd1_Cat == "Pos" & All$Ebf1_Cat == "Pos" & All$Htr4_Cat == "Pos"), "Category"] <- "Drd1+/Htr4+/Ebf1+"

All %>% count(Category)


library(ggplot2)

p1 <- ggplot(data = subset(All,subset=(Category == "Drd1+ only")),aes(x = Difference_x,y = Difference_y)) +
  geom_point(colour = "maroon",alpha = 0.5,shape = 16,size = 3) +
  theme_classic() +
  coord_fixed(1) +
  #NoGrid() +
  xlim(c(-6000,6000)) +
  ylim(c(-3000,3000)) +
  ggtitle("Drd1+/Ebf1+") +
  xlab("Distance from\nX Midline (uM)") +
  ylab("Distance from\nY Midline (uM)") +
  labs(x = "",
       y = "")

p2 <- ggplot(data = subset(All,subset=(Category == "Drd1+/Ebf1+")),aes(x = Difference_x,y = Difference_y)) +
  geom_point(colour = "darkseagreen",alpha = 0.5,shape = 16,size = 3) +
  theme_classic() +
  coord_fixed(1) +
  #NoGrid() +
  xlim(c(-6000,6000)) +
  ylim(c(-3000,3000)) +
  ggtitle("Drd1+/Ebf1+") +
  xlab("Distance from\nX Midline (uM)") +
  ylab("Distance from\nY Midline (uM)") +
  labs(x = "",
       y = "")


p3 <- ggplot(data = subset(All,subset=(Category == "Drd1+/Htr4+")),aes(x = Difference_x,y = Difference_y)) +
  geom_point(colour = "cornflowerblue",alpha = 0.5,shape = 16,size = 3) +
  theme_classic() +
  coord_fixed(1) +
  #NoGrid() +
  xlim(c(-6000,6000)) +
  ylim(c(-3000,3000)) +
  ggtitle("Drd1+/Htr4+") +
  xlab("Distance from\nX Midline (uM)") +
  ylab("Distance from\nY Midline (uM)") +
  labs(x = "",
       y = "")


#p5 <- ggplot(data = subset(All,subset=(Category == "Drd1+/Htr4+/Ebf1+")),aes(x = Difference_x,y = Difference_y)) +
# geom_point(colour = "grey59",alpha = 0.5,shape = 16,size = 3) +
#theme_classic() +
#coord_fixed(1) +
##xlim(c(-6000,6000)) +
#ylim(c(-3000,3000)) +
#ggtitle("Drd1+/Htr4+/Ebf1+") +
#xlab("Distance from\nX Midline (uM)") +
#ylab("Distance from\nY Midline (uM)") +
#labs(x = "",
#    y = "")

p6 <- ggplot(data = subset(All,subset=(Category == "Drd1+/Htr4+/Ebf1+")),aes(x=Difference_x,y=Difference_y)) +
  geom_point(colour = "grey59",alpha=0.5,shape=16,size=3) +
  theme_classic() +
  coord_fixed(1) +
  xlim(c(-6000,6000)) +
  ylim(c(-3000,3000)) +
  ggtitle("Drd1+/Htr4+/Ebf1+") +
  geom_segment(aes(x = 1000, y = 750, xend = 2000, yend = 750)) +
  annotate("text",x = 1500,y = 800,label = "1000uM") +
  geom_segment(aes(x = -1000, y = 750, xend = -1500, yend = 750)) +
  annotate("text",x = -1250,y = 800,label = "500uM")

p4 <- ggplot(data = All,aes(x = Difference_x,y = Difference_y)) +
  #geom_point(data = subset(Drd1_All,subset=(Category == "Drd1+ Only")), colour = "limegreen",alpha = 0.5,shape = 16,size = 3) +
  geom_point(data = subset(All,subset=(Category == "Drd1+ only")), colour = "maroon",alpha = 0.5,shape = 16,size = 3) +
  geom_point(data = subset(All,subset=(Category == "Drd1+/Ebf1+")), colour = "darkseagreen",alpha = 0.5,shape = 16,size = 3) +
  geom_point(data = subset(All,subset=(Category == "Drd1+/Htr4+")), colour = "cornflowerblue",alpha = 0.5,shape = 16,size = 3) +
  #geom_point(data = subset(All,subset=(Category == "Drd1+/Htr4+/Ebf1+")), colour = "grey59",alpha = 0.5,shape = 16,size = 3) +
  theme_classic() +
  coord_fixed(1) +
  #NoGrid() +
  xlim(c(-6000,6000)) +
  ylim(c(-3000,3000)) +
  ggtitle("Merge") +
  xlab("Distance from\nX Midline (uM)") +
  ylab("Distance from\nY Midline (uM)") +
  labs(x = "",
       y = "")

ggsave(plot = p1,
       filename = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/Drd1posonly.pdf",
       height = 5.87,
       width = 16)

ggsave(plot = p6,
       filename = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/Drd1pos_Ebf1pos_Htr4pos_withSCALE.pdf",
       height = 5.87,
       width = 16)

ggsave(plot = p5,
       filename = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/Drd1pos_Ebf1pos_Htr4pos_withaxes1.pdf",
       height = 5.87,
       width = 16)

ggsave(plot = p2,
       filename = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/Drd1Pos_Ebf1Pos_withAxes1.pdf",
       height = 5.87,
       width = 16)
ggsave(plot = p3,
       filename = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/Drd1Pos_Htr4Pos_withAxes1.pdf",
       height = 5.87,
       width = 16)
ggsave(plot = p4,
       filename = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/Merge_withAxes1.pdf",
       height = 5.87,
       width = 16)


ggsave(plot = p2 + theme(plot.title = element_text(hjust = 0.5),
                         axis.line=element_blank(),
                         axis.text.x=element_blank(),
                         axis.text.y=element_blank(),
                         axis.ticks=element_blank(),
                         axis.title.x=element_blank(),
                         axis.title.y=element_blank()),
       filename = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/Drd1Pos_Ebf1Pos_withoutAxes1.png",
       height = 5.87,
       width = 16)
ggsave(plot = p3 + theme(plot.title = element_text(hjust = 0.5),
                         axis.line=element_blank(),
                         axis.text.x=element_blank(),
                         axis.text.y=element_blank(),
                         axis.ticks=element_blank(),
                         axis.title.x=element_blank(),
                         axis.title.y=element_blank()),
       filename = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/Drd1Pos_Htr4Pos_withoutAxes1.pdf",
       height = 5.87,
       width = 16)
ggsave(plot = p4 + theme(plot.title = element_text(hjust = 0.5),
                         axis.line=element_blank(),
                         axis.text.x=element_blank(),
                         axis.text.y=element_blank(),
                         axis.ticks=element_blank(),
                         axis.title.x=element_blank(),
                         axis.title.y=element_blank()),
       filename = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/Merge_withoutAxes1.pdf",
       height = 5.87,
       width = 16)

ggsave(plot = p1 + theme(plot.title = element_text(hjust = 0.5),
                         axis.line=element_blank(),
                         axis.text.x=element_blank(),
                         axis.text.y=element_blank(),
                         axis.ticks=element_blank(),
                         axis.title.x=element_blank(),
                         axis.title.y=element_blank()),
       filename = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/Drd1only_withoutaxes1.pdf",
       height = 5.87,
       width = 16)

#SessionInfo 20230103
R version 4.2.2 (2022-10-31)
Platform: aarch64-apple-darwin20 (64-bit)
Running under: macOS Ventura 13.1

Matrix products: default
LAPACK: /Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/lib/libRlapack.dylib

locale:
  [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
  [1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
  [1] ggplot2_3.4.0 scales_1.2.1  dplyr_1.0.10 

loaded via a namespace (and not attached):
  [1] pillar_1.8.1     compiler_4.2.2   bslib_0.4.2      jquerylib_0.1.4  tools_4.2.2      digest_0.6.31    jsonlite_1.8.4   evaluate_0.19    lifecycle_1.0.3 
[10] tibble_3.1.8     gtable_0.3.1     pkgconfig_2.0.3  rlang_1.0.6      cli_3.5.0        rstudioapi_0.14  yaml_2.3.6       xfun_0.36        fastmap_1.1.0   
[19] withr_2.5.0      stringr_1.5.0    knitr_1.41       generics_0.1.3   vctrs_0.5.1      sass_0.4.4       grid_4.2.2       tidyselect_1.2.0 glue_1.6.2      
[28] R6_2.5.1         fansi_1.0.3      rmarkdown_2.19   farver_2.1.1     magrittr_2.0.3   htmltools_0.5.4  rsconnect_0.8.28 colorspace_2.0-3 labeling_0.4.2  
[37] utf8_1.2.2       stringi_1.7.8    munsell_0.5.0    cachem_1.0.6    
