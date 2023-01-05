library(dplyr)

#NAc Whole

####G1F####
#read in DAPI ROIs
G1F_ROIs       <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G1F_NAc_Dapi_Results.csv",header = TRUE)
G1F_ROIs$Index <- 1:nrow(G1F_ROIs)
#Ebf1
G1F_Ebf1 <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G1F_NAc_Ebf1_Results.csv",header = TRUE)
G1F_Ebf1$Sample <- "G1F"
G1F_Ebf1$Category <- ifelse(G1F_Ebf1$Max >70,"Pos","Neg")
G1F_Ebf1 <- select(G1F_Ebf1, c("Category","Mean", "Sample"))
#Drd1
G1F_Drd1 <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G1F_NAc_Drd1_Results.csv",header = TRUE)
G1F_Drd1$Sample <- "G1F"
G1F_Drd1$Category <- ifelse(G1F_Drd1$Mean >20,"Pos","Neg")
G1F_Drd1 <- select(G1F_Drd1, c("Category","Mean", "Sample"))
#Htr4
G1F_Htr4 <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G1F_NAc_Htr4_Results.csv",header = TRUE)
G1F_Htr4$Category <- ifelse(G1F_Htr4$Max >35,"Pos","Neg")
G1F_Htr4$Sample <- "G1F"
G1F_Htr4 <- select(G1F_Htr4, c("Category","Mean", "Sample"))

#Categories
G1F_All_Data <- data.frame(Ebf1               = G1F_Ebf1$Mean,
                           Ebf1_Category      = G1F_Ebf1$Category,
                           Drd1          = G1F_Drd1$Mean,
                           Drd1_Category = G1F_Drd1$Category,
                           Htr4             = G1F_Htr4$Mean,
                           Htr4_Category    = G1F_Htr4$Category, 
                           Sex = "Female",
                           Category         = NA)

G1F_All_Data$Category <- ifelse(G1F_All_Data$Htr4_Category == "Pos" & G1F_All_Data$Drd1_Category == "Neg" & G1F_All_Data$Ebf1_Category == "Pos",
                                "Drd1- Ebf1+ Htr4+",
                                ifelse(G1F_All_Data$Htr4_Category == "Neg" & G1F_All_Data$Drd1_Category == "Pos" & G1F_All_Data$Ebf1_Category == "Pos",
                                       "Drd1+ Ebf1+ Htr4-",
                                       ifelse(G1F_All_Data$Htr4_Category == "Pos" & G1F_All_Data$Drd1_Category == "Pos" & G1F_All_Data$Ebf1_Category == "Pos",
                                              "Drd1+ Ebf1+ Htr4+",
                                              ifelse(G1F_All_Data$Htr4_Category == "Neg" & G1F_All_Data$Drd1_Category == "Pos" & G1F_All_Data$Ebf1_Category == "Neg",
                                                     "Drd1+ Ebf1- Htr4-",
                                                     ifelse(G1F_All_Data$Htr4_Category == "Pos" & G1F_All_Data$Drd1_Category == "Neg" & G1F_All_Data$Ebf1_Category == "Neg",
                                                            "Drd1- Ebf1- Htr4+",
                                                            ifelse(G1F_All_Data$Htr4_Category == "Neg" & G1F_All_Data$Drd1_Category == "Neg" & G1F_All_Data$Ebf1_Category == "Pos",
                                                                   "Drd1- Ebf1+ Htr4-",
                                                                   ifelse(G1F_All_Data$Htr4_Category == "Pos" & G1F_All_Data$Drd1_Category == "Pos" & G1F_All_Data$Ebf1_Category == "Neg",
                                                                          "Drd1+ Ebf1- Htr4+", "Other")))))))

G1F_All_Data %>% count(Category)


#G2F
#Read in the DAPI ROI positions
G2F_ROIs       <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G2F_NAc_Dapi_Results.csv",header = TRUE)
G2F_ROIs$Index <- 1:nrow(G2F_ROIs)
#Ebf1
G2F_Ebf1 <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G2F_NAc_Ebf1_Results.csv",header = TRUE)
G2F_Ebf1$Sample <- "G2F"
G2F_Ebf1$Category <- ifelse(G2F_Ebf1$Max >30,"Pos","Neg")
G2F_Ebf1 <- select(G2F_Ebf1, c("Category","Mean", "Sample"))
#Drd1
G2F_Drd1 <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G2F_NAc_Drd1_Results.csv",header = TRUE)
G2F_Drd1$Sample <- "G2F"
G2F_Drd1$Category <- ifelse(G2F_Drd1$Mean >6,"Pos","Neg")
G2F_Drd1 <- select(G2F_Drd1, c("Category","Mean", "Sample"))
#Htr4
G2F_Htr4 <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G2F_NAc_Htr4_Results.csv",header = TRUE)
G2F_Htr4$Category <- ifelse(G2F_Htr4$Max >40,"Pos","Neg")
G2F_Htr4$Sample <- "G2F"
G2F_Htr4 <- select(G2F_Htr4, c("Category","Mean", "Sample"))

#Categories
G2F_All_Data <- data.frame(Ebf1               = G2F_Ebf1$Mean,
                           Ebf1_Category      = G2F_Ebf1$Category,
                           Drd1          = G2F_Drd1$Mean,
                           Drd1_Category = G2F_Drd1$Category,
                           Htr4             = G2F_Htr4$Mean,
                           Htr4_Category    = G2F_Htr4$Category, 
                           Sex = "Female",
                           Category         = NA)

G2F_All_Data$Category <- ifelse(G2F_All_Data$Htr4_Category == "Pos" & G2F_All_Data$Drd1_Category == "Neg" & G2F_All_Data$Ebf1_Category == "Pos",
                                "Drd1- Ebf1+ Htr4+",
                                ifelse(G2F_All_Data$Htr4_Category == "Neg" & G2F_All_Data$Drd1_Category == "Pos" & G2F_All_Data$Ebf1_Category == "Pos",
                                       "Drd1+ Ebf1+ Htr4-",
                                       ifelse(G2F_All_Data$Htr4_Category == "Pos" & G2F_All_Data$Drd1_Category == "Pos" & G2F_All_Data$Ebf1_Category == "Pos",
                                              "Drd1+ Ebf1+ Htr4+",
                                              ifelse(G2F_All_Data$Htr4_Category == "Neg" & G2F_All_Data$Drd1_Category == "Pos" & G2F_All_Data$Ebf1_Category == "Neg",
                                                     "Drd1+ Ebf1- Htr4-",
                                                     ifelse(G2F_All_Data$Htr4_Category == "Pos" & G2F_All_Data$Drd1_Category == "Neg" & G2F_All_Data$Ebf1_Category == "Neg",
                                                            "Drd1- Ebf1- Htr4+",
                                                            ifelse(G2F_All_Data$Htr4_Category == "Neg" & G2F_All_Data$Drd1_Category == "Neg" & G2F_All_Data$Ebf1_Category == "Pos",
                                                                   "Drd1- Ebf1+ Htr4-",
                                                                   ifelse(G2F_All_Data$Htr4_Category == "Pos" & G2F_All_Data$Drd1_Category == "Pos" & G2F_All_Data$Ebf1_Category == "Neg",
                                                                          "Drd1+ Ebf1- Htr4+", "Other")))))))

G2F_All_Data %>% count(Category)


#G3F
#Read in the DAPI ROI positions
G3F_ROIs       <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G3F_NAc_Dapi_Results.csv",header = TRUE)
G3F_ROIs$Index <- 1:nrow(G3F_ROIs)
#Ebf1
G3F_Ebf1 <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G3F_NAc_Ebf1_Results.csv",header = TRUE)
G3F_Ebf1$Sample <- "G3F"
G3F_Ebf1$Category <- ifelse(G3F_Ebf1$Max >30,"Pos","Neg")
G3F_Ebf1 <- select(G3F_Ebf1, c("Category","Mean", "Sample"))
#Drd1
G3F_Drd1 <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G3F_NAc_Drd1_Results.csv",header = TRUE)
G3F_Drd1$Sample <- "G3F"
G3F_Drd1$Category <- ifelse(G3F_Drd1$Mean >6,"Pos","Neg")
G3F_Drd1 <- select(G3F_Drd1, c("Category","Mean", "Sample"))
#Htr4
G3F_Htr4 <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G3F_NAc_Htr4_Results.csv",header = TRUE)
G3F_Htr4$Category <- ifelse(G3F_Htr4$Max >60,"Pos","Neg")
G3F_Htr4$Sample <- "G3F"
G3F_Htr4 <- select(G3F_Htr4, c("Category","Mean", "Sample"))

#Categories
G3F_All_Data <- data.frame(Ebf1               = G3F_Ebf1$Mean,
                           Ebf1_Category      = G3F_Ebf1$Category,
                           Drd1          = G3F_Drd1$Mean,
                           Drd1_Category = G3F_Drd1$Category,
                           Htr4             = G3F_Htr4$Mean,
                           Htr4_Category    = G3F_Htr4$Category, 
                           Sex = "Female",
                           Category         = NA)

G3F_All_Data$Category <- ifelse(G3F_All_Data$Htr4_Category == "Pos" & G3F_All_Data$Drd1_Category == "Neg" & G3F_All_Data$Ebf1_Category == "Pos",
                                "Drd1- Ebf1+ Htr4+",
                                ifelse(G3F_All_Data$Htr4_Category == "Neg" & G3F_All_Data$Drd1_Category == "Pos" & G3F_All_Data$Ebf1_Category == "Pos",
                                       "Drd1+ Ebf1+ Htr4-",
                                       ifelse(G3F_All_Data$Htr4_Category == "Pos" & G3F_All_Data$Drd1_Category == "Pos" & G3F_All_Data$Ebf1_Category == "Pos",
                                              "Drd1+ Ebf1+ Htr4+",
                                              ifelse(G3F_All_Data$Htr4_Category == "Neg" & G3F_All_Data$Drd1_Category == "Pos" & G3F_All_Data$Ebf1_Category == "Neg",
                                                     "Drd1+ Ebf1- Htr4-",
                                                     ifelse(G3F_All_Data$Htr4_Category == "Pos" & G3F_All_Data$Drd1_Category == "Neg" & G3F_All_Data$Ebf1_Category == "Neg",
                                                            "Drd1- Ebf1- Htr4+",
                                                            ifelse(G3F_All_Data$Htr4_Category == "Neg" & G3F_All_Data$Drd1_Category == "Neg" & G3F_All_Data$Ebf1_Category == "Pos",
                                                                   "Drd1- Ebf1+ Htr4-",
                                                                   ifelse(G3F_All_Data$Htr4_Category == "Pos" & G3F_All_Data$Drd1_Category == "Pos" & G3F_All_Data$Ebf1_Category == "Neg",
                                                                          "Drd1+ Ebf1- Htr4+", "Other")))))))

G3F_All_Data %>% count(Category)


####G4M####
G4M_ROIs       <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G4M_NAc_Dapi_Results.csv",header = TRUE)
G4M_ROIs$Index <- 1:nrow(G4M_ROIs)
#Ebf1
G4M_Ebf1 <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G4M_NAc_Ebf1_Results.csv",header = TRUE)
G4M_Ebf1$Sample <- "G4M"
G4M_Ebf1$Category <- ifelse(G4M_Ebf1$Max >30,"Pos","Neg")
G4M_Ebf1 <- select(G4M_Ebf1, c("Category","Mean", "Sample"))
#Drd1
G4M_Drd1 <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G4M_NAc_Drd1_Results.csv",header = TRUE)
G4M_Drd1$Sample <- "G4M"
G4M_Drd1$Category <- ifelse(G4M_Drd1$Mean >6,"Pos","Neg")
G4M_Drd1 <- select(G4M_Drd1, c("Category","Mean", "Sample"))
#Htr4
G4M_Htr4 <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G4M_NAc_Htr4_Results.csv",header = TRUE)
G4M_Htr4$Category <- ifelse(G4M_Htr4$Max >40,"Pos","Neg")
G4M_Htr4$Sample <- "G4M"
G4M_Htr4 <- select(G4M_Htr4, c("Category","Mean", "Sample"))

#Categories
G4M_All_Data <- data.frame(Ebf1               = G4M_Ebf1$Mean,
                           Ebf1_Category      = G4M_Ebf1$Category,
                           Drd1          = G4M_Drd1$Mean,
                           Drd1_Category = G4M_Drd1$Category,
                           Htr4             = G4M_Htr4$Mean,
                           Htr4_Category    = G4M_Htr4$Category, 
                           Sex = "Female",
                           Category         = NA)

G4M_All_Data$Category <- ifelse(G4M_All_Data$Htr4_Category == "Pos" & G4M_All_Data$Drd1_Category == "Neg" & G4M_All_Data$Ebf1_Category == "Pos",
                                "Drd1- Ebf1+ Htr4+",
                                ifelse(G4M_All_Data$Htr4_Category == "Neg" & G4M_All_Data$Drd1_Category == "Pos" & G4M_All_Data$Ebf1_Category == "Pos",
                                       "Drd1+ Ebf1+ Htr4-",
                                       ifelse(G4M_All_Data$Htr4_Category == "Pos" & G4M_All_Data$Drd1_Category == "Pos" & G4M_All_Data$Ebf1_Category == "Pos",
                                              "Drd1+ Ebf1+ Htr4+",
                                              ifelse(G4M_All_Data$Htr4_Category == "Neg" & G4M_All_Data$Drd1_Category == "Pos" & G4M_All_Data$Ebf1_Category == "Neg",
                                                     "Drd1+ Ebf1- Htr4-",
                                                     ifelse(G4M_All_Data$Htr4_Category == "Pos" & G4M_All_Data$Drd1_Category == "Neg" & G4M_All_Data$Ebf1_Category == "Neg",
                                                            "Drd1- Ebf1- Htr4+",
                                                            ifelse(G4M_All_Data$Htr4_Category == "Neg" & G4M_All_Data$Drd1_Category == "Neg" & G4M_All_Data$Ebf1_Category == "Pos",
                                                                   "Drd1- Ebf1+ Htr4-",
                                                                   ifelse(G4M_All_Data$Htr4_Category == "Pos" & G4M_All_Data$Drd1_Category == "Pos" & G4M_All_Data$Ebf1_Category == "Neg",
                                                                          "Drd1+ Ebf1- Htr4+", "Other")))))))

G4M_All_Data %>% count(Category)


####G5M####
G5M_ROIs       <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G5M_NAc_Dapi_Results.csv",header = TRUE)
G5M_ROIs$Index <- 1:nrow(G5M_ROIs)
#Ebf1
G5M_Ebf1 <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G5M_NAc_Ebf1_Results.csv",header = TRUE)
G5M_Ebf1$Sample <- "G5M"
G5M_Ebf1$Category <- ifelse(G5M_Ebf1$Max >20,"Pos","Neg")
G5M_Ebf1 <- select(G5M_Ebf1, c("Category","Mean", "Sample"))
#Drd1
G5M_Drd1 <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G5M_NAc_Drd1_Results.csv",header = TRUE)
G5M_Drd1$Sample <- "G5M"
G5M_Drd1$Category <- ifelse(G5M_Drd1$Mean >5,"Pos","Neg")
G5M_Drd1 <- select(G5M_Drd1, c("Category","Mean", "Sample"))
#Htr4
G5M_Htr4 <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G5M_NAc_Htr4_Results.csv",header = TRUE)
G5M_Htr4$Category <- ifelse(G5M_Htr4$Max >50,"Pos","Neg")
G5M_Htr4$Sample <- "G5M"
G5M_Htr4 <- select(G5M_Htr4, c("Category","Mean", "Sample"))

#Categories
G5M_All_Data <- data.frame(Ebf1               = G5M_Ebf1$Mean,
                           Ebf1_Category      = G5M_Ebf1$Category,
                           Drd1          = G5M_Drd1$Mean,
                           Drd1_Category = G5M_Drd1$Category,
                           Htr4             = G5M_Htr4$Mean,
                           Htr4_Category    = G5M_Htr4$Category, 
                           Sex = "Female",
                           Category         = NA)

G5M_All_Data$Category <- ifelse(G5M_All_Data$Htr4_Category == "Pos" & G5M_All_Data$Drd1_Category == "Neg" & G5M_All_Data$Ebf1_Category == "Pos",
                                "Drd1- Ebf1+ Htr4+",
                                ifelse(G5M_All_Data$Htr4_Category == "Neg" & G5M_All_Data$Drd1_Category == "Pos" & G5M_All_Data$Ebf1_Category == "Pos",
                                       "Drd1+ Ebf1+ Htr4-",
                                       ifelse(G5M_All_Data$Htr4_Category == "Pos" & G5M_All_Data$Drd1_Category == "Pos" & G5M_All_Data$Ebf1_Category == "Pos",
                                              "Drd1+ Ebf1+ Htr4+",
                                              ifelse(G5M_All_Data$Htr4_Category == "Neg" & G5M_All_Data$Drd1_Category == "Pos" & G5M_All_Data$Ebf1_Category == "Neg",
                                                     "Drd1+ Ebf1- Htr4-",
                                                     ifelse(G5M_All_Data$Htr4_Category == "Pos" & G5M_All_Data$Drd1_Category == "Neg" & G5M_All_Data$Ebf1_Category == "Neg",
                                                            "Drd1- Ebf1- Htr4+",
                                                            ifelse(G5M_All_Data$Htr4_Category == "Neg" & G5M_All_Data$Drd1_Category == "Neg" & G5M_All_Data$Ebf1_Category == "Pos",
                                                                   "Drd1- Ebf1+ Htr4-",
                                                                   ifelse(G5M_All_Data$Htr4_Category == "Pos" & G5M_All_Data$Drd1_Category == "Pos" & G5M_All_Data$Ebf1_Category == "Neg",
                                                                          "Drd1+ Ebf1- Htr4+", "Other")))))))

G5M_All_Data %>% count(Category)

####G6M####
G6M_ROIs       <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G6M_NAc_Dapi_Results.csv",header = TRUE)
G6M_ROIs$Index <- 1:nrow(G6M_ROIs)
#Ebf1
G6M_Ebf1 <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G6M_NAc_Ebf1_Results.csv",header = TRUE)
G6M_Ebf1$Sample <- "G6M"
G6M_Ebf1$Category <- ifelse(G6M_Ebf1$Max >30,"Pos","Neg")
G6M_Ebf1 <- select(G6M_Ebf1, c("Category","Mean", "Sample"))
#Drd1
G6M_Drd1 <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G6M_NAc_Drd1_Results.csv",header = TRUE)
G6M_Drd1$Sample <- "G6M"
G6M_Drd1$Category <- ifelse(G6M_Drd1$Mean >3.5,"Pos","Neg")
G6M_Drd1 <- select(G6M_Drd1, c("Category","Mean", "Sample"))
#Htr4
G6M_Htr4 <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G6M_NAc_Htr4_Results.csv",header = TRUE)
G6M_Htr4$Category <- ifelse(G6M_Htr4$Max >30,"Pos","Neg")
G6M_Htr4$Sample <- "G6M"
G6M_Htr4 <- select(G6M_Htr4, c("Category","Mean", "Sample"))

#Categories
G6M_All_Data <- data.frame(Ebf1               = G6M_Ebf1$Mean,
                           Ebf1_Category      = G6M_Ebf1$Category,
                           Drd1          = G6M_Drd1$Mean,
                           Drd1_Category = G6M_Drd1$Category,
                           Htr4             = G6M_Htr4$Mean,
                           Htr4_Category    = G6M_Htr4$Category, 
                           Sex = "Female",
                           Category         = NA)

G6M_All_Data$Category <- ifelse(G6M_All_Data$Htr4_Category == "Pos" & G6M_All_Data$Drd1_Category == "Neg" & G6M_All_Data$Ebf1_Category == "Pos",
                                "Drd1- Ebf1+ Htr4+",
                                ifelse(G6M_All_Data$Htr4_Category == "Neg" & G6M_All_Data$Drd1_Category == "Pos" & G6M_All_Data$Ebf1_Category == "Pos",
                                       "Drd1+ Ebf1+ Htr4-",
                                       ifelse(G6M_All_Data$Htr4_Category == "Pos" & G6M_All_Data$Drd1_Category == "Pos" & G6M_All_Data$Ebf1_Category == "Pos",
                                              "Drd1+ Ebf1+ Htr4+",
                                              ifelse(G6M_All_Data$Htr4_Category == "Neg" & G6M_All_Data$Drd1_Category == "Pos" & G6M_All_Data$Ebf1_Category == "Neg",
                                                     "Drd1+ Ebf1- Htr4-",
                                                     ifelse(G6M_All_Data$Htr4_Category == "Pos" & G6M_All_Data$Drd1_Category == "Neg" & G6M_All_Data$Ebf1_Category == "Neg",
                                                            "Drd1- Ebf1- Htr4+",
                                                            ifelse(G6M_All_Data$Htr4_Category == "Neg" & G6M_All_Data$Drd1_Category == "Neg" & G6M_All_Data$Ebf1_Category == "Pos",
                                                                   "Drd1- Ebf1+ Htr4-",
                                                                   ifelse(G6M_All_Data$Htr4_Category == "Pos" & G6M_All_Data$Drd1_Category == "Pos" & G6M_All_Data$Ebf1_Category == "Neg",
                                                                          "Drd1+ Ebf1- Htr4+", "Other")))))))

G6M_All_Data %>% count(Category)




#NAc CORE

####G1F####
#Read in the DAPI ROI positions
G1FCore_ROIs       <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G1F_NAcC_Dapi_Results.csv",header = TRUE)
G1FCore_ROIs$Index <- 1:nrow(G1FCore_ROIs)
#Ebf1
G1FCore_Ebf1 <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G1F_NAcC_Ebf1_Results.csv",header = TRUE)
G1FCore_Ebf1$Sample <- "G1F"
G1FCore_Ebf1$Category <- ifelse(G1FCore_Ebf1$Max >30,"Pos","Neg")
G1FCore_Ebf1 <- select(G1FCore_Ebf1, c("Category","Mean", "Sample"))
#Drd1
G1FCore_Drd1 <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G1F_NAcC_Drd1_Results.csv",header = TRUE)
G1FCore_Drd1$Sample <- "G1F"
G1FCore_Drd1$Category <- ifelse(G1FCore_Drd1$Mean >5,"Pos","Neg")
G1FCore_Drd1 <- select(G1FCore_Drd1, c("Category","Mean", "Sample"))
#Htr4
G1FCore_Htr4 <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G1F_NAcC_Htr4_Results.csv",header = TRUE)
G1FCore_Htr4$Category <- ifelse(G1FCore_Htr4$Max >35,"Pos","Neg")
G1FCore_Htr4$Sample <- "G1F"
G1FCore_Htr4 <- select(G1FCore_Htr4, c("Category","Mean", "Sample"))

#Categories
G1FCore_All_Data <- data.frame(Ebf1               = G1FCore_Ebf1$Mean,
                               Ebf1_Category      = G1FCore_Ebf1$Category,
                               Drd1          = G1FCore_Drd1$Mean,
                               Drd1_Category = G1FCore_Drd1$Category,
                               Htr4             = G1FCore_Htr4$Mean,
                               Htr4_Category    = G1FCore_Htr4$Category, 
                               Sex = "Female",
                               Category         = NA)

G1FCore_All_Data$Category <- ifelse(G1FCore_All_Data$Htr4_Category == "Pos" & G1FCore_All_Data$Drd1_Category == "Neg" & G1FCore_All_Data$Ebf1_Category == "Pos",
                                    "Drd1- Ebf1+ Htr4+",
                                    ifelse(G1FCore_All_Data$Htr4_Category == "Neg" & G1FCore_All_Data$Drd1_Category == "Pos" & G1FCore_All_Data$Ebf1_Category == "Pos",
                                           "Drd1+ Ebf1+ Htr4-",
                                           ifelse(G1FCore_All_Data$Htr4_Category == "Pos" & G1FCore_All_Data$Drd1_Category == "Pos" & G1FCore_All_Data$Ebf1_Category == "Pos",
                                                  "Drd1+ Ebf1+ Htr4+",
                                                  ifelse(G1FCore_All_Data$Htr4_Category == "Neg" & G1FCore_All_Data$Drd1_Category == "Pos" & G1FCore_All_Data$Ebf1_Category == "Neg",
                                                         "Drd1+ Ebf1- Htr4-",
                                                         ifelse(G1FCore_All_Data$Htr4_Category == "Pos" & G1FCore_All_Data$Drd1_Category == "Neg" & G1FCore_All_Data$Ebf1_Category == "Neg",
                                                                "Drd1- Ebf1- Htr4+",
                                                                ifelse(G1FCore_All_Data$Htr4_Category == "Neg" & G1FCore_All_Data$Drd1_Category == "Neg" & G1FCore_All_Data$Ebf1_Category == "Pos",
                                                                       "Drd1- Ebf1+ Htr4-",
                                                                       ifelse(G1FCore_All_Data$Htr4_Category == "Pos" & G1FCore_All_Data$Drd1_Category == "Pos" & G1FCore_All_Data$Ebf1_Category == "Neg",
                                                                              "Drd1+ Ebf1- Htr4+", "Other")))))))

G1FCore_All_Data %>% count(Category)
######G2F#####
G2FCore_ROIs       <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G2F_NAcC_Dapi_Results.csv",header = TRUE)
G2FCore_ROIs$Index <- 1:nrow(G2FCore_ROIs)
#Ebf1
G2FCore_Ebf1 <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G2F_NAcC_Ebf1_Results.csv",header = TRUE)
G2FCore_Ebf1$Sample <- "G2F"
G2FCore_Ebf1$Category <- ifelse(G2FCore_Ebf1$Max >20,"Pos","Neg")
G2FCore_Ebf1 <- select(G2FCore_Ebf1, c("Category","Mean", "Sample"))
#Drd1
G2FCore_Drd1 <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G2F_NAcC_Drd1_Results.csv",header = TRUE)
G2FCore_Drd1$Sample <- "G2F"
G2FCore_Drd1$Category <- ifelse(G2FCore_Drd1$Mean >3.5,"Pos","Neg")
G2FCore_Drd1 <- select(G2FCore_Drd1, c("Category","Mean", "Sample"))
#Htr4
G2FCore_Htr4 <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G2F_NAcC_Htr4_Results.csv",header = TRUE)
G2FCore_Htr4$Category <- ifelse(G2FCore_Htr4$Max >30,"Pos","Neg")
G2FCore_Htr4$Sample <- "G2F"
G2FCore_Htr4 <- select(G2FCore_Htr4, c("Category","Mean", "Sample"))

#Categories
G2FCore_All_Data <- data.frame(Ebf1               = G2FCore_Ebf1$Mean,
                               Ebf1_Category      = G2FCore_Ebf1$Category,
                               Drd1          = G2FCore_Drd1$Mean,
                               Drd1_Category = G2FCore_Drd1$Category,
                               Htr4             = G2FCore_Htr4$Mean,
                               Htr4_Category    = G2FCore_Htr4$Category, 
                               Sex = "Female",
                               Category         = NA)

G2FCore_All_Data$Category <- ifelse(G2FCore_All_Data$Htr4_Category == "Pos" & G2FCore_All_Data$Drd1_Category == "Neg" & G2FCore_All_Data$Ebf1_Category == "Pos",
                                    "Drd1- Ebf1+ Htr4+",
                                    ifelse(G2FCore_All_Data$Htr4_Category == "Neg" & G2FCore_All_Data$Drd1_Category == "Pos" & G2FCore_All_Data$Ebf1_Category == "Pos",
                                           "Drd1+ Ebf1+ Htr4-",
                                           ifelse(G2FCore_All_Data$Htr4_Category == "Pos" & G2FCore_All_Data$Drd1_Category == "Pos" & G2FCore_All_Data$Ebf1_Category == "Pos",
                                                  "Drd1+ Ebf1+ Htr4+",
                                                  ifelse(G2FCore_All_Data$Htr4_Category == "Neg" & G2FCore_All_Data$Drd1_Category == "Pos" & G2FCore_All_Data$Ebf1_Category == "Neg",
                                                         "Drd1+ Ebf1- Htr4-",
                                                         ifelse(G2FCore_All_Data$Htr4_Category == "Pos" & G2FCore_All_Data$Drd1_Category == "Neg" & G2FCore_All_Data$Ebf1_Category == "Neg",
                                                                "Drd1- Ebf1- Htr4+",
                                                                ifelse(G2FCore_All_Data$Htr4_Category == "Neg" & G2FCore_All_Data$Drd1_Category == "Neg" & G2FCore_All_Data$Ebf1_Category == "Pos",
                                                                       "Drd1- Ebf1+ Htr4-",
                                                                       ifelse(G2FCore_All_Data$Htr4_Category == "Pos" & G2FCore_All_Data$Drd1_Category == "Pos" & G2FCore_All_Data$Ebf1_Category == "Neg",
                                                                              "Drd1+ Ebf1- Htr4+", "Other")))))))

G2FCore_All_Data %>% count(Category)

###G3F#####
#Read in the DAPI ROI positions
G3FCore_ROIs       <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G3F_NAcC_dapi_Results.csv",header = TRUE)
G3FCore_ROIs$Index <- 1:nrow(G3FCore_ROIs)
#Ebf1
G3FCore_Ebf1 <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G3F_NAcC_Ebf1_Results.csv",header = TRUE)
G3FCore_Ebf1$Sample <- "G3F"
G3FCore_Ebf1$Category <- ifelse(G3FCore_Ebf1$Max >25,"Pos","Neg")
G3FCore_Ebf1 <- select(G3FCore_Ebf1, c("Category","Mean", "Sample"))
#Drd1
G3FCore_Drd1 <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G3F_NAcC_Drd1_Results.csv",header = TRUE)
G3FCore_Drd1$Sample <- "G3F"
G3FCore_Drd1$Category <- ifelse(G3FCore_Drd1$Mean >3,"Pos","Neg")
G3FCore_Drd1 <- select(G3FCore_Drd1, c("Category","Mean", "Sample"))
#Htr4
G3FCore_Htr4 <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G3F_NAcC_Htr4_Results.csv",header = TRUE)
G3FCore_Htr4$Category <- ifelse(G3FCore_Htr4$Max >30,"Pos","Neg")
G3FCore_Htr4$Sample <- "G3F"
G3FCore_Htr4 <- select(G3FCore_Htr4, c("Category","Mean", "Sample"))

#Categories
G3FCore_All_Data <- data.frame(Ebf1               = G3FCore_Ebf1$Mean,
                               Ebf1_Category      = G3FCore_Ebf1$Category,
                               Drd1          = G3FCore_Drd1$Mean,
                               Drd1_Category = G3FCore_Drd1$Category,
                               Htr4             = G3FCore_Htr4$Mean,
                               Htr4_Category    = G3FCore_Htr4$Category, 
                               Sex = "Female",
                               Category         = NA)

G3FCore_All_Data$Category <- ifelse(G3FCore_All_Data$Htr4_Category == "Pos" & G3FCore_All_Data$Drd1_Category == "Neg" & G3FCore_All_Data$Ebf1_Category == "Pos",
                                    "Drd1- Ebf1+ Htr4+",
                                    ifelse(G3FCore_All_Data$Htr4_Category == "Neg" & G3FCore_All_Data$Drd1_Category == "Pos" & G3FCore_All_Data$Ebf1_Category == "Pos",
                                           "Drd1+ Ebf1+ Htr4-",
                                           ifelse(G3FCore_All_Data$Htr4_Category == "Pos" & G3FCore_All_Data$Drd1_Category == "Pos" & G3FCore_All_Data$Ebf1_Category == "Pos",
                                                  "Drd1+ Ebf1+ Htr4+",
                                                  ifelse(G3FCore_All_Data$Htr4_Category == "Neg" & G3FCore_All_Data$Drd1_Category == "Pos" & G3FCore_All_Data$Ebf1_Category == "Neg",
                                                         "Drd1+ Ebf1- Htr4-",
                                                         ifelse(G3FCore_All_Data$Htr4_Category == "Pos" & G3FCore_All_Data$Drd1_Category == "Neg" & G3FCore_All_Data$Ebf1_Category == "Neg",
                                                                "Drd1- Ebf1- Htr4+",
                                                                ifelse(G3FCore_All_Data$Htr4_Category == "Neg" & G3FCore_All_Data$Drd1_Category == "Neg" & G3FCore_All_Data$Ebf1_Category == "Pos",
                                                                       "Drd1- Ebf1+ Htr4-",
                                                                       ifelse(G3FCore_All_Data$Htr4_Category == "Pos" & G3FCore_All_Data$Drd1_Category == "Pos" & G3FCore_All_Data$Ebf1_Category == "Neg",
                                                                              "Drd1+ Ebf1- Htr4+", "Other")))))))

G3FCore_All_Data %>% count(Category)

####G4M####
G4MCore_ROIs       <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G4M_NAcC_Dapi_Results.csv",header = TRUE)
G4MCore_ROIs$Index <- 1:nrow(G4MCore_ROIs)
#Ebf1
G4MCore_Ebf1 <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G4M_NAcC_Ebf1_Results.csv",header = TRUE)
G4MCore_Ebf1$Sample <- "G4M"
G4MCore_Ebf1$Category <- ifelse(G4MCore_Ebf1$Max >25,"Pos","Neg")
G4MCore_Ebf1 <- select(G4MCore_Ebf1, c("Category","Mean", "Sample"))
#Drd1
G4MCore_Drd1 <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G4M_NAcC_Drd1_Results.csv",header = TRUE)
G4MCore_Drd1$Sample <- "G4M"
G4MCore_Drd1$Category <- ifelse(G4MCore_Drd1$Mean >4,"Pos","Neg")
G4MCore_Drd1 <- select(G4MCore_Drd1, c("Category","Mean", "Sample"))
#Htr4
G4MCore_Htr4 <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G4M_NAcC_Htr4_Results.csv",header = TRUE)
G4MCore_Htr4$Category <- ifelse(G4MCore_Htr4$Max >30,"Pos","Neg")
G4MCore_Htr4$Sample <- "G4M"
G4MCore_Htr4 <- select(G4MCore_Htr4, c("Category","Mean", "Sample"))

#Categories
G4MCore_All_Data <- data.frame(Ebf1               = G4MCore_Ebf1$Mean,
                               Ebf1_Category      = G4MCore_Ebf1$Category,
                               Drd1          = G4MCore_Drd1$Mean,
                               Drd1_Category = G4MCore_Drd1$Category,
                               Htr4             = G4MCore_Htr4$Mean,
                               Htr4_Category    = G4MCore_Htr4$Category, 
                               Sex = "Female",
                               Category         = NA)

G4MCore_All_Data$Category <- ifelse(G4MCore_All_Data$Htr4_Category == "Pos" & G4MCore_All_Data$Drd1_Category == "Neg" & G4MCore_All_Data$Ebf1_Category == "Pos",
                                    "Drd1- Ebf1+ Htr4+",
                                    ifelse(G4MCore_All_Data$Htr4_Category == "Neg" & G4MCore_All_Data$Drd1_Category == "Pos" & G4MCore_All_Data$Ebf1_Category == "Pos",
                                           "Drd1+ Ebf1+ Htr4-",
                                           ifelse(G4MCore_All_Data$Htr4_Category == "Pos" & G4MCore_All_Data$Drd1_Category == "Pos" & G4MCore_All_Data$Ebf1_Category == "Pos",
                                                  "Drd1+ Ebf1+ Htr4+",
                                                  ifelse(G4MCore_All_Data$Htr4_Category == "Neg" & G4MCore_All_Data$Drd1_Category == "Pos" & G4MCore_All_Data$Ebf1_Category == "Neg",
                                                         "Drd1+ Ebf1- Htr4-",
                                                         ifelse(G4MCore_All_Data$Htr4_Category == "Pos" & G4MCore_All_Data$Drd1_Category == "Neg" & G4MCore_All_Data$Ebf1_Category == "Neg",
                                                                "Drd1- Ebf1- Htr4+",
                                                                ifelse(G4MCore_All_Data$Htr4_Category == "Neg" & G4MCore_All_Data$Drd1_Category == "Neg" & G4MCore_All_Data$Ebf1_Category == "Pos",
                                                                       "Drd1- Ebf1+ Htr4-",
                                                                       ifelse(G4MCore_All_Data$Htr4_Category == "Pos" & G4MCore_All_Data$Drd1_Category == "Pos" & G4MCore_All_Data$Ebf1_Category == "Neg",
                                                                              "Drd1+ Ebf1- Htr4+", "Other")))))))

G4MCore_All_Data %>% count(Category)

####G5M####
G5MCore_ROIs       <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G5M_NAcC_Dapi_Results.csv",header = TRUE)
G5MCore_ROIs$Index <- 1:nrow(G5MCore_ROIs)
#Ebf1
G5MCore_Ebf1 <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G5M_NAcC_Ebf1_Results.csv",header = TRUE)
G5MCore_Ebf1$Sample <- "G5M"
G5MCore_Ebf1$Category <- ifelse(G5MCore_Ebf1$Max >25,"Pos","Neg")
G5MCore_Ebf1 <- select(G5MCore_Ebf1, c("Category","Mean", "Sample"))
#Drd1
G5MCore_Drd1 <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G5M_NAcC_Drd1_Results.csv",header = TRUE)
G5MCore_Drd1$Sample <- "G5M"
G5MCore_Drd1$Category <- ifelse(G5MCore_Drd1$Mean >5,"Pos","Neg")
G5MCore_Drd1 <- select(G5MCore_Drd1, c("Category","Mean", "Sample"))
#Htr4
G5MCore_Htr4 <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G5M_NAcC_Htr4_Results.csv",header = TRUE)
G5MCore_Htr4$Category <- ifelse(G5MCore_Htr4$Max >50,"Pos","Neg")
G5MCore_Htr4$Sample <- "G5M"
G5MCore_Htr4 <- select(G5MCore_Htr4, c("Category","Mean", "Sample"))

#Categories
G5MCore_All_Data <- data.frame(Ebf1               = G5MCore_Ebf1$Mean,
                               Ebf1_Category      = G5MCore_Ebf1$Category,
                               Drd1          = G5MCore_Drd1$Mean,
                               Drd1_Category = G5MCore_Drd1$Category,
                               Htr4             = G5MCore_Htr4$Mean,
                               Htr4_Category    = G5MCore_Htr4$Category, 
                               Sex = "Female",
                               Category         = NA)

G5MCore_All_Data$Category <- ifelse(G5MCore_All_Data$Htr4_Category == "Pos" & G5MCore_All_Data$Drd1_Category == "Neg" & G5MCore_All_Data$Ebf1_Category == "Pos",
                                    "Drd1- Ebf1+ Htr4+",
                                    ifelse(G5MCore_All_Data$Htr4_Category == "Neg" & G5MCore_All_Data$Drd1_Category == "Pos" & G5MCore_All_Data$Ebf1_Category == "Pos",
                                           "Drd1+ Ebf1+ Htr4-",
                                           ifelse(G5MCore_All_Data$Htr4_Category == "Pos" & G5MCore_All_Data$Drd1_Category == "Pos" & G5MCore_All_Data$Ebf1_Category == "Pos",
                                                  "Drd1+ Ebf1+ Htr4+",
                                                  ifelse(G5MCore_All_Data$Htr4_Category == "Neg" & G5MCore_All_Data$Drd1_Category == "Pos" & G5MCore_All_Data$Ebf1_Category == "Neg",
                                                         "Drd1+ Ebf1- Htr4-",
                                                         ifelse(G5MCore_All_Data$Htr4_Category == "Pos" & G5MCore_All_Data$Drd1_Category == "Neg" & G5MCore_All_Data$Ebf1_Category == "Neg",
                                                                "Drd1- Ebf1- Htr4+",
                                                                ifelse(G5MCore_All_Data$Htr4_Category == "Neg" & G5MCore_All_Data$Drd1_Category == "Neg" & G5MCore_All_Data$Ebf1_Category == "Pos",
                                                                       "Drd1- Ebf1+ Htr4-",
                                                                       ifelse(G5MCore_All_Data$Htr4_Category == "Pos" & G5MCore_All_Data$Drd1_Category == "Pos" & G5MCore_All_Data$Ebf1_Category == "Neg",
                                                                              "Drd1+ Ebf1- Htr4+", "Other")))))))

G5MCore_All_Data %>% count(Category)

####G6M####
G6MCore_ROIs       <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G6M_NAcC_Dapi_Results.csv",header = TRUE)
G6MCore_ROIs$Index <- 1:nrow(G6MCore_ROIs)
#Ebf1
G6MCore_Ebf1 <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G6M_NAcC_Ebf1_Results.csv",header = TRUE)
G6MCore_Ebf1$Sample <- "G6M"
G6MCore_Ebf1$Category <- ifelse(G6MCore_Ebf1$Max >25,"Pos","Neg")
G6MCore_Ebf1 <- select(G6MCore_Ebf1, c("Category","Mean", "Sample"))
#Drd1
G6MCore_Drd1 <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G6M_NAcC_Drd1_Results.csv",header = TRUE)
G6MCore_Drd1$Sample <- "G6M"
G6MCore_Drd1$Category <- ifelse(G6MCore_Drd1$Mean >4,"Pos","Neg")
G6MCore_Drd1 <- select(G6MCore_Drd1, c("Category","Mean", "Sample"))
#Htr4
G6MCore_Htr4 <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G6M_NAcC_Htr4_Results.csv",header = TRUE)
G6MCore_Htr4$Category <- ifelse(G6MCore_Htr4$Max >35,"Pos","Neg")
G6MCore_Htr4$Sample <- "G6M"
G6MCore_Htr4 <- select(G6MCore_Htr4, c("Category","Mean", "Sample"))

#Categories
G6MCore_All_Data <- data.frame(Ebf1               = G6MCore_Ebf1$Mean,
                               Ebf1_Category      = G6MCore_Ebf1$Category,
                               Drd1          = G6MCore_Drd1$Mean,
                               Drd1_Category = G6MCore_Drd1$Category,
                               Htr4             = G6MCore_Htr4$Mean,
                               Htr4_Category    = G6MCore_Htr4$Category, 
                               Sex = "Female",
                               Category         = NA)

G6MCore_All_Data$Category <- ifelse(G6MCore_All_Data$Htr4_Category == "Pos" & G6MCore_All_Data$Drd1_Category == "Neg" & G6MCore_All_Data$Ebf1_Category == "Pos",
                                    "Drd1- Ebf1+ Htr4+",
                                    ifelse(G6MCore_All_Data$Htr4_Category == "Neg" & G6MCore_All_Data$Drd1_Category == "Pos" & G6MCore_All_Data$Ebf1_Category == "Pos",
                                           "Drd1+ Ebf1+ Htr4-",
                                           ifelse(G6MCore_All_Data$Htr4_Category == "Pos" & G6MCore_All_Data$Drd1_Category == "Pos" & G6MCore_All_Data$Ebf1_Category == "Pos",
                                                  "Drd1+ Ebf1+ Htr4+",
                                                  ifelse(G6MCore_All_Data$Htr4_Category == "Neg" & G6MCore_All_Data$Drd1_Category == "Pos" & G6MCore_All_Data$Ebf1_Category == "Neg",
                                                         "Drd1+ Ebf1- Htr4-",
                                                         ifelse(G6MCore_All_Data$Htr4_Category == "Pos" & G6MCore_All_Data$Drd1_Category == "Neg" & G6MCore_All_Data$Ebf1_Category == "Neg",
                                                                "Drd1- Ebf1- Htr4+",
                                                                ifelse(G6MCore_All_Data$Htr4_Category == "Neg" & G6MCore_All_Data$Drd1_Category == "Neg" & G6MCore_All_Data$Ebf1_Category == "Pos",
                                                                       "Drd1- Ebf1+ Htr4-",
                                                                       ifelse(G6MCore_All_Data$Htr4_Category == "Pos" & G6MCore_All_Data$Drd1_Category == "Pos" & G6MCore_All_Data$Ebf1_Category == "Neg",
                                                                              "Drd1+ Ebf1- Htr4+", "Other")))))))

G6MCore_All_Data %>% count(Category)




#NAc SHELL

####G1F####
#read in DAPI rois
G1FShell_ROIs       <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G1F_NAcS_Dapi_Results.csv",header = TRUE)
G1FShell_ROIs$Index <- 1:nrow(G1FShell_ROIs)
#Ebf1
G1FShell_Ebf1 <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G1F_NAcS_Ebf1_Results.csv",header = TRUE)
G1FShell_Ebf1$Sample <- "G1F"
G1FShell_Ebf1$Category <- ifelse(G1FShell_Ebf1$Max >25,"Pos","Neg")
G1FShell_Ebf1 <- select(G1FShell_Ebf1, c("Category","Mean", "Sample"))
#Drd1
G1FShell_Drd1 <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G1F_NAcS_Drd1_Results.csv",header = TRUE)
G1FShell_Drd1$Sample <- "G1F"
G1FShell_Drd1$Category <- ifelse(G1FShell_Drd1$Mean >6,"Pos","Neg")
G1FShell_Drd1 <- select(G1FShell_Drd1, c("Category","Mean", "Sample"))
#Htr4
G1FShell_Htr4 <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G1F_NAcS_Htr4_Results.csv",header = TRUE)
G1FShell_Htr4$Category <- ifelse(G1FShell_Htr4$Max >35,"Pos","Neg")
G1FShell_Htr4$Sample <- "G1F"
G1FShell_Htr4 <- select(G1FShell_Htr4, c("Category","Mean", "Sample"))

#Categories
G1FShell_All_Data <- data.frame(Ebf1               = G1FShell_Ebf1$Mean,
                                Ebf1_Category      = G1FShell_Ebf1$Category,
                                Drd1          = G1FShell_Drd1$Mean,
                                Drd1_Category = G1FShell_Drd1$Category,
                                Htr4             = G1FShell_Htr4$Mean,
                                Htr4_Category    = G1FShell_Htr4$Category, 
                                Sex = "Female",
                                Category         = NA)

G1FShell_All_Data$Category <- ifelse(G1FShell_All_Data$Htr4_Category == "Pos" & G1FShell_All_Data$Drd1_Category == "Neg" & G1FShell_All_Data$Ebf1_Category == "Pos",
                                     "Drd1- Ebf1+ Htr4+",
                                     ifelse(G1FShell_All_Data$Htr4_Category == "Neg" & G1FShell_All_Data$Drd1_Category == "Pos" & G1FShell_All_Data$Ebf1_Category == "Pos",
                                            "Drd1+ Ebf1+ Htr4-",
                                            ifelse(G1FShell_All_Data$Htr4_Category == "Pos" & G1FShell_All_Data$Drd1_Category == "Pos" & G1FShell_All_Data$Ebf1_Category == "Pos",
                                                   "Drd1+ Ebf1+ Htr4+",
                                                   ifelse(G1FShell_All_Data$Htr4_Category == "Neg" & G1FShell_All_Data$Drd1_Category == "Pos" & G1FShell_All_Data$Ebf1_Category == "Neg",
                                                          "Drd1+ Ebf1- Htr4-",
                                                          ifelse(G1FShell_All_Data$Htr4_Category == "Pos" & G1FShell_All_Data$Drd1_Category == "Neg" & G1FShell_All_Data$Ebf1_Category == "Neg",
                                                                 "Drd1- Ebf1- Htr4+",
                                                                 ifelse(G1FShell_All_Data$Htr4_Category == "Neg" & G1FShell_All_Data$Drd1_Category == "Neg" & G1FShell_All_Data$Ebf1_Category == "Pos",
                                                                        "Drd1- Ebf1+ Htr4-",
                                                                        ifelse(G1FShell_All_Data$Htr4_Category == "Pos" & G1FShell_All_Data$Drd1_Category == "Pos" & G1FShell_All_Data$Ebf1_Category == "Neg",
                                                                               "Drd1+ Ebf1- Htr4+", "Other")))))))

G1FShell_All_Data %>% count(Category)

#G2F
#Read in the DAPI ROI positions
G2FShell_ROIs       <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G2F_NAcS_Dapi_Results.csv",header = TRUE)
G2FShell_ROIs$Index <- 1:nrow(G2FShell_ROIs)
#Ebf1
G2FShell_Ebf1 <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G2F_NAcS_Ebf1_Results.csv",header = TRUE)
G2FShell_Ebf1$Sample <- "G2F"
G2FShell_Ebf1$Category <- ifelse(G2FShell_Ebf1$Max >25,"Pos","Neg")
G2FShell_Ebf1 <- select(G2FShell_Ebf1, c("Category","Mean", "Sample"))
#Drd1
G2FShell_Drd1 <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G2F_NAcS_Drd1_Results.csv",header = TRUE)
G2FShell_Drd1$Sample <- "G2F"
G2FShell_Drd1$Category <- ifelse(G2FShell_Drd1$Mean >5,"Pos","Neg")
G2FShell_Drd1 <- select(G2FShell_Drd1, c("Category","Mean", "Sample"))
#Htr4
G2FShell_Htr4 <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G2F_NAcS_Htr4_Results.csv",header = TRUE)
G2FShell_Htr4$Category <- ifelse(G2FShell_Htr4$Max >35,"Pos","Neg")
G2FShell_Htr4$Sample <- "G2F"
G2FShell_Htr4 <- select(G2FShell_Htr4, c("Category","Mean", "Sample"))

#Categories
G2FShell_All_Data <- data.frame(Ebf1               = G2FShell_Ebf1$Mean,
                                Ebf1_Category      = G2FShell_Ebf1$Category,
                                Drd1          = G2FShell_Drd1$Mean,
                                Drd1_Category = G2FShell_Drd1$Category,
                                Htr4             = G2FShell_Htr4$Mean,
                                Htr4_Category    = G2FShell_Htr4$Category, 
                                Sex = "Female",
                                Category         = NA)

G2FShell_All_Data$Category <- ifelse(G2FShell_All_Data$Htr4_Category == "Pos" & G2FShell_All_Data$Drd1_Category == "Neg" & G2FShell_All_Data$Ebf1_Category == "Pos",
                                     "Drd1- Ebf1+ Htr4+",
                                     ifelse(G2FShell_All_Data$Htr4_Category == "Neg" & G2FShell_All_Data$Drd1_Category == "Pos" & G2FShell_All_Data$Ebf1_Category == "Pos",
                                            "Drd1+ Ebf1+ Htr4-",
                                            ifelse(G2FShell_All_Data$Htr4_Category == "Pos" & G2FShell_All_Data$Drd1_Category == "Pos" & G2FShell_All_Data$Ebf1_Category == "Pos",
                                                   "Drd1+ Ebf1+ Htr4+",
                                                   ifelse(G2FShell_All_Data$Htr4_Category == "Neg" & G2FShell_All_Data$Drd1_Category == "Pos" & G2FShell_All_Data$Ebf1_Category == "Neg",
                                                          "Drd1+ Ebf1- Htr4-",
                                                          ifelse(G2FShell_All_Data$Htr4_Category == "Pos" & G2FShell_All_Data$Drd1_Category == "Neg" & G2FShell_All_Data$Ebf1_Category == "Neg",
                                                                 "Drd1- Ebf1- Htr4+",
                                                                 ifelse(G2FShell_All_Data$Htr4_Category == "Neg" & G2FShell_All_Data$Drd1_Category == "Neg" & G2FShell_All_Data$Ebf1_Category == "Pos",
                                                                        "Drd1- Ebf1+ Htr4-",
                                                                        ifelse(G2FShell_All_Data$Htr4_Category == "Pos" & G2FShell_All_Data$Drd1_Category == "Pos" & G2FShell_All_Data$Ebf1_Category == "Neg",
                                                                               "Drd1+ Ebf1- Htr4+", "Other")))))))

G2FShell_All_Data %>% count(Category)


#G3F
#Read in the DAPI ROI positions
G3FShell_ROIs       <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G3F_NAcS_Dapi_Results.csv",header = TRUE)
G3FShell_ROIs$Index <- 1:nrow(G3FShell_ROIs)
#Ebf1
G3FShell_Ebf1 <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G3F_NAcS_Ebf1_Results.csv",header = TRUE)
G3FShell_Ebf1$Sample <- "G3F"
G3FShell_Ebf1$Category <- ifelse(G3FShell_Ebf1$Max >25,"Pos","Neg")
G3FShell_Ebf1 <- select(G3FShell_Ebf1, c("Category","Mean", "Sample"))
#Drd1
G3FShell_Drd1 <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G3F_NAcS_Drd1_Results.csv",header = TRUE)
G3FShell_Drd1$Sample <- "G3F"
G3FShell_Drd1$Category <- ifelse(G3FShell_Drd1$Mean >4,"Pos","Neg")
G3FShell_Drd1 <- select(G3FShell_Drd1, c("Category","Mean", "Sample"))
#Htr4
G3FShell_Htr4 <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G3F_NAcS_Htr4_Results.csv",header = TRUE)
G3FShell_Htr4$Category <- ifelse(G3FShell_Htr4$Max >35,"Pos","Neg")
G3FShell_Htr4$Sample <- "G3F"
G3FShell_Htr4 <- select(G3FShell_Htr4, c("Category","Mean", "Sample"))

#Categories
G3FShell_All_Data <- data.frame(Ebf1               = G3FShell_Ebf1$Mean,
                                Ebf1_Category      = G3FShell_Ebf1$Category,
                                Drd1          = G3FShell_Drd1$Mean,
                                Drd1_Category = G3FShell_Drd1$Category,
                                Htr4             = G3FShell_Htr4$Mean,
                                Htr4_Category    = G3FShell_Htr4$Category, 
                                Sex = "Female",
                                Category         = NA)

G3FShell_All_Data$Category <- ifelse(G3FShell_All_Data$Htr4_Category == "Pos" & G3FShell_All_Data$Drd1_Category == "Neg" & G3FShell_All_Data$Ebf1_Category == "Pos",
                                     "Drd1- Ebf1+ Htr4+",
                                     ifelse(G3FShell_All_Data$Htr4_Category == "Neg" & G3FShell_All_Data$Drd1_Category == "Pos" & G3FShell_All_Data$Ebf1_Category == "Pos",
                                            "Drd1+ Ebf1+ Htr4-",
                                            ifelse(G3FShell_All_Data$Htr4_Category == "Pos" & G3FShell_All_Data$Drd1_Category == "Pos" & G3FShell_All_Data$Ebf1_Category == "Pos",
                                                   "Drd1+ Ebf1+ Htr4+",
                                                   ifelse(G3FShell_All_Data$Htr4_Category == "Neg" & G3FShell_All_Data$Drd1_Category == "Pos" & G3FShell_All_Data$Ebf1_Category == "Neg",
                                                          "Drd1+ Ebf1- Htr4-",
                                                          ifelse(G3FShell_All_Data$Htr4_Category == "Pos" & G3FShell_All_Data$Drd1_Category == "Neg" & G3FShell_All_Data$Ebf1_Category == "Neg",
                                                                 "Drd1- Ebf1- Htr4+",
                                                                 ifelse(G3FShell_All_Data$Htr4_Category == "Neg" & G3FShell_All_Data$Drd1_Category == "Neg" & G3FShell_All_Data$Ebf1_Category == "Pos",
                                                                        "Drd1- Ebf1+ Htr4-",
                                                                        ifelse(G3FShell_All_Data$Htr4_Category == "Pos" & G3FShell_All_Data$Drd1_Category == "Pos" & G3FShell_All_Data$Ebf1_Category == "Neg",
                                                                               "Drd1+ Ebf1- Htr4+", "Other")))))))

G3FShell_All_Data %>% count(Category)


####G4M####
#read in DAPI ROIs
G4MShell_ROIs       <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G4M_NAcS_Dapi_Results.csv",header = TRUE)
G4MShell_ROIs$Index <- 1:nrow(G4MShell_ROIs)
#Ebf1
G4MShell_Ebf1 <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G4M_NAcS_Ebf1_Results.csv",header = TRUE)
G4MShell_Ebf1$Sample <- "G4M"
G4MShell_Ebf1$Category <- ifelse(G4MShell_Ebf1$Max >25,"Pos","Neg")
G4MShell_Ebf1 <- select(G4MShell_Ebf1, c("Category","Mean", "Sample"))
#Drd1
G4MShell_Drd1 <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G4M_NAcS_Drd1_Results.csv",header = TRUE)
G4MShell_Drd1$Sample <- "G4M"
G4MShell_Drd1$Category <- ifelse(G4MShell_Drd1$Mean >4,"Pos","Neg")
G4MShell_Drd1 <- select(G4MShell_Drd1, c("Category","Mean", "Sample"))
#Htr4
G4MShell_Htr4 <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G4M_NAcS_Htr4_Results.csv",header = TRUE)
G4MShell_Htr4$Category <- ifelse(G4MShell_Htr4$Max >35,"Pos","Neg")
G4MShell_Htr4$Sample <- "G4M"
G4MShell_Htr4 <- select(G4MShell_Htr4, c("Category","Mean", "Sample"))

#Categories
G4MShell_All_Data <- data.frame(Ebf1               = G4MShell_Ebf1$Mean,
                                Ebf1_Category      = G4MShell_Ebf1$Category,
                                Drd1          = G4MShell_Drd1$Mean,
                                Drd1_Category = G4MShell_Drd1$Category,
                                Htr4             = G4MShell_Htr4$Mean,
                                Htr4_Category    = G4MShell_Htr4$Category, 
                                Sex = "Female",
                                Category         = NA)

G4MShell_All_Data$Category <- ifelse(G4MShell_All_Data$Htr4_Category == "Pos" & G4MShell_All_Data$Drd1_Category == "Neg" & G4MShell_All_Data$Ebf1_Category == "Pos",
                                     "Drd1- Ebf1+ Htr4+",
                                     ifelse(G4MShell_All_Data$Htr4_Category == "Neg" & G4MShell_All_Data$Drd1_Category == "Pos" & G4MShell_All_Data$Ebf1_Category == "Pos",
                                            "Drd1+ Ebf1+ Htr4-",
                                            ifelse(G4MShell_All_Data$Htr4_Category == "Pos" & G4MShell_All_Data$Drd1_Category == "Pos" & G4MShell_All_Data$Ebf1_Category == "Pos",
                                                   "Drd1+ Ebf1+ Htr4+",
                                                   ifelse(G4MShell_All_Data$Htr4_Category == "Neg" & G4MShell_All_Data$Drd1_Category == "Pos" & G4MShell_All_Data$Ebf1_Category == "Neg",
                                                          "Drd1+ Ebf1- Htr4-",
                                                          ifelse(G4MShell_All_Data$Htr4_Category == "Pos" & G4MShell_All_Data$Drd1_Category == "Neg" & G4MShell_All_Data$Ebf1_Category == "Neg",
                                                                 "Drd1- Ebf1- Htr4+",
                                                                 ifelse(G4MShell_All_Data$Htr4_Category == "Neg" & G4MShell_All_Data$Drd1_Category == "Neg" & G4MShell_All_Data$Ebf1_Category == "Pos",
                                                                        "Drd1- Ebf1+ Htr4-",
                                                                        ifelse(G4MShell_All_Data$Htr4_Category == "Pos" & G4MShell_All_Data$Drd1_Category == "Pos" & G4MShell_All_Data$Ebf1_Category == "Neg",
                                                                               "Drd1+ Ebf1- Htr4+", "Other")))))))

G4MShell_All_Data %>% count(Category)


####G5M####
#Read in DAPI ROIs
G5MShell_ROIs       <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G5M_NAcS_Dapi_Results.csv",header = TRUE)
G5MShell_ROIs$Index <- 1:nrow(G5MShell_ROIs)
#Ebf1
G5MShell_Ebf1 <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G5M_NAcS_Ebf1_Results.csv",header = TRUE)
G5MShell_Ebf1$Sample <- "G5M"
G5MShell_Ebf1$Category <- ifelse(G5MShell_Ebf1$Max >25,"Pos","Neg")
G5MShell_Ebf1 <- select(G5MShell_Ebf1, c("Category","Mean", "Sample"))
#Drd1
G5MShell_Drd1 <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G5M_NAcS_Drd1_Results.csv",header = TRUE)
G5MShell_Drd1$Sample <- "G5M"
G5MShell_Drd1$Category <- ifelse(G5MShell_Drd1$Mean >5,"Pos","Neg")
G5MShell_Drd1 <- select(G5MShell_Drd1, c("Category","Mean", "Sample"))
#Htr4
G5MShell_Htr4 <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G5M_NAcS_Htr4_Results.csv",header = TRUE)
G5MShell_Htr4$Category <- ifelse(G5MShell_Htr4$Max >30,"Pos","Neg")
G5MShell_Htr4$Sample <- "G5M"
G5MShell_Htr4 <- select(G5MShell_Htr4, c("Category","Mean", "Sample"))

#Categories
G5MShell_All_Data <- data.frame(Ebf1               = G5MShell_Ebf1$Mean,
                                Ebf1_Category      = G5MShell_Ebf1$Category,
                                Drd1          = G5MShell_Drd1$Mean,
                                Drd1_Category = G5MShell_Drd1$Category,
                                Htr4             = G5MShell_Htr4$Mean,
                                Htr4_Category    = G5MShell_Htr4$Category, 
                                Sex = "Female",
                                Category         = NA)

G5MShell_All_Data$Category <- ifelse(G5MShell_All_Data$Htr4_Category == "Pos" & G5MShell_All_Data$Drd1_Category == "Neg" & G5MShell_All_Data$Ebf1_Category == "Pos",
                                     "Drd1- Ebf1+ Htr4+",
                                     ifelse(G5MShell_All_Data$Htr4_Category == "Neg" & G5MShell_All_Data$Drd1_Category == "Pos" & G5MShell_All_Data$Ebf1_Category == "Pos",
                                            "Drd1+ Ebf1+ Htr4-",
                                            ifelse(G5MShell_All_Data$Htr4_Category == "Pos" & G5MShell_All_Data$Drd1_Category == "Pos" & G5MShell_All_Data$Ebf1_Category == "Pos",
                                                   "Drd1+ Ebf1+ Htr4+",
                                                   ifelse(G5MShell_All_Data$Htr4_Category == "Neg" & G5MShell_All_Data$Drd1_Category == "Pos" & G5MShell_All_Data$Ebf1_Category == "Neg",
                                                          "Drd1+ Ebf1- Htr4-",
                                                          ifelse(G5MShell_All_Data$Htr4_Category == "Pos" & G5MShell_All_Data$Drd1_Category == "Neg" & G5MShell_All_Data$Ebf1_Category == "Neg",
                                                                 "Drd1- Ebf1- Htr4+",
                                                                 ifelse(G5MShell_All_Data$Htr4_Category == "Neg" & G5MShell_All_Data$Drd1_Category == "Neg" & G5MShell_All_Data$Ebf1_Category == "Pos",
                                                                        "Drd1- Ebf1+ Htr4-",
                                                                        ifelse(G5MShell_All_Data$Htr4_Category == "Pos" & G5MShell_All_Data$Drd1_Category == "Pos" & G5MShell_All_Data$Ebf1_Category == "Neg",
                                                                               "Drd1+ Ebf1- Htr4+", "Other")))))))

G5MShell_All_Data %>% count(Category)

####G6M####
G6MShell_ROIs       <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G6M_NAcS_Dapi_Results.csv",header = TRUE)
G6MShell_ROIs$Index <- 1:nrow(G6MShell_ROIs)
#Ebf1
G6MShell_Ebf1 <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G6M_NAcS_Ebf1_Results.csv",header = TRUE)
G6MShell_Ebf1$Sample <- "G6M"
G6MShell_Ebf1$Category <- ifelse(G6MShell_Ebf1$Max >35,"Pos","Neg")
G6MShell_Ebf1 <- select(G6MShell_Ebf1, c("Category","Mean", "Sample"))
#Drd1
G6MShell_Drd1 <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G6M_NAcS_Drd1_Results.csv",header = TRUE)
G6MShell_Drd1$Sample <- "G6M"
G6MShell_Drd1$Category <- ifelse(G6MShell_Drd1$Mean >7,"Pos","Neg")
G6MShell_Drd1 <- select(G6MShell_Drd1, c("Category","Mean", "Sample"))
#Htr4
G6MShell_Htr4 <- read.csv(file = "/Volumes/Day-Lab$/Dalton/RNAscope/df000005_RP.JT.NAc/Figure/Images for analysis/G6M_NAcS_Htr4_Results.csv",header = TRUE)
G6MShell_Htr4$Category <- ifelse(G6MShell_Htr4$Max >25,"Pos","Neg")
G6MShell_Htr4$Sample <- "G6M"
G6MShell_Htr4 <- select(G6MShell_Htr4, c("Category","Mean", "Sample"))

#Categories
G6MShell_All_Data <- data.frame(Ebf1               = G6MShell_Ebf1$Mean,
                                Ebf1_Category      = G6MShell_Ebf1$Category,
                                Drd1          = G6MShell_Drd1$Mean,
                                Drd1_Category = G6MShell_Drd1$Category,
                                Htr4             = G6MShell_Htr4$Mean,
                                Htr4_Category    = G6MShell_Htr4$Category, 
                                Sex = "Female",
                                Category         = NA)

G6MShell_All_Data$Category <- ifelse(G6MShell_All_Data$Htr4_Category == "Pos" & G6MShell_All_Data$Drd1_Category == "Neg" & G6MShell_All_Data$Ebf1_Category == "Pos",
                                     "Drd1- Ebf1+ Htr4+",
                                     ifelse(G6MShell_All_Data$Htr4_Category == "Neg" & G6MShell_All_Data$Drd1_Category == "Pos" & G6MShell_All_Data$Ebf1_Category == "Pos",
                                            "Drd1+ Ebf1+ Htr4-",
                                            ifelse(G6MShell_All_Data$Htr4_Category == "Pos" & G6MShell_All_Data$Drd1_Category == "Pos" & G6MShell_All_Data$Ebf1_Category == "Pos",
                                                   "Drd1+ Ebf1+ Htr4+",
                                                   ifelse(G6MShell_All_Data$Htr4_Category == "Neg" & G6MShell_All_Data$Drd1_Category == "Pos" & G6MShell_All_Data$Ebf1_Category == "Neg",
                                                          "Drd1+ Ebf1- Htr4-",
                                                          ifelse(G6MShell_All_Data$Htr4_Category == "Pos" & G6MShell_All_Data$Drd1_Category == "Neg" & G6MShell_All_Data$Ebf1_Category == "Neg",
                                                                 "Drd1- Ebf1- Htr4+",
                                                                 ifelse(G6MShell_All_Data$Htr4_Category == "Neg" & G6MShell_All_Data$Drd1_Category == "Neg" & G6MShell_All_Data$Ebf1_Category == "Pos",
                                                                        "Drd1- Ebf1+ Htr4-",
                                                                        ifelse(G6MShell_All_Data$Htr4_Category == "Pos" & G6MShell_All_Data$Drd1_Category == "Pos" & G6MShell_All_Data$Ebf1_Category == "Neg",
                                                                               "Drd1+ Ebf1- Htr4+", "Other")))))))

G6MShell_All_Data %>% count(Category)


#20230103 Session Info
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
