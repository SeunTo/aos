#A
library(readxl)
library(dplyr)
library(forcats)
setwd("C:/Users/u1552345/Dropbox/Ecology Analysis/Ecology_Analysis")
PWB <- read_excel("popWB.xlsx")

F4Funnel <- Grp %>% 
  left_join(PWB, by = "Countries")

FinalFinal <- F4Funnel %>%
  left_join(Cat25, by = "Countries")

FinalFinal$Catastrophic10 <- as.numeric(FinalFinal$Catastrophic10)
FinalFinal$pop <- as.numeric(FinalFinal$pop)
FinalFinal$VHI_CHE <- as.numeric(FinalFinal$VHI_CHE)
FinalFinal$SHI_CHE <- as.numeric(FinalFinal$SHI_CHE)
FinalFinal$ExT_CHE <- as.numeric(FinalFinal$ExT_CHE)
FinalFinal$HDI <- as.numeric(FinalFinal$HDI)

FinalFun <- na.omit(F4Funnel)

FinalFun <- FinalFun[-c(21,35), ]
  

if(!require("berryFunctions")) install.packages("berryFunctions")
library(tidyverse)
dataFunnel <- FinalFinal %>%
  select(Countries, Catastrophic10, pop) %>%
  mutate(event = Catastrophic10/1000000, number = pop/1000000)


r <- dataFunnel$event
n <- dataFunnel$number
Name <- dataFunnel$Countries


funnelPlot(x=r, n=n, labels=Name,
           at=list(cex=0.5, col="black"), 
           ap=list(cex=0.7, col="black"), 
           am=list(lwd = 3, col="black"), 
           a2=list(lty=9, lwd = 3,  col="blue"), 
           a3=list(lty=5, lwd = 3,  col="red"), 
           ylab = "Percentage of people that spent > 10 OOP",
           xlab = "Population")


dev.copy2pdf(file="Figure5.pdf")
dev.off()


# B
data(medpar)
library(FunnelPlotR)
library(COUNT)
library(ggplot2)
FinalFinal$Countries<-factor(FinalFinal$Countries)

FinalFinal <- na.omit(FinalFinal)

mod <- glm(Catastrophic10 ~ GDP + HDI + CHE_GDP + ExT_CHE + GGHED_CHE + SHI_CHE + VHI_CHE + pop, family="poisson", data= FinalFinal)
summary(mod)
FinalFinal$prds <- predict(mod, type="response")
a<-funnel_plot(numerator=FinalFinal$Catastrophic10, denominator=FinalFinal$prds, group = FinalFinal$Countries, 
               title = 'Catastropic Expenditure Funnel plot', Poisson_limits = TRUE,
               OD_adjust = FALSE,label_outliers = TRUE, return_elements = "plot")
a



