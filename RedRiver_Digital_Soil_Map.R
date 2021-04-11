options(shiny.maxRequestSize = 300*1024^5) 

library(raster)
library(ggplot2)
library(dplyr)
library(rlang)
library(plotly)
library(hrbrthemes)
library(pracma)
library(ptw)
library(shiny)
library(shinydashboard)
library(leaflet)
# File import
rs1 <- list.files("D:/UGA-project/Output/Soil_Class/FAW_20", pattern = ".tif$", full.names = T)
#Stack rasters
rs1 <- stack(rs1)
#Reprojection
projection(rs1) <- crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
# Fix threshold for moist and dry FAW values
rc <- function(x){ifelse(x <=0.269,0, ifelse(x > 0.269,1, NA))}
# Assign threshold values to each pixel
sc <- calc(rs1, fun=rc)
#Function to calculate the maximum sum of consecutive moist days <90
f1 <- function(x) {
  max(sum(x[with(rle(x == 1), rep(lengths * values < 90, lengths))]))
}
#Maximum consecutive moist days at pixel level
max_cons_moist<-calc(sc, f1)

#sum moist days >=90
sum_moist1 <- calc(sc, function(x,na.rm) sum(x) >=90)
#sum moist days >90
sum_moist2 <- calc(sc, function(x,na.rm) sum(x) >90)
#sum moist days >183
sum_moist3 <- calc(sc, function(x,na.rm) sum(x) >183)
#sum moist days >183
sum_moist4 <- calc(sc, function(x,na.rm) sum(x) >=183)
#Function to calculate the maximum sum of consecutive dry days >183
f2 <- function(x) {
  max(sum(x[with(rle(x == 0), rep(lengths * values >183, lengths))]))
}
#Maximum consecutive dry days at pixel level
max_cons_dry<-calc(sc, f2)

#To find out the sum dry days >183 at pixel level
sum_dry2 <- calc(sc, function(x,na.rm) sum(x)>183)
#To find out the sum dry days >=90 at pixel level
sum_dry1 <- calc(sc, function(x,na.rm) sum(x)>=90)
# Import files
rs2 <- list.files("D:/UGA-project/PRISM_TMAX_4KM_DAILY_TIF", pattern = ".tif$", full.names = T)
# Stack files
rs2 <- stack(rs2)
#Define spatial reference
rs2_prj <- projectRaster(rs2, rs1)
#Mask to area of interest
rs2_prj <- mask(rs2_prj, rs1)
#Import files
rs3 <- list.files("D:/UGA-project/PRISM_TMIN_4KM_DAILY_TIF", pattern = ".tif$", full.names = T)
#Stack rasters
rs3 <- stack(rs3)
#Project rasters
rs3_prj <- projectRaster(rs3, rs1) 
#Mask to area of interest
rs3_prj <- mask(rs3_prj, rs1)

#Mean maximum and minimum temperatures
stm <- overlay(rs2_prj,rs3_prj, fun=function(x,y) {(((x+y)/2) + 1)})
#Filters stm layers more than 5 
hot.5.1 <- calc(stm, function(x,na.rm) x >5)
#Filters stm less than 6 
hot.6.1 <- calc(stm, function(x,na.rm) x <6)
#Filters stm less than and equals 6
hot.6.2 <- calc(stm, function(x,na.rm) x >=6)
#Filters stm more than 8
hot.8.1 <- calc(stm, function(x,na.rm) x >8)
#Filters stm less than 22
hot.22.1 <- calc(stm, function(x,na.rm) x <22)
#Filters stm more than and equals 22
hot.22.2 <- calc(stm, function(x,na.rm) x >=22)

#A raster in the stack is indexed using the double bracket and 
#the raster values are indexed using a single bracket

#Stacks maximum temperature values for summer (June, July, August) 
SUMT <- stack(rs2_prj[[153:215]])
#Mean summer temperature
SUMT <- calc(SUMT, mean)

#Stacks minimum temperature for winter (December, January February) 
WINT <- stack(rs3_prj[[c(1:59,335:366)]])
#Mean winter temperature
WINT <- calc(WINT, mean)

# Difference between summer and winter mean temperature
DIFF <- (SUMT - WINT)
DIFF1 <- calc(DIFF, function(x,na.rm) x >6)
DIFF2 <- calc(DIFF, function(x,na.rm) x >=6)

#Four months following summer solstice (June 21 to Oct 20)
SUMM_SL_4 <- stack(rs1[[c(174:296)]])
SUMM_SL_4 <- calc(SUMM_SL_4, fun =mean)

#Function to calculate the maximum sum of consecutive dry days 
#following four months of summer solstice <45
f3 <- function(x) {
  max(sum(x[with(rle(x == 0), rep(lengths * values <45, lengths))]))
}
#Maximum sum of consecutive dry days following four months of summer solstice <45 
max_cons_SUMMSL4_dry<-calc(SUMM_SL_4, f3)

#Stacks four months following winter solstice (Sept 22 to Jan 21)
WINT_SL_4 <- stack(rs1[[267:366, 1:21]])
#Assign values at pixel level
WINT_SL_4 <- calc(WINT_SL_4, fun = mean)

#Function to calculate the maximum sum of consecutive moist days 
#following four months of winter solstice >=45
f4 <- function(x) {
  max(sum(x[with(rle(x == 0), rep(lengths * values >=45, lengths))]))
}
#Maximum consecutive moist days for four months following winter solstice
max_cons_WINTSL4_moist<-calc(WINT_SL_4, f4)

#Soil class Identification
Aridic <-((sum_dry1*hot.5.1) + (max_cons_moist* hot.8.1))
Udic <- ((sum_moist1 * hot.22.1 * DIFF1) + (max_cons_SUMMSL4_dry))
Ustic <- ((sum_dry2 * hot.22.2 * DIFF1 * sum_moist3) + 
            (sum_dry2 * hot.22.1 * DIFF2) + (sum_moist4 * hot.5.1) +
            (max_cons_WINTSL4_moist * max_cons_SUMMSL4_dry))
Xeric <- ((max_cons_SUMMSL4_dry * max_cons_WINTSL4_moist) +  
            (sum_moist2 * hot.5.1) + (sum_moist2 * hot.8.1))
#Mean soil classes
ras <- mean(Aridic+Udic+Ustic+Xeric)
#plot(ras, col = terrain.colors(5))
#plot(ras,col=rainbow(4))
plot(ras, legend = FALSE, #col = terrain.colors(5),
     col = c("#0C2C84", "#41B6C4", "gold1", "red"), 
     axes = TRUE
)
legend("bottomleft", adj = c(0, 0.6),
       legend = c("Aridic", "Udic", "Ustic", "Xeric"),
       fill = c("#0C2C84", "#41B6C4", "gold1", "red"),
       #fill = terrain.colors(5),
       border = TRUE,
       bty = "o") # turn off legend border)
