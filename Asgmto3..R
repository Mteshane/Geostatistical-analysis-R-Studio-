# Install necessary packages: You only need to run this part once
install.packages(c("sp", "gstat", "lattice", "automap"))

# Load the required libraries
#_________________________________________________________
library(sp)
library(gstat)
library(lattice)
library(automap)
options(prompt="> ", continue="+ ", digits=3, width=70, show.signif.stars=T)
rm(list=ls()) 

# Install necessary packages: You only need to run this part once "sp", "gstat", "lattice", "automap"


# Load the required libraries
#_________________________________________________________

# get help for the gstat package:
#______________


# Check for example data sets in the sp package
#______________

data(package="gstat")
# Description of the data from the package documentation
data(meuse)
?meuse

# Load the data into R and explore the data

# The gstat package assumes that data are projected, 
# i.e. they should not be provided as latitude/longitude # and have to tell it which vectors are the coordinates

# The command below specifies the corrdinates and changes the dataframe 
# to a "SpatialPointDataFrame", which is necessary for variogram modelling

coordinates(meuse) <- ~ x + y    
proj4string(meuse) <- CRS("+init=epsg:28992")

# Explore the projected data
class(meuse)
str(meuse)




# Read the interpolation grid data (set of regularly-spaced points that covers 
# the study area and specify locations were predictions will be done)
data(meuse.grid)
class(meuse.grid)


coordinates(meuse.grid) <-  ~ x + y   
proj4string(meuse.grid) <- CRS("+init=epsg:28992")
class(meuse.grid)

# Plot the interpolation grid to see what you have
plot(meuse.grid, main = "Meuse Interpolation Grid with Control Points")
points(meuse,pch=10)
# Promote the grid data 
gridded(meuse.grid) <- TRUE
class(meuse.grid)


# Read in the meuse river outline

data(meuse.riv)
meuse.sr = SpatialPolygons(list(Polygons(list(Polygon(meuse.riv)),"meuse.riv")))
meuse.lt = list("sp.polygons", meuse.sr, fill = "grey")


#QUESTION01

# Summary Statistics
summary(meuse$cadmium) 
summary(log(meuse$cadmium))

# Steam and leaf plot
stem(meuse$cadmium) 
stem(log(meuse$cadmium))

# Histogram and Q-Q Plots
par(mfrow=c(3,2)) # 6 figures per page arranged in 3 rows and 2 columns
hist(meuse$cadmium, n=20, main = "Histogram of Cadmium (ppm)")
hist(log(meuse$cadmium), n=20, main = "Histogram of Log-Cadmium  (ppm)")

boxplot(meuse$cadmium, main = "Boxplot of Cadmium (ppm)")
boxplot(log(meuse$cadmium), main = "Boxplot of Log-Cadmium (ppm)")

qqnorm(meuse$cadmium,  main = "Q-Q Plot of Cadmium  (ppm)")
qqnorm(log(meuse$cadmium),  main = "Q-Q Plot of Log-Cadmium (ppm)")

par(mfrow=c(1,1)) # Reset to default plotting of 1 figure per page

# Question 2

# Obtain a bubble plot of the cadmium data with the dots scaled to concentration
bubble(meuse, "cadmium", col = "green" , main = "cadmium concentrations (ppm)")

# Obtain a dot plot and add some context to the plot
sp.theme(TRUE)
spplot(meuse, "cadmium", key.space = "right", col.regions=bpy.colors(),
       main = "Cadmium concentrations (ppm)",
       scales = list(draw = TRUE),# Add a reference by showing the coordinate system
       sp.layout= list("sp.polygons", meuse.sr, fill = "lightblue")
       # Add geographic reference (Meuse river boundaries)
)
# Question 3
# Map of distance to river: 
meuse.grid$sqrtdist = sqrt(meuse.grid$dist)
spplot(meuse.grid["sqrtdist"], col.regions = bpy.colors() , 
       sp.layout = list("sp.points", meuse, col = 3, cex=.5), 
       main = "Distance to river")

# Cadmium (ppm) vs. distance to river
xyplot(log(cadmium)~sqrt(dist), as.data.frame(meuse),
       main="Scatterplot of Cadmium vs. distance")


# Fit the regression model
cadmium.lm <- lm(log(cadmium)~sqrt(dist), meuse)

# Get a summary of the regression model
summary(cadmium.lm)

# Get diagnostic plots
layout(matrix(1:4, ncol=2))
plot(cadmium.lm, add.smooth = FALSE)
layout(1)

# Get Predicted Values and Standard Error of Fit for all locations on the grid
meuse.grid$lzn.fit <- predict(cadmium.lm, meuse.grid)
meuse.grid$se.fit <- predict(cadmium.lm, meuse.grid, se.fit=TRUE)$se.fit

# Plot the predicted values
spplot(meuse.grid, "lzn.fit", sp.layout = meuse.lt,
       main = "Log(Cadmium) - ppm: Regression Interpolation \n Predicted values")

# Plot the Standard Error of fit
spplot(meuse.grid, "se.fit", sp.layout = meuse.lt,
       main = "Log(Cadmium) - ppm: Regression Interpolation \n Standard Error of fit")

# Q4
meuse.grid$tr1 = krige(log(cadmium) ~ 1, meuse, meuse.grid, degree = 1)$var1.pred
spplot(meuse.grid, c("tr1", "tr2", "tr3"), sp.layout = meuse.lt,
       main = "Log(cadmium) - ppm \n Trend Surface Interpolation")
str(meuse.grid)
# Trend surface up to degree 2http://127.0.0.1:29233/graphics/3260c39c-1555-4b7d-82e8-0cbabff3fd4f.png
meuse.grid$tr2 = krige(log(cadmium) ~ 1, meuse, meuse.grid, degree = 2)$var1.pred

# Trend surface up to degree 3
meuse.grid$tr3 = krige(log(cadmium) ~ 1, meuse, meuse.grid, degree = 3)$var1.pred

# my own code Plotting degree 2 and 3 trend surfaces

spplot(meuse.grid, c("tr2", "tr3"), sp.layout = meuse.lt,
       main = c("Degree 2 Trend Surface", "Degree 3 Trend Surface"))

#Q5
# Run the IDW interpolation
# Check the effect of power (p) on the IDW interpolation 
# by changing the value of "idp", to the following values (1, 2.5, 5, 10)

meuse.grid$idwp05 = idw(log(cadmium)  ~ 1, meuse, meuse.grid, idp = 0.5)$var1.pred
# Plot the outputs and the control points 
spplot(meuse.grid, c("idwp05", "idwp1", "idwp2.5","idwp5", "idwp10"), 
       sp.layout = list("sp.points", meuse, col = 3, cex=.5),   
       main = "Log(Cadmium) - ppm , IDW Interpolation ")

#Q6
install.packages(c("sp", "gstat", "spatstat"))

# Creating lagged scatterplots
hscat(log(cadmium) ~ 1, meuse, (0:9) * 100)

# Creating a variogram cloud
plot(variogram(log(cadmium) ~ 1, meuse, cloud = TRUE))

# Creating a sample variogram (binned variogram) plot
plot(variogram(log(cadmium) ~ 1, meuse))

# Creating variograms at four different angles
plot(variogram(log(cadmium) ~ 1, meuse, alpha = c(0, 45, 90, 135)))

# Creating a variogram plot with overridden cutoff and interval width values
plot(variogram(log(cadmium) ~ 1, meuse, cutoff = 1000, width = 50))

# Creating a variogram with specified interval boundaries for the distance vector
variogram(log(cadmium) ~ 1, meuse, boundaries = c(0, 50, 100, seq(250, 1500, 250)))

# Creating a variogram
v <- variogram(log(cadmium) ~ 1, meuse)

# Plotting the variogram
plot(v)

# Fitting a spherical variogram model with specific initial parameter values
v.fit <- fit.variogram(v, vgm(1, "Sph", 800, 1))
# Print the fitted model parameters
print(v.fit)

# Performing partial fitting of variogram coefficients
fit.variogram(v, vgm(1, "Sph", 800, 0.06), fit.sills = c(FALSE, TRUE))

# Performing REML fitting of a variogram model
fit.variogram.reml(log(cadmium) ~ 1, meuse, model = vgm(0.6, "Sph", 800, 0.06))

# Calculating variograms for different directions
v.dir <- variogram(log(cadmium) ~ 1, meuse, alpha = (0:3) * 45)

# Fitting an anisotropic spherical variogram model
v.anis <- vgm(.6, "Sph", 1600, .05, anis = c(45, 0.3))

# Plotting experimental variograms and the fitted anisotropic model
plot(v.dir, v.anis)

# Creating a variogram map
plot(variogram(log(cadmium) ~ 1, meuse, map = TRUE, cutoff = 1000, width = 100))

# Performing simple kriging with a specified 'beta'
lz.sk <- krige(log(cadmium) ~ 1, meuse, meuse.grid, v.fit, beta = 5.9)

install.packages("ggplot2")

library(sp)
library(ggplot2)

# Create data frames with interpolated values
lz.sk <- krige(log(cadmium)~1, meuse, meuse.grid, v.fit, beta = 5.9)

sk.df <- data.frame(x = meuse.grid$x, y = meuse.grid$y, predicted = lz.sk$predict)
ok.df <- data.frame(x = meuse.grid$x, y = meuse.grid$y, predicted = lz.ok$predict)

# Convert data frames to SpatialPointsDataFrames
sk.points <- SpatialPointsDataFrame(coords = sk.df[, c("x", "y")], data = sk.df)
ok.points <- SpatialPointsDataFrame(coords = ok.df[, c("x", "y")], data = ok.df)


# Assuming 'meuse' is your training data and 'meuse.grid' is your grid data
# Assuming 'v.fit' is a variogram model
lz.sk <- krige(log(cadmium) ~ 1, meuse, meuse.grid, v.fit, beta = 5.9)

# Print out the structure of lz.sk to verify it's as expected
print(str(lz.sk))

# Check if predictions can be obtained from lz.sk
predictions <- predict(lz.sk)

# Print out the structure of predictions to see if it's valid
print(str(predictions))

# Create the data frame with predictions
sk.df <- data.frame(x = meuse.grid$x, y = meuse.grid$y, predicted = predictions$var1.pred)

# Print out the structure of sk.df to ensure it's created as expected
print(str(sk.df))
# Plot the maps
par(mfrow = c(1, 2))
plot(meuse, pch = 19, col = "gray", main = "Measured Values")
plot(sk.points, pch = 19, col = "blue", add = TRUE, cex = 0.7, main = "Simple Kriging")
plot(ok.points, pch = 19, col = "red", add = TRUE, cex = 0.7, main = "Ordinary Kriging")
legend("topright", legend = c("Measured", "Simple Kriging", "Ordinary Kriging"), 
       pch = 19, col = c("gray", "blue", "red"), bty = "n")
