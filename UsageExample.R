
# setwd("Set this to the directory that contains the downloaded file")


source(file="SS4TrtSelSurv.R")

# Parameter settings
normalmean = 0.5;
normalsd = sqrt(1/12)

#  Parameter for the number of MC runs for each estimation (fixed)
# The k-values for input
ExpPropOfDeaths = .75;
NumberOfBootstraps = 200;
mydesign = "strategy"
k1 = 0.25; k2 = 0.75;
k3 = 0.75; k4 = 0.25;
t0=5;

start1 = Sys.time()
temp1 <- Calculate.sample.size(targetwidth=0.10,k1,k2,k3,k4,t0,
                        mydesign=mydesign,ExpPropOfDeaths=ExpPropOfDeaths,
                        NumberOfBootstraps=NumberOfBootstraps)
end1 = Sys.time()

print(temp1);

plot(as.numeric(temp1[[3]]),as.numeric(temp1[[4]]),xlab="Sample Size", ylab="Inverse squared Width")


