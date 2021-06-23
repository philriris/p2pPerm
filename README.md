# p2pPerm
## A function for retrieving metrics for resistance-resilience from _rcarbon_ `SpdPermTest` objects.

`p2pPerm()` is an adaptation of Edinborough et al.'s post-hoc statistical test for demographic events in written and oral history (https://doi.org/10.1073/pnas.1713012114). The original function - `p2pTest()` in _rcarbon_ - is for use with objects of class `SpdModelTest`. This function extends the principle to `SpdPermTest` objects, with similar options for summed probability distribution (`type = spd`) or Rate of Change (`type = roc`) versions. Focal marks can also be manually specified with the `focalm` parameter. It assumes the `permTest()` function has been called with `raw = TRUE`. 

The function also has options to specify `type = res` or `type = ros`. Following, principally, Nimmo et al. (2015; https://doi.org/10.1016/j.tree.2015.07.008), Cantarello et al. (2017; https://doi.org/10.1002/ece3.3491), and Van Meerbeek et al. (2021; https://doi.org/10.1111/1365-2745.13651), this will perform post-hoc tests for **resistance** and **resilience** on marks of an `SpdPermTest` object over a user-defined interval. Informally, these two metrics are defined as the ability to absorb disturbances and "bounce back" following disturbances, respectively. They are normalised relative to the value of the SPD, or rate of change, at the start of the interval of interest. In an analogous fashion to `p2pTest()`, it employs the simulation envelope produced by `permTest()` to perform a two-sided test of significance for both metrics. By default, the `p3` parameter is the minimum of the SPD or RoC curve between `p1` and `p2` but it can be manually specified.

The function outputs a list containing the value and p-value of both metrics, as well as the 'lag'. This is the minimum of the SPD or RoC curve in the interval between the start and end dates, and is reported as the number of years between the date of the minimum and start date. 

## Usage

```R
devtools::source_url("https://raw.githubusercontent.com/philriris/p2pPerm/main/R/p2pPerm.R")

library(rcarbon)

# Prepare data
data(emedyd) # Load data for Eastern Mediterranean Younger Dryas
caldates <- calibrate(x=emedyd$CRA, errors=emedyd$Error, normalised=FALSE) # Calibrate
bins <- binPrep(sites=emedyd$SiteName, ages=emedyd$CRA, h=200) # Bin

# Mark permutation test
ydperm <- permTest(caldates, marks=emedyd$Region, bins=bins, nsim=1000, runm=100, backsight=100, 
timeRange=c(13400,11000), datenormalised=FALSE, raw=TRUE) # Time range roughly +/- 500 years of Younger Dryas

# not run
# plot(ydperm, focalm=1)
# plot(ydperm, focalm=2)
# plot(ydperm, focalm=3)

par(mfrow=c(3,1))
slevant <- p2pPerm(x=ydperm,p1=13000,p2=12500, type="ros", focalm = "1", plot=TRUE) # Rate of Change version
nlevant <- p2pPerm(x=ydperm,p1=13000,p2=12500, type="ros", focalm = "2", plot=TRUE) 
anatolia <- p2pPerm(x=ydperm,p1=13000,p2=12500, type="ros", focalm = "3", plot=TRUE) 
dev.off()
```

![image](https://user-images.githubusercontent.com/37066355/122895977-86ed5300-d340-11eb-9591-6345675a14cc.png)

```
plot(y=slevant$resistance,x=slevant$resilience,  xlim=c(-1,1), ylim=c(-1,1), col="green", pch=19, cex=2,
     ylab="Resistance", xlab="Resilience")
points(y=nlevant$resistance,x=nlevant$resilience,  xlim=c(-1,1), ylim=c(-1,1), col="purple", pch=19, cex=2)
points(y=anatolia$resistance,x=anatolia$resilience, xlim=c(-1,1), ylim=c(-1,1), col="darkorange", pch=19, cex=2)
```

![image](https://user-images.githubusercontent.com/37066355/122896120-a97f6c00-d340-11eb-94e0-61ac5ba9853b.png)


