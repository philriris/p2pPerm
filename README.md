# p2pPerm
## A function for retrieving metrics for resistance-resilience from _rcarbon_ `SpdPermTest` objects.

`p2pPerm()` is an adaptation of Edinborough et al.'s post-hoc statistical test for demographic events in written and oral history (https://doi.org/10.1073/pnas.1713012114). The original function - `p2pTest()` in _rcarbon_ - is for use with objects of class `SpdModelTest`. This function extends the principle to `SpdPermTest` objects, with similar options for summed probability distribution (`type = spd`) or Rate of Change (`type = roc`) versions. Focal marks can also be manually specified with the `focalm` parameter. It assumes the `permTest()` function has been called with `raw = TRUE`. 

The function also has options to specify `type = res` or `type = ros`. Following, principally, Nimmo et al. (2015; https://doi.org/10.1016/j.tree.2015.07.008), Cantarello et al. (2017; https://doi.org/10.1002/ece3.3491), and Van Meerbeek et al. (2021; https://doi.org/10.1111/1365-2745.13651). This will perform post-hoc tests for **resistance** and **resilience** on marks of an `SpdPermTest` object over a user-defined interval. Informally, these two metrics are defined as the ability to absorb disturbances and "bounce back" following disturbances, respectively. They are normalised relative to the value of the SPD, or rate of change, at the start of the interval of interest. In an analogous fashion to `p2pTest()`, it employs the simulation envelope produced by `permTest()` to perform a two-sided test of significance for both metrics.

The function outputs a list containing the value and p-value of both metrics, as well as the 'lag'. This is the minimum of the SPD or RoC curve in the interval between the start and end dates, and is reported as the number of years between the date of the minimum and start date. 
