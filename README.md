---
output: html_document
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/",
  out.width = "100%" 
)


```

<center>FORSEARCH </center>


   There are a large number of computer programs that apply a statistical method to a database.  The user assumes that the data follow a model and knows that the method is appropriate for that model. There are far fewer programs that check whether all the observations follow that model or, if some do not clearly follow it, what impact these outliers have on the estimation of unknown parameters of the model, both fixed and random.  The forsearch package evaluates the compliance of each observation to a model that could be analyzed with the lm() function or the glm() function of the stats package or the lme() function of the nlme package. 

   One approach to this goal would be to start with the complete analysis and to remove each observation in turn and reanalyze, looking for large changes in the estimates.  Then remove pairs of observations and repeat the process.  Clearly, the number of pairs, triplets, foursomes, etc quickly becomes too large to manage.  The forsearch approach is to begin with a minimal number of observations, and to increase the number by one until all observations are included.  For example, to estimate p parameters of a linear regression model, only p observations are needed to permit the mathematics to run to completion.  Naturally, no estimate of variation is possible, but this is not needed at this stage.  Atkinson and Riani (2000) refer to this as Step 1. 

In Step 2, the number of observations in the set is increased by 1 until all observations are included.  The first set of observations is chosen in such a way that there is little chance that it includes an outlier.  Subsequent sets are defined by adding the observation that least disturbs the estimates. This process causes the most outlying observation(s) to be added at the end of Step 2. In some stages of Step 2, one or more of the observations already in the set would be dropped out to make way for a better fitting set, but the set would still be incremented in size by 1.

In forsearch, Step 1 is accomplished by sampling all the observations to determine which ones should be chosen for the initial set. The size of this sample is set by the user (initial.sample) in forsearch_lm(), forsearch_glm(), and forsearch_lme().

During Step 2, the output of each run of a forsearch function is saved.  For example, fixed parameter estimates and the estimate of underlying variation are saved. For lme(), the other variance estimates are also saved. Other functions of the forsearch package are used to plot the changes in these statistics over the stages of Step 2. Outliers create characteristic plots.  The different plots indicate where the outlier(s) have an impact on the subsequent formal analysis of the data and where they do not. 

The plot functions can be configured to use Cairo graphics, if desired, and to export a Windows metafile, if desired.

If the observations that best fit the model are known (for example, from previous runs of a forsearch function, they can be manually entered in a call to these functions, thereby skipping Step 1 and saving some execution time. 

In this version of the package (2.0.0), we address data that would be analyzed as a linear fixed effects model, a linear mixed effects model, or a generalized linear model.  For the latter, acceptable families are binomial, gamma, and Poisson.  There are no limitations on the link functions of these families.

Some of the plotting functions can be used on any forsearch object and some are specific to the function used to generate the object.

## Installation

You can install the development version of forsearch as follows:

```
install.packages("forsearch")
```

## Examples

```
These are three basic examples which show you how to solve common problems:

Assume that the database (DB) consists of 97 observations (rows) and that there are 6 columns--the response 
y and 5 independent variables, x1 through x5. Also assume that there may be interaction between x1 and x2 
that will be evaluated. So there are not 5, but 7 fixed parameters: (intercept), x1, x2, x1:x2, x3, x4, 
and x5. We want to know how any outlier(s) would affect subsequent estimates of these parameters and the 
underlying variation. We want to assess each observation's leverage, and Cook's distance (a measure of 
impact on combinations of the parameters). 

The first thing is to add a column of observation numbers.  This must be an explicit column and it must be the 
first column. This is the case for all forsearch functions. 

Some graphical procedures restrict the number of plot lines in the interest of clarity. It is possible to specify 
which parameters will be displayed on each graph.  In this case, we might want to plot the intercept, x1, x2 and 
x1:x2 on one graph and the intercept, x3, x4 and x5 on a second graph.   

```

```
* Augment the data frame with Observation numbers

    Observation <- 1:97
    DB.O <- data.frame(Observation,DB)

* Produce the observation statistics

    DBfor.1 <- forsearch_lm(formula=y~x1*x2+x3+x4+x5,data=DB.O,initial.sample=100)

* Display the resulting file

    showme(DBfor.1)

* Get parameter display index

    identifyFixedCoeffs(DB.O)

* Plot the statistics

    plotdiag.residuals(data=DBfor.1)
    plotdiag.params.fixed(data=DBfor.1,coeff.codenums=c(1,2,3,4))
    plotdiag.params.fixed(data=DBfor.1,coeff.codenums=c(1,5,6,7))
    plotdiag.s2(data=DBfor.1)
    plotdiag.leverage(data=DBfor.1)
    plotdiag.Cook(data=DBfor.1)
    plotdiag.tstats(data=DBfor.1)
```

```

In the second example we consider a grouped data situation (Pinheiro and Bates, 2000). Three varieties of alfalfa 
were cut on four dates in a 3x4 pattern. The experimental units were assigned in 6 blocks, each subdivided into 
4 plots. Both the variety and the block effects will be treated as random. Assume that a previous analysis 
indicates that there is little if any effect of date of cutting.    

```

```

* Import the dataset

    data("nlme::Alfalfa")

* Remove the existing grouping structure

    Alf <- Alfalfa

* Augment the data frame with Observation number

    Observation <- 1:72
    Alf.O <- data.frame(Observation,Alf)

* Produce the observation statistics

    Alfalfa.O.formle<-forsearch_lme(fixed=Yield~1,data=Alf.O,random= ~1 | Block/Variety,
      formula=Yield~1|Block/Variety,response.column=5,initial.sample=100,robs=2)

* Display the resulting file

    showmelme(Alfalfa.O.formle)

* Get parameter display index

    identifyCoeffs(Alf.O)

* Plot the statistics

    plotdiag.residuals(data=Alfalfa.O.formle)
    plotdiag.params.fixed(data=Alfalfa.O.formle,coeff.codenums=NULL)
    plotdiag.params.random(data=Alfalfa.O.formle,coeff.codenums=c(1,5,6,7))
    plotdiag.leverage(data=Alfalfa.O.formle)
    plotdiag.Cook(data=Alfalfa.O.formle)
    plotdiag.tstats(data=Alfalfa.O.formle)

```

```
The third example is a set of 67 records of British train accidents (Atkinson and Riani, 2000). The number of 
deaths in the accident is related to the type of rolling stock, the year and month of the accident, and the 
amount of traffic on the railway system at the time. A single variable (ym) is first constructed from the 
year and month.

```

```
* Augment the data frame with Observation number

    Observation <- 1:67
    train.O <- data.frame(Observation, train)

* Produce the observation statistics, using 100 samples each with 5 observations

    train.O.for<-forsearch_glm(initial.sample=100,formula=Deaths~log(Traffic) + ym, family=poisson("log"),
          data=train.O, skip.step1= NULL, cobs=5)

* Display the resulting file

    showmegl(trainO.for)

* Get parameter display index

    identifyCoeffs(train.O)

* Plot the statistics

    plotdiag.deviances(data=train.O.for,null=FALSE)
    plotdiag.deviances(data=train.O.for,null=TRUE)
    plotdiag.params.fixed(data=train.O.for,coeff.codenums=NULL)
    plotdiag.deviance.residuals(data=train.O.for)
    plotdiag.phihatx(data=train.O.for)
    plotdiag.leverage(data=train.O.for)
    plotdiag.AICX(data=train.O.for)
    plotdiag.tstats(data=train.O.formle)

```

```
The table below summarizes the use of plotdiag functions:

```

```
------------------- | -------------------------------- | --- | --- | ---
plotdiag.xxx	    | plots data named 	               | lm	 | lme | glm  
------------------- | -------------------------------- | --- | --- | ---  
AICX                | AIC			                   |     |     | X  
Cook	            | Modified Cook distance	       | X   | X   |  
deviance.residuals	| Deviance residuals and augments  | 	 |     | X  
deviances	        | Residual deviance, Null deviance |	 |	   | X  
leverage	        | Leverage	                       | X	 | X   | X  
params.fixed	    | Fixed parameter estimates	       | X	 | X   | X  
params.random	    | Random parameter estimates       |	 | X   |  
phihatx	            | PhiHat			               |     |     | X  
residuals	        | Standardized residuals	       | X	 | X   |  
s2	                | s^2	                           | X	 |	   |  
tstats	            | t statistics	                   | X	 | X   | X  
```

