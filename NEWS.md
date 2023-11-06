11/6/2023 forsearch version 4.2.0 
==================================
* Major changes: None
  
* Moderate changes:
  + Added plotdiag.ANOX2 and applied it to 
      proportionality test output of forsearch_cph, and test of null
      hypothesis in forsearch_lm and forsearch_lme
  + Increased scope of models that can by handled by forsearch_lme by
      revising forsearch_lme and bStep2. Removed need for bStep1
  + Added vector of Step 1 observation numbers to all forsearch_xxx functions    

* Minor changes and bug fixes:
  + Added proportionality test (cox.zph) to forsearch_cph
  + Added test of null hypothesis (anova) to forsearch_lm and forsearch_lme
  + Added proportionality test to showme for forsearch_cph
  + Added null hypothesis test to showme for forsearch_lm and forsearch_lme
  + Correct failure to define response.colnum in forsearch_lme that prevented
     skip of Step1
  + Changed argument "addline" in plotdiag.Cook to single value default
  
9/12/2023 forsearch version 4.1.0 
=================================
* Major changes: None

* Moderate changes: None

* Minor changes and bug fixes:   
  + Replaced argument diagnose=T/F with begin.diagnose numeric indicator   
  + Removed argument diagnose and diagnostic code in several small functions     
  + Added examples in all forsearch_xxx functions   
  + Cleaned up documentation wording   

8/8/2023 forsearch version 4.0.0 
================================
* Major changes:
  + Added function forsearch_cph for Cox proportional hazards
  + Added function cStep2 for step 2 of Cox proportional hazard function
  + Added functions plotdiag.Wald, plotdiag.lrt, and plotdiag.loglik

* Moderate changes:
   + Combined showme, showmegl, and showmelme into showme function  
   + Added forsearch_cph to showme function

* Minor changes and bug fixes:
   + Correct display of ANOVA structure in forsearch_glm for family=Gamma
   + Correct indexing in aStep1
   + Correct function ID in plotdiag.leverage
   + Set default for addline in plotdiag.s2
   + Added message in variablelist regarding single observations in some factor levels
   + Correct coding error in picksome function

2/17/2023 forsearch version 3.2.0 
=================================
* Major changes: None

* Moderate changes:
   + Added guidance to prevent linearly dependent X'X matrices in forsearch_glm
   + Added capability of forsearch_glm to handle ANOVA and ANCOVA models

* Minor changes and bug fixes:
   + Correct coding error in picksome function

1/4/2023 forsearch version 3.0.1 
================================
* Major changes: None

* Moderate changes: None

* Minor changes and bug fixes:
   + Fixed response column identifier in forsearch_lme, aStep1, aStep2, bStep1 and bStep2
   + Fixed need for handling NaN in certain data for forsearch_glm   
   + Clarified help for forsearch_glm for two response column numbers 

12/14/2022 forsearch version 3.0.0 
==================================
* Major changes:
   + Rewrote forsearch_lme to clarify code for generation of Steps 1 and 2
   + Extracted Step 1 and Step 2 increases in sample size to separate functions aStep1 and aStep2
      from forsearch_lm and forsearch_lme
   + Added functions aStep1 and aStep2 for linear model functions
   + Added functions bStep1 and bStep2 for linear mixed effects models

* Moderate changes:
   + Added plotdiag.all function to create drafts of all applicable graphs
   + Added functions variablelist() and picksome() to subdivide observations into their factor levels
   + Added capability of forsearch_lm, forsearch_glm and forsearch_lme to process factor variables 
      and interactions with regression variables
   + Added capability of forsearch_lm to calculate leverage of observations over subsets
   + Added examples for forsearch_lm, forsearch_glm, and forsearch_lme

* Minor changes and bug fixes:
   + Changed README to add plotdiag.all, variablelist, and picksome functions
   + Changed some plotting functions to allow no fit of a line to the data
   + Added option for all forsearch_xxx functions to blind initial display of data analysis structure, 
      anticipating future developments 
 

4/21/2022 forsearch version 2.3.0
=================================
* Major changes: None

* Moderate changes:
   + Correct subgroup selection procedure in forsearch_lme
   + Calculate root mean square of random coefficients instead of the coefficients themselves
     in forsearch_lme   
   + Correct plotdiag.params.random to correspond with changed calculation of RMS
   + Add plotdiag.fit3 to plot AIC, BIC, and log likelihood for forsearch_lme objects

* Minor changes and bug fixes:
   + Print summary of assumed analysis in forsearch_lm and in forsearch_lme to permit 
     formulation check
   + Modify README and documentation of functions

3/24/2022 forsearch version 2.2.0
=================================
* Moderate changes:
   + Clean up README and help files
   + Modify forsearch_glm to handle observations of 0 or N (binomial)

3/17/2022 forsearch version 2.1.0
=================================
* Moderate changes:
   + Correct the omission of README file

3/16/2022 forsearch version 2.0.0
=================================
* Major changes:
   + Added forsearch_glm() function
   + Added plotdiag.deviance.residuals() function
   + Added plotdiag.phihatx() function
   + Added plotdiag.AICX() function
   + Added plotdiag.deviances() function

* Moderate changes:
   + Added the collection of t statistics in forsearch_lm() and forsearch_lme() for subsequent plotting 
   + Added plotdiag.tstats() function to plot t statistics
   + Added showmegl() function
   + Changed some existing plotdiag functions to accommodate forsearch_glm()

* Minor changes and bug fixes:
   + Added t statistics to list of summary prints in showme() and showmelme()
   + Modified search.history() to account for case when all current observations replaced by new ones
   + Enabled user to enter caption on all graphs
   + Corrected 2nd example in README file
   + Updated README file
   + Removed extra squaring of response in plotdiag.residuals() function


1/27/2022 forsearch version 1.0.0
=================================
First version of this package


