1/7/2025 forsearch version 6.4.0
================================
* Major changes: None
 
* Moderate changes:
  + Added vignette on Study Credibility  
  + Added arguments to forsearch_nls to liberalize controls
  + Added function logist3 to calculate 3-parameter logistic 
  
* Minor changes and bug fixes:
  + Removed phases in forsearch_nls in favor of sections
  + Returned to use of stats::predict for all forsearch
    functions
  + Updated all examples  

11/15/2024 forsearch version 6.3.0
==================================
* Major changes: None
 
* Moderate changes: 
  + Added candprep support function
  + Optimized step 1 of every forsearch_foo function
  + Identified observation set within step 1 as yielding
    lowest median squared error WITHIN each factor
    subset
  + Changed to 2-phase procedure in step 2 to ensure all
    parameters can be estimated at every stage in step 2
  + Added arguments to forsearch_lme to liberalize controls
  + Changed from predict to predict.lme, predict.glm and
    predict.coxph in xStep1 and xStep2 
 
* Minor changes and bug fixes: 
  + Modified showme to limit output for forsearch_lme
  + Added graphs for forsearch_lme
  + Corrected use of deviances in forsearch_glm
  + Removed forsearch_nls for corrections and improvements

7/13/2024 forsearch version 6.2.0
=================================
* Major changes: None
 
* Moderate changes: 
  + Revised aStep1 to improve step 1 for _lm
  + Added bStep1 and cStep1 to improve step 1 for _lme and _cph

* Minor changes and bug fixes: 
  + Added option to disturb event times slightly to avoid duplicates in forsearch_cph
    observations
  + Changed message in plotdiag.ANOX2 if proportionality not tested in
    proportional hazard data
  + Minor text changes to vignettes
  + Corrections to examples to correspond to function changes

6/22/2024 forsearch version 6.1.0 
=================================
* Major changes: None  
 
* Moderate changes:
  + Enabled models with constructed variables (those that appear in the formula
    but not in the database)

* Minor changes and bug fixes: None 
 
3/30/2024 forsearch version 6.0.0 
=================================
* Major changes:  
  + Added functions forsearch_nls, eStep1, and eStep2 for nonlinear statistics
  
* Moderate changes: None

* Minor changes and bug fixes:
  + Added code to function showme for forsearch_nls objects and to function
      plotdiag.allgraphs for forsearch_nls graphics 
  + Modified vignette regarding where we get observations to cover nonlinear
      models

2/16/2024 forsearch version 5.1.0 
=================================
* Major changes: None 

* Moderate changes: 
  + Completely revised flow of forsearch_glm to better accommodate models
       containing independent variables as factors 
  + Added vignette regarding exploration of the search history

* Minor changes and bug fixes:
  + Fixed proportionality test bug (na's) in plotdiag.ANOX2 for Cox regression
  + Changed use of n.obs.per.level in picksome to be additive
  + Remove unused argument of bStep2
  + Added blind.label to arguments of plotdiag.allgraphs to allow
       user to request labeling of graph to reflect blinded database
  
1/13/2024 forsearch version 5.0.0 
=================================
* Major changes: 
  + Completely revised flow of forsearch functions _lm, _lme and _cph
       to better accommodate models containing independent variables 
       as factors
  
* Moderate changes: 
  + Added vignette describing strategy of the forsearch procedure
  + Added ability to skip test of proportionality in forsearch_cph

* Minor changes and bug fixes:
  + Modified diagnostic printouts for forsearch_cph
  + Removed option in variablelist and in picksome for identifying function
  + Integrated diagnostic printouts of aStep1, aStep2 and bStep2 
       with those of calling functions
  + Added requirement that forsearch_cph have at least 1 factor variable
  + Fixed bug preventing use of multiple observations in Step 1.
  + Cleaned up examples 

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


