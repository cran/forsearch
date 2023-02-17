
2/14/2023 forsearch version 3.1.0 
=================================
     Major changes:
None
     Moderate changes:
Added guidance to prevent linearly dependent X'X matrices in forsearch_glm
Added capability of forsearch_glm to handle ANOVA and ANCOVA models

     Minor changes and bug fixes:
Correct coding error in picksome function

1/4/2023 forsearch version 3.0.1 
================================
     Major changes:
None
     Moderate changes:
None
     Minor changes and bug fixes:
Fixed response column identifier in forsearch_lme, aStep1, aStep2, bStep1 and bStep2
Fixed need for handling NaN is certain data for forsearch_glm   
Clarified help for forsearch_glm for two response column numbers 

12/14/2022 forsearch version 3.0.0 
==================================
     Major changes:
Rewrote forsearch_lme to clarify code for generation of Steps 1 and 2
Extracted Step 1 and Step 2 increases in sample size to separate functions aStep1 and aStep2
   from forsearch_lm and forsearch_lme
Added functions aStep1 and aStep2 for linear model functions
Added functions bStep1 and bStep2 for linear mixed effects models

     Moderate changes:
Added plotdiag.all function to create drafts of all applicable graphs
Added functions variablelist() and picksome() to subdivide observations into their factor levels
Added capability of forsearch_lm, forsearch_glm and forsearch_lme to process factor variables 
   and interactions with regression variables
Added capability of forsearch_lm to calculate leverage of observations over subsets
Added examples for forsearch_lm, forsearch_glm, and forsearch_lme

     Minor changes and bug fixes:
Changed README to add plotdiag.all, variablelist, and picksome functions
Changed some plotting functions to allow no fit of a line to the data
Added option for all forsearch_xxx functions to blind initial display of data analysis structure, anticipating future developments 
 

4/21/2022 forsearch version 2.3.0
=================================
     Major changes:
None

     Moderate changes:
Correct subgroup selection procedure in forsearch_lme
Calculate root mean square of random coefficients instead of the coefficients themselves
     in forsearch_lme   
Correct plotdiag.params.random to correspond with changed calculation of RMS
Add plotdiag.fit3 to plot AIC, BIC, and log likelihood for forsearch_lme objects

     Minor changes and bug fixes:
Print summary of assumed analysis in forsearch_lm and in forsearch_lme to permit 
     formulation check
Modify README and documentation of functions



3/24/2022 forsearch version 2.2.0
=================================
Clean up README and help files
Modify forsearch_glm to handle observations of 0 or N (binomial)


3/17/2022 forsearch version 2.1.0
=================================
Correct the omission of README file


3/16/2022 forsearch version 2.0.0
=================================
    Major changes:
Added forsearch_glm() function
Added plotdiag.deviance.residuals() function
Added plotdiag.phihatx() function
Added plotdiag.AICX() function
Added plotdiag.deviances() function

    Moderate changes:
Added the collection of t statistics in forsearch_lm() and forsearch_lme() for subsequent plotting 
Added plotdiag.tstats() function to plot t statistics
Added showmegl() function
Changed some existing plotdiag functions to accommodate forsearch_glm()

    Minor changes and bug fixes:
Added t statistics to list of summary prints in showme() and showmelme()
Modified search.history() to account for case when all current observations replaced by new ones
Enabled user to enter caption on all graphs
Corrected 2nd example in README file
Updated README file
Removed extra squaring of response in plotdiag.residuals() function


1/27/2022 forsearch version 1.0.0
=================================
First version of this package


