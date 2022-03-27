
4/17/2022 forsearch version 2.2.0
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


