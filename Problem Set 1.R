
# Instructions for Problem Set 1 ------------------------------------------

# Using the ‘twindata’ dataset packaged with OpenMx, evaluate the role of genetic and
# environmental factors in explaining variation in body mass index in young males (MZ: zyg=2;
#                                                                                  DZ: zyg=4).
# Provide a summary of results in tables (goodness-of-fit statistics and estimates) and write a
# summary paragraph(s) with key results, addressing model assumptions.
# Maximum 1 page including tables, due Tuesday 9/1/20 at 9am to hmaes@vcu.edu.
#
#go to https://hermine-maes.squarespace.com/#/one/ for sample code to run the various twin models


# Univariate Saturated Model -------------------------------------------------------------------


library(OpenMx)

head(twinData)
0.78/2 #most likely will use ADE model 
