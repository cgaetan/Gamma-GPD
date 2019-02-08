# Gamma-GPD-process
---
Hierarchical space-time modeling of asymptotically independent exceedances
---

This is the repository for the code related to the paper 

Bacro, J.-N., Gaetan, C., Opitz, T., Toulemonde, G. (2019) *Hierarchical space-time modeling of asymptotically independent exceedances*, Journal of the American Statistical Association, to appear

Data 
----


Data set in the paper consists of hourly observations at 50 rainfall stations for
the years 1993 to 2014 in the mountainous areas of the Cevennes and the
Alps (France). Extreme rainfall events usually occur during fall season;
therefore, our application considers only observations from the
September to November months. Elevations at gauge sites range from 140
to 1000 meters. Data have been extracted from the database of the French
Met Office, Météo France.

The data cannot freely distributed, but access can be requested to Météo France for academic purposes. 

Data can be downloaded from the following web interface of
[Metéo France](https://publitheque.meteo.fr/okapi/accueil/okapiWebPubli/index.jsp)
site.

To obtain the stations used in the application study in our manuscript,
the following identifiers of measurement sites must be requested:

\"26047001\" \"26050001\" \"26100001\" \"26113003\" \"26115001\"
\"26124001\" \"26142001\" \"26167001\" \"26179001\" \"26202001\"
\"26258001\" \"26264002\" \"26313001\" \"26321003\" \"26348001\"
\"26361001\" \"84125003\"



Third-party code 
----



To efficiently calculate the ellipse intersection area, we follow

Hughes, G. B., and Chraibi, M. (2012), *Calculating ellipse overlap
areas*, Computing and Visualization in Science, 15, 291--301.*

and we use selected part of the code in http://github.com/chraibi/EEOver, namely 

Roots3And4.c, solvers-2.c, solvers.h, config.h, program_constants.h



Instructions for Use 
----

Code has been tested on a workstation running a Linux OS


1\) Install the [GSL](https://www.gnu.org/software/gsl/)  library. 

For instance under a Linux Debian
distribution, open a terminal and type

sudo apt-get install libgsl0ldbl libgsl-dev

For Mac users, please follow the instructions [here](http://macappstore.org/gsl/) 

2\) Compile the source files and then link all specified object files into a shared object

For instance under  Linux OS open a terminal and type

R CMD SHLIB spt-gamma.c zsolve_quartic.c solvers-2.c Roots3And4.c
-L/usr/lib/x86\_64-linux-gnu -lgsl -lgslcblas



3\) file simulation-gamma.R contains the code for a simulated example. To run this 

R CMD BATCH simulation-gamma.R
