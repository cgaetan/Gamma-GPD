# Gamma-GPD

This is the repository for the code related to the paper 

Jean-Noël Bacro, Carlo Gaetan, Thomas Opitz & Gwladys Toulemonde (2020) Hierarchical Space-Time Modeling of Asymptotically Independent Exceedances With an Application to Precipitation Data, Journal of the American Statistical Association, 115:530, 555-569, DOI: 10.1080/01621459.2019.1617152



Data 
----


Data set in the paper consists of hourly observations at 50 rainfall stations for
the years 1993 to 2014 in the mountainous areas of the Cevennes and the
Alps (France). Extreme rainfall events usually occur during fall season;
therefore, our application considers only observations from the
September to November months. Elevations at gauge sites range from 140
to 1000 meters. Data have been extracted from the database of the French
Met Office, Météo France.

The data cannot be freely distributed, but access can be requested to Météo France for academic purposes. 

Data can be downloaded from the following web interface of
[Météo France](https://publitheque.meteo.fr/okapi/accueil/okapiWebPubli/index.jsp)
site.

To obtain the stations used in the application study in our manuscript,
the following identifiers of measurement sites must be requested:

 7011004    7025001    7032003    7117001    7131001  7159001    7172002    7228001    7286004    7330004  
 7334003    26002003   26047001   26050001   26115001 26124001   26176001   26177001   26202001   26292002 
 26298001   26313001   26361001   30068001   30087002 30132003   30273001   30299001   34008001   34028003 
 34030002   34186001   34209002   34217001   34317001 34344001   43096001   43150001   43234005   43268005 
 48027003   48069001   48095004   48171001   84036001 84082001   84085003   84094001   84107002   84150001 



Third-party code 
----

We follow

Hughes, G. B., and Chraibi, M. (2012), *Calculating ellipse overlap
areas*, Computing and Visualization in Science, 15, 291--301.

in order to efficiently calculate the ellipse intersection area.

We used a selected part of the code from [their repository](http://github.com/chraibi/EEOver), namely 

Roots3And4.c, solvers-2.c, solvers.h, zsolve_quartic.c, config.h, program_constants.h



Instructions for use 
----

Code has been tested on a workstation running a Linux OS


1\) Install the [GSL](https://www.gnu.org/software/gsl/)  library. 

For instance under a Linux Debian
distribution, open a terminal and type

sudo apt-get install libgsl0ldbl libgsl-dev

For Mac users, please follow the instructions [here](http://macappstore.org/gsl/) 

2\) Compile the source files and then link all specified object files into a shared object

For instance under  Linux OS open a terminal and type

R CMD SHLIB spt-gamma.c zsolve_quartic.c solvers-2.c Roots3And4.c -lgsl -lgslcblas



3\) file simulation-gamma.R contains the code for a simulated example. To run this 

R CMD BATCH simulation-gamma.R
