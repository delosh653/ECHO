# Welcome to ECHO!

Note: currently in beta testing.

ECHO (Extended Circadian Harmonic Oscillations) is an R-powered application designed to find and identify circadian rhythms from your data using extended harmonic oscillators. To read more about our inital work on this project, see here (https://dl.acm.org/citation.cfm?id=3107420&CFID=1006798692&CFTOKEN=88816139): Hannah De los Santos, Emily J. Collins, Jennifer M. Hurley, and Kristin P. Bennett. 2017. Circadian Rhythms in Neurospora Exhibit Biologically Relevant Driven and Damped Harmonic Oscillations. In Proceedings of the 8th ACM International Conference on Bioinformatics, Computational Biology,and Health Informatics (ACM-BCB '17). ACM, New York, NY, USA, 455-463. DOI: https://doi.org/10.1145/3107411.3107420 

We also can also run another common circadian rhythm detection system, JTK_CYCLE. If you use their detection system, please cite them here (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3119870/): Hughes ME, Hogenesch JB, Kornacker K. JTK_CYCLE: an efficient non-parametric algorithm for detecting rhythmic components in genome-scale datasets. Journal of biological rhythms. 2010;25(5):372-380. doi:10.1177/0748730410379711.

# Use and First-Time Set-Up Instructions

Thank you for downloading ECHO (Extended Circadian Harmonic Oscillations)! ECHO is is an app designed to find and identify circadian rhythms from your data, then visualize the results. This guide will lead you in first time set-up and use. Pictures have been provided for ease of use, using Windows 7, in the files ECHO README.docx and ECHO_README.pdf, found above. A double asterisk indicates the step has an explanation below, and a tilde indicates the step is first-time set up only.

Last updated: 11/17/17 (ECHO version 0.3)

Steps: 
1.	** ~ Download Firefox (https://www.mozilla.org/en-US/firefox/new/) or Chrome (https://www.google.com/chrome/browser/desktop/index.html) and make it your default browser.
2.	~ Download R, if you do not already have it: https://www.r-project.org/
3.	~ Download RStudio, if you do not already have it (RStudio Desktop is sufficient): https://www.rstudio.com/products/rstudio/download/
4.	Open RStudio.
5.	~ Copy and paste the following text into the console window (bottom right window of the RStudio Session), then press enter:

install.packages("rstudioapi")
install.packages("shiny")
install.packages("ggplot2")
install.packages("VennDiagram")
install.packages("reshape2")
install.packages("minpack.lm")
install.packages("doParallel")
install.packages("foreach")
install.packages("iterators")
install.packages("doSNOW")
install.packages("colorRamps")
install.packages("fields")

This will install these packages (a set of functions that this application uses) onto your computer. This may ask for your input, so just say yes to the questions asked. If you run into errors saying “yes,” just say no instead. Note: this may take some time.

6.	Open app.R, which should be included in the .zip file you downloaded and also contained this README. It should open in the top left window of your RStudio session.

7.	In the top right corner of the app.R window, you should see the button, “Run App”. Click on the small downwards arrow next to it and choose “Run External”. 

8.	Now click “Run App”. This should open the ECHO application in your now default browser window (either Firefox or Chrome). The picture below is a representation in Firefox.

9.	Have fun!

** Why do I have to install either Firefox or Chrome, you ask? Why not Internet Explorer, or some other browser? Well, it is known there are problems downloading files when viewing shiny apps in Internet Explorer, so we definitely want to avoid that. However, I have not tested this app in browsers like Microsoft Edge, Safari, etc. If you can verify that these work, please let me know at delosh@rpi.edu.

# Contact Information and Bug Reporting

As you may have noticed, this is still in beta testing! Therefore, we would love to hear your feedback on the program. For general feedback, email delosh@rpi.edu with the subject line "ECHO Feedback".

If you run into any errors, please email delosh@rpi.edu with the following (subject line: ECHO Error): 
- a short desciption of your problem
- ECHO version number 
- your dataset/file(s) 
- your exact settings for the run (a screenshot will do) 
- your exact error from the console window (a screenshot will do) 

Contact:
Hannah De los Santos /
email: delosh@rpi.edu /
Rensselaer Polytechnic Institute
