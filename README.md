# Welcome to ECHO!

<p align="center">
<img src="ECHO Shiny App/www/wc1_Neurospora_Replicates_Unsmoothed.PNG" width="500" />
</p>

This is the first step in the PAICE (Pipeline for Amplitude Integration of Circadian Exploration) Suite! This suite of tools provides high-throughput applications for circadian, ultradian, and infradian rhythms. The second step, ENCORE, can be found [here](https://github.com/delosh653/ENCORE), and can be run after using the application.

## README Outline

* Overview
* Use and First-Time Set-Up Instructions
* ECHO Features
* ECHO R Package
* Contact Information and Bug Reporting
* FAQ

## Overview

Note: currently in beta testing. Since this app is currently in development, please check here for updates, especially for an updated ECHO citation.

ECHO (Extended Circadian Harmonic Oscillators) is an R-powered application designed to find and visualize circadian rhythms from your data using extended harmonic oscillators. To read more about this work and cite us, see [*ECHO: an Application for Detection and Analysis of Oscillators Identifies Metabolic Regulation on Genome-Wide Circadian Output*](https://www.biorxiv.org/content/10.1101/690941v1) by H. De los Santos et al. (2019) To read more about the inital work on this project, see [*Circadian Rhythms in Neurospora Exhibit Biologically Relevant Driven and Damped Harmonic Oscillations*](https://dl.acm.org/citation.cfm?id=3107420&CFID=1006798692&CFTOKEN=88816139) by H. De los Santos et al. (2017) We also can also run another common circadian rhythm detection system, JTK_CYCLE. If you use their detection system, please cite them: [*JTK_CYCLE: an efficient non-parametric algorithm for detecting rhythmic components in genome-scale datasets*](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3119870/) by  ME Hughes et al. (2011)

All images created by ECHO using data from [*Analysis of clock-regulated genes in Neurospora reveals widespread posttranscriptional control of metabolic potential*](https://www.ncbi.nlm.nih.gov/pubmed/25362047) by J. Hurley, et al. (2014)

## Use and First-Time Set-Up Instructions

Thank you for downloading ECHO! ECHO is is an app designed to find and identify circadian rhythms from your data, then visualize the results. This guide will lead you in first time set-up and use. Pictures have been provided for ease of use, using Windows 7, in the files ECHO README.docx and ECHO_README.pdf, found above. A double asterisk indicates the step has an explanation below, and a tilde indicates the step is first-time set up only.

Steps: 
1.	** ~ Download [Firefox](https://www.mozilla.org/en-US/firefox/new/) or [Chrome](https://www.google.com/chrome/browser/desktop/index.html) and make it your default browser.
2.	~ [Download R](https://www.r-project.org/), if you do not already have it. 
3.	~ [Download RStudio](https://www.rstudio.com/products/rstudio/download/), if you do not already have it (RStudio Desktop is sufficient).
4.	Open RStudio.
5.	~ Copy and paste the following text into the console window (bottom left window of the RStudio Session), then press enter:

```r
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
install.packages("boot")
```

This will install these packages (a set of functions that this application uses) onto your computer. This may ask for your input, so just say yes to the questions asked. If you run into errors saying “yes,” just say no instead. Note: this may take some time.

6.	Open app.R, which should be included in the .zip file you downloaded and also contained README.pdf and README.docx. It should open in the top left window of your RStudio session.

7.	In the top right corner of the app.R window, you should see the button, “Run App”. Click on the small downwards arrow next to it and choose “Run External”. 

8.	Now click “Run App”. This should open the ECHO application in your now default browser window (either Firefox or Chrome). The picture below is a representation in Firefox.

9.	Have fun!

** Why do I have to install either Firefox or Chrome, you ask? Why not Internet Explorer, or some other browser? Well, it is known there are problems downloading files when viewing shiny apps in Internet Explorer, so we definitely want to avoid that. However, I have not tested this app in browsers like Microsoft Edge, Safari, etc. If you can verify that these work, please let me know at delosh@rpi.edu.

## ECHO Features

ECHO's interface is divided into two sections: **Finding Rhythms** and **Visualizing Results**.

Within the **Finding Rhythms** tab, you can upload your data (.csv) and enter its information, such as time point range, resolution (in hours), and amount and type of replicates. You can then choose from a variety of preprocessing steps including smoothing, removing unexpressed genes, and removing linear baselines. As mentioned above, you can also choose to run JTK_CYCLE from this tab. You can then download your results as both a .csv (for viewing) and a .RData (for visualizations).

In the **Visualizing Results** tab, simply upload the .RData file from your results and choose from several visualization and gene subset exploration options. You can explore subsets of data under the "Gene List" tab and sort by the various output parameters, such as Period or P-Value. You can also choose from a host of automatically-generated visualizations, including Venn diagrams, heat maps, gene expression plots (with or without replicates visualized), and parameter density graphs (examples displayed below).

<p align="center">
<img src="ECHO Shiny App/www/venn_diagram_by_adj_include_overdamped_Neurospora_Replicates_Unsmoothed.PNG" width="200" /> <img src="ECHO Shiny App/www/heat_map_Neurospora_Replicates_Smoothed_ECHO.PNG" width="200" /> <img src="ECHO Shiny App/www/wc1_Neurospora_Replicates_Unsmoothed.PNG" width="300" /> <img src="ECHO Shiny App/www/wc1_gene_expression_wo_rep_Neurospora_Replicates_Unsmoothed.PNG" width="300" /> <img src="ECHO Shiny App/www/forcing_coefficient_density_Neurospora_Replicates_Smoothed_ECHO.PNG" width="300" />
</p>

## ECHO R Package

ECHO's methodology is now available as an R package on CRAN! To download and use ECHO as a package, enter the following in the R console:

```r
install.packages("echo.find")
library(echo.find)
```

With this, you then have access to finding rhythms in data with one function, echo_find(). For more information on how to use this package and its functionality, check out the [echo.find vignette](https://cran.r-project.org/web/packages/echo.find/vignettes/echo-vignette.html).

Note that using this package requires knowledge of coding in R. If you having no coding knowledge, we recommend that you download and use the app as directed above. Also note that this version of ECHO does not take advantage of parallelism that the ECHO app does and therefore takes longer to run (only using one core, rather than using all cores except one). Further, there is no console output to show a progress bar of how long the output will take. If you would prefer built in parallelism and a progress bar, we highly recommend that you use the app. However, we have thought of several workarounds -- if interested, feel free to contact us with the "Feedback" form below.

## Contact Information and Bug Reporting

As you may have noticed, this is still in beta testing! Therefore, we would love to hear your feedback on the program. For general feedback, email delosh@rpi.edu with the subject line "ECHO Feedback".

If you run into any errors, please email delosh@rpi.edu with the following (subject line: "ECHO Error"): 
- a short desciption of your problem
- ECHO version number 
- your dataset/file(s) (this may be a sample of at least 50% of the data)
- your exact settings for the run (a screenshot will do) 
- your exact error from the console window (a screenshot will do)

However, please read the FAQ below before sending error reports.

Contact:
Hannah De los Santos /
email: delosh@rpi.edu /
Rensselaer Polytechnic Institute

## FAQ

**Q:** My dataset isn't working for some reason and I'm not using the latest ECHO version! Why?

**A:** Please update to the current ECHO version (you can see this at the top of the github repository by checking out what the latest commit was), since this may have been corrected in an update. If this problem persists, feel free to send another email!

---

**Q:** Does ECHO work with 24-hour time course length data to determine circadian rhythms?

**A:** Yes, it does! However, I would not categorize any rhythms convincingly "Damped", "Harmonic", or "Forced", since you need more than one cycle to determine this.

---

**Q:** My data has starting points/ending points/resolution of less than an hour, or a fractional amount of hours! How do I run this through ECHO?

**A:** If you have resolution of less than an hour, please enter the fraction into the box, in the form: numerator/denominator. For example, if my resolution was every 10 minutes (or 6 times every hour), I would enter: 1/6. This fractional form extends to starting and ending points as well. You must enter the fraction, NOT a mixed number. For example, if my starting time was 16 hours, 10 minutes, my starting time would be: 97/6. (This stems from the following calculation: (1/6 x 16) + (1/6))

---

**Q:** I get the following warning and errors when I try to run my data through ECHO:

```r
Warning in read.table(file = file, header = header, sep = sep, quote = quote,  :
  incomplete final line found by readTableHeader on '....csv'

Warning: Error in +: non-conformable arrays
  76: avg_all_rep [damped_oscillator_master.R#157]
  72: observeEventHandler [/Users/delosh/Desktop/ECHO-master/ECHO Shiny App/app.R#442]
```

**A:** Double check that your data is actually a .csv, i.e. actually comma-delimited. Just because it has ".csv" as a file extension may not necessarily mean that it's comma-delimited! You can check this by opening your file in a text editor, such as Notepad. If your data is comma-delimited, you should see commas between your entries. If it is not, you may see spaces or tabs. In this case, please save your data as a proper comma-delimited file.

---

**Q:** I get the following warnings and errors when I try to visualize my results using the .RData file:
```r
Warning in circ_us & circ_jtk:
  longer object length is not a multiple of shorter object length
...
Warning: Error in data.frame: arguments imply differing number of rows
```
**A:** Double check that you do not have any duplicate gene names in your dataset. If you do, differentiate these names somehow. JTK_CYCLE automatically removes duplicate names while ECHO doesn't, leading to this error.

---

**Q:** I get the following error when I try to run my dataset using ECHO:
```r
Warning: Error in {: task 3 failed - "missing value where TRUE/FALSE needed"
```
**A:** Check your data csv to make sure that your first row is data labels, your first column is only character strings (such as "ABCD"), and your numerical data only contains numbers and blank spaces. In R, there is a simple way to check this (enter the following code in your console window):
```r
install.packages("readr") # if you already have downloaded this package, put a # before install
library(readr)
genes <- read_csv("YOUR FILE EXTENSION HERE (place between quotes, must use forward slashes)")
```
Here is an example of a file extension: "C:/Users/delosh/Documents/example.csv" . 
If your data is correctly structured, you should see the following (replace "Gene Symbol" with whatever your gene name column is titled):
```r
Parsed with column specification:
cols(
  .default = col_double(),
  `Gene Symbol` = col_character()
)
See spec(...) for full column specifications.
```
If you have problems in your data (in this error case, nonnumerical data in your numerical data), you will get the following error (or similar):
```r
Parsed with column specification:
cols(
  .default = col_double(),
  `Gene Symbol` = col_character()
)
See spec(...) for full column specifications.
Warning: 72 parsing failures.
 row  col expected  actual
1250 16.1 a double #DIV/0!
1250 16.2 a double #DIV/0!
1250 18.1 a double #DIV/0!
.... .... ........ .......
See problems(...) for more details.
```
In this case, nonnumerical data was intermixed with numerical data (#DIV/0! was where a number or a blank cell should have been). This user would have to go back and change these to the appropriate number or blank space.

