# Circadian Rhythms GUI
# by Hannah De los Santos

# libraries and resources to load once ----

library(rstudioapi)
library(shiny)
library(ggplot2)
library(VennDiagram)
library(reshape2)
library(minpack.lm)
library(doParallel)
library(foreach)
library(iterators)
library(doSNOW)
library(colorRamps)
library(fields)

#https://stackoverflow.com/questions/3452086/getting-path-of-an-r-script/35842176#35842176
# set working directory - only works in RStudio (with rstudioapi)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# clear workspace (sorry!)
rm(list=ls())

options(shiny.maxRequestSize=80*1024^2) # enlarge the max file acceptance size

# https://stackoverflow.com/questions/13445435/ggplot2-aes-string-fails-to-handle-names-starting-with-numbers-or-containing-s
# ggName -> changes a string so it is enclosed in back-ticks.
#   This can be used to make column names that have spaces (blanks)
#   or non-letter characters acceptable to ggplot2.
#   This version of the function is vectorized with sapply.
ggname <- function(x) {
  if (class(x) != "character") {
    return(x)
  }
  y <- sapply(x, function(s) {
    if (!grepl("^`", s)) {
      s <- paste("`", s, sep="", collapse="")
    }
    if (!grepl("`$", s)) {
      s <- paste(s, "`", sep="", collapse="")
    }
  }
  )
  y 
}

# UI ----

ui <- fluidPage(
  tags$head( # solid black line for aesthetics
    tags$style(HTML("hr {border-top: 1px solid #b8babc;}"))
  ),
  
  headerPanel("ECHO: Finding and Visualizing Circadian Rhythms"), # title  
  
  tabsetPanel( 
    # ui for finding rhythms (ECHO) ----
    tabPanel("Find Rhythms",fluid=TRUE,
             sidebarLayout(
               sidebarPanel( # sidebar with options
                 div(class="header", checked=NA,
                     ("Required fields, with the exception of example data, are marked with *. Required fields in all cases are marked with **.")),
                 
                 div(style="display: inline-block; vertical align: center; width: 350px;",
                     fileInput("data_file", "Upload Data File (.csv) *", accept=c(".csv"))),
                 
                 div(style="display: inline-block; vertical align: top; width: 5px;",
                     actionButton("file_help", icon("question", lib="font-awesome"))),
                 uiOutput("Help_file"),
                 
                 checkboxInput("use_example", "Use Example?", value = FALSE, width = NULL),
                 
                 textInput("project", "Project Title **"),
                 
                 div(style="display: inline-block; width: 80px;",
                     textInput("begin", "Start Time *")),
                 
                 div(style="display: inline-block; width: 80px;",
                     textInput("end", "End Time *")),
                 
                 div(style="display: inline-block; width: 90px;",
                     textInput("resol", "Resolution *")),
                 
                 div(style="display: inline-block; vertical-align: top; width: 5px;",
                     actionButton("time_help", icon("question", lib="font-awesome"))),
                 uiOutput("Help_time"),
                 
                 div(style="display: inline-block;",
                     radioButtons("tied", "Type of replicates? *", choices = NULL, selected = NULL,
                                  inline = TRUE, width = NULL, choiceNames = c("None","Paired","Unpaired"), choiceValues = c("none","TRUE","FALSE"))), 
                 # one replicate defaults to paired replicates (smoothed without average)
                 
                 div(style="display: inline-block; vertical-align:top;  width: 20px;",
                     actionButton("tying_help", icon("question", lib="font-awesome"))),
                 uiOutput("Help_tying"),
                 
                 div(style="display: inline-block; width: 180px;",
                     numericInput("num_reps_res", "Number of replicates? *", 1, min = 1)),
                 
                 div(class="header", checked=NA,
                     tags$b("Looking for rhythms between:")),
                 
                 div(style="display: inline-block; width: 70px;",
                     textInput("low", "(lower:)")),
                 
                 div(style="display: inline-block; width: 70px;",
                     textInput("high", "(upper:)")),
                 
                 div(style="display: inline-block; vertical-align:top;  width: 20px;",
                     actionButton("limit_help", icon("question", lib="font-awesome"))),
                 uiOutput("Help_limit"),
                 
                 div(style="display: inline-block;",
                     checkboxInput("smooth", "Smooth data?", value = FALSE, width = NULL)),
                 div(style="display: inline-block; vertical-align:top;  width: 20px;",
                     actionButton("smooth_help", icon("question", lib="font-awesome"))),
                 uiOutput("Help_smooth"),
                 
                 div(style="display: inline-block;",
                     checkboxInput("rem_unexpr", "Remove unexpressed genes?", value = FALSE, width = NULL)),
                 div(style="display: inline-block; width: 95px;",
                     numericInput("rem_unexpr_amt", "Threshold %?", value = 70,min = 0, max = 100, step = 1, width = NULL)),
                 div(style="display: inline-block; vertical-align:top;  width: 20px;",
                     actionButton("remove_help", icon("question", lib="font-awesome"))),tags$br(),
                 uiOutput("Help_remove"),
                 
                 div(style="display: inline-block;",
                     checkboxInput("is_normal", "Normalize data?", value = FALSE, width = NULL)),
                 div(style="display: inline-block; vertical-align:top;  width: 20px;",
                     actionButton("normal_help", icon("question", lib="font-awesome"))),tags$br(),
                 uiOutput("Help_normal"),
                 
                 div(style="display: inline-block;",
                     checkboxInput("is_de_linear_trend", "Remove linear trends?", value = FALSE, width = NULL)),
                 div(style="display: inline-block; vertical-align:top;  width: 20px;",
                     actionButton("de_linear_trend_help", icon("question", lib="font-awesome"))),tags$br(),
                 uiOutput("Help_de_linear_trend"),
                 
                 checkboxInput("run_jtk", "Run JTK_CYCLE as well?", value = FALSE, width = NULL),
                 
                 # action buttons for running echo code and downloading results
                 actionButton("find_rhythms", "Find Rhythms!"),tags$br(),
                 hr(),
                 downloadButton("downloadECHO", "Download ECHO CSV"),
                 downloadButton("downloadJTK", "Download JTK CSV"),
                 downloadButton("downloadRData", "Download RData")
               )
               ,
               mainPanel(tabsetPanel(
                 # picture
                 # instructions
                 tabPanel("About Finding Rhythms",fluidRow(verbatimTextOutput("finish"),
                   tags$div(class="header", checked = NA,
                            list(
                              HTML('<center><img src="wc1_Neurospora_Replicates_Unsmoothed.PNG" style="width:500px"></center><br>'),
                              tags$p(tags$b("  Welcome to ECHO!")),
                              tags$p("ECHO (Extended Circadian Harmonic Oscillations) is an app designed to find and identify circadian rhythms from your data.
                                     Just upload a .csv of your data and enter its parameters to get started.
                                     For explanations for what specific terms mean, click the '?' buttons to their right.
                                     Recommendations are also made in the '?' text.
                                     Please note: in order to download results, you must open this app in the browser
                                     window before starting to run."),
                              tags$p("If you'd just like to try this out or get a good idea of what data format is necessary,
                                     just check 'Run Example'. This will run the example .csv that came with your download of
                                     ECHO. This example data is fabricated expression data from 2 to 48 hours with 2 hour
                                     resolution and 3 replicates. Random missing data is also included."),
                              tags$p("When you run your data, a progress bar will display in the bottom left corner showing the stage of progress (ECHO, JTK_CYCLE (if run), finish). For ECHO, another progress bar will display in the console window to show directly how far one is in the ECHO progress. Upon finishing, the results will display above. You can then download the ECHO results (.csv), 
                                     the results for use in the visualization tab (.RData), and the JTK_CYCLE results 
                                     (.csv) (if run)."),
                              tags$p(HTML("If you are using the results from this app or want to learn about its methods, please <a href='https://dl.acm.org/citation.cfm?id=3107420&CFID=826084181&CFTOKEN=52238765'>cite us</a>. Additionally, if you are using JTK_CYCLE results, please <a href='https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3119870/'>cite them</a>.")),
                              ("If you run into any errors, please email delosh@rpi.edu with the following (subject 
                               line: ECHO Error):"),tags$br(),
                              "- a short desciption of your problem" ,tags$br(),
                              "- ECHO version number",tags$br(),
                              "- your dataset/file(s)",tags$br(),
                              "- your exact settings for the run (a screenshot will do)",tags$br(),
                              "- your exact error from the console window (a screenshot will do)",tags$br(),tags$br(),
                              "All images created by ECHO using data from:",tags$br(),
                              "Hurley, J. et al. 2014. PNAS. 111 (48) 16995-17002. Analysis of clock-regulated genes in Neurospora reveals widespread posttranscriptional control of metabolic potential. doi:10.1073/pnas.1418963111 ",
                              tags$br(),tags$br(),
                              tags$p("ECHO Version 1.9")
                              ))
                              )),
                 
                 tabPanel("Interpreting Your Result Files",
                          tags$div(class="header",checked=NA,
                                   list(tags$p("Once you run ECHO (and possibly JTK_CYCLE), there are three results files that you have the possibility of downloading: and ECHO CSV, a JTK_CYCLE CSV, and an .RData file. Below, you can find an explanation of the contents of each of the CSV file. If a column has an asterisk next to it, you can check our published paper for more in-depth information."),
                                        
                                     ("ECHO CSV has the following columns:"),tags$br(),
                                     ("- Gene Name: The name of each expression, as entered in the original data file."),tags$br(),
                                     ("- Convergence: Whether or not ECHO's fitting method converged for each expression. This column will also contain markers for 'No Deviation', which means a constant gene, or 'Unexpressed', if removing unexpressed genes was selected during the run."),tags$br(),
                                     ("- Iterations: Amount of iterations for ECHO's fitting method."),tags$br(),
                                     ("- Forcing.Coefficient*: Parameter which states the amount of amplitude change over time in the system."),tags$br(),
                                     ("- Oscillation Type*: States the expression's category based on forcing coefficient (driven, damped, harmonic, overexpressed, repressed)."),tags$br(),
                                     ("- Amplitude*: Parameter describing initial amplitude of expression."),tags$br(),
                                     ("- Radian.Frequency*: Parameter describing frequency of oscillations, in radians."),tags$br(),
                                     ("- Period*: States the time for one complete oscillation, assumed to be in hours."),tags$br(),
                                     ("- Phase Shift*: Parameter describing the amount the oscillator is shifted, in radians."),tags$br(),
                                     ("- Hours Shifted: Desribes the amount the oscillator is shifted in hours, calculated from phase shift and fitted period. This is the time of the first peak of the oscillation, relative to starting time."),tags$br(),
                                     ("- Equilibrium Value*: Parameter describing the center of the oscillator, i.e. the line the expression oscillates around."),tags$br(),
                                     ("- Tau*: Calculated to determine p-values using Kendall's Tau."),tags$br(),
                                     ("- P-Value*: Significance of ECHO fit, unadjusted."),tags$br(),
                                     ("- BH Adj P-Value*: Significance of ECHO fit, adjusted using the Benjamini-Hochberg criterion. Corrects for multiple hypothesis testing."),tags$br(),
                                     ("- BY Adj P-Value*: Significance of ECHO fit, adjusted using the Benhamini-Yekutieli criterion (more stringent). Corrects for multiple hypothesis testing."),tags$br(),
                                     ("- Original TPX.R: Your original data for time point (TP) X, and replicate R."),tags$br(),
                                     ("- Fitted TPX.R: ECHO's fitted data for time point (TP) X, and replicate R."),tags$br(),tags$br(),
                                     
                                     ("JTK_CYCLE CSV has the following columns:"),tags$br(),
                                     ("- Gene Name: The name of each expression, as entered in the original data file."),tags$br(),
                                     ("- BY.Q: Significance of JTK fit (ADJ.P), adjusted using the Benhamini-Yekutieli criterion (more stringent). Corrects for multiple hypothesis testing."),tags$br(),
                                     ("- BH.Q: Significance of JTK fit (ADJ.P), adjusted using the Benjamini-Hochberg criterion. Corrects for multiple hypothesis testing."),tags$br(),
                                     ("- ADJ.P: Significance of JTK fit, adjusted internally."),tags$br(),
                                     ("- PER: Period of one complete oscillation."),tags$br(),
                                     ("- LAG: Hours delayed for expression."),tags$br(),
                                     ("- AMP: Amplitude."),tags$br(),
                                     ("- Original data, complete with original names, then follows."),tags$br(),tags$br(),
                                     
                                     tags$p("The .RData file contains a series of R objects that are necessary for the automatic visualizations on the next tab. These objects include the ECHO and JTK output and user input information.")
                                   )))
                 # tabPanel("Gene List")
                 #tabPanel("Gene Clustering",plotOutput("plot_genes", height = "650px")),
                 #tabPanel("Silhouette Plot",plotOutput("plot_sil", height = "650px"))
                                     )
                                     )
                 )
             )
    ,
    # ui for visualization ----
    tabPanel("Visualization", fluid=TRUE, sidebarLayout(
      sidebarPanel(#sidebar for decisions
        div(class="header", checked=NA,
            tags$b("Required fields in all cases are marked with *.")),
        
        fileInput("file","Upload Results File (.RData) *", accept=c(".RData", ".Rds")),
        
        div(style="display: inline-block;",
            checkboxInput("no_restrict_gamma", "Include repressed/overexpressed genes?", value = FALSE, width = NULL)),
        div(style="display: inline-block; vertical-align:top;  width: 20px;",
            actionButton("restrict_help", icon("question", lib="font-awesome"))),
        uiOutput("Help_restrict"),
        
        selectInput("viz","Choose Visualization Method",c("Venn Diagram","Gene Expression Graph","Gene Expression with Replicates","Heat Map","Parameter Density Graph (PDG)")),
        
        textInput("gene_name","Enter Gene to View (for Gene Expression with/without Replicates)"),
        
        selectInput("pval_cat","Enter P-Value Adjustment to View (for PDG, Heat Maps, Gene Lists)",c("BH Adj P-Value", "BY Adj P-Value", "P-Value")),
        
        numericInput("pval_cutoff","Enter P-Value Significance Cutoff (for PDG, Heat Maps, Gene Lists)",min = 0, max = 1, step = .01, value = 0.05),
        
        selectInput("coeff","Choose Coefficient to View (for PDG)",c("Forcing.Coefficient","Amplitude","Radian.Frequency","Period","Phase Shift","Hours Shifted", "Equilibrium Value", "Tau", "P-Value", "BH Adj P-Value", "BY Adj P-Value")),
        
        div(style="display: inline-block;",
            selectInput("subset_look","Subset of Data to View (for PDG and Gene Lists)",c("None","ECHO","JTK_CYCLE","Driven","Damped","Harmonic","Repressed","Overexpressed","Nonconverged","Nonstarter","No Deviation","Both ECHO and JTK","Only JTK","Only ECHO","Neither ECHO and JTK"))),
        
        div(style="display: inline-block; vertical-align:top;  width: 20px;",
            actionButton("subset_help", icon("question", lib="font-awesome"))),
        uiOutput("Help_subset"),
        
        div(style="display: inline-block;",
            selectInput("heat_subset_look","Subset of Data to View (for Heat Maps)",c("None","ECHO","JTK_CYCLE","Driven","Damped","Harmonic","Repressed","Overexpressed","Both ECHO and JTK","Only JTK","Only ECHO","Neither ECHO and JTK"))),
        div(style="display: inline-block; vertical-align:top;  width: 20px;",
            actionButton("heat_subset_help", icon("question", lib="font-awesome"))),
        uiOutput("Help_heat_subset"),
        
        div(style="display: inline-block;",
            textInput("heat_subset_rep","Enter Replicate to View (for Heat Maps)", value = "all")),
        div(style="display: inline-block; vertical-align:top;  width: 20px;",
            actionButton("heat_subset_rep_help", icon("question", lib="font-awesome"))),
        uiOutput("Help_heat_subset_rep"),
        
        # action buttons for visualization and downloading results
        actionButton("go","Update"),
        hr(),
        downloadButton('downloadData', 'Download Gene List'),
        downloadButton('downloadPlot', 'Download Plot *')
      ),
      # where the output happens
      mainPanel(
        tabsetPanel(
          tabPanel("Visualization",fluidRow(plotOutput("plot_viz", height = "650px"),
                                            verbatimTextOutput("text"))),
          tabPanel("Gene List", dataTableOutput("table")),
          tabPanel("About Visualizing Results",
                   tags$div(class="header", checked = NA,
                            list(
                              tags$p("Once you have run your data, you can visualize and explore your results using the .RData
                                     file available for download after your run. (*: Note that .PNG download is not available for Venn Diagrams and Heat Maps, but right-click copy and paste works well.) To update visualizations or gene lists, select options from each of the drop down menus to the left and press Update. There are five types of visualizations/
                                     explorations currently available:"),
                              HTML('<center>'),tags$b("Gene Lists:"),HTML('</center>'),
                              tags$p("A subset of ECHO parameter results for exploration, which can be sorted by any parameter.These subsets are specified by the P-Value Type, P-Value Cutoff, and
                                     Subset of Data to View."),
                              
                              HTML('<center><img src="venn_diagram_by_adj_include_overdamped_Neurospora_Replicates_Unsmoothed.PNG" style="width:200px"></center><br>'),
                              HTML('<center>'),tags$b("Venn Diagram:"),HTML('</center>'),
                              tags$p("A Venn Diagram of ECHO results compared to JTK_CYCLE results (if available). A text 
                                     summary also appears below the Venn Diagram. If JTK_CYCLE results are not available,
                                     a text summary of ECHO results appear below a blank plot. In addition, user inputs are displayed below the text summary."),
                              
                              HTML('<center><img src="heat_map_Neurospora_Replicates_Unsmoothed.PNG" style="width:200px"></center><br>'),
                              HTML('<center>'),tags$b("Heat Map:"),HTML('</center>'),
                              tags$p("A heat map of the expression data for indicated subset of results (not including constant, unexpressed, or overly noisy data). Each expression appears on the y axis and time on the x axis. Expressions are sorted from top to bottom by phase and normalized to be in the range [-1,1]. If 'all' is entered for heat maps, an average of all expressions for replicates is shown. Otherwise, a replicate number is entered to display the expressions for a specific replicate."),
                              
                              HTML('<center><img src="wc1_gene_expression_wo_rep_Neurospora_Replicates_Unsmoothed.PNG" style="width:300px"></center><br>'),
                              HTML('<center>'),tags$b("Gene Expressions:"),HTML('</center>'),
                              tags$p("A plot of original data and fitted data. A summary of the parameters for that gene appears
                                     below the plot. If one replicate, the original data is plotted as
                                     a blank line and the fitted data is plotted as a green line. If multiple replicates,
                                     the original data is plotted with shading between the max and min at each data point and
                                     the fitted data is plotted as a black line."),
                              
                              HTML('<center><img src="wc1_Neurospora_Replicates_Unsmoothed.PNG" style="width:300px"></center><br>'),
                              HTML('<center>'),tags$b("Gene Expressions with Replicates:"),HTML('</center>'),
                              tags$p("If multiple replicates, the original data is plotted with shading between the max and min
                                     at each data point, as well as the original data for each replicate (up to 8 replicates).
                                     The fitted data is plotted as a black line. A summary of the parameters for that gene 
                                     appears below the plot. Replicates are colored in the following manner (1 to 8): red, blue, green, yellow, purple, pink, orange, magenta."),
                              
                              HTML('<center><img src="gamma_density_Neurospora_Replicates_Unsmoothed.PNG" style="width:300px"></center><br>'),
                              HTML('<center>'),tags$b("Parameter Density Graph (PDG):"),HTML('</center>'),
                              tags$p("A density plot (a smoothed histogram, or frequency graph) for specific parameters for a
                                     subset of the data. These subsets are specified by the P-Value Type, P-Value Cutoff, and
                                     Subset of Data to View."))
                              )
                              )
        )
      )
    )
    ) 
    )
)


# Server ----

server <- function(input,output){ # aka the code behind the results
  # all functionality in find rhythms tab
  {
    # help text ----
    
    output$Help_file=renderUI({ # csv help
      if(input$file_help%%2){
        helpText("The data must be in a .csv file with the following format: first row is column labels, first column has gene labels/names, and all other columns have expression data. This expression data must be ordered by time point then by replicate, and must have evenly spaced time points. Any missing data must have cells left blank.")
      }
      else{
        return()
      }
    })
    output$Help_time=renderUI({ # time inputs help
      if(input$time_help%%2){
        helpText("These numbers indicate the beginning, end, and resolution (difference between time points) of the time points in your data. For example, data that begins at 2 hours, ends at 48 hours, with a resolution of 2 hours.")
      }
      else{
        return()
      }
    })
    output$Help_smooth=renderUI({ # data smoothing help
      if(input$smooth_help%%2){
        helpText("Indicates whether data should be smoothed, which also depends on type of data. If checked, the data is weighted smoothed over a rolling window of 3 points, with each of the points having weights of 1,2,1 respectively. If paired data, each replicate is smoothed independently. If unpaired data, each time point is smoothed by itself, centered, and the average expression per time point on either side. Smoothed data will be returned in output files. Note: this will increase running time.")
      }
      else{
        return()
      }
    })
    output$Help_limit=renderUI({ # upper and lower limits for rhythms help
      if(input$limit_help%%2){
        helpText("Upper and lower limits to look for rhythms. For example, one could look for hours between 20 and 26. Note 1: if running, JTK_CYCLE, these must be multiples of resolution. Note 2: if limits left blank, rhythms of any length within timecourse will be considered, and JTK_CYCLE cannot run.")
      }
      else{
        return()
      }
    })
    output$Help_remove=renderUI({ #remove unexpressed genes help
      if(input$remove_help%%2){
        helpText("If checked, marks and does not fit genes that have low detection/expression in the dataset. A gene is considered adequately expressed for modeling if a nonzero value is available for at least the specified threshold percentage of the total time points. A 70% threshold is recommended. By default, genes with a zero value at all time points are marked and not fit to allow appropriate calculation of multiple hypothesis corrections. Recommended.")
      }
      else{
        return()
      }
    })
    output$Help_tying=renderUI({
      if((input$tying_help)%%2){ # paired or unpaired help
        helpText("Indicates whether replicates are paired (replicates are related within a time course) or unpaired (replicates are not related to any specific time course).")
      }
      else{
        return()
      }
    })
    output$Help_restrict=renderUI({
      if(input$restrict_help%%2){ # forcing coefficient help
        helpText("If checked, repressed (strong decreases in amplitude over time) and overexpressed (strong increases in amplitude over time) genes will be included in circadian calculations. Not recommended.")
      }
      else{
        return()
      }
    })
    output$Help_subset=renderUI({
      if(input$subset_help%%2){ # forcing coefficient help
        helpText("Identifies which subset of data to view for specified visualizations. For more information regarding what terms such as 'Both ECHO and JTK' and 'Only JTK' mean, please view Venn Diagram summary. Note 1: Forcing coefficient values are restricted by ECHO significance speficied. Note 2: Some subsets are not available if JTK_CYCLE results are not available. ")
      }
      else{
        return()
      }
    })
    output$Help_normal=renderUI({
      if(input$normal_help%%2){ # forcing coefficient help
        helpText("Normalizes data by row using the normal distribution (subtract each row by row mean and divide by row standard deviation). Normalized data is returned in results, rather than original data. Normalized data will be returned in output files. Recommended for un-normalized data. If you have normalized your data in any way, this option is strongly not recommended.")
      }
      else{
        return()
      }
    })
    output$Help_heat_subset=renderUI({
      if(input$heat_subset_help%%2){ # forcing coefficient help
        helpText("Identifies which subset of data to view for heat map visualizations. For more information regarding what terms such as 'Both ECHO and JTK' and 'Only JTK' mean, please view Venn Diagram summary. Note 1: Forcing coefficient values are restricted by ECHO significance speficied. Note 2: Some subsets are not available if JTK_CYCLE results are not available. ")
      }
      else{
        return()
      }
    })
    output$Help_heat_subset_rep=renderUI({
      if(input$heat_subset_rep_help%%2){ # forcing coefficient help
        helpText("Identifies which replicate of data to view for heat map visualizations. If 'all' is entered, average of all replicates will be displayed. Otherwise, enter number of replicate to display (1 through total number of replicates). Note: visualization of specific replicates only recommended for paired data.")
      }
      else{
        return()
      }
    })
    output$Help_de_linear_trend=renderUI({
      if(input$de_linear_trend_help%%2){ # forcing coefficient help
        helpText("If checked, removes linear trend, or baseline, from data. For paired replicates, the linear trend is computed and removed from each replicate. For unpaired data, the linear trend is computed and removed from all replicates together.")
      }
      else{
        return()
      }
    })
    
    
    # run echo with maybe jtk ----
    observeEvent(input$find_rhythms, {
      withProgress(message = 'Finding Rhythms!', value = 0, {
        
        if (input$run_jtk && input$high != "" && input$low != ""){
          incProgress(1/3,detail = paste("Running ECHO. Started on:",Sys.time()))
        }
        else{
          incProgress(1/2,detail = paste("Running ECHO. Started on:",Sys.time()))
        }
        source("damped_oscillator_master.R",local = TRUE) # load the package
        
        # echo run ----
        # upload data
        print("Beginning to run ECHO...")
        print(paste("Started on:",Sys.time()))
        if (input$use_example){
          # load example datasets and images
          genes <- read.csv("DataExample.csv", header = TRUE, stringsAsFactors = FALSE)
          
          # creating times sequence used for the genes
          num_reps <- 3
          begin <- 2 # beginning
          end <- 48 # end
          resol <- 2 # resolution of data
          timen <- seq(begin,end,resol) # the times for cicadian rhythms
          tied <- FALSE # not paired
        }
        else{ # the user's own dataset
          genes <- read.csv(input$data_file$datapath, header = TRUE)
          num_reps <- as.numeric(input$num_reps_res) # number of replicates
          
          # creating times sequence used for the genes
          begin <- as.numeric(sapply(input$begin, function(x) eval(parse(text=x)))) # beginning
          end <- as.numeric(sapply(input$end, function(x) eval(parse(text=x)))) # end
          resol <- as.numeric(sapply(input$resol, function(x) eval(parse(text=x)))) # resolution of data
          timen <- seq(begin,end,resol) # the times for cicadian rhythms
          if (input$tied=="none"){ # one replicate, default to true paired-ness
            tied <- TRUE
          } else {# more than one replicate
            tied <- as.logical(input$tied) # the type of replicate
          }
          
        }
        # if genes only has one row, add another constant row to alleviate automatic vectorization
        if (nrow(genes)==1){
          # creating a constant row and adding it to genes
          add_row <- data.frame(matrix(0L, 1, ncol(genes)))
          add_row[1,1] <- "not considering"
          colnames(add_row) <- colnames(genes)
          genes <- rbind(genes,add_row)
          add_one <- TRUE # marker for appropriate displays for progress bar
        }
        else{
          add_one <- FALSE # marker for appropriate displays for progress bar
        }
        
        # run all genes:
        
        rem_unexpr <- input$rem_unexpr # indicator for removing unexpressed genes
        rem_unexpr_amt <- (input$rem_unexpr_amt)/100 # threshold for removing unexpressed genes, converted to a decimal
        # if yes, check for genes that are unexpressed before preprocessing
        if (rem_unexpr){
          rem_unexpr_vect <- genes_unexpressed_all(rem_unexpr_amt)
        } else{
          rem_unexpr_vect <- rep(FALSE,nrow(genes))
        }
        
        # normalize and store original data
        if (input$is_normal){
          norm_list <- normalize_all()
          genes <- norm_list$dat
        }
        
        # remove baseline
        if (input$is_de_linear_trend){
          genes <- de_linear_trend_all(timen,num_reps,tied)
        }
        
        # getting average data, for more than one replicate
        avg_genes <- avg_all_rep(num_reps)
        
        # figuring out whether smoothing is wanted
        if (!input$smooth){ # no smoothing
          is_smooth <- FALSE
          is_weighted <- FALSE
        }
        else{ # yes smoothing, weighted or unweighted
          is_smooth <- TRUE
          is_weighted <- TRUE # only offer weighted smoothing
          
          # create smoothed matrix, reassign genes
          if(is_smooth){ # smooth the data, if requested
            if (tied){ # if paired replicates
              genes <- smoothing_all_tied(is_weighted, num_reps)
            }
            else{ # if unpaired replicates
              genes <- smoothing_all_untied(is_weighted, num_reps)
            }
          }
        }
        
        if (is_smooth){ # average should correspond to smoothed data
          # not unsmoothed data
          avg_genes <- avg_all_rep(num_reps)
        }
        
        
        # figuring out whether a range is wanted, adjusting accordingly
        if (input$low ==""){ # empty low input, adjust to time series
          if (resol >= 1){
            low <- 2*pi/resol
            low_end <- resol
          }
          else{ # if the begining is >=0, smallest period available is 1
            low <- 2*pi/1
            low_end <- 1
          }
        } else{ # there is a low input
          low_input <- as.numeric(sapply(input$low, function(x) eval(parse(text=x))))
          low <- 2*pi/low_input
          low_end <- low_input
        }
        if (input$high ==""){ # empty high input, adjust to time series
          high <- 2*pi/(resol*length(timen))
          high_end <- (resol*length(timen))
        } else{ # there is a high input
          high_input <- as.numeric(sapply(input$high, function(x) eval(parse(text=x))))
          high <- 2*pi/high_input
          high_end <- high_input
        }
        
        start.time <- Sys.time() # begin counting time
        
        # if more than one replicate or requested, an exact distribution is needed
        #if (num_reps > 1 ){ 
          # create exact distribution for pvalues
          jtklist <- jtkdist(length(timen), reps = num_reps) 
          jtk.alt <- list() # preallocate pvalue distribution for missing data
        #}
        
        # prepare for parallelism
        cores = detectCores() # dectect how many processors
        cl <- makeCluster(cores[1]-1) # not to overload your computer, need one for OS
        registerDoSNOW(cl)
        
        # making a progress bar
        if (add_one){
          print(paste("Percentage finished out of",nrow(genes)-1,"expression:"))
        } else {
          print(paste("Percentage finished out of",nrow(genes),"expressions:"))
        }
        pb <- txtProgressBar(max = nrow(genes), style = 3)
        progress <- function(n) setTxtProgressBar(pb, n)
        opts <- list(progress = progress)
        
        # where we put the result
        total_results <- foreach (i=1:nrow(genes), .combine = rbind, .packages='minpack.lm',.options.snow = opts) %dopar% {
          calculate_param(i, timen, resol, num_reps, tied = tied, is_smooth = is_smooth, is_weighted = is_weighted,low = low,high = high,rem_unexpr = rem_unexpr, rem_unexpr_amt = rem_unexpr_amt, jtklist)
        }
        close(pb)
        
        stopCluster(cl) # stop using the clusters
        
        # renaming columns of the final results
        colnames(total_results) <- c("Gene Name","Convergence","Iterations","Forcing.Coefficient","Oscillation Type","Amplitude","Radian.Frequency","Period","Phase Shift","Hours Shifted","Equilibrium Value", "Tau", "P-Value", paste(rep("Original TP",length(rep(timen, each = num_reps))),rep(timen, each = num_reps),rep(".",length(rep(timen, each = num_reps))),rep(c(1:num_reps), length(timen)),sep=""), paste(rep("Fitted TP",length(timen)),timen,sep=""))
        
        # remove the fake row I added if there is only one gene
        if (add_one){
          total_results <- total_results[-nrow(total_results),]
        }
        
        adjusted_p_val_us <- p.adjust(unlist(total_results$`P-Value`), method = "BH") # benjamini-hochberg adjust p-values
        total_results <- cbind(total_results[,c(1:13)],`BH Adj P-Value` = adjusted_p_val_us, total_results[,c(14:ncol(total_results))]) # assign to data frame
        
        # adding the benjamini-hochberg-yekutieli p-value adjustment
        total_results <- cbind(total_results[,c(1:14)],`BY Adj P-Value` = p.adjust(unlist(total_results$`P-Value`), method = "BY"), total_results[,c(15:ncol(total_results))])
        
        
        # time measured - output
        end.time <- Sys.time()
        time.taken.echo <- difftime(end.time,start.time,units = "mins")
        
        # assigning everything to the global environment in order to save them
        total_results <<- total_results
        timen <<- timen
        num_reps <<- num_reps
        low <<- low
        high <<- high
        low_end <<- low_end
        high_end <<- high_end
        
        print(paste("Ended on:",Sys.time())) # show ending time in console window
        
        user_input <<- list("ECHO_end_date" = Sys.time(),
                            "begin"=begin,
                            "end"=end,
                            "resol"=resol,
                            "tied"=tied,
                            "is_smooth"=input$smooth,
                            "rem_unexpr"=rem_unexpr,
                            "rem_unexpr_amt"=rem_unexpr_amt,
                            "is_normal"=input$is_normal,
                            "is_de_linear_trend"=input$is_de_linear_trend,
                            "run_jtk"=input$run_jtk,
                            "v_num"=1.9) # VERSION NUMBER
        
        # jtk run -----
        
        # if possible or desired, run JTK_CYCLE
        if ((input$run_jtk && input$high != "" && input$low != "") || (input$run_jtk && input$use_example)){
          incProgress(1/3,detail = paste("Running JTK. Started on:",Sys.time()))
          print("Beginning to run JTK...")
          print(paste("Started on:",Sys.time()))
          
          source("JTK_CYCLEv3.1.R", local = TRUE) # load JTK_CYCLE script
          
          start.time <- Sys.time() # begin counting time
          options(stringsAsFactors=FALSE)
          # if smoothed, the smoothed data appears in the total_results dataframe
          if(input$smooth){
            data <- total_results[,c(16:(15+(length(timen)*num_reps)))]
            data <- cbind(genes[,1],data)
            colnames(data)[1] <- "Gene.Name"
          } else{ # otherwise, just use the inputted data
            data <- genes 
          }
          # remove unexpressed genes if desired
          if (rem_unexpr){
            # logical of which entries to remove
            rem_choices <- !(total_results$Convergence=="Unexpressed" | total_results$Convergence=="No Deviation")
            rem_choices[is.na(rem_choices)] <- TRUE
            
            data <- data[rem_choices,]
          } else { # always remove constant genes
            # logical of which entries to remove
            rem_choices <- !(total_results$Convergence=="No Deviation")
            rem_choices[is.na(rem_choices)] <- TRUE
            
            data <- data[rem_choices,]
          }
          annot <- data.frame("Gene.Name" = data[,1],stringsAsFactors = FALSE)
          
          data <- data[,-1]
          jtkdist(length(timen), num_reps) # total time points, replicates per time point
          
          # looking for rhythms between high-low hours (i.e. between low/resol and high/resol time points per cycle).
          if (!input$use_example){
            periods <- ((as.numeric(input$low)/resol):(as.numeric(input$high)/resol))
          } else {
            periods <- ((20/resol):(26/resol))
            
          }
          jtk.init(periods,resol)  # resol is the number of hours between time points
          
          { # run jtk
            res <- apply(data,1,function(z) {
              jtkx(z)
              c(JTK.ADJP,JTK.PERIOD,JTK.LAG,JTK.AMP)
            })
            res <- as.data.frame(t(res))
            bhq <- p.adjust(unlist(res[,1]),"BH")
            byq <- p.adjust(unlist(res[,1]),"BY")
            res <- cbind(byq,bhq,res)
            colnames(res) <- c("BY.Q","BH.Q","ADJ.P","PER","LAG","AMP")
            JTK_results <- cbind(annot,res,data)
            JTK_results <- JTK_results[order(res$ADJ.P,-res$AMP),]
          }
          # now to add unexamined genes (unexpressed, constant, etc.)
          missing_genes <- setdiff(total_results[,1],JTK_results[,1])
          if (length(missing_genes) > 0){
            blank <- as.data.frame(matrix(NA,length(missing_genes),ncol(JTK_results)))
            blank[,1] <- missing_genes
            colnames(blank) <- colnames(JTK_results)
            JTK_results <- rbind(JTK_results,blank)
          }
          
          end.time <- Sys.time() # end counting time
          time.taken.jtk <- difftime(end.time,start.time,units = "mins")
          
          JTK_results <<- JTK_results
          
          print(paste("Ended on:",Sys.time())) # show ending time in console window
        }
        
        # download results ----
        
        output$finish <- renderPrint({cat("Done!\n")
          cat(paste("ECHO Time:",time.taken.echo,"mins\n"))
          if ((input$run_jtk && input$high != "" && input$low != "") | (input$run_jtk && input$use_example)){
            cat(paste("JTK Time:",time.taken.jtk,"mins"))
          }
        })
        
        # function to download gene list subset as a csv
        # output$downloadZip <- downloadHandler(
        #   filename = function() { paste(input$project,'.zip', sep='') },
        #   content = function(file) {
        #     #go to a temp dir to avoid permission issues
        #     owd <- setwd(tempdir())
        #     on.exit(setwd(owd))
        #     files <- NULL;
        #     
        #     #write each output to a csv file, save the name
        #     fileName <- paste('ECHO_',input$project,'.csv', sep='')
        #     write.csv(total_results,fileName, row.names = F)
        #     files <- c(fileName,files)
        #     if (input$run_jtk && input$high != "" && input$low != ""){
        #       fileName <- paste('JTK_',input$project,'.csv', sep='')
        #       write.csv(JTK_results,fileName,sep = ';', row.names = F)
        #       files <- c(fileName,files)
        #     }
        #     fileName <- paste(input$project,'.RData', sep='')
        #     save.image(fileName)
        #     files <- c(fileName,files)
        #     
        #     #create the zip file
        #     zip(file,files)
        #   }
        # )
        
        output$downloadECHO <- downloadHandler( # ECHO results
          filename = function() { paste('ECHO_',input$project,'.csv', sep='') },
          content = function(file) {
            fileName <- paste('ECHO_',input$project,'.csv', sep='')
            write.csv(total_results,file, row.names = F)
          }
        )
        if ((input$run_jtk && input$high != "" && input$low != "") | (input$run_jtk && input$use_example)){
          output$downloadJTK <- downloadHandler( # JTK_CYCLE results
            filename = function() { paste('JTK_',input$project,'.csv', sep='') },
            content = function(file) {
              fileName <- paste('JTK_',input$project,'.csv', sep='')
              write.csv(JTK_results,file, row.names = F)
            }
          )
        }
        output$downloadRData <- downloadHandler( # Visualization results
          filename = function() { paste(input$project,'.RData', sep='') },
          content = function(file) {
            if ((input$run_jtk && input$high != "" && input$low != "") | (input$run_jtk && input$use_example)){
              save(file=file,list=c("JTK_results","total_results","num_reps","low","high","timen","user_input","high_end","low_end"))
            } else {
              save(file=file,list=c("total_results","num_reps","low","high","timen","user_input","high_end","low_end"),file)
            }
            
          }
        )
        
        if (input$run_jtk && input$high != "" && input$low != ""){
          incProgress(1/3,detail = paste("Finished!",Sys.time()))
        }
        else{
          incProgress(1/2,detail = paste("Finished!",Sys.time()))
        }
      })
      
    })
    
  }
  
  
  # all functionality in vizualization tab
  # output depends on the "visualize!" button on the UI
  observeEvent(input$go, {
    load(input$file$datapath) # load .RData file -- contains evenronment produced from running our oscillator code
    
    # call upon user input
    is_jtk <- user_input$run_jtk
    
    #low_end <- floor(2*pi/low)
    #high_end <- ceiling(2*pi/high)
    
    circ_us <- (total_results$Period<=high_end & total_results$Period >=low_end & total_results[,input$pval_cat]<as.numeric(input$pval_cutoff)) # circadian genes
    circ_us[is.na(circ_us)] <- FALSE # na's are false
    if (!input$no_restrict_gamma){
      circ_us[total_results$`Oscillation Type` == "Repressed" | total_results$`Oscillation Type` == "Overexpressed"] <- FALSE # restricting gamma for circadian genes
    }
    
    
    if (is_jtk){
      # getting JTK circadian genes
      JTK_results[,1] <- factor(JTK_results[,1], levels = unique(total_results$`Gene Name`)) # match JTK rows to our rows
      JTK_results <- JTK_results[order(JTK_results[,1]),] # reorder JTK rows to be the same as our result's rows
      
      if (input$pval_cat == "P-Value"){
        circ_jtk <- (JTK_results$ADJ.P<as.numeric(input$pval_cutoff) & JTK_results$PER >=low_end & JTK_results$PER <= high_end) #jtk's circadian genes
      }
      else if(input$pval_cat == "BH Adj P-Value"){
        circ_jtk <- (JTK_results$BH.Q<as.numeric(input$pval_cutoff) & JTK_results$PER >=low_end & JTK_results$PER <= high_end) #jtk's circadian genes
      }
      else{
        circ_jtk <- (JTK_results$BY.Q<as.numeric(input$pval_cutoff) & JTK_results$PER >=low_end & JTK_results$PER <= high_end) #jtk's circadian genes
      }
      circ_jtk[is.na(circ_jtk)] <- FALSE # na's are false
      # confusion matrix of circadian genes
      #           jtk
      #        yes     no
      # us yes both    ours
      #    no  theirs  diff
      both <- (circ_us & circ_jtk)
      theirs <- (!circ_us & circ_jtk)
      ours <- (circ_us & !circ_jtk)
      diff <- (!circ_us & !circ_jtk)
    }
    else{
      circ_jtk <- NA
      both <- NA
      theirs <- NA
      ours <- NA
      diff <- NA
    }
    
    # logicals
    nas_found <- (is.na(total_results$Forcing.Coefficient)) # amount of nas
    nonconv <- ((total_results$Convergence==0)) # amount of results that didn't converge
    nodev <- (total_results$Convergence=="No Deviation") # results with no standard deviation deviation
    nodev[is.na(nodev)] <- FALSE # na's are false
    noexpr <- (total_results$Convergence=="Unexpressed") # results with no standard deviation deviation
    noexpr[is.na(noexpr)] <- FALSE # na's are false
    # only showing significant, by our measures
    driven <- total_results$`Oscillation Type` == "Driven" & circ_us
    damped <- total_results$`Oscillation Type` == "Damped" & circ_us
    harmonic <- total_results$`Oscillation Type` == "Harmonic" & circ_us
    
    # adjust significance for over damped/driven genes
    circ_over <-  (total_results$Period<=high_end & total_results$Period >=low_end & total_results[,input$pval_cat]<as.numeric(input$pval_cutoff)) # circadian genes
    circ_over[is.na(circ_over)] <- FALSE # na's are false
    overexpressed <- total_results$`Oscillation Type` == "Overexpressed" & circ_over
    repressed <- total_results$`Oscillation Type` == "Repressed" & circ_over
    
    # gene list output based on specified subset
    { 
      if (!is_jtk && (input$subset_look == "JTK" || input$subset_look == "Both ECHO and JTK" || input$subset_look == "Only JTK" || input$subset_look == "Only ECHO" || input$subset_look == "Neither ECHO and JTK")){
        df<- data.frame()
      }
      else if(input$subset_look == "None"){
        if (is_jtk){
          df<-as.data.frame(cbind(total_results[,1:15],JTK_results[,2:4]))
        } else {
          df<-as.data.frame(total_results[,1:15])
        }
      }
      else if(input$subset_look == "JTK_CYCLE"){
        if(is_jtk){
          df <-as.data.frame(cbind(total_results[circ_jtk,1:15],JTK_results[circ_jtk,2:4]))
        } else {
          df<-as.data.frame(total_results[circ_jtk,1:15])
        }
        
      }
      else if(input$subset_look == "ECHO"){
        if(is_jtk){
          df<-as.data.frame(cbind(total_results[circ_us,1:15],JTK_results[circ_us,2:4]))
        } else {
          df<-as.data.frame(total_results[circ_us,1:15])
        }
      }
      else if(input$subset_look == "Damped"){
        if(is_jtk){
          df<-as.data.frame(cbind(total_results[damped,1:15],JTK_results[damped,2:4]))
        } else {
          df<-as.data.frame(total_results[damped,1:15])
        }
      }
      else if(input$subset_look == "Driven"){
        if(is_jtk){
          df<-as.data.frame(cbind(total_results[driven,1:15],JTK_results[driven,2:4]))
        } else {
          df<-as.data.frame(total_results[driven,1:15])
        }
      }
      else if(input$subset_look == "Harmonic"){
        if(is_jtk){
          df<-as.data.frame(cbind(total_results[harmonic,1:15],JTK_results[harmonic,2:4]))
        } else {
          df<-as.data.frame(total_results[harmonic,1:15])
        }
      }
      else if(input$subset_look == "Repressed"){
        if(is_jtk){
          df<-as.data.frame(cbind(total_results[repressed,1:15],JTK_results[repressed,2:4]))
        } else {
          df<-as.data.frame(total_results[repressed,1:15])
        }
        
      }
      else if(input$subset_look == "Overexpressed"){
        if(is_jtk){
          df<-as.data.frame(cbind(total_results[overexpressed,1:15],JTK_results[overexpressed,2:4]))
        } else {
          df<-as.data.frame(total_results[overexpressed,1:15])
        }
      }
      else if(input$subset_look == "Nonconverged"){
        if(is_jtk){
          df<-as.data.frame(cbind(total_results[nonconv,1:15],JTK_results[nonconv,2:4]))
        } else {
          df<-as.data.frame(total_results[nonconv,1:15])
        }
      }
      else if(input$subset_look == "Nonstarter"){
        if(is_jtk){
          df<-as.data.frame(cbind(total_results[nas_found-nodev,1:15],JTK_results[nas_found-nodev,2:4]))
        } else {
          df<-as.data.frame(total_results[nas_found-nodev,1:15])
        }
      }
      else if(input$subset_look == "No Deviation"){
        if(is_jtk){
          df<-as.data.frame(cbind(total_results[nodev,1:15],JTK_results[nodev,2:4]))
        } else {
          df<-as.data.frame(total_results[nodev,1:15])
        }
      }
      else if(input$subset_look == "Both ECHO and JTK"){
        if(is_jtk){
          df<-as.data.frame(cbind(total_results[both,1:15],JTK_results[both,2:4]))
        } else {
          df<-as.data.frame(total_results[both,1:15])
        }
      }
      else if(input$subset_look == "Only JTK"){
        if(is_jtk){
          df<-as.data.frame(cbind(total_results[theirs,1:15],JTK_results[theirs,2:4]))
        } else {
          df<-as.data.frame(total_results[theirs,1:15])
        }
      }
      else if(input$subset_look == "Only ECHO"){
        if(is_jtk){
          df<-as.data.frame(cbind(total_results[ours,1:15],JTK_results[ours,2:4]))
        } else {
          df<-as.data.frame(total_results[ours,1:15])
        }
      }
      else if(input$subset_look == "Neither ECHO and JTK"){
        if(is_jtk){
          df<-as.data.frame(cbind(total_results[diff,1:15],JTK_results[diff,2:4]))
        } else {
          df<-as.data.frame(total_results[diff,1:15])
        }
      }
      
      output$table <- renderDataTable({
        df
      })
    }
    
    # function to download gene list subset as a csv
    output$downloadData <- downloadHandler(
      filename = function() { paste(input$subset_look, '_subset.csv', sep='') },
      content = function(file) {
        write.csv(df, file, row.names = FALSE)
      }
    )
    
    if (num_reps > 1){ # create a data frame for gene visualizations
      if(sum(total_results$`Gene Name`==input$gene_name)!=0){
        rep_genes <- total_results[total_results$'Gene Name'==input$gene_name,16:(15+(length(timen)*num_reps))]
        
        ribbon.df <- data.frame(matrix(ncol = 4+num_reps, nrow = length(timen)))
        colnames(ribbon.df) <- c("Times","Fit","Min","Max", paste(rep("Rep",num_reps),c(1:num_reps), sep=".")) # assigning column names
        ribbon.df$Times <- timen
        ribbon.df$Fit <- t(total_results[total_results$'Gene Name'==input$gene_name,c((16+(length(timen)*num_reps)):ncol(total_results))]) # assigning the fit
        ribbon.df$Min <- sapply(seq(1,ncol(rep_genes), by = num_reps), function(x) min(unlist(rep_genes[,c(x:(num_reps-1+x))]), na.rm = TRUE)) # getting min values of replicates
        ribbon.df$Max <- sapply(seq(1,ncol(rep_genes), by = num_reps), function(x) max(unlist(rep_genes[,c(x:(num_reps-1+x))]), na.rm = TRUE)) # getting max values of replicates
        for (i in 1:num_reps){ # assign each of the replicates
          ribbon.df[,4+i] <- t(rep_genes[,seq(i,ncol(rep_genes),by=num_reps)])
        }
      }
    }
    
    
    # Visualization Ouptuts ----
    
    if(input$viz == "Gene Expression Graph"){
      if(sum(total_results$`Gene Name`==input$gene_name)==0){ # if the gene name inputted isn't in set
        output$text <- renderPrint({cat("")}) # no text
        plot_viz<-ggplot() # no visualization
      }
      else{
        # summary in the UI is all of the statistics vital to the gene (name, p-value, etc.)
        output$text <- renderPrint({
          cat(paste("Gene Name:",total_results$`Gene Name`[total_results$`Gene Name`==input$gene_name],"\n"))
          cat(paste("Convergence:", total_results$Convergence[total_results$`Gene Name`==input$gene_name],"\n"))
          cat(paste("Iterations:",total_results$Iterations[total_results$`Gene Name`==input$gene_name],"\n"))
          cat(paste("Forcing Coefficient:", total_results$Forcing.Coefficient[total_results$`Gene Name`==input$gene_name],"\n"))
          cat(paste("Oscillation Type:",total_results$`Oscillation Type`[total_results$`Gene Name`==input$gene_name],"\n"))
          cat(paste("Amplitude:", total_results$Amplitude[total_results$`Gene Name`==input$gene_name],"\n"))
          cat(paste("Radian.Frequency:",total_results$Radian.Frequency[total_results$`Gene Name`==input$gene_name],"\n"))
          cat(paste("Period:",total_results$Period[total_results$`Gene Name`==input$gene_name],"\n"))
          cat(paste("Phase Shift:",total_results$`Phase Shift`[total_results$`Gene Name`==input$gene_name],"\n"))
          cat(paste("Hours Shifted:",total_results$`Hours Shifted`[total_results$`Gene Name`==input$gene_name],"\n"))
          cat(paste("P-Value:",total_results$`P-Value`[total_results$`Gene Name`==input$gene_name],"\n"))
          cat(paste("BH Adj P-Value:",total_results$`BH Adj P-Value`[total_results$`Gene Name`==input$gene_name],"\n"))
          cat(paste("BY Adj P-Value:",total_results$`BY Adj P-Value`[total_results$`Gene Name`==input$gene_name],"\n"))
        })
        if (num_reps == 1){ # single replicate
          # generate a graph of gene expression with original data and fitted data
          
          # getting the total results: original and fitted values
          data.m <- data.frame(matrix(0,length(timen),3))
          colnames(data.m) <- c("Original","Fit","Times")
          data.m$Original <- as.numeric(total_results[total_results$`Gene Name`==input$gene_name,c(16:(length(timen)+15))])
          data.m$Fit <- as.numeric(total_results[total_results$`Gene Name`==input$gene_name,-c(1:(length(timen)+15))])
          data.m$Times <- timen
          
          # create gene expression plot
          col_vect <- c("Original"="black","Fit"="green")
          plot_viz <- ggplot(data = data.m,aes(x=Times))+
            geom_line(aes(y=Original,colour="Original"))+
            geom_line(aes(y=Fit,colour="Fit"))+
            scale_color_manual("",values=col_vect)+
            ggtitle(paste(input$gene_name))+
            theme(text= element_text(size = 20),plot.title = element_text(hjust = .5),
                  legend.position = "bottom",legend.direction = "horizontal")+
            labs(x="Hours",y="Expression")
          
        }
        else{ # multiple replicates
          col_vect <- c("Original"="gray","Fit"="black")
          #Plot the fit line and min/max shading
          plot_viz<-ggplot(data = ribbon.df,aes(x=Times))+ # declare the dataframe and main variables
            geom_ribbon(aes(x=Times, ymax=Max, ymin=Min, colour = "Original"), fill = "gray", alpha = .5)+ # create shading
            geom_line(aes(y=Fit,colour="Fit"))+ # plot fit line
            scale_color_manual("",values=col_vect)+
            scale_fill_manual("",values=col_vect)+
            ggtitle(paste(input$gene_name))+
            theme(text= element_text(size = 20),plot.title = element_text(hjust = .5),
                  legend.position = "bottom",legend.direction = "horizontal")+
            labs(x="Hours", y="Expression") #Labels for axes
        }
      }
    }
    else if(input$viz=="Gene Expression with Replicates"){
      #Plot the fit line and min/max shading
      #Plot the fit line with shading and replicate expressions
      
      if(num_reps == 1 ||sum(total_results$`Gene Name`==input$gene_name)==0){ # if the gene name inputted isn't in set
        output$text <- renderPrint({cat(paste(""))}) # no text
        plot_viz<-ggplot() # no plot
      }
      else{
        output$text <- renderPrint({  # summary in the UI is all of the statistics vital to the gene (name, p-value, etc.)
          cat(paste("Gene Name:",total_results$`Gene Name`[total_results$`Gene Name`==input$gene_name],"\n"))
          cat(paste("Convergence:", total_results$Convergence[total_results$`Gene Name`==input$gene_name],"\n"))
          cat(paste("Iterations:",total_results$Iterations[total_results$`Gene Name`==input$gene_name],"\n"))
          cat(paste("Forcing Coefficient:", total_results$Forcing.Coefficient[total_results$`Gene Name`==input$gene_name],"\n"))
          cat(paste("Oscillation Type:",total_results$`Oscillation Type`[total_results$`Gene Name`==input$gene_name],"\n"))
          cat(paste("Amplitude:", total_results$Amplitude[total_results$`Gene Name`==input$gene_name],"\n"))
          cat(paste("Radian.Frequency:",total_results$Radian.Frequency[total_results$`Gene Name`==input$gene_name],"\n"))
          cat(paste("Period:",total_results$Period[total_results$`Gene Name`==input$gene_name],"\n"))
          cat(paste("Phase Shift:",total_results$`Phase Shift`[total_results$`Gene Name`==input$gene_name],"\n"))
          cat(paste("Hours Shifted:",total_results$`Hours Shifted`[total_results$`Gene Name`==input$gene_name],"\n"))
          cat(paste("P-Value:",total_results$`P-Value`[total_results$`Gene Name`==input$gene_name],"\n"))
          cat(paste("BH Adj P-Value:",total_results$`BH Adj P-Value`[total_results$`Gene Name`==input$gene_name],"\n"))
          cat(paste("BY Adj P-Value:",total_results$`BY Adj P-Value`[total_results$`Gene Name`==input$gene_name],"\n"))
        })
        
        if(num_reps > 8){ # cannot visualize more than 8 replicates
          plot_viz<-ggplot()
        }
        else{
          color_bar <- c("Rep. 1"="red","Rep. 2"="blue","Rep. 3"="green",
                         "Rep. 4"="yellow","Rep. 5"="purple","Rep. 6"="pink",
                         "Rep. 7"="orange","Rep. 8"="magenta","Fit"="black",
                         "Original"="grey")
          
          plot_viz<-ggplot(data = ribbon.df,aes(x=Times))+ # declare the dataframe and main variables
            geom_ribbon(aes(x=Times, ymax=Max, ymin=Min, colour="Original"),
                        fill = "gray", alpha = 0.5)+ # create shading
            geom_line(aes(y=Fit,colour="Fit"))+ # fitted values
            ggtitle(paste(input$gene_name))+ # gene name is title
            scale_color_manual("",values=color_bar)+
            scale_fill_manual("",values=color_bar)+
            theme(text= element_text(size = 20),plot.title = element_text(hjust = .5),
                  legend.position = "bottom",legend.direction = "horizontal")+
            labs(x="Hours", y="Expression") #Label for axes
          
          # add specific replicate lines 
          for (i in 1:num_reps){
            plot_viz <- plot_viz + geom_line(data = ribbon.df,aes_string(x="Times",y=paste("Rep",i,sep = ".")), colour=color_bar[i], alpha=0.6)
          }
          
        }
        
      }
    }
    else if(input$viz == "Parameter Density Graph (PDG)"){
      { # each of these if statements create density graphs based on the input coefficients
        # and the subset specified
        # the summary is a 5-number summary of the coefficient
        
        # you can't have JTK and look at JTK_specific subsets
        if (!is_jtk && (input$subset_look == "JTK" || input$subset_look == "Both ECHO and JTK" || input$subset_look == "Only JTK" || input$subset_look == "Only ECHO" || input$subset_look == "Neither ECHO and JTK")){
          output$text <- renderPrint({
            "N/A"
          })
          plot_viz<-ggplot()
        }
        else if(input$subset_look == "None"){
          output$text <- renderPrint({
            summary(total_results[,input$coeff])
          })
          
          plot_viz<-ggplot(total_results[!is.na(total_results$Forcing.Coefficient),], aes_string(ggname(input$coeff))) +
            geom_density()+
            ggtitle(paste(input$subset_look,": ",input$coeff, sep = ""))+
            theme(text= element_text(size = 20),plot.title = element_text(hjust = .5))
          
        }
        else if(input$subset_look == "JTK_CYCLE"){
          output$text <- renderPrint({
            summary(total_results[circ_jtk,input$coeff])
          })
          
          plot_viz<-ggplot(total_results[circ_jtk,], aes_string(ggname(input$coeff))) +
            geom_density()+
            ggtitle(paste(input$subset_look,": ",input$coeff, sep = ""))+
            theme(text= element_text(size = 20),plot.title = element_text(hjust = .5))
          
        }
        else if(input$subset_look == "ECHO"){
          output$text <- renderPrint({
            summary(total_results[circ_us,input$coeff])
          })
          
          plot_viz<-ggplot(total_results[circ_us,], aes_string(ggname(input$coeff))) +
            geom_density()+
            ggtitle(paste(input$subset_look,": ",input$coeff, sep = ""))+
            theme(text= element_text(size = 20),plot.title = element_text(hjust = .5))
          
        }
        else if(input$subset_look == "Damped"){
          output$text <- renderPrint({
            summary(total_results[damped,input$coeff])
          })
          
          plot_viz<-ggplot(total_results[damped,], aes_string(ggname(input$coeff))) +
            geom_density()+
            ggtitle(paste(input$subset_look,": ",input$coeff, sep = ""))+
            theme(text= element_text(size = 20),plot.title = element_text(hjust = .5))
          
        }
        else if(input$subset_look == "Driven"){
          output$text <- renderPrint({
            summary(total_results[driven,input$coeff])
          })
          
          plot_viz<-ggplot(total_results[driven,], aes_string(ggname(input$coeff))) +
            geom_density()+
            ggtitle(paste(input$subset_look,": ",input$coeff, sep = ""))+
            theme(text= element_text(size = 20),plot.title = element_text(hjust = .5))
          
        }
        else if(input$subset_look == "Harmonic"){
          output$text <- renderPrint({
            summary(total_results[harmonic,input$coeff])
          })
          
          plot_viz<-ggplot(total_results[harmonic,], aes_string(ggname(input$coeff))) +
            geom_density()+
            ggtitle(paste(input$subset_look,": ",input$coeff, sep = ""))+
            theme(text= element_text(size = 20),plot.title = element_text(hjust = .5))
          
        }
        else if(input$subset_look == "Overexpressed"){
          output$text <- renderPrint({
            summary(total_results[overexpressed,input$coeff])
          })
          
          plot_viz<-ggplot(total_results[overexpressed,], aes_string(ggname(input$coeff))) +
            geom_density()+
            ggtitle(paste(input$subset_look,": ",input$coeff, sep = ""))+
            theme(text= element_text(size = 20),plot.title = element_text(hjust = .5))
          
        }
        else if(input$subset_look == "Repressed"){
          output$text <- renderPrint({
            summary(total_results[repressed,input$coeff])
          })
          
          plot_viz<-ggplot(total_results[repressed,], aes_string(ggname(input$coeff))) +
            geom_density()+
            ggtitle(paste(input$subset_look,": ",input$coeff, sep = ""))+
            theme(text= element_text(size = 20),plot.title = element_text(hjust = .5))
          
        }
        else if(input$subset_look == "Nonconverged"){
          output$text <- renderPrint({
            summary(total_results[nonconv,input$coeff])
          })
          
          plot_viz<-ggplot(total_results[nonconv & !is.na(total_results$Forcing.Coefficient),], aes_string(ggname(input$coeff))) +
            geom_density()+
            ggtitle(paste(input$subset_look,": ",input$coeff, sep = ""))+
            theme(text= element_text(size = 20),plot.title = element_text(hjust = .5))
        }
        else if(input$subset_look == "Nonstarter"){
          output$text <- renderPrint({
            summary(total_results[nas_found-nodev,input$coeff])
          })
          
          plot_viz<- ggplot(total_results[(nas_found-nodev) & !is.na(total_results$Forcing.Coefficient),], aes_string(ggname(input$coeff))) +
            geom_density()+
            ggtitle(paste(input$subset_look,": ",input$coeff, sep = ""))+
            theme(text= element_text(size = 20),plot.title = element_text(hjust = .5))
        }
        else if(input$subset_look == "No Deviation"){
          output$text <- renderPrint({
            summary(total_results[nodev,input$coeff])
          })
          
          ggplot(total_results[nodev & !is.na(total_results$Forcing.Coefficient),], aes_string(ggname(input$coeff))) +
            geom_density()+
            ggtitle(paste(input$subset_look,": ",input$coeff, sep = ""))+
            theme(text= element_text(size = 20),plot.title = element_text(hjust = .5))
          
        }
        else if(input$subset_look == "Both ECHO and JTK"){
          output$text <- renderPrint({
            summary(total_results[both,input$coeff])
          })
          
          plot_viz<-  ggplot(total_results[both,], aes_string(ggname(input$coeff))) +
            geom_density()+
            ggtitle(paste(input$subset_look,": ",input$coeff, sep = ""))+
            theme(text= element_text(size = 20),plot.title = element_text(hjust = .5))
        }
        else if(input$subset_look == "Only JTK"){
          output$text <- renderPrint({
            summary(total_results[theirs,input$coeff])
          })
          
          plot_viz<- ggplot(total_results[theirs & !is.na(total_results$Forcing.Coefficient),], aes_string(ggname(input$coeff))) +
            geom_density()+
            ggtitle(paste(input$subset_look,": ",input$coeff, sep = ""))+
            theme(text= element_text(size = 20),plot.title = element_text(hjust = .5))
        }
        else if(input$subset_look == "Only ECHO"){
          output$text <- renderPrint({
            summary(total_results[ours,input$coeff])
          })
          
          plot_viz<-ggplot(total_results[ours,], aes_string(ggname(input$coeff))) +
            geom_density()+
            ggtitle(paste(input$subset_look,": ",input$coeff, sep = ""))+
            theme(text= element_text(size = 20),plot.title = element_text(hjust = .5))
        }
        else if(input$subset_look == "Neither ECHO and JTK"){
          output$text <- renderPrint({
            summary(total_results[diff,input$coeff])
          })
          
          plot_viz<-ggplot(total_results[diff & !is.na(total_results$Forcing.Coefficient),], aes_string(ggname(input$coeff))) +
            geom_density()+
            ggtitle(paste(input$subset_look,": ",input$coeff, sep = ""))+
            theme(text= element_text(size = 20),plot.title = element_text(hjust = .5))
        }
      }
    }
    else if(input$viz == "Venn Diagram"){
      # generates a venn diagram of JTK output and extended harmonic oscillator output
      # summary is an overall summary of how extended harmonic oscillators did, both internally and compared to JTK
      plot_viz <- ggplot()
      
      if (!is_jtk){
        output$text <- renderPrint({ 
          # generate summary of outputs given by ECHO
          cat(paste("Nonconverged:",sum(nonconv,na.rm=TRUE),"\n"))
          cat(paste("Nonstarter (i.e., too much noise):", sum(nas_found-nodev,na.rm=TRUE),"\n"))
          cat(paste("Unexpressed:", sum(noexpr,na.rm=TRUE),"\n"))
          cat(paste("No Deviation (constant):",sum(nodev,na.rm=TRUE),"\n"))
          cat(paste("Circadian (ECHO):", sum(circ_us),"\n"))
          cat(paste("  Damped:", sum(damped & circ_us,na.rm=TRUE),"\n"))
          cat(paste("  Driven:", sum(driven & circ_us,na.rm=TRUE),"\n"))
          cat(paste("  Harmonic:", sum(harmonic & circ_us,na.rm=TRUE),"\n"))
          cat(paste("  Overexpressed:", sum(overexpressed & circ_us,na.rm=TRUE),"\n"))
          cat(paste("  Repressed:", sum(repressed & circ_us,na.rm=TRUE),"\n"))
          
          # also put user inputs
          cat("\n")
          cat("User Inputs:\n")
          cat(paste("ECHO End Date and Time: ",user_input$ECHO_end_date,"\n"))
          cat(paste("Begin: ",user_input$begin,"\n"))
          cat(paste("End: ",user_input$end,"\n"))
          cat(paste("Resolution: ",user_input$resol,"\n"))
          cat(paste("Number of Replicates: ",num_reps,"\n"))
          cat(paste("Seeking Rhythms, Lower End: ",low_end,"\n"))
          cat(paste("Seeking Rhythms, Higher End: ",high_end,"\n"))
          cat(paste("Paired Replicates: ",user_input$tied,"\n"))
          cat(paste("Smoothing?: ",user_input$is_smooth,"\n"))
          cat(paste("Remove unexpressed genes?: ",user_input$rem_unexpr,"\n"))
          cat(paste("Remove unexpressed genes, percentage: ",user_input$rem_unexpr_amt,"\n"))
          cat(paste("Normalize data?: ",user_input$is_normal,"\n"))
          cat(paste("Remove linear trend?: ",user_input$is_de_linear_trend,"\n"))
          cat(paste("Run JTK?: ",user_input$run_jtk,"\n"))
          cat(paste("ECHO Version No.: ",user_input$v_num,"\n"))
        })
        plot_viz <- renderPlot({plot_viz <- ggplot()}) # no venn diagram
      }
      else{ # generate summary of outputs given by ECHO and JTK
        output$text <- renderPrint({ # generate summary
          cat(paste("Nonconverged:",sum(nonconv,na.rm=TRUE),"\n"))
          cat(paste("Nonstarter (i.e., too much noise):", sum(nas_found-nodev,na.rm=TRUE),"\n"))
          cat(paste("Unexpressed:", sum(noexpr,na.rm=TRUE),"\n"))
          cat(paste("No Deviation (constant):",sum(nodev,na.rm=TRUE),"\n"))
          cat(paste("Circadian (ECHO):", sum(circ_us),"\n"))
          cat(paste("  Damped:", sum(damped & circ_us,na.rm=TRUE),"\n"))
          cat(paste("  Driven:", sum(driven & circ_us,na.rm=TRUE),"\n"))
          cat(paste("  Harmonic:", sum(harmonic & circ_us,na.rm=TRUE),"\n"))
          cat(paste("  Overexpressed:", sum(overexpressed & circ_us,na.rm=TRUE),"\n"))
          cat(paste("  Repressed:", sum(repressed & circ_us,na.rm=TRUE),"\n"))
          cat(paste("Circadian (JTK_CYCLE):", sum(circ_jtk),"\n"))
          cat("Confusion Matrix of Circadian Genes
                                JTK_CYCLE
                       yes                  no
              ECHO yes Both ECHO and JTK    Only ECHO
                   no  Only JTK             Neither ECHO nor JTK\n")
          cat(paste("Both ECHO and JTK:",sum(both,na.rm=TRUE),"\n"))
          cat(paste("Only JTK:",sum(theirs,na.rm=TRUE),"\n"))
          cat(paste("Only ECHO:",sum(ours,na.rm=TRUE),"\n"))
          cat(paste("Neither ECHO nor JTK:",sum(diff,na.rm=TRUE),"\n"))
          
          # also put user inputs
          cat("\n")
          cat("User Inputs:\n")
          cat(paste("ECHO End Date and Time: ",user_input$ECHO_end_date,"\n"))
          cat(paste("Begin: ",user_input$begin,"\n"))
          cat(paste("End: ",user_input$end,"\n"))
          cat(paste("Resolution: ",user_input$resol,"\n"))
          cat(paste("Number of Replicates: ",num_reps,"\n"))
          cat(paste("Seeking Rhythms, Lower End: ",low_end,"\n"))
          cat(paste("Seeking Rhythms, Higher End: ",high_end,"\n"))
          cat(paste("Paired Replicates: ",user_input$tied,"\n"))
          cat(paste("Smoothing?: ",user_input$is_smooth,"\n"))
          cat(paste("Remove unexpressed genes?: ",user_input$rem_unexpr,"\n"))
          cat(paste("Remove unexpressed genes, percentage: ",user_input$rem_unexpr_amt,"\n"))
          cat(paste("Normalize data?: ",user_input$is_normal,"\n"))
          cat(paste("Remove linear trend?: ",user_input$is_de_linear_trend,"\n"))
          cat(paste("Run JTK?: ",user_input$run_jtk,"\n"))
          cat(paste("ECHO Version No.: ",user_input$v_num,"\n"))
        })
        
        output$plot_viz <- renderPlot({
          # makes a venn diagram of overlap between JTK and ECHO
          grid.newpage()
          draw.pairwise.venn(sum(circ_us),sum(circ_jtk),sum(circ_us & circ_jtk), category = c("ECHO","JTK_CYCLE"),lty = rep("blank",2), fill = c("light blue","pink"), alpha = rep(.5,2), cat.pos = c(0, 0), cat.dist=rep(.025,2),cat.fontfamily = rep("sans",2),fontfamily ="sans",cex = rep(3, 3), cat.cex = rep(3, 2))
        })
      }
    }
    else if(input$viz == "Heat Map"){ # to create a heat map
      plot_viz <- ggplot()
      
      start.time <- Sys.time() # begin counting time
      
      if (!is_jtk && (input$heat_subset_look == "JTK" || input$heat_subset_look == "Both ECHO and JTK" || input$heat_subset_look == "Only JTK" || input$heat_subset_look == "Only ECHO" || input$heat_subset_look == "Neither ECHO and JTK" || input$heat_subset_rep != "all" && is.na(as.numeric(input$heat_subset_rep)) || (input$heat_subset_rep != "all" && as.numeric(input$heat_subset_rep) > num_reps) || input$heat_subset_rep != "all" && as.numeric(input$heat_subset_rep) <= 0)){
        output$text <- renderPrint({
          "N/A"
        })
        plot_viz<-ggplot()
        total_results_na <- data.frame()
      }
      else if(input$heat_subset_look == "None"){
        total_results_na <-  total_results[!is.na(total_results$Amplitude),]; # filter out the rows where the paremeters are NA
        
      }
      else if(input$heat_subset_look == "JTK_CYCLE"){
        total_results_na <-  total_results[circ_jtk,]; # filter out the rows where the paremeters are NA
      }
      else if(input$heat_subset_look == "ECHO"){
        # filter out the rows where the paremeters are NA and our characteristic
        total_results_na <-  total_results[circ_us,]; 
        total_results_na <- total_results_na[!is.na(total_results_na$Amplitude),] 
        
      }
      else if(input$heat_subset_look == "Damped"){
        # filter out the rows where the paremeters are NA and our characteristic
        total_results_na <-  total_results[damped,]; 
        total_results_na <- total_results_na[!is.na(total_results_na$Amplitude),] 
        
      }
      else if(input$heat_subset_look == "Driven"){
        # filter out the rows where the paremeters are NA and our characteristic
        total_results_na <-  total_results[driven,]; 
        total_results_na <- total_results_na[!is.na(total_results_na$Amplitude),] 
      }
      else if(input$heat_subset_look == "Harmonic"){
        # filter out the rows where the paremeters are NA and our characteristic
        total_results_na <-  total_results[harmonic,]; 
        total_results_na <- total_results_na[!is.na(total_results_na$Amplitude),] 
      }
      else if(input$heat_subset_look == "Overexpressed"){
        # filter out the rows where the paremeters are NA and our characteristic
        total_results_na <-  total_results[overexpressed,]; 
        total_results_na <- total_results_na[!is.na(total_results_na$Amplitude),] 
      }
      else if(input$heat_subset_look == "Repressed"){
        # filter out the rows where the paremeters are NA and our characteristic
        total_results_na <-  total_results[repressed,]; 
        total_results_na <- total_results_na[!is.na(total_results_na$Amplitude),] 
      }
      else if(input$heat_subset_look == "Both ECHO and JTK"){
        # filter out the rows where the paremeters are NA and our characteristic
        total_results_na <-  total_results[both,]; 
        total_results_na <- total_results_na[!is.na(total_results_na$Amplitude),] 
      }
      else if(input$heat_subset_look == "Only JTK"){
        total_results_na <-  total_results[theirs,]; # filter out the rows where the paremeters are NA
      }
      else if(input$heat_subset_look == "Only ECHO"){
        # filter out the rows where the paremeters are NA and our characteristic
        total_results_na <-  total_results[ours,]; 
        total_results_na <- total_results_na[!is.na(total_results_na$Amplitude),] 
      }
      else if(input$heat_subset_look == "Neither ECHO and JTK"){
        # filter out the rows where the paremeters are NA and our characteristic
        total_results_na <-  total_results[diff,]; 
        total_results_na <- total_results_na[!is.na(total_results_na$Amplitude),] 
      }
      #total_results_na <-  total_results[!is.na(total_results$Amplitude),]; # filter out the rows where the paremeters are NA
      
      # don't let it output if jtk requested without JTK_results
      if (!(!is_jtk && (input$heat_subset_look == "JTK" || input$heat_subset_look == "Both ECHO and JTK" || input$heat_subset_look == "Only JTK" || input$heat_subset_look == "Only ECHO" || input$heat_subset_look == "Neither ECHO and JTK"))){
        
        if (input$heat_subset_look != "JTK_CYCLE" & input$heat_subset_look != "Only JTK"){
          phase <- total_results_na$`Phase Shift`
          amp <- total_results_na$Amplitude
          
          #loop through every row in replace negative amplitudes with positive ones
          # adjust the phase shift accordingly
          for (i in 1:length(phase)) {
            if(amp[i]<0){
              amp[i] <- -1*amp[i]
              #total_results_na[i,16:(15+length(timen)*num_reps)] <- -1*total_results_na[i,16:(15+length(timen)*num_reps)]
              phase[i] <- phase[i]+pi
            }
            if(phase[i]>2*pi){
              phase[i] <- phase[i]-2*pi
            }
            if(phase[i]<0){
              phase[i] <- phase[i]+2*pi
            }
          }
        } else {
          if(input$heat_subset_look == "JTK_CYCLE"){
            phase <- -1*JTK_results$LAG[circ_jtk]
            amp <- JTK_results$AMP[circ_jtk]
          } else {
            phase <- -1*JTK_results$LAG[theirs]
            amp <- JTK_results$AMP[theirs]
          }
          
        }
        
        
        #get matrix of just the relative expression over time
        hm_mat <- as.matrix(total_results_na[,16:(15+length(timen)*num_reps)])
        
        if (input$heat_subset_rep == "all"){ # make an average of replicates
          #if there are replicates, average the relative expression for each replicate
          mtx_reps <- list() # to store actual matrix
          mtx_count <- list() # to store how many are NA
          for (i in 1:num_reps){
            mtx_reps[[i]] <- hm_mat[, seq(i,ncol(hm_mat), by=num_reps)]
            mtx_count[[i]] <- is.na(mtx_reps[[i]])
            mtx_reps[[i]][is.na(mtx_reps[[i]])] <- 0
          }
          repmtx <- matrix(0L,ncol = length(timen),nrow = nrow(hm_mat))+num_reps # to store how many we should divide by
          hm_mat <- matrix(0L,ncol = length(timen),nrow = nrow(hm_mat)) # to store the final result
          for (i in 1:num_reps){
            hm_mat <- hm_mat + mtx_reps[[i]] # sum the replicates
            repmtx <- repmtx - mtx_count[[i]] # how many replicates are available for each time point
          }
          repmtx[repmtx==0] <- NA # to avoid division by 0 and induce NAs if there are no time points available
          hm_mat <- hm_mat/repmtx
          
        } else { # display only one time point
          # subset the heat map
          rep_looking_at <- as.numeric(input$heat_subset_rep)
          hm_mat <- hm_mat[,seq(rep_looking_at,ncol(hm_mat),by=num_reps)]
        }
        
        # center rows around mean
        # vector of row means
        all_row_mean <- rowMeans(hm_mat, na.rm = TRUE)
        hm_mat <- hm_mat - all_row_mean
        
        #normalize each row to be between -1 and 1
        for (i in 1:length(phase)){
          # gene_mean <-mean(as.matrix(hm_mat[i,]),na.rm=TRUE)
          # for (j in 1:length(timen)){
          #   hm_mat[i,j] <- (hm_mat[i,j]-gene_mean)
          # }
          
          gene_max <- max(abs((hm_mat[i,])),na.rm = TRUE)
          hm_mat[i,] <- hm_mat[i,]/gene_max 
          # for (j in 1:length(timen)){
          #   hm_mat[i,j] <- (hm_mat[i,j]/gene_max)
          # }
        }
        
        #sort by phase shift
        hm_mat <- hm_mat[order(phase),]
        
        output$plot_viz <- renderPlot({
          par(mar = c(1,2,1,2))
          image.plot(t(hm_mat),col = blue2yellow(256),xlab = "Hours",ylab = "Expressions",axes=FALSE,lab.breaks=NULL)
          
        })
        
        # time taken to produce heat map
        end.time <- Sys.time()
        time.taken.heat <- end.time - start.time
        output$text <- renderPrint({cat(paste("Render Time:",time.taken.heat,"seconds\n"))
          cat("Expressions are each row and time points are each column.\n")
          cat(paste("Number of expressions:",nrow(hm_mat)))})
      }
    }
    
    
    #output plot
    if (input$viz != "Venn Diagram" && input$viz != "Heat Map"){
      output$plot_viz<- suppressWarnings(renderPlot({
        suppressWarnings(plot_viz)
      }))
    }
    
    
    #function to download png of plot
    output$downloadPlot <- downloadHandler(
      filename = function() { paste(input$viz, '.png', sep='') },
      content = function(file) 
      {
        png(file)
        print(plot_viz)
        dev.off()
      },
      contentType='image/png'
    )
  })
}

# run ShinyApp ----

shinyApp(ui, server)