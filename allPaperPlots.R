#!/usr/bin/env Rscript
suppressMessages({
	library(tidyverse)
	library(ggplot2)
	library(cowplot)
	library(ggbeeswarm)
	library(ggsignif)
	library(rstatix)
	library(data.table)
	library(see)
	library(dunn.test)
	library(ROCit)
	library(survival)
	library(survminer) ## Use install_github("rschauner/survminer") version
	library(parallel)
	library(glmnet)
	library(emmeans)})

###Configuration
################
###Plotting config
##Dimensions from Cancer Research
columnWidthInch <- 3.125
pageWidthInch <- columnWidthInch*2+3/16
minTextSize <- 6
maxTextSize <- 12
minLineWidth <- 0.5

#ThisTheme variables
fontSize <- minTextSize
tinyFontSize <- minTextSize
smallFontSize <- minTextSize
largeFontSize <- round(fontSize * 1.25)
gridLabelFontSize <- largeFontSize

thisTheme <- theme_cowplot(
  #font_size = round(minTextSize*14/11) #makes sure tiny is >minTextSize
  font_size = fontSize,
  rel_tiny = tinyFontSize/fontSize,
  rel_small = smallFontSize/fontSize,
  rel_large = largeFontSize/fontSize
)

thisThemeFun <- function(x){return(thisTheme)}
theme_set(thisTheme)

#Derived font sizing
ommnibusTextSize <- fontSize/.pt
ggsignifTextSize <- fontSize/.pt

#Plot element sizing
pointRangeSize <- 0.5/.pt
pointRangeLineWidth <- 1.5/.pt
forestPointSize <- 6/.pt
geomLineWidth <- 0.5/.pt
quasirandomSize <- 1/.pt

##Plot sizing
aspectRatio.plot <- 1

#ggsignif
ggsignifTextVjust <- -0.3
ggsignifTipLength <- 0.025
ggsignifStepIncrease <- 0.17

#Other plotting parameters
dodgeWidth <- 0.75

##Colors
pointRangeColor='#980043'

#Cross-sectional
crossColors=c("#003C6799","#8F770099","#A7303099")
names(crossColors)=c("Pure", "Synchronous", "Invasive")
crossNames=c("Pure DCIS","Synchronous DCIS","Synchronous IDC")
names(crossNames)=c("Pure","Synchronous","Invasive")
crossNamesSplit=c("Pure\nDCIS","Synchronous\nDCIS","Synchronous\nIDC")
names(crossNamesSplit)=c("Pure","Synchronous","Invasive")

#Longitudinal
longColors=c("#7AA6DC99","#EFC00099","#CD534C99")
names(longColors)=c(0,1,2)
longNames=c("Nonrecurrents","Recurrents","Progressors")
names(longNames)=c(0,1,2)
longShortNames=c("Nonrec.","Rec.","Prog.")
names(longShortNames)=c(0,1,2)

###Execution parameters
minLassoFreq <- 0.9
nLasso <- 100
theSeed <- 20
alpha <- 0.05
alphaPlot <- 0.1
RNGkind("L'Ecuyer-CMRG") #Non-standard random number generator that makes mclapply is reproducible

###Data Modification
estimatedMarginC1 <- 0.5
estimatedMarginC2 <- 2.5
estimatedMarginC3 <- 10

###FUNCTIONS
############

#' Generates a latex table from a data.table
#' 
#' @param x data.table
#' @param outfile outputfile name
#' @param precission number of significant digits to print if -4 < exponent < precission, otherwise exponential notation (precission for sprintf g)
#' @param alpha print rows with p.val <= alpha
#' @param pvalColumn name of p.val column 
#' @param caption title
#' @param footer table footnote
#' @param supplementary Add S to the table number
textTable=function(x,outfile="outfile.tex",precission=2, alpha=0.05, pvalColumn="p.val",caption=NULL,footer=NULL,tableNumber=1, supplementary=T) {
  require(xtable)
  colClasses=sapply(x,class)
  display=c("s",ifelse(colClasses=="factor" | colClasses=="character","s","g"))
  
  bold <- function(x) {paste('{\\textbf{',sanitize(x),'}}', sep ='')}
  
  sink(outfile)
  #Environment
  cat("\\documentclass[convert={density=300,outext=.png},varwidth=\\maxdimen]{standalone}
\t\\usepackage{booktabs}
\t\\usepackage[flushleft]{threeparttable} %flushleft here plus [para] in the tablenotes for notes without items
\t\\usepackage[labelfont=bf]{caption} %bold Table X.")
  if(supplementary) {cat("\n\t\\renewcommand{\\thetable}{S\\arabic{table}}")}
  #Staring document
  cat("\n\t\\begin{document}")
  #Setting table counter
  cat(sprintf("\n\t\\setcounter{table}{%d}",tableNumber-1))
  #Floating environment with threeparttable to have table label and footer
  cat("\n\t\\begin{table}
\t\t\\begin{threeparttable}[b]")
  #Add title if needed
  if(!is.null(caption)){
    cat(sprintf("\n\t\t\t\\caption{\\textbf{%s}}",caption))
  }
  #Table body
  print(xtable(x[get(pvalColumn)<=alpha,][order(get(pvalColumn))],
               digits=precission,
               display=display),
        include.rownames=F,
        include.colnames=T,
        sanitize.colnames.function=bold,
        booktabs=T,
        floating=F, #I need to add specific options to threeparttable, so I need to add table manually
  )
  #Table footer
  if(!is.null(footer)){
    cat("\n\t\t\t\\begin{tablenotes}[para]
\t\t\t\t\\footnotesize")
    cat(sprintf("\n\t\t\t\t\\item %s",footer))
    cat("\n\t\t\t\\end{tablenotes}")
  }
  cat("\n\t\t\\end{threeparttable}
\t\\end{table}
\\end{document}")
  sink()
}

#' Returns the optimal ROC threshold given true and predicted values
#' 
#' @param py predicted values
#' @param y true binary classification values
#' @param criterium criterium to pick the threshold. One of minTopLeft (smallest distance to the top-left corner of the ROC curve), J (maximum Youden's J), and absJ (abs(J))
#' @return threshold to split py in 0 and ones
getROCThreshold <- function(py, y, criterium=c("minTopLeft", "J", "absJ")) {
  require(ROCit)
  criterium <- match.arg(criterium)
  roc <- rocit(py,y)
  
  th <- NULL
  if(criterium=="minTopLeft"){
    minTopLeft <- (1-roc$TPR)^2 + roc$FPR^2
    th <- roc$Cutoff[which.min(minTopLeft)]
  }else{
    J <- roc$TPR - roc$FPR
    if(criterium=="J"){
      th <- roc$Cutoff[which.max(J)]
    }else{
      th <- roc$Cutoff[which.max(abs(J))]
    }
  }
  return(th)
}

#' Returns the optimal classification threshold for a univariate binomial regression model
#' 
#' @param model glm or cox univariate regression
#' @param criterium criterium to pick the threshold using getROCThreshold
#' @return threshold for the X variable 
getThreshold <- function(model,criterium="J"){
  py <- predict(model, type="link")
  linkThreshold <- getROCThreshold(py,model$y,criterium)
  xThreshold <- (linkThreshold - as.numeric(coef(model)[1]))/as.numeric(coef(model)[2])
  
  #Equivalent but in probability space (type="response"), for which then we need to calculate the logit
  #py <- predict(model, type="response")
  #pThreshold <- getROCThreshold(py,model$y,criterium)
  #xThreshold=(log(pThreshold/(1-pThreshold)) - as.numeric(coef(model)[1]))/as.numeric(coef(model)[2])
  
  return(xThreshold)
}

#' Returns a factor indicating, for each sample, if it is in the high or low risk in cox regression
#' It uses the Relative Risk with sample reference and threshold given by either the mean (RR:1) or selected under a ROC criterium
#' 
#' @param x coxph model
#' @param criterium RR or one of getROCThreshold criteria
#' @return factor indicating "low" or "high" according to their risk
getStrata=function(model, y = NULL, criterium = "RR"){
  if(!is(model,"coxph")) stop("This function only works with coxph models")
  py <- predict(model,type = "risk",reference="sample")
  
  if (criterium == "RR") {
    th <- 1
  } else {
    if(is.null(y)) stop("This function needs the true values to find the best threshold under the ROC")
    th <- getROCThreshold(py, y, criterium)
  }
  
  return(as.factor(ifelse(py < th,"low","high")))
}

#' Returns the summary of an univariate cox regression
#' 
#' @param x variable name
#' @param theData data.table with the data
#' @param scaleVar boolean that indicates if the variable x will be scaled or not
#' @param timeToEventVar name of the variable (column in theData) that indicates the time to event
#' @param eventVar name of the variable (column in theData) that indicates if the event was reached or not
#' @return data.table with the columns Variable, p.val, C, HR, HRlowCI, HRupperCI, scaled, outcome
getCox=function(x,theData,scaleVar=T,timeToEventVar,eventVar){
  varname=x
  if(scaleVar==T){
    var=paste0("scale(",varname,")")
  } else {
    var=varname
  }
  model=summary(
    coxph(
      as.formula(paste0("Surv(",timeToEventVar,",",eventVar,") ~ ",var)),
      data=theData))
  concordance=model$concordance["C"]
  p.val=model$waldtest["pvalue"]
  hr=model$coefficients[which(colnames(model$coefficients)=="exp(coef)")]
  lowerCI=model$conf.int[which(colnames(model$conf.int)=="lower .95")]
  upperCI=model$conf.int[which(colnames(model$conf.int)=="upper .95")]
  return(data.table(Variable=varname, p.val=p.val,C=concordance, HR=hr, HRlowCI=lowerCI, HRupperCI=upperCI, scaled=scaleVar, outcome=eventVar))
}

#' Returns the summary of an univariate binomial regression
#' 
#' @param x variable name
#' @param theData data.table with the data
#' @param scaleVar boolean that indicates if the variable x will be scaled or not
#' @param eventVar name of the variable (column in theData) that indicates if the event was reached or not
#' @return data.table with the columns Variable, p.val, coef, aic, scaled, outcome
getLogisticModel=function(x,theData,scaleVar=T,eventVar){
  varname <- x
  if(scaleVar==T){
    var <- paste0("scale(",varname,")")
  } else {
    var <- varname
  }
  model <- summary(
    glm(
      as.formula(paste0(eventVar," ~ ",var)),
      data=theData,
      family = "binomial"))
  p.val <- model$coefficients[var,"Pr(>|z|)"]
  coefficient <- model$coefficients[var,"Estimate"]
  aic <- model$aic
  return(data.table(Variable=varname, p.val=p.val, coef=coefficient, aic=aic, outcome=eventVar))
}

#' Returns a Kaplan-Meier plot using ggsurvplot
#' 
#' @param kmModel 
#' @param threshold threshold value used to generate the two strata
#' @param labelConstant a string constant to add to the threshold for the label, e.g., threshold's units
#' @param ylab Y-axis label
#' @param xlab X-axis label
#' @param colors 
#' @param title Plot title
#' @return ggsurvplot
getKMPlot=function(kmModel,threshold,labelConstant,ylab = "Survival probability",xlab = "Time (month)", colors=NULL,title=NULL){
  require(cowplot)
  #Making strata labels
  kmLabels=paste(sep=" ",c("≥","<"),threshold,labelConstant)
  kmColors=NULL
  
  #Finding the strata order
  if(length(names(kmModel$strata))!=2) stop("The KM model must have only 2 strata")
  highPos=grep("high",names(kmModel$strata),value=F)
  if(length(highPos)!=1) stop ('Error selecting the high label strata, one level should include "high" in its name')
  
  #Trying to make strata labels and colors safe, ggsurvplot does not use names for labels like ggplot, only positional order!
  if(is.null(colors) | length(colors)!=2) stop("Specify only two colors, one per strata, in the high low order")
  
  if(highPos==2){
    kmLabels=kmLabels[c(2,1)]
    kmColors=colors[c(2,1)]
  } else {
    kmColors=colors
  }
  
  kmPlot=ggsurvplot(kmModel,
                    surv.median.line = "hv",
                    conf.int = TRUE,
                    pval=T,
                    #Labels
                    legend = 'none',
                    xlab = xlab,
                    ylab = ylab, 
                    #Hardcoded dangerous
                    legend.labs = kmLabels,
                    palette = kmColors,
                    ###Risk table config
                    risk.table=T,
                    tables.theme = theme_cleantable(),
                    tables.height = 0.2,
                    pval.size = smallFontSize/.pt,
                    risk.table.fontsize = fontSize/.pt,
                    ggtheme=thisThemeFun())
  
  if(!is.null(title)) kmPlot$plot = kmPlot$plot + labs(title=title)
  
  kmPlot$table <- kmPlot$table + labs(
    title    = NULL,
  )
  return(kmPlot)
}

#' Modification of ggforest for better format control 
myggforest <- function(model, data = NULL,
                       main = "Hazard ratio", 
                       cpositions=c(0.02, 0.22, 0.4),
                       fontsize = 11,
                       relFontSizeAnno = 0.75,
                       pointSize = 2,
                       refLabel = "reference", 
                       sortTerms = T,
                       termDict = NULL,
                       noDigits=2) {
  require(broom)
  require(grid)
  require(gridExtra)
  require(grDevices)
  require(stats)
  require(data.table)
  
  conf.high <- conf.low <- estimate <- NULL
  stopifnot(inherits(model, "coxph"))
  
  # get data and variables/terms from cox model
  data  <- survminer:::.get_data(model, data = data)
  terms <- attr(model$terms, "dataClasses")[-1]
  # removed as requested in #388
  #  terms <- terms[intersect(names(terms),
  #    gsub(rownames(anova(model))[-1], pattern = "`", replacement = ""))]
  
  # use broom to get some required statistics
  coef <- as.data.frame(tidy(model, conf.int = TRUE))
  gmodel <- glance(model)
  
  # extract statistics for every variable
  allTerms <- lapply(seq_along(terms), function(i){
    var <- names(terms)[i]
    if (terms[i] %in% c("factor", "character")) {
      adf <- as.data.frame(table(data[, var]))
      cbind(var = var, adf, pos = 1:nrow(adf))
    }
    else if (terms[i] == "numeric") {
      data.frame(var = var, Var1 = "", Freq = nrow(data),
                 pos = 1)
    }
    else {
      vars = grep(paste0("^", var, "*."), coef$term, value=TRUE)
      data.frame(var = vars, Var1 = "", Freq = nrow(data),
                 pos = seq_along(vars))
    }
  })
  allTermsDF <- do.call(rbind, allTerms)
  colnames(allTermsDF) <- c("var", "level", "N", "pos")
  inds <- apply(allTermsDF[,1:2], 1, paste0, collapse="")
  
  # use broom again to get remaining required statistics
  rownames(coef) <- gsub(coef$term, pattern = "`", replacement = "")
  toShow <- cbind(allTermsDF, coef[inds,])[,c("var", "level", "N", "p.value", "estimate", "conf.low", "conf.high", "pos")]
  toShowExp <- toShow[,5:7]
  toShowExp[is.na(toShowExp)] <- 0
  toShowExp <- format(exp(toShowExp), digits=noDigits)
  toShowExpClean <- data.frame(toShow,
                               pvalue = signif(toShow[,4],noDigits+1),
                               toShowExp)
  toShowExpClean$stars <- paste0(round(toShowExpClean$p.value, noDigits+1), " ",
                                 ifelse(toShowExpClean$p.value < 0.05, "*",""),
                                 ifelse(toShowExpClean$p.value < 0.01, "*",""),
                                 ifelse(toShowExpClean$p.value < 0.001, "*",""))
  toShowExpClean$ci <- paste0("(",toShowExpClean[,"conf.low.1"]," - ",toShowExpClean[,"conf.high.1"],")")
  toShowExpClean$estimate.1[is.na(toShowExpClean$estimate)] = refLabel
  toShowExpClean$stars[which(toShowExpClean$p.value < 0.001)] = "<0.001 ***"
  toShowExpClean$stars[is.na(toShowExpClean$estimate)] = ""
  toShowExpClean$ci[is.na(toShowExpClean$estimate)] = ""
  toShowExpClean$estimate[is.na(toShowExpClean$estimate)] = 0
  toShowExpClean$var = as.character(toShowExpClean$var)
  toShowExpClean$var[duplicated(toShowExpClean$var)] = ""
  # make label strings:
  # toShowExpClean$N <- paste0("(N=",toShowExpClean$N,")")
  
  #flip order
  toShowExpClean <- toShowExpClean[nrow(toShowExpClean):1, ]
  
  rangeb <- range(toShowExpClean$conf.low, toShowExpClean$conf.high, na.rm = TRUE)
  breaks <- axisTicks(rangeb/2, log = TRUE, nint = 7)
  rangeplot <- rangeb
  # make plot twice as wide as needed to create space for annotations
  rangeplot[1] <- rangeplot[1] - diff(rangeb)
  # increase white space on right for p-vals:
  rangeplot[2] <- rangeplot[2] + .15 * diff(rangeb)
  
  width <- diff(rangeplot)
  # y-coordinates for labels:
  y_variable <- rangeplot[1] +  cpositions[1] * width
  y_nlevel <- rangeplot[1]  +  cpositions[2] * width
  y_cistring <- rangeplot[1]  +  cpositions[3] * width
  y_stars <- rangeb[2]
  x_annotate <- seq_len(nrow(toShowExpClean))
  
  # geom_text relFontSizeAnno is in mm (https://github.com/tidyverse/ggplot2/issues/1828)
  annot_size_mm <- relFontSizeAnno *
    as.numeric(convertX(unit(fontsize, "pt"), "mm"))
  
  # Rename terms
  if(!is.null(termDict)){
    dtTerm <- data.table(newVar=termDict,var=names(termDict))
    setDT(toShowExpClean)
    toShowExpClean <- merge(dtTerm,toShowExpClean,by = "var",all.y = T)
    toShowExpClean[,`:=`(oldVar = var)]
    toShowExpClean[,`:=`(var=newVar)]
    toShowExpClean[is.na(newVar),`:=`(var=oldVar)]
  }
  
  
  # Sort using the effect size
  if(sortTerms){
    toShowExpClean <- toShowExpClean[order(abs(toShowExpClean$estimate)),]
  }
  
  p <- ggplot(toShowExpClean, aes(seq_along(var), exp(estimate))) +
    geom_rect(aes(xmin = seq_along(var) - .5, xmax = seq_along(var) + .5,
                  ymin = exp(rangeplot[1]), ymax = exp(rangeplot[2]),
                  fill = ordered(seq_along(var) %% 2 + 1))) +
    scale_fill_manual(values = c("#FFFFFF33", "#00000033"), guide = "none") +
    geom_point(pch = 15, size = pointSize) +
    geom_errorbar(aes(ymin = exp(conf.low), ymax = exp(conf.high)), width = 0.15) +
    geom_hline(yintercept = 1, linetype = 3) +
    coord_flip(ylim = exp(rangeplot)) +
    ggtitle(main) +
    scale_y_log10(
      name = "",
      labels = sprintf("%g", breaks),
      expand = c(0.02, 0.02),
      breaks = breaks) +
    theme_light(base_size = fontSize) +
    theme(panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          legend.position = "none",
          panel.border=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    theme(plot.title = theme_get()$plot.title) + #WARNING HARDCODED
    xlab("") +
    annotate(geom = "text", x = x_annotate, y = exp(y_variable),
             label = toShowExpClean$var, fontface = "bold", hjust = 0,
             size = annot_size_mm) +
    annotate(geom = "text", x = x_annotate, y = exp(y_nlevel), hjust = 0,
             label = toShowExpClean$level, vjust = -0.1, size = annot_size_mm) +
    # annotate(geom = "text", x = x_annotate, y = exp(y_nlevel),
    #          label = toShowExpClean$N, fontface = "italic", hjust = 0,
    #          vjust = ifelse(toShowExpClean$level == "", .5, 1.1),
    #          size = annot_size_mm) +
    annotate(geom = "text", x = x_annotate, y = exp(y_cistring),
             label = toShowExpClean$estimate.1, size = annot_size_mm,
             vjust = ifelse(toShowExpClean$estimate.1 == "reference", .5, -0.1)) +
    annotate(geom = "text", x = x_annotate, y = exp(y_cistring),
             label = toShowExpClean$ci, size = annot_size_mm,
             vjust = 1.1,  fontface = "italic") +
    annotate(geom = "text", x = x_annotate, y = exp(y_stars),
             label = toShowExpClean$stars, size = annot_size_mm,
             hjust = -0.2,  fontface = "italic") +
    annotate(geom = "text", x = 0.5, y = exp(y_variable),
             label = paste0("n: ",unique(toShowExpClean$N),"; # Events: ", gmodel$nevent, "\nLog-Rank p: ",
                            format.pval(gmodel$p.value.log, eps = .001),
                            "; CI: ", round(gmodel$concordance,2)),
             size = annot_size_mm, hjust = 0, vjust = 1.2,  fontface = "italic")
  # annotate(geom = "text", x = 0.5, y = exp(y_variable),
  #          label = paste0("# Events: ", gmodel$nevent, "; Global p-value (Log-Rank): ",
  #                         format.pval(gmodel$p.value.log, eps = ".001"), " \nAIC: ", round(gmodel$AIC,2),
  #                         "; Concordance Index: ", round(gmodel$concordance,2)),
  #          size = annot_size_mm, hjust = 0, vjust = 1.2,  fontface = "italic")
  # switch off clipping for p-vals, bottom annotation:
  gt <- ggplot_gtable(ggplot_build(p))
  gt$layout$clip[gt$layout$name == "panel"] <- "off"
  # grid.draw(gt)
  # invisible(p)
  ggpubr::as_ggplot(gt)
}

#' Returns a forest plot with the regressors of a Cox model, and a Kaplan-Meier plot (using ggsurvplot) stratifying the patients in two groups 
#' 
#' @param theseDataPatients data.frame that contains the column patient with the list of patients to subset theseData
#' @param coxModel cox model that contains the variables to include from theseData
#' @param timeToEventVar name fo the column in theseData that contains the time to event
#' @param eventVar name of the column in theseData that contains the binary variable with the outcome
#' @param theseColors of the two strata/patient groups
#' @param kmYlab y-lab of the km plot
#' @param kmTitle title of the km plot
#' @param forestTitle Title of the forest plot 
#' @param termDict vector of human-readable terms named with the original terms for re-labeling
#' @param rel_height Relative height of the survival table
#' @return list of a forest plot and a km plot
getForestAndKMCox <- function(theseDataPatients, coxModel, timeToEventVar, eventVar, theseColors, 
                              kmTitle = "KM plot", 
                              kmYlab = "Event-free Survival", 
                              forestTitle = "    Proportional Hazards Model",
                              termDict = NULL,
                              rel_height = 0.2){
  criterium <- "J"
  grouping <- getStrata(coxModel,theseDataPatients[,get(eventVar)], criterium)
  th <- getROCThreshold(predict(coxModel,type = "risk",reference="sample"), theseDataPatients[,get(eventVar)], criterium)
  
  km <- surv_fit(Surv(get(timeToEventVar), get(eventVar)) ~ grouping, data = theseDataPatients)
  
  medianTable <- surv_median(km)
  medianDiff <-  abs(diff(medianTable$median))
  
  if(is.na(medianDiff)){
    setDT(medianTable)
    strataWoM <- medianTable[is.na(median),unique(strata)]
    medianDiff <- paste0(">",max(theseDataPatients[grouping==gsub(pattern = "grouping=",replacement = "",x = strataWoM),TimeToEvent_mo]) - medianTable[strata!=strataWoM,median])
  }
  
  print(paste0("KM delta median time to event: ", medianDiff,", medians: ", paste0(collapse=",",medianTable$median)))
  
  kmPlot <- getKMPlot(km, threshold = round(th,digits = 1), labelConstant = "RR", ylab = kmYlab, colors = theseColors, title = kmTitle)
  forestPlot <- myggforest(coxModel,
                           data = theseDataPatients,
                           main = forestTitle,
                           fontsize = largeFontSize,
                           relFontSizeAnno = fontSize/largeFontSize,
                           pointSize = forestPointSize,
                           cpositions = c(0.02, 0.26, 0.4),
                           noDigits = 2,
                           termDict = termDict)
  
  plot1 <- forestPlot
  
  plot2 <- plot_grid(kmPlot$plot + theme(plot.margin = unit(c(0.01, 0, 0, 0.04), "npc")),
                     kmPlot$table + theme(plot.margin = unit(c(0, 0, 0, 0.04), "npc")),
                     ncol=1,
                     rel_heights = c(1,rel_height),
                     align = 'v') #cowplot to the rescue!
  
  return(list(plot1,plot2))
}

#' Returns a forest plot with the regressors of a Cox model, and a Kaplan-Meier plot (using ggsurvplot) stratifying the patients in two groups
#' It extends getForestAndKMCox to do the patient selection and calculate the model. It also removes variables selected by the LASSO that generate problems in the cox regression (uncommon)
#' 
#' @param imputedData the imputed patient data (here all variables may be imputed and scaled)
#' @param originalData the patient data without imputation to select patients based on missing data in the variables of interest. Variables can be missing or imputed on purpose not to take them into consideration 
#' @param freqs LASSO frequencies
#' @param minFreq minimum LASSO frequency for a variable to be included in the Cox model
#' @param timeVar name of the column in theseData that contains the time to event
#' @param eventVar name of the column in theseData that contains the binary variable with the outcome
#' @param colors of the two strata/patient groups
#' @param kmTitle title of the km plot
#' @param forestTitle Title of the forest plot 
#' @param termDict vector of human-readable terms named with the original terms for re-labeling
#' @param rel_height Relative height of the survival table
#' @param sigOnly Logical. Only include variables that are significant in the plot 
#' @param originalSDs Named vector with the SDs of the imputedData before scaling to print re-scaled Hazard Ratios (not in plots)
#' @return list of a forest plot and a km plot
getForestAndKMCoxFromData <- function(imputedData,originalData,freqs,minFreq,eventVar,colors,kmTitle,forestTitle,termDict,rel_height,timeVar="TimeToEvent_mo",sigOnly = F,originalSDs=NULL){
  removeFreq <- NULL

  repeat{
    freqs <- freqs[!var%in%removeFreq]
    selectedData <- imputedData[getFullDataPatients(originalData,freqs,minFreq)]
    thisFormula <- getFormulaLassoFreq(freqs,paste0("Surv(",timeVar,",",eventVar,")"),minFreq)
    if(!is.null(thisFormula)){
      removeFreq <- tryCatch(
        {
          thisModel <- coxph(thisFormula,selectedData)
          NULL
        },
        error = function(cond){
          stop(cond)
        },
        warning = function(cond){
          thisModel <- coxph(thisFormula,selectedData)
          removeCoeff <- as.numeric(gsub(pattern = ".*before variable  ([0-9]*).*",replacement = "\\1",x = cond))
          warning(paste0("Removing problematic parameter ",names(thisModel$coefficients[removeCoeff])))
          names(thisModel$coefficients[removeCoeff])
        })
      
    } else {
      return(NULL)
    }
    if(is.null(removeFreq)){break}
  }
  
  if(sigOnly) {
    thisModel <- getModelSignificantOnlyModel(thisModel,selectedData)
  }
  
  if(!is.null(originalSDs)){
    suppressWarnings(
    print(rbindlist(list(melt(data.table(t(rescaleHR(exp(thisModel$coefficients),originalSDs[names(thisModel$coefficients)]))),value.name = "HR")[,`:=`(type="Rescaled")]
              ,melt(data.table(t(exp(thisModel$coefficients))),value.name = "HR")[,`:=`(type="Scaled")])))
    )
  }
  
  thisPlot <- getForestAndKMCox(theseDataPatients = selectedData,
                                coxModel = thisModel,
                                timeToEventVar = timeVar,
                                eventVar = eventVar,
                                theseColors = colors,
                                kmTitle = kmTitle,
                                forestTitle = forestTitle,
                                termDict = termDict,
                                rel_height = rel_height)
  
  return(thisPlot)
}


#' Returns a Kaplan-Meier plot (using ggsurvplot) from a univariate Cox regression by stratifying the patients in two groups 
#' 
#' @param theseData the patient data
#' @param thesePatients list of patients to subset theseData
#' @param coxModel cox model that contains the variables to include from theseData
#' @param timeToEventVar name fo the column in theseData that contains the time to event
#' @param eventVar name of the column in theseData that contains the binary variable with the outcome
#' @param theseColors of the two strata/patient groups
#' @param ylab y-lab of the km plot
#' @param xname name for the threshold variable in the plot
#' @param title title of the km plot
#' @return list of a forest plot and a km plot
getKMPlotUnivariateCox <- function(theseData, thesePatients, coxModel, timeToEventVar, eventVar, theseColors, 
                                   title = "KM plot", 
                                   ylab = "Event-free Survival",
                                   xname = NULL){
  
  if(length(coxModel$coefficients) > 1)
    stop("This function is intended to work with univariate Cox regressions only")
  
  if(is.null(xname)) {
    xname <- names(coxModel$coefficients)
  }
  
  criterium <- "J"
  grouping <- getStrata(coxModel,theseData[thesePatients,get(eventVar)], criterium)
  
  th <- getROCThreshold(predict(coxModel,type = "risk",reference="sample"), theseData[thesePatients,get(eventVar)], criterium)
  
  #There is probably a better way of doing this, but this works
  xScaleTh <- optimize(function(x){
    return(abs(th - predict(coxModel,type = "risk",reference="sample",newdata = setnames(data.table(c(x)),names(coxModel$coefficients)))))},
    lower = min(theseData[thesePatients,get(names(coxModel$coefficients))]), 
    upper = max(theseData[thesePatients,get(names(coxModel$coefficients))]))
  
  km <- surv_fit(Surv(get(timeToEventVar), get(eventVar)) ~ grouping, data = theseData[thesePatients,])
  
  medianTable <- surv_median(km)
  medianDiff <-  abs(diff(medianTable$median))
  
  if(is.na(medianDiff)){
    setDT(medianTable)
    strataWoM <- medianTable[is.na(median),unique(strata)]
    medianDiff <- paste0("> ", max(theseData[thesePatients,][grouping==gsub(pattern = "grouping=",replacement = "",x = strataWoM),TimeToEvent_mo]) - medianTable[strata!=strataWoM,median])
  }
  
  print(paste0("KM delta median time to event: ", medianDiff,", medians: ", paste0(collapse=",",medianTable$median)))
  
  kmPlot <- getKMPlot(km, threshold = round(xScaleTh$minimum,digits = 1), labelConstant = xname, ylab = ylab, colors = theseColors, title = title)
  
  # finalPlot <- plot_grid(kmPlot$plot + theme(plot.margin = unit(c(0.01, 0, 0, 0.04), "npc")),
  #                    kmPlot$table + theme(plot.margin = unit(c(0, 0, 0, 0.04), "npc")),
  #                    ncol=1,
  #                    rel_heights = c(1,0.25),
  #                    align = 'v') #cowplot to the rescue!
  # 
  return(kmPlot)
}

#' Generates a data.table with a summary of the frequency at which each dependent variable was selected using LASSO
#' It runs in parallel using mc.cores cores
#' 
#' @param seed random number generator seed
#' @param x x matrix
#' @param y response y
#' @param lambda LASSO's lambda value to select
#' @param ... Passed to cv.glmnet, include here the family, etc.
#' @return list with a data.table with variables and coefficients and an error element
getLassoFreq=function(seed=theSeed,n=100,x,y,lambda=c("lambda.1se","lambda.min"),...){
  set.seed(seed)
  lambda=match.arg(lambda)
  theColnames=c("(Intercept)",colnames(x))
  theTable=data.table(t(vector(mode="numeric",length=length(theColnames))))
  setnames(theTable,theColnames)
  
  #the loop/lapply here
  theResults=rbindlist(mclapply(seq(1,n),mc.cores=mc.cores,FUN=function(i){
    fitModel=cv.glmnet(x=as.matrix(x),y=y,...)
    theTables=coef(fitModel,s=lambda)
    if(class(theTables)=="list"){
      theTables=theTables[[1]]
    }
    toAdd1=rownames(theTables)[theTables[,1]!=0]
    data.table(i=i, var=toAdd1)
  }))
  
  return(theResults[,.(pModel=.N/n),by=var][order(-pModel)])
}

#' Generates a formula including all selected variables depending on their frequency after running getLassoFreq
#' 
#' @param x data.table resulting from getLassoFreq
#' @param yname name of the y variable 
#' @param minLassoFreq minimum frequency a variable must have to be included in the formula
#' @return formula to use in a different lm 
getFormulaLassoFreq=function(x,yname,minLassoFreq=0.1){
  passVars <- x[var!="(Intercept)" & pModel >= minLassoFreq,var]
  if(length(passVars)>0){
    return(as.formula(paste0(yname," ~ ",paste(collapse=" + ",x[var!="(Intercept)" & pModel >= minLassoFreq,var]))))
  } else {
    return(NULL)
  }
}

#' Returns a list of patients that have full data for the variables selected from freqData that are present in patientData
#' 
#' @param patientData data.table 
#' @param freqData data.table resulting from getLassoFreq
#' @param minfreq minimum frequency a variable must have to be included in the formula
#' @param freqName name of the column containing frequency information in freqData
#' @param varName name of the column containing 
#' @return formula to use in a different lm 
getFullDataPatients <- function(patientData, freqData, minfreq, freqName = "pModel", varName = "var", patientName = "patient") {
  selectedVars <- freqData[get(freqName)>=minfreq,get(varName)]
  finalVars <- intersect(selectedVars, colnames(patientData))
  if(length(finalVars) < 1) {
    stop("No overlap between variables selected in freqData and patientData")#, so all patients in patientData will be returned")
  }
  if(length(finalVars) != length(selectedVars)){
    warning(sprintf("Not all selected variables were found in the patient data. Missing variables: %s",paste(collapse = ", ",selectedVars[!selectedVars %in% finalVars])))
  }
  return(patientData[,.(patient = get(patientName),complete=all(!is.na(.SD))),.SDcols=finalVars,by=.I][complete==T,patient])
}

#' Reformulates a coxph keeping significant regressors only
#' 
#' @param x coxph model
#' @param data data to run the new model. If null, original data
#' @return list with a data.table with variables and coefficients and an error element
getModelSignificantOnlyModel=function(x,data=NULL){
  if(is.null(x)) return(NULL)
  if(!is(x,"coxph")) stop("This function only works with coxph models")
  pvals=summary(x)$coefficients[,"Pr(>|z|)"]
  rhs = paste(collapse=" + ", names(pvals[pvals<=alpha]))
  newformula=copy(x$formula)
  newformula=reformulate(rhs,newformula[[2]])
  if(is.null(data)) data=x$call$data
  return(do.call("coxph",list(formula = newformula, data = data)))
}

#' Subsets a table of LASSO frequencies given their Cox significance
#' 
#' @param x coxph model
#' @param freqs LASSO frequencies
#' @return freqs of only significant variables
getFreqsSignificantOnlyModel=function(x,freqs){
  if(is.null(x)) return(NULL)
  if(!is(x,"coxph")) stop("This function only works with coxph models")
  pvals <- summary(x)$coefficients[,"Pr(>|z|)"]
  vars <- pvals[pvals<=alpha]
  if(nrow(freqs[var %in% names(vars),])==0) return(NULL)
  return(freqs[var %in% names(vars),])
}

#' Re-scales a Hazard Ratio from an scaled variable to the original scale
#' @param hr hazard ratio (exp(b))
#' @param originalSD sd of the original (unscaled) variable
#' @param perXOriginal number to divide the original scale by to re-calculate new units (e.g., per 10 SNVs)
#' @return rescaled HR
rescaleHR <- function(hr,originalSD,perXOriginal=1){
  exp(log(hr)*perXOriginal/originalSD)
}

##########

###Data parsing
###############
configFile <- paste(sep="/",Sys.getenv("ManuscriptScripts_DCISRecurrenceVsProgression"),"configFile")
if(!file.exists(configFile)){
    stop("Configuration file configFile not found. Edit configFile.in to adapt it to your system, save it as configFile, and export the ManuscriptScripts_DCISRecurrenceVsProgression environment variable before running this script")
}
source(configFile)

###Cross-sectional
##SNV
aim1GeneticData=fread(file = aim1GeneticDataFile)
aim1GeneticData[,`:=`(Type=factor(Type,levels=c("Pure","Synchronous","Invasive")))]

##Phenotypic
#Intensity
aim1IntScoreData <- readRDS(file = aim1IntScoreData)
aim1IntScoreData[,Type:=factor(aim1IntScoreData$Kind,levels = c("Pure DCIS","Synchronous DCIS"),labels = c("Pure","Synchronous"))]
aim1IntScoreData[,Kind:=NULL]
setnames(aim1IntScoreData,old="DCIS",new="Patient")

#Supplementary: 1 data point per patient at random
aim1IntRNDScoreData  <- readRDS(file = aim1IntRNDScoreDataFile)
aim1IntRNDScoreData[,Type:=factor(aim1IntRNDScoreData$Kind,levels = c("Pure DCIS","Synchronous DCIS"),labels = c("Pure","Synchronous"))]
aim1IntRNDScoreData[,Kind:=NULL]
setnames(aim1IntRNDScoreData,old="DCIS",new="Patient")

#Diversity
aim1DivergenceDataWide = readRDS(file = aim1DivergenceDataFile)
aim1DivergenceScores = c("EMD_HER2", "EMD_GLUT1", "EMD_ER", "EMD_FOXP3", "CDI_CA9" , "CDI_ER")
aim1DivergenceIDCols = c("Kind","DCIS")
aim1DivergenceDataLong=melt(aim1DivergenceDataWide,measure.vars = aim1DivergenceScores,id.vars = aim1DivergenceIDCols,variable.name = "Category",value.name = "Score")[!is.na(Kind),]
aim1DivergenceDataLong[,`:=`(Type=factor(aim1DivergenceDataLong$Kind,levels = c("Pure DCIS","Synchronous DCIS"),labels = c("Pure","Synchronous")),Marker=gsub("[^_]*_","",Category),Test=gsub("([^_]*).*","\\1",Category))]
aim1DivergenceDataLong[,Kind:=NULL]
aim1DivergenceDataLong[Test=="CDI",`:=`(Score=1-Score)] ##So that we are measuring Divergence/heterogeneity instead of similarity

#Supplementary: all markers
aim1EMDistData <- readRDS(file = aim1EMDistDataFile)
setnames(aim1EMDistData,old=c("DCIS","Kind","EMDist"),new=c("Patient","Type","Score"))
aim1EMDistData[,`:=`(Type=factor(Type,levels = c("Pure DCIS","Synchronous DCIS"),labels = c("Pure","Synchronous")))]

aim1CDIScoreData <- readRDS(file = aim1CDIScoreDataFile)
setnames(aim1CDIScoreData,old=c("DCIS","Kind","cdiScore"),new=c("Patient","Type","Score"))
aim1CDIScoreData[,`:=`(Type=factor(Type,levels = c("Pure DCIS","Synchronous DCIS"),labels = c("Pure","Synchronous")),Score=1-Score)]

###Longitudinal
##SNV
aim4SNVData=fread(file = aim4SNVDataFile)
aim4SNVData[,`:=`(Cohort=factor(Cohort))]
setnames(aim4SNVData,old=c("Sample","SNV","Divergence"),new=c("patient","SNVBurden","SNVDivergence"))
setkey(aim4SNVData,"patient")
#CNA
aim4CNADataByPatient = readRDS(file = aim4CNADataByPatientFile)[,.(patient,Cohort,CNADivergence=1-Similarity,MeanAlteredP,AlteredP)]
aim4CNADataBySample = readRDS(file = aim4CNADataBySampleFile)[,.(patient,Cohort,Sample,AlteredP)]

##Phenotypic
aim4ERData = readRDS(file = aim4ERDataFile)
aim4ERData[,`:=`(OriginalCohort=Cohort,Cohort=factor(Cohort,levels=longNames,labels = names(longNames)))]
setnames(aim4ERData,old="intScore1",new="intScore")

aim4GLUTData = readRDS(file = aim4GLUTDataFile)
aim4GLUTData[,`:=`(OriginalCohort=Cohort,Cohort=factor(Cohort,levels=longNames,labels = names(longNames)))]

aim4ERPosERData = readRDS(file = aim4ERPosERDataFile)
aim4ERPosERData[,`:=`(OriginalCohort=Cohort,Cohort=factor(Cohort,levels=longNames,labels = names(longNames)))]
setnames(aim4ERPosERData,old="intScore1",new="intScore")

aim4ERPosGLUTData = readRDS(file = aim4ERPosGLUTDataFile)
aim4ERPosGLUTData[,`:=`(OriginalCohort=Cohort,Cohort=factor(Cohort,levels=longNames,labels = names(longNames)))]

aim4PhenotypicHeterData = readRDS(file = aim4PhenotypicHeterDataFile)

#Extracting only the phenotypic data from this data file
aim4PhenotypicHeterData = aim4PhenotypicHeterData[,.(patient=Patient,EMD_GLUT1,MIS_GLUT1,MIS_ER,MSE_ER,CDI_ER)]
aim4PhenotypicHeterData <- aim4PhenotypicHeterData[,`:=`(toRemove=all(is.na(.SD))),.SDcols=c("MIS_GLUT1","MIS_ER"),by=.I][toRemove==F,!"toRemove"]
setkey(aim4PhenotypicHeterData,"patient")

##Clinical
theClinicalData=fread(file = theClinicalDataFile)
theClinicalData=theClinicalData[,1:42] #WARNING HARDCODED: Only DCIS-related variables
theClinicalData[,`:=`(patient=gsub("N?RZ?$","",PatientID))]
setkey(theClinicalData,patient)

#Fixing some variables
#Logic
logicalCols=c("HispanicOrLatino","Censored","DCISWasAtSurgicalMargin","Necrosis","MicroCalcs","DCISMicrocalc","BenignMicrocalc","ER","PR","HER2_IHC","RadiationDCIS","HormonalTherapyDCIS","TamoxifenDCIS")
theClinicalData[,(logicalCols):=lapply(.SD,as.logical),.SDcols=logicalCols]
#Factor
factorCols=c("MenopausalStatus","NuclearGrade","Race","RaceN","Cohort","LastContactStatus","SurgicalProcedure","AxillaryDissection","Laterality","Quadrant","FinalMargin","NecrosisType","MicrocalcType","HormonalTherapyType")
theClinicalData[,(factorCols):=lapply(.SD,as.factor),.SDcols=factorCols]
#Numeric
numericCols=c("AgeAtDCISDiagnosis","Height_cm","Weight_kg","DCISDiagnosis_y","TimeToEvent_mo","LastContact_y","SurgicalProcedure_y","NNodesExaminedN","EstimatedMargin","DCISSize_cm")
theClinicalData[,(numericCols):=lapply(.SD,as.numeric),.SDcols=numericCols]
#Calculating some variables
#Events
theClinicalData[,`:=`(EventRecurrence=!Censored)]
theClinicalData[,`:=`(EventProgression=EventRecurrence)]
theClinicalData[Cohort==1,`:=`(EventProgression=F)]
#BMI
theClinicalData[,`:=`(BMI=Weight_kg/((Height_cm/100)^2))]
#MenopausalStatus
theClinicalData[,`:=`(MenopausalStatusCompleted=MenopausalStatus)]
theClinicalData[is.na(MenopausalStatus),`:=`(MenopausalStatusCompleted=ifelse(AgeAtDCISDiagnosis>=55,"PostMenopausal","PreMenopausal"))]
#DCISTreatment
theClinicalData[,`:=`(DCISTreatment=as.character(SurgicalProcedure))]
theClinicalData[SurgicalProcedure=="Lumpectomy" & RadiationDCIS==T,`:=`(DCISTreatment="LumpectomyRadiation")]
theClinicalData[SurgicalProcedure=="Lumpectomy" & RadiationDCIS==F,`:=`(DCISTreatment="LumpectomyOnly")]
theClinicalData[DCISTreatment=="Lumpectomy",`:=`(DCISTreatment=NA)]
theClinicalData[,`:=`(DCISTreatment=as.factor(DCISTreatment))]
#EstimatedMargin recalculation
theClinicalData[,`:=`(Margin_mm_Fixed=as.numeric(Margin_mm))]
theClinicalData[grep("[<->]",Margin_mm),.(patient,Margin_mm)] #WARNING HARDCODED: These cases were coded as <2.
theClinicalData[patient=="DCIS-317",`:=`(Margin_mm_Fixed=1)] #WARNING HARDCODED: This case was coded as <2. We are fixing it as 1. Making it NA does not change any results qualitatively
theClinicalData[patient=="DCIS-372",`:=`(Margin_mm_Fixed=1)] #WARNING HARDCODED: This case was coded as <2. We are fixing it as 1. Making it NA does not change any results qualitatively
theClinicalData[,`:=`(EstimatedMargin=NULL)][,`:=`(EstimatedMargin=Margin_mm_Fixed)][is.na(EstimatedMargin) & !is.na(FinalMargin),`:=`(EstimatedMargin=sapply(FinalMargin,FUN=function(x){if(x=="Final margin less than 1mm"){estimatedMarginC1}else if (x=="Final margin between 1-10mm"){estimatedMarginC2}else{estimatedMarginC3}}))]
#Adequate margin >=2mm. Margin_mm has priority because it was re-measured (better/most updated)
theClinicalData[,`:=`(AdequateMargin=ifelse(Margin_mm_Fixed>=2,T,F))][is.na(AdequateMargin) & FinalMargin!="Final margin between 1-10mm",`:=`(AdequateMargin=ifelse(FinalMargin=="Final margin greater than 10mm",T,F))]


#All variables that are clean enough to be used for LASSO, if extremely correlated, only one and the others commented out. If recoded version of a raw covariate, annotated but not included
clinicalVariablesThatMayBeConsidered=c("AgeAtDCISDiagnosis", 
                                       "MenopausalStatusCompleted", #recoded MenopausalStatus to use the Age at Diagnosis to reduce the missing data (threshold = 55 years, only for those without this information)
                                       "RaceN",#recoded Race
                                       "Height_cm", 
                                       "Weight_kg",
                                       "BMI",
                                       "AxillaryDissection",
                                       "NNodesExaminedN",# recoded NNodesExamined 
                                       "DCISWasAtSurgicalMargin",
                                       "EstimatedMargin", #recoded "FinalMargin" using estimatedMarginCX variables + Margin_mm
                                       "DCISSize_cm",
                                       "NuclearGrade",
                                       "NecrosisType",
                                       #"Necrosis",
                                       "MicrocalcType",
                                       #"MicroCalcs",
                                       #"DCISMicroCalcs",
                                       #"BenignMicroCalcs",
                                       "ER",
                                       "PR",
                                       "DCISTreatment",
                                       #"SurgicalProcedure", #Recoded as DCISTreatment
                                       #"RadiationDCIS", #Recoded as DCISTreatment
                                       "HormonalTherapyDCIS",
                                       #"TamoxifenDCIS",
                                       #"HormonalTherapyType",
                                       #To discard below this line, reason as a comment
                                       "DCISDiagnosis_y", #Discarded because we expect it to be associated with the outcome, since the study finished at a given time
                                       #"SurgicalProcedure_y", #chosen DCISDiagnosis_y, r=0.99 between the two
                                       "HispanicOrLatino", #too few to take into account (7 only) 
                                       "Laterality", #In principle not important
                                       "Quadrant" #~41% of missing data
)
clinicalVariablesForNow=clinicalVariablesThatMayBeConsidered[1:(length(clinicalVariablesThatMayBeConsidered)-4)]

###WARNING: Deleting patient information from these data.tables to make sure we are using the most updated patient information from clinicalData
#####################################################################################
for (theseData in list(aim4ERData,aim4GLUTData,aim4ERPosERData,aim4ERPosGLUTData)){
  theseData[,`:=`(Cohort=NULL,Recurrence=NULL,Progression=NULL,OriginalCohort=NULL)]
}
aim4ERData <- merge(aim4ERData,theClinicalData[,.(Patient=patient,Cohort=Cohort,OriginalCohort=longNames[Cohort])],by="Patient")
aim4GLUTData <- merge(aim4GLUTData,theClinicalData[,.(Patient=patient,Cohort=Cohort,OriginalCohort=longNames[Cohort])],by="Patient")
aim4ERPosERData <- merge(aim4ERPosERData,theClinicalData[,.(Patient=patient,Cohort=Cohort,OriginalCohort=longNames[Cohort])],by="Patient")
aim4ERPosGLUTData <- merge(aim4ERPosGLUTData,theClinicalData[,.(Patient=patient,Cohort=Cohort,OriginalCohort=longNames[Cohort])],by="Patient")

##Combining the datsetes for multivariable modeling
allAim4Data=merge(merge(aim4SNVData[,.(patient,SNVBurden,SNVDivergence)],aim4CNADataByPatient[,.(patient,CNADivergence,MeanAlteredP,AlteredP)],by="patient",all=T),aim4PhenotypicHeterData,by="patient",all=T)
allAim4DataWClinical <- theClinicalData[allAim4Data]
allAim4Data <- merge(allAim4Data,allAim4DataWClinical[,.(patient,Cohort,EventRecurrence,EventProgression)])

humanReadableHotEncodedVariables=setNames(c('SNV Burden',
                                            "SNV Divergence",
                                            "Height",
                                            "Surgical Margin",
                                            "ER+",
                                            "Lumpectomy Only",
                                            "GLUT1 Intensity",
                                            "GLUT1 Divergence",
                                            "DCIS microcalcs",
                                            "Race: Other",
                                            "CNA Divergence",
                                            "Race: Black",
                                            "Adequate Margin",
                                            "No microcalcs"),
                                          c("SNVBurden",
                                            "SNVDivergence",
                                            "Height_cm",
                                            "EstimatedMargin",
                                            "ERTRUE",
                                            "DCISTreatmentLumpectomyOnly",
                                            "MIS_GLUT1",
                                            "EMD_GLUT1",
                                            "MicrocalcTypeDCIS",
                                            "RaceNOther",
                                            "CNADivergence",
                                            "RaceNBlack",
                                            "AdequateMarginTRUE",
                                            "MicrocalcTypeNone")
)

###FIGURES###
#############

#Figure 2
#########
dataFig1=dcast(aim1GeneticData,value.var="SNV",Patient~Type)[,.(Patient,diffSNVs=Invasive-Synchronous)][!is.na(diffSNVs)][aim1GeneticData[,.SD,keyby=Patient],]

{
  ylimits=c(0,dataFig1[,max(SNV)]*1.1)#WARNING HARDCODED with padding for geom_signif
  Ns=dataFig1[,.N,keyby=Type][,`:=`(label=paste(sep="\n",crossNamesSplit[as.character(Type)],N))][]
  thisLabels=Ns[,label]
  names(thisLabels)=Ns[,Type]
  
  fig2A=ggplot(dataFig1[Type!="Invasive",],aes(x=Type,y=SNV,fill=Type)) +
    geom_violin(color=NA) +
    #geom_boxplot(width=0.1,outlier.shape = NA,show.legend=F) +
    geom_quasirandom(alpha=1,size=quasirandomSize,shape=20,aes(color=Type),show.legend = F) +
    stat_summary(fun.min=function(z) { quantile(z,0.25) },fun.max=function(z) { quantile(z,0.75) },fun=median,geom="pointrange",color=pointRangeColor,show.legend=F,linewidth=pointRangeLineWidth,size=pointRangeSize) +
    # geom_signif(comparisons=list(c("Pure","Synchronous")),
    #             test="wilcox.test",
    #             map_signif_level = function(p){sprintf("%.1g", p)},
    #             textsize = ggsignifTextSize,
    #             vjust = ggsignifTextVjust) +
    scale_y_continuous(name="# SNVs per patient") +#,limits=ylimits) +
    scale_x_discrete(name="",labels=thisLabels) +
    scale_fill_manual(values=crossColors)+
    scale_color_manual(values=crossColors)+
    theme(legend.position="none") +
    labs(title = "DCIS SNV Burden")
  
  # Test for the manuscript, without the 4 highly-altered cases the differences are not significant anymore
  # outliers <- dataFig1[Type!="Pure",.(Type=factor(Type),SNV,Patient,diffSNVs)][!is.na(diffSNVs),] %>% filter(Type=="Invasive") %>% arrange(desc(SNV)) %>% top_n(4) %>% pull(Patient)
  # pvalWithoutOutliers <- dataFig1[Type!="Pure",.(Type=factor(Type),SNV,Patient,diffSNVs)][!is.na(diffSNVs),] %>% anti_join(outliers,by = c("Patient")) %>%
  #   sign_test(SNV ~ Type) %>%
  #   add_significance() %>%
  #   mutate(p = sprintf("%.1g", p)) %>%
  #   pull(p)
  
  pval=dataFig1[Type!="Pure",.(Type=factor(Type),SNV,Patient,diffSNVs)][!is.na(diffSNVs),] %>%
    sign_test(SNV ~ Type) %>%
    add_significance() %>%
    mutate(p = sprintf("%.1g", p)) %>%
    pull(p)
  
  ylimits=c(0,dataFig1[Type=="Invasive",max(SNV)]*1.1)#WARNING HARDCODED with padding for geom_signif
  Ns=dataFig1[,.SD,key=Patient][dataFig1[Type=="Invasive",.(Patient),],][,.N,keyby=Type][,`:=`(label=paste(sep="\n",crossNamesSplit[as.character(Type)],N))][]
  thisLabels=Ns[,label]
  names(thisLabels)=Ns[,Type]
  
  fig2B=ggplot(dataFig1[,.SD,key=Patient][dataFig1[Type=="Invasive",.(Patient),],],aes(x=Type,y=SNV,fill=Type)) +
    geom_violinhalf(show.legend=F,flip=c(1,0),color=NA) +
    geom_line(aes(group=Patient),size = geomLineWidth) +
    #geom_line(aes(group=Patient,color=diffSNVs),linewidth=1) +
    #geom_boxplot(width=0.1,outlier.shape = NA,show.legend=F) +
    geom_point(alpha=1,size=quasirandomSize,shape=20,show.legend=F,) +
    stat_summary(fun.min=function(z) { quantile(z,0.25) },fun.max=function(z) { quantile(z,0.75) },fun=median,geom="pointrange",color=pointRangeColor,show.legend=F,linewidth=pointRangeLineWidth,size=pointRangeSize) +
    geom_signif(comparisons=list(c("Synchronous","Invasive")),
                annotations=c(pval),
                textsize = ggsignifTextSize,
                vjust = ggsignifTextVjust) +
    scale_y_continuous(name="# SNVs per stage",limits=ylimits) +
    scale_x_discrete(name="",labels=thisLabels[c(2,3,1)]) +
    scale_fill_manual(values=crossColors)+
    #scale_color_viridis_c(name="Difference",option="magma") +
    #scale_color_gradient2(name="Difference",low=lowIntensintyColor,high=highIntensityColor) +
    labs(title = "DCIS vs. IDC SNV Burden")
  
  fig2B
  
  #Figure 2C
  dataFig3=dcast(aim1GeneticData,value.var="Divergence",Patient~Type)[,.(Patient,diffDivergences=Invasive-Synchronous)][!is.na(diffDivergences)][aim1GeneticData[,.SD,keyby=Patient],]
  
  ylimits=c(0,aim1GeneticData[Type!="Invasive",max(Divergence,na.rm = T)] * 1.1)#WARNING HARDCODED with padding for geom_signif
  Ns=dataFig3[Type!="Invasive",][!is.na(Divergence),.N,keyby=Type][,`:=`(label=paste(sep="\n",crossNamesSplit[as.character(Type)],N))][]
  thisLabels=Ns[,label]
  names(thisLabels)=Ns[,Type]
  
  fig2C=ggplot(dataFig3[Type!="Invasive",],aes(x=Type,y=Divergence,fill=Type)) +
    geom_violin(color=NA) +
    geom_quasirandom(alpha=1,size=quasirandomSize,shape=20,aes(color=Type)) +
    stat_summary(fun.min=function(z) { quantile(z,0.25) },fun.max=function(z) { quantile(z,0.75) },fun=median,geom="pointrange",color=pointRangeColor,show.legend=F,linewidth=pointRangeLineWidth,size=pointRangeSize) +
    geom_signif(comparisons=list(c("Pure","Synchronous")),
                test="wilcox.test",
                map_signif_level = function(p){sprintf("%.1g", p)},
                textsize = ggsignifTextSize,
                vjust = ggsignifTextVjust) +
    scale_y_continuous(name="Divergence (%)",limits=ylimits) +
    scale_x_discrete(name="",labels=thisLabels) +
    scale_fill_manual(values=crossColors)+
    scale_color_manual(values=crossColors)+
    theme(legend.position="none") +
    labs(title = "DCIS SNV Divergence")
  
  fig2=plot_grid(fig2A,fig2B,fig2C, labels=c("A","B","C"),nrow = 1,label_size = gridLabelFontSize)
  
  ncol <- 3
  nrow <- 1
  save_plot(fig2,file=paste0(outDir,"fig2.pdf"),base_width = pageWidthInch, base_height = pageWidthInch*nrow/ncol/aspectRatio.plot)
  save_plot(fig2,file=paste0(outDir,"fig2.png"),base_width = pageWidthInch, base_height = pageWidthInch*nrow/ncol/aspectRatio.plot)
}

# Supplementary figure
{
  pval=dataFig3[Type!="Pure",.(Type=factor(Type),Divergence,Patient,diffDivergences)][!is.na(diffDivergences),] %>%
    sign_test(Divergence ~ Type) %>%
    add_significance() %>%
    mutate(p = sprintf("%.1g", p)) %>%
    pull(p)
  
  ylimits=c(0,aim1GeneticData[Type!="Pure",max(Divergence,na.rm = T)] * 1.1)#WARNING HARDCODED with padding for geom_signif
  Ns=dataFig3[Type!="Pure",][!is.na(Divergence),.N,keyby=Type][,`:=`(label=paste(sep="\n",crossNamesSplit[as.character(Type)],N))][]
  thisLabels=Ns[,label]
  names(thisLabels)=Ns[,Type]
  
  figSIDCDivergence=ggplot(dataFig3[Type!="Pure",],aes(x=Type,y=Divergence,fill=Type)) +
    geom_violinhalf(show.legend=F,flip=c(1,0),color=NA) +
    geom_line(aes(group=Patient), size = geomLineWidth) +
    #geom_line(aes(group=Patient,color=diffDivergences),linewidth=1) +
    geom_point(alpha=1,size=quasirandomSize,shape=20,show.legend=F) +
    stat_summary(fun.min=function(z) { quantile(z,0.25) },fun.max=function(z) { quantile(z,0.75) },fun=median,geom="pointrange",color=pointRangeColor,show.legend=F,linewidth=pointRangeLineWidth,size=pointRangeSize) +
    geom_signif(comparisons=list(c("Synchronous","Invasive")),
                annotations=c(pval),
                textsize = ggsignifTextSize,
                vjust = ggsignifTextVjust) +
    scale_y_continuous(name="SNV Divergence (%)",limits=ylimits) +
    scale_x_discrete(name="",labels=thisLabels[c(2,3,1)]) +
    scale_fill_manual(values=crossColors)+
    #scale_color_viridis_c(name="Difference",option="magma") +
    labs(title = "Within vs. Between SNV Divergence")
  
  figSIDCDivergence
  
  ncol <- 1
  nrow <- 1
  save_plot(figSIDCDivergence,file=paste0(outDir,"figSIDCDivergence.pdf"),base_width = columnWidthInch, base_height = columnWidthInch*nrow/ncol/aspectRatio.plot)
  save_plot(figSIDCDivergence,file=paste0(outDir,"figSIDCDivergence.png"),base_width = columnWidthInch, base_height = columnWidthInch*nrow/ncol/aspectRatio.plot)
}
#Figure 3
#########

#Mod params
ggsignifPosModFig3 <- 1.3
topYlimModFig3 <- 2
pointRangeScaleFig3 <- 1 #Just in case we need to scale the size of pointRange point and range for very large plots (1 = unused)
#
{
  pwc = aim1IntScoreData %>%
    group_by(Marker) %>%
    #pairwise_wilcox_test( ##We can't correct for multiple tests using this function, because it is not adjusting per comparison (group by Marker) but per pairwise comparison (we only have one pair, Pure vs. Synchronous)
    wilcox_test(
      paired=F,
      as.formula(paste("intScore", "~", "Type"))) %>% 
    mutate(p.adj.signif=NULL) %>% 
    adjust_pvalue(method="holm") %>%
    arrange(p)
  
  ypos=aim1IntScoreData %>%
    group_by(Marker) %>%
    summarize(ypos=max(intScore))
  
  #To add the x_positions I need to do something kind of sketchy, first re-sorting the factor (Marker), which is fine, but then sorting the calculated xmin and xmax (which are not in the proper order, because they were calculated internally with the original data, saved in the attribute args)
  #Alternatively, I calculate the xpos here for any Marker, and then use the position to merge by row, this works better in my opinion
  xpos=pwc %>%
    add_x_position(x="Marker", dodge = dodgeWidth) %>%
    mutate(x = sort(x)) %>%
    mutate(xmin = sort(xmin)) %>%
    mutate(xmax = sort(xmax)) %>%
    select(x,xmin,xmax)
  
  sortedMarkers=pwc$Marker
  
  ylimits=c(min(aim1IntScoreData[intScore>0,intScore],na.rm = T),max(aim1IntScoreData$intScore,na.rm = T) * topYlimModFig3)
  thisLabels=as.character(pwc %>% mutate(label=paste0(Marker,"\n",n1,"  ",n2)) %>%pull(label))
  names(thisLabels)=as.character(pwc %>% pull(Marker))
  
  thisPwc = pwc %>%
    mutate(Marker=factor(Marker,levels=sortedMarkers))%>%
    drop_na() %>%
    inner_join(ypos,by="Marker") %>%
    mutate(x=row_number()) %>%
    inner_join(xpos,by="x") %>%
    filter(p<=alphaPlot) %>%
    mutate(ypos=ypos*ggsignifPosModFig3)
  
  thisLabels=as.character(pwc %>% mutate(label=paste0(Marker,"\n",n1,"  ",n2)) %>%pull(label))
  names(thisLabels)=as.character(pwc %>% pull(Marker))
  
  fig3AFull=ggplot(aim1IntScoreData[,c(.SD,.(sortedMarker=factor(Marker,levels = sortedMarkers)))],aes(x=sortedMarker,y=intScore,fill=Type))+
    geom_violin(position = position_dodge(width=dodgeWidth),color=NA) +
    geom_quasirandom(alpha=1,size=quasirandomSize,shape=20,dodge.width = dodgeWidth,aes(color=Type),show.legend = F) +
    stat_summary(fun.min=function(z) { quantile(z,0.25) },fun.max=function(z) { quantile(z,0.75) },fun=median,geom="pointrange",color=pointRangeColor,show.legend=F,linewidth=pointRangeLineWidth*pointRangeScaleFig3,size=pointRangeSize*pointRangeScaleFig3,position = position_dodge(width=dodgeWidth)) +
    geom_signif(y_position = log10(thisPwc$ypos),
                annotation=sprintf("%0.1g",thisPwc$p),
                xmin=thisPwc$xmin,
                xmax=thisPwc$xmax,
                tip_length = ggsignifTipLength,
                textsize = ggsignifTextSize,
                vjust = ggsignifTextVjust) + 
    scale_y_log10(name="Mean Intensity Score (MIS)",limits=ylimits) +
    scale_x_discrete(name="",labels=thisLabels) +
    scale_fill_manual(values=crossColors)+
    scale_color_manual(values=crossColors)+
    theme(legend.position="none") + ##WARNING HARDCODED
    labs(title="Phenotypic intensities")
  
  #Fig3A, only significant ones (without multiple-test correction!)
  thisSortedMarkers=sortedMarkers[pwc$p<=alpha]
  
  thisPwc = pwc %>% 
    mutate(Marker=factor(Marker,levels=thisSortedMarkers))%>%
    drop_na() %>%
    inner_join(ypos,by="Marker") %>%
    mutate(x=row_number()) %>%
    inner_join(xpos,by="x") %>%
    filter(p<=alphaPlot) %>%
    mutate(ypos=ypos*ggsignifPosModFig3)
  
  fig3A=ggplot(aim1IntScoreData[,c(.SD,.(sortedMarker=factor(Marker,levels = sortedMarkers)))][Marker%in%pwc[pwc$p<=alpha,]$Marker,],aes(x=sortedMarker,y=intScore,fill=Type))+
    geom_violin(position = position_dodge(width=dodgeWidth),color=NA) +
    geom_quasirandom(alpha=1,size=quasirandomSize,shape=20,dodge.width = dodgeWidth,aes(color=Type)) +
    stat_summary(fun.min=function(z) { quantile(z,0.25) },fun.max=function(z) { quantile(z,0.75) },fun=median,geom="pointrange",color=pointRangeColor,show.legend=F,linewidth=pointRangeLineWidth*pointRangeScaleFig3,size=pointRangeSize*pointRangeScaleFig3,position = position_dodge(width=dodgeWidth)) +
    geom_signif(y_position = log10(thisPwc$ypos),
                annotation=sprintf("%0.1g",thisPwc$p),
                xmin=thisPwc$xmin,
                xmax=thisPwc$xmax,
                textsize = ggsignifTextSize,
                vjust = ggsignifTextVjust) +
    scale_y_log10(name="Mean Intensity Score (MIS)",limits=ylimits) +
    scale_x_discrete(name="",labels=thisLabels) +
    scale_fill_manual(values=crossColors,labels=crossNamesSplit)+
    scale_color_manual(values=crossColors,labels=crossNamesSplit)+
    theme(legend.position = "none") ##WARNING HARDCODED
  
  #Fig3B, only non-significant ones (without multiple-test correction!)
  thisSortedMarkers=sortedMarkers[pwc$p>alpha]
  
  thisPwc = pwc %>% 
    mutate(Marker=factor(Marker,levels=thisSortedMarkers))%>%
    drop_na() %>%
    inner_join(ypos,by="Marker") %>%
    mutate(x=row_number()) %>%
    inner_join(xpos,by="x") %>%
    filter(p<=alphaPlot) %>%
    mutate(ypos=ypos*ggsignifPosModFig3)
  
  fig3B=ggplot(aim1IntScoreData[,c(.SD,.(sortedMarker=factor(Marker,levels = sortedMarkers)))][Marker%in%pwc[pwc$p>=alpha,]$Marker,],aes(x=sortedMarker,y=intScore,fill=Type))+
    geom_violin(position = position_dodge(width=dodgeWidth),color=NA) +
    geom_quasirandom(alpha=1,size=quasirandomSize,shape=20,dodge.width = dodgeWidth,aes(color=Type)) +
    stat_summary(fun.min=function(z) { quantile(z,0.25) },fun.max=function(z) { quantile(z,0.75) },fun=median,geom="pointrange",color=pointRangeColor,show.legend=F,linewidth=pointRangeLineWidth*pointRangeScaleFig3,size=pointRangeSize*pointRangeScaleFig3,position = position_dodge(width=dodgeWidth)) +
    geom_signif(y_position = log10(thisPwc$ypos),
                annotation=sprintf("%0.1g",thisPwc$p),
                xmin=thisPwc$xmin,
                xmax=thisPwc$xmax,
                textsize = ggsignifTextSize,
                vjust = ggsignifTextVjust) +
    scale_y_log10(name=NULL, limits=ylimits) +
    scale_x_discrete(name="",labels=thisLabels) +
    scale_fill_manual(values=crossColors,labels=crossNamesSplit)+
    scale_color_manual(values=crossColors,labels=crossNamesSplit)+
    guides(y="none")# +
  # theme(legend.text = element_text(size=9),
  #       legend.title = element_text(size=11),
  #       axis.text.x = element_text(size=rel(0.6)),) ##WARNING HARDCODED
  
  fig3B
  
  fig3R1=plot_grid(fig3A,fig3B,labels=c("A","B"),rel_widths = c(1,3),align = 'h',label_size = gridLabelFontSize)
  fig3R1
  
  #Figure 3R2
  ###########
  
  #Figure 3R2C
  pwc = aim1DivergenceDataLong %>%
    filter(Test=="EMD") %>%
    group_by(Marker) %>%
    #pairwise_wilcox_test( ##We can't correct for multiple tests using this function, because it is not adjusting per comparison (group by Marker) but per pairwise comparison (we only have one pair, Pure vs. Synchronous)
    wilcox_test(
      paired=F,
      as.formula(paste("Score", "~", "Type"))) %>% 
    adjust_pvalue(method="holm") %>%
    arrange(p)
  
  ypos=aim1DivergenceDataLong %>%
    filter(Test=="EMD") %>%
    group_by(Marker) %>%
    summarize(ypos=max(Score,na.rm=T))
  
  xpos=pwc %>%
    add_x_position(x="Marker", dodge = dodgeWidth) %>%
    mutate(x = sort(x)) %>%
    mutate(xmin = sort(xmin)) %>%
    mutate(xmax = sort(xmax)) %>%
    select(x,xmin,xmax)
  
  sortedMarkers=pwc$Marker
  
  thisPwc = pwc %>% 
    mutate(Marker=factor(Marker,levels=sortedMarkers))%>%
    drop_na() %>%
    inner_join(ypos,by="Marker") %>%
    mutate(x=row_number()) %>%
    inner_join(xpos,by="x") %>%
    filter(p<=alphaPlot) %>%
    mutate(ypos=ypos*ggsignifPosModFig3)
  
  ylimits=c(min(aim1DivergenceDataLong[Test=="EMD" & Score > 0,Score],na.rm=T),max(aim1DivergenceDataLong[Test=="EMD" & Score > 0,Score],na.rm=T) * topYlimModFig3)
  thisLabels=as.character(pwc %>% mutate(label=paste0(Marker,"\n",n1,"  ",n2)) %>%pull(label)) #WARNING HARDCODED
  names(thisLabels)=as.character(pwc %>% pull(Marker))
  
  fig3R2C=ggplot(aim1DivergenceDataLong[Test=="EMD",c(.SD,.(sortedMarker=factor(Marker,levels = sortedMarkers)))],aes(x=sortedMarker,y=Score,fill=Type))+
    geom_violin(position = position_dodge(width=dodgeWidth),color=NA) +
    geom_quasirandom(alpha=1,size=quasirandomSize,shape=20,dodge.width = dodgeWidth,aes(color=Type)) +
    stat_summary(fun.min=function(z) { quantile(z,0.25) },fun.max=function(z) { quantile(z,0.75) },fun=median,geom="pointrange",color=pointRangeColor,show.legend=F,linewidth=pointRangeLineWidth*pointRangeScaleFig3,size=pointRangeSize*pointRangeScaleFig3,position = position_dodge(width=dodgeWidth)) +
    geom_signif(y_position = log10(thisPwc$ypos),
                annotation=sprintf("%0.1g",thisPwc$p),
                xmin=thisPwc$xmin,
                xmax=thisPwc$xmax,
                tip_length = ggsignifTipLength,
                textsize = ggsignifTextSize,
                vjust = ggsignifTextVjust) + 
    scale_y_log10(name="Between-sample Divergence (EMD)",limits=ylimits) +
    scale_x_discrete(name="",labels=thisLabels) +
    scale_fill_manual(values=crossColors)+
    scale_color_manual(values=crossColors)+
    theme(legend.position="none") +##WARNING HARDCODED
    labs(title="Between-sample Phenotypic Divergence")
  
  #Figure 3R2D
  pwc = aim1DivergenceDataLong %>%
    filter(Test=="CDI") %>%
    group_by(Marker) %>%
    #pairwise_wilcox_test( ##We can't correct for multiple tests using this function, because it is not adjusting per comparison (group by Marker) but per pairwise comparison (we only have one pair, Pure vs. Synchronous)
    wilcox_test(
      paired=F,
      as.formula(paste("Score", "~", "Type"))) %>% 
    adjust_pvalue(method="holm") %>%
    arrange(p)
  
  ypos=aim1DivergenceDataLong %>%
    filter(Test=="CDI") %>%
    group_by(Marker) %>%
    summarize(ypos=max(Score,na.rm=T))
  
  xpos=pwc %>%
    add_x_position(x="Marker", dodge = dodgeWidth) %>%
    mutate(x = sort(x)) %>%
    mutate(xmin = sort(xmin)) %>%
    mutate(xmax = sort(xmax)) %>%
    select(x,xmin,xmax)
  
  sortedMarkers=pwc$Marker
  
  thisPwc = pwc %>% 
    mutate(Marker=factor(Marker,levels=sortedMarkers))%>%
    drop_na() %>%
    inner_join(ypos,by="Marker") %>%
    mutate(x=row_number()) %>%
    inner_join(xpos,by="x") %>%
    filter(p<=alphaPlot) %>%
    mutate(ypos=ypos*ggsignifPosModFig3)
  
  ylimits=c(min(aim1DivergenceDataLong[Test=="CDI" & Score > 0,Score],na.rm=T),max(aim1DivergenceDataLong[Test=="CDI" & Score > 0,Score],na.rm=T) * topYlimModFig3)
  thisLabels=as.character(pwc %>% mutate(label=paste0(Marker,"\n",n1,"  ",n2)) %>%pull(label)) #WARNING HARDCODED
  names(thisLabels)=as.character(pwc %>% pull(Marker))
  
  fig3R2D=ggplot(aim1DivergenceDataLong[Test=="CDI",c(.SD,.(sortedMarker=factor(Marker,levels = sortedMarkers)))],aes(x=sortedMarker,y=Score,fill=Type))+
    geom_violin(position = position_dodge(width=dodgeWidth),color=NA) +
    geom_quasirandom(alpha=1,size=quasirandomSize,shape=20,dodge.width = dodgeWidth,aes(color=Type)) +
    stat_summary(fun.min=function(z) { quantile(z,0.25) },fun.max=function(z) { quantile(z,0.75) },fun=median,geom="pointrange",color=pointRangeColor,show.legend=F,linewidth=pointRangeLineWidth*pointRangeScaleFig3,size=pointRangeSize*pointRangeScaleFig3,position = position_dodge(width=dodgeWidth)) +
    geom_signif(y_position = log10(thisPwc$ypos),
                annotation=sprintf("%0.1g",thisPwc$p),
                xmin=thisPwc$xmin,
                xmax=thisPwc$xmax,
                tip_length = ggsignifTipLength,
                textsize = ggsignifTextSize,
                vjust = ggsignifTextVjust) + 
    scale_y_log10(name="Within-sample divergence (CDI)", limits=ylimits) +
    scale_x_discrete(name="",labels=thisLabels) +
    scale_fill_manual(values=crossColors,labels=crossNamesSplit) +
    scale_color_manual(values=crossColors,labels=crossNamesSplit) +
    labs(title="Within-sample Phenotypic Divergence")
  
  fig3R2D
  
  fig3R2=plot_grid(fig3R2C,fig3R2D,labels=c("C","D"),label_size = gridLabelFontSize)
  fig3R2
  
  fig3R1 <- plot_grid(fig3AFull,labels=c("A"),label_size = gridLabelFontSize)
  fig3R2 <- plot_grid(fig3R2C,fig3R2D,labels=c("B","C"),label_size = gridLabelFontSize)
  fig3 <- plot_grid(fig3R1,fig3R2,nrow = 2,label_size = gridLabelFontSize)
  ncol <- 3
  nrow <- 2
  save_plot(fig3,file=paste0(outDir,"fig3.pdf"),base_width = pageWidthInch, base_height = pageWidthInch*nrow/ncol/aspectRatio.plot)
  save_plot(fig3,file=paste0(outDir,"fig3.png"),base_width = pageWidthInch, base_height = pageWidthInch*nrow/ncol/aspectRatio.plot)
}

# Supplementary figures
{
  pwc = aim1EMDistData %>%
    group_by(Marker) %>%
    #pairwise_wilcox_test( ##We can't correct for multiple tests using this function, because it is not adjusting per comparison (group by Marker) but per pairwise comparison (we only have one pair, Pure vs. Synchronous)
    wilcox_test(
      paired=F,
      as.formula(paste("Score", "~", "Type"))) %>% 
    mutate(p.adj.signif=NULL) %>% 
    adjust_pvalue(method="holm") %>%
    arrange(p)
  
  ypos=aim1EMDistData %>%
    group_by(Marker) %>%
    summarize(ypos=max(Score))
  
  #To add the x_positions I need to do something kind of sketchy, first re-sorting the factor (Marker), which is fine, but then sorting the calculated xmin and xmax (which are not in the proper order, because they were calculated internally with the original data, saved in the attribute args)
  #Alternatively, I calculate the xpos here for any Marker, and then use the position to merge by row, this works better in my opinion
  xpos=pwc %>%
    add_x_position(x="Marker", dodge = dodgeWidth) %>%
    mutate(x = sort(x)) %>%
    mutate(xmin = sort(xmin)) %>%
    mutate(xmax = sort(xmax)) %>%
    select(x,xmin,xmax)
  
  sortedMarkers=pwc$Marker
  
  ylimits=c(min(aim1EMDistData[Score>0,Score],na.rm = T),max(aim1EMDistData$Score,na.rm = T) * topYlimModFig3)
  thisLabels=as.character(pwc %>% mutate(label=paste0(Marker,"\n",n1,"  ",n2)) %>%pull(label))
  names(thisLabels)=as.character(pwc %>% pull(Marker))
  
  thisPwc = pwc %>%
    mutate(Marker=factor(Marker,levels=sortedMarkers))%>%
    drop_na() %>%
    inner_join(ypos,by="Marker") %>%
    mutate(x=row_number()) %>%
    inner_join(xpos,by="x") %>%
    filter(p<=alphaPlot) %>%
    mutate(ypos=ypos*ggsignifPosModFig3)
  
  thisLabels=as.character(pwc %>% mutate(label=paste0(Marker,"\n",n1,"  ",n2)) %>%pull(label))
  names(thisLabels)=as.character(pwc %>% pull(Marker))
  
  fig3BFull=ggplot(aim1EMDistData[,c(.SD,.(sortedMarker=factor(Marker,levels = sortedMarkers)))],aes(x=sortedMarker,y=Score,fill=Type))+
    geom_violin(position = position_dodge(width=dodgeWidth),color=NA) +
    geom_quasirandom(alpha=1,size=quasirandomSize,shape=20,dodge.width = dodgeWidth,aes(color=Type),show.legend = F) +
    stat_summary(fun.min=function(z) { quantile(z,0.25) },fun.max=function(z) { quantile(z,0.75) },fun=median,geom="pointrange",color=pointRangeColor,show.legend=F,linewidth=pointRangeLineWidth*pointRangeScaleFig3,size=pointRangeSize*pointRangeScaleFig3,position = position_dodge(width=dodgeWidth)) +
    geom_signif(y_position = log10(thisPwc$ypos),
                annotation=sprintf("%0.1g",thisPwc$p),
                xmin=thisPwc$xmin,
                xmax=thisPwc$xmax,
                tip_length = ggsignifTipLength,
                textsize = ggsignifTextSize,
                vjust = ggsignifTextVjust) + 
    scale_y_log10(name="Between-sample Divergence (EMD)",limits=ylimits) +
    scale_x_discrete(name="",labels=thisLabels) +
    scale_fill_manual(values=crossColors)+
    scale_color_manual(values=crossColors)+
    theme(legend.position="none") + ##WARNING HARDCODED
    labs(title="Between-sample Phenotypic Divergence")
  
  ncol <- 2
  nrow <- 1
  save_plot(fig3BFull,file=paste0(outDir,"fig3BFull.pdf"),base_width = pageWidthInch, base_height = pageWidthInch*nrow/ncol/aspectRatio.plot)
  save_plot(fig3BFull,file=paste0(outDir,"fig3BFull.png"),base_width = pageWidthInch, base_height = pageWidthInch*nrow/ncol/aspectRatio.plot)
  
  pwc = aim1CDIScoreData %>%
    group_by(Marker) %>%
    #pairwise_wilcox_test( ##We can't correct for multiple tests using this function, because it is not adjusting per comparison (group by Marker) but per pairwise comparison (we only have one pair, Pure vs. Synchronous)
    wilcox_test(
      paired=F,
      as.formula(paste("Score", "~", "Type"))) %>% 
    mutate(p.adj.signif=NULL) %>% 
    adjust_pvalue(method="holm") %>%
    arrange(p)
  
  ypos=aim1CDIScoreData %>%
    group_by(Marker) %>%
    summarize(ypos=max(Score))
  
  #To add the x_positions I need to do something kind of sketchy, first re-sorting the factor (Marker), which is fine, but then sorting the calculated xmin and xmax (which are not in the proper order, because they were calculated internally with the original data, saved in the attribute args)
  #Alternatively, I calculate the xpos here for any Marker, and then use the position to merge by row, this works better in my opinion
  xpos=pwc %>%
    add_x_position(x="Marker", dodge = dodgeWidth) %>%
    mutate(x = sort(x)) %>%
    mutate(xmin = sort(xmin)) %>%
    mutate(xmax = sort(xmax)) %>%
    select(x,xmin,xmax)
  
  sortedMarkers=pwc$Marker
  
  ylimits=c(min(aim1CDIScoreData[Score>0,Score],na.rm = T),max(aim1CDIScoreData$Score,na.rm = T) * topYlimModFig3)
  thisLabels=as.character(pwc %>% mutate(label=paste0(Marker,"\n",n1,"  ",n2)) %>%pull(label))
  names(thisLabels)=as.character(pwc %>% pull(Marker))
  
  thisPwc = pwc %>%
    mutate(Marker=factor(Marker,levels=sortedMarkers))%>%
    drop_na() %>%
    inner_join(ypos,by="Marker") %>%
    mutate(x=row_number()) %>%
    inner_join(xpos,by="x") %>%
    filter(p<=alphaPlot) %>%
    mutate(ypos=ypos*ggsignifPosModFig3)
  
  thisLabels=as.character(pwc %>% mutate(label=paste0(Marker,"\n",n1,"  ",n2)) %>%pull(label))
  names(thisLabels)=as.character(pwc %>% pull(Marker))
  
  fig3CFull=ggplot(aim1CDIScoreData[,c(.SD,.(sortedMarker=factor(Marker,levels = sortedMarkers)))],aes(x=sortedMarker,y=Score,fill=Type))+
    geom_violin(position = position_dodge(width=dodgeWidth),color=NA) +
    geom_quasirandom(alpha=1,size=quasirandomSize,shape=20,dodge.width = dodgeWidth,aes(color=Type),show.legend = F) +
    stat_summary(fun.min=function(z) { quantile(z,0.25) },fun.max=function(z) { quantile(z,0.75) },fun=median,geom="pointrange",color=pointRangeColor,show.legend=F,linewidth=pointRangeLineWidth*pointRangeScaleFig3,size=pointRangeSize*pointRangeScaleFig3,position = position_dodge(width=dodgeWidth)) +
    geom_signif(y_position = log10(thisPwc$ypos),
                annotation=sprintf("%0.1g",thisPwc$p),
                xmin=thisPwc$xmin,
                xmax=thisPwc$xmax,
                tip_length = ggsignifTipLength,
                textsize = ggsignifTextSize,
                vjust = ggsignifTextVjust) + 
    scale_y_log10(name="Within-sample Divergence (CDI)",limits=ylimits) +
    scale_x_discrete(name="",labels=thisLabels) +
    scale_fill_manual(values=crossColors)+
    scale_color_manual(values=crossColors)+
    theme(legend.position="none") + ##WARNING HARDCODED
    labs(title="Within-sample Phenotypic Divergence")
  
  ncol <- 2
  nrow <- 1
  save_plot(fig3CFull,file=paste0(outDir,"fig3CFull.pdf"),base_width = pageWidthInch, base_height = pageWidthInch*nrow/ncol/aspectRatio.plot)
  save_plot(fig3CFull,file=paste0(outDir,"fig3CFull.png"),base_width = pageWidthInch, base_height = pageWidthInch*nrow/ncol/aspectRatio.plot)
  
  pwc = aim1IntRNDScoreData %>%
    group_by(Marker) %>%
    #pairwise_wilcox_test( ##We can't correct for multiple tests using this function, because it is not adjusting per comparison (group by Marker) but per pairwise comparison (we only have one pair, Pure vs. Synchronous)
    wilcox_test(
      paired=F,
      as.formula(paste("intRNDScore", "~", "Type"))) %>% 
    mutate(p.adj.signif=NULL) %>% 
    adjust_pvalue(method="holm") %>%
    arrange(p)
  
  ypos=aim1IntRNDScoreData %>%
    group_by(Marker) %>%
    summarize(ypos=max(intRNDScore))
  
  #To add the x_positions I need to do something kind of sketchy, first re-sorting the factor (Marker), which is fine, but then sorting the calculated xmin and xmax (which are not in the proper order, because they were calculated internally with the original data, saved in the attribute args)
  #Alternatively, I calculate the xpos here for any Marker, and then use the position to merge by row, this works better in my opinion
  xpos=pwc %>%
    add_x_position(x="Marker", dodge = dodgeWidth) %>%
    mutate(x = sort(x)) %>%
    mutate(xmin = sort(xmin)) %>%
    mutate(xmax = sort(xmax)) %>%
    select(x,xmin,xmax)
  
  sortedMarkers=pwc$Marker
  
  ylimits=c(min(aim1IntRNDScoreData[intRNDScore>0,intRNDScore],na.rm = T),max(aim1IntRNDScoreData$intRNDScore,na.rm = T) * topYlimModFig3)
  thisLabels=as.character(pwc %>% mutate(label=paste0(Marker,"\n",n1,"  ",n2)) %>%pull(label))
  names(thisLabels)=as.character(pwc %>% pull(Marker))
  
  thisPwc = pwc %>%
    mutate(Marker=factor(Marker,levels=sortedMarkers))%>%
    drop_na() %>%
    inner_join(ypos,by="Marker") %>%
    mutate(x=row_number()) %>%
    inner_join(xpos,by="x") %>%
    filter(p<=alphaPlot) %>%
    mutate(ypos=ypos*ggsignifPosModFig3)
  
  thisLabels=as.character(pwc %>% mutate(label=paste0(Marker,"\n",n1,"  ",n2)) %>%pull(label))
  names(thisLabels)=as.character(pwc %>% pull(Marker))
  
  fig3AFullRND=ggplot(aim1IntRNDScoreData[,c(.SD,.(sortedMarker=factor(Marker,levels = sortedMarkers)))],aes(x=sortedMarker,y=intRNDScore,fill=Type))+
    geom_violin(position = position_dodge(width=dodgeWidth),color=NA) +
    geom_quasirandom(alpha=1,size=quasirandomSize,shape=20,dodge.width = dodgeWidth,aes(color=Type),show.legend = F) +
    stat_summary(fun.min=function(z) { quantile(z,0.25) },fun.max=function(z) { quantile(z,0.75) },fun=median,geom="pointrange",color=pointRangeColor,show.legend=F,linewidth=pointRangeLineWidth*pointRangeScaleFig3,size=pointRangeSize*pointRangeScaleFig3,position = position_dodge(width=dodgeWidth)) +
    geom_signif(y_position = log10(thisPwc$ypos),
                annotation=sprintf("%0.1g",thisPwc$p),
                xmin=thisPwc$xmin,
                xmax=thisPwc$xmax,
                tip_length = ggsignifTipLength,
                textsize = ggsignifTextSize,
                vjust = ggsignifTextVjust) + 
    scale_y_log10(name="Intensity",limits=ylimits) +
    scale_x_discrete(name="",labels=thisLabels) +
    scale_fill_manual(values=crossColors)+
    scale_color_manual(values=crossColors)+
    theme(legend.position="none") + ##WARNING HARDCODED
    labs(title="Phenotypic intensities")
  
  ncol <- 2
  nrow <- 1
  save_plot(fig3AFullRND,file=paste0(outDir,"fig3AFullRND.pdf"),base_width = pageWidthInch, base_height = pageWidthInch*nrow/ncol/aspectRatio.plot)
  save_plot(fig3AFullRND,file=paste0(outDir,"fig3AFullRND.png"),base_width = pageWidthInch, base_height = pageWidthInch*nrow/ncol/aspectRatio.plot)
}


#Longitudinal

#Figure 4
#########

#Mod params
topYlimModFig4 <- 1.5
topYlimModFig41S <- 1.25 
topYlimModFig40S <- 1

{
  #Fig4R1A
  #SNV
  kruskal.test(SNVBurden~Cohort,data=allAim4Data)
  dunnTest=dunn.test(x=allAim4Data[,SNVBurden],g=allAim4Data[,Cohort],kw = T,method = "hs", altp = T)
  
  Ns=allAim4Data[!is.na(SNVBurden),.N,keyby=Cohort][,`:=`(label=paste(sep="\n",longShortNames[as.character(Cohort)],N))][]
  
  thisLabels=Ns[,label]
  names(thisLabels)=Ns[,Cohort]
  ylimits=c(min(allAim4Data[,SNVBurden],na.rm=T),max(allAim4Data[,SNVBurden],na.rm=T) * topYlimModFig4)
  
  fig4A=ggplot(allAim4Data,aes(x=Cohort,group=Cohort,y=SNVBurden,fill=Cohort))+
    geom_violin(aes(fill=Cohort),color=NA)+
    geom_quasirandom(alpha=1,size=quasirandomSize,shape=20,aes(color=Cohort)) +
    geom_signif(comparisons = list(c("1","2"),c("0","2")),annotation = c(sprintf("%.1g",dunnTest$altP.adjusted[which(dunnTest$comparisons=="1 - 2")]),
                                                                         sprintf("%.1g",dunnTest$altP.adjusted[which(dunnTest$comparisons=="0 - 2")])), 
                step_increase = ggsignifStepIncrease,
                tip_length = ggsignifTipLength,
                textsize = ggsignifTextSize,
                vjust = ggsignifTextVjust)+
    stat_summary(fun.min=function(z) { quantile(z,0.25) },fun.max=function(z) { quantile(z,0.75) },fun=median,geom="pointrange",color=pointRangeColor,show.legend=F,linewidth=pointRangeLineWidth,size=pointRangeSize,position = position_dodge(width=dodgeWidth)) +
    #stat_summary(geom = "text",fun.data=myMedianAnnotation,fun.args=list(yjust=0,label="ñ: ",ndigits=2),size=5,hjust=-0.3) +
    annotate("text",label=sprintf("Omnibus p = %.1g",kruskal.test(SNVBurden~Cohort,data=allAim4Data)$p.value), y= Inf, x =Inf,vjust=1, hjust=1, size = ommnibusTextSize) +
    scale_y_continuous(name="# SNVs", limits = ylimits) +
    scale_x_discrete(name="",labels=thisLabels) +
    scale_fill_manual(values = longColors)+
    scale_color_manual(values = longColors)+
    theme(legend.position="none") +
    labs(title="SNV Burden")
  
  fig4A
  
  #Fig4R1B
  #CNA
  mixedEffectsAnovaSQRTalteredP=aov(data=aim4CNADataBySample,sqrt(AlteredP)~Cohort+Error(patient))
  summary(mixedEffectsAnovaSQRTalteredP)[[1]]
  eMeansResults=as.data.frame(summary(emmeans(mixedEffectsAnovaSQRTalteredP,pairwise~Cohort))$contrasts)
  
  Ns=aim4CNADataBySample[,.N,keyby=Cohort][,`:=`(label=paste(sep="\n",longShortNames[as.character(Cohort)],N))][]
  thisLabels=Ns[,label]
  names(thisLabels)=Ns[,Cohort]
  ylimits=c(min(aim4CNADataBySample[,AlteredP*100],na.rm=T),max(aim4CNADataBySample[,AlteredP*100],na.rm=T) * topYlimModFig4)
  
  fig4B=ggplot(aim4CNADataBySample,aes(x=Cohort,group=Cohort,y=AlteredP*100,fill=Cohort))+
    geom_violin(aes(fill=Cohort),color=NA)+
    geom_quasirandom(alpha=1,size=quasirandomSize,shape=20,aes(color=Cohort)) +
    #geom_signif(comparisons = list(c("0","1"),c("1","2"),c("0","2")),annotation = c(round(eMeansResults[eMeansResults[,1]=="Cohort0 - Cohort1","p.value"],digits=3),
    #                                                                                round(eMeansResults[eMeansResults[,1]=="Cohort1 - Cohort2","p.value"],digits=3),
    #                                                                                round(eMeansResults[eMeansResults[,1]=="Cohort0 - Cohort2","p.value"],digits=3)), step_increase = 0.05,tip_length = 0.025,)+
    geom_signif(comparisons = list(c("1","2"),
                                   c("0","2")),
                annotation = c(sprintf("%.1g",eMeansResults[eMeansResults[,1]=="Cohort1 - Cohort2","p.value"]),
                               sprintf("%.1g",eMeansResults[eMeansResults[,1]=="Cohort0 - Cohort2","p.value"])), 
                step_increase = ggsignifStepIncrease,
                tip_length = ggsignifTipLength,
                textsize = ggsignifTextSize,
                vjust = ggsignifTextVjust)+
    stat_summary(fun.min=function(z) { quantile(z,0.25) },fun.max=function(z) { quantile(z,0.75) },fun=median,geom="pointrange",color=pointRangeColor,linewidth=pointRangeLineWidth,size=pointRangeSize,position = position_dodge(width=dodgeWidth)) +
    #stat_summary(geom = "text",fun.data=myMedianAnnotation,fun.args=list(yjust=0,label="ñ: ",ndigits=2),size=5,hjust=-0.3) +
    #annotate("text",label=paste0("Mixed-effects Anova SQRT-transformed p= ",round(summary(mixedEffectsAnovaSQRTalteredP)[[1]][[1]][["Pr(>F)"]][1],digits = 3)), y= Inf, x =Inf,vjust=1, hjust=1) +
    annotate("text",label=sprintf("Omnibus p = %.1g", summary(mixedEffectsAnovaSQRTalteredP)[[1]][[1]][["Pr(>F)"]][1]), y= Inf, x =Inf,vjust=1, hjust=1, size = ommnibusTextSize) +
    scale_y_continuous(name="% altered genome", limits = ylimits) +
    scale_x_discrete(name="",labels=thisLabels) +
    scale_fill_manual(values = longColors)+
    scale_color_manual(values = longColors)+
    theme(legend.position="none") +
    labs(title="CNA Burden")
  
  fig4B
  
  #Fig4R2C
  #SNV
  dunnTest=dunn.test(x=aim4SNVData$SNVDivergence,g=aim4SNVData$Cohort,kw = T,method = "hs", altp = T)
  
  Ns=aim4SNVData[!is.na(SNVDivergence),.N,keyby=Cohort][,`:=`(label=paste(sep="\n",longShortNames[as.character(Cohort)],N))][]
  thisLabels=Ns[,label]
  names(thisLabels)=Ns[,Cohort]
  ylimits=c(min(aim4SNVData[,SNVDivergence],na.rm=T),max(aim4SNVData[,SNVDivergence],na.rm=T) * topYlimModFig40S)
  
  divergenceSNVPlot=ggplot(aim4SNVData,aes(x=Cohort,y=SNVDivergence,fill=Cohort)) +
    geom_violin(aes(fill=Cohort),color=NA)+
    geom_quasirandom(alpha=1,size=quasirandomSize,shape=20,aes(color=Cohort),show.legend = F) +
    stat_summary(fun.min=function(z) { quantile(z,0.25) },fun.max=function(z) { quantile(z,0.75) },fun=median,geom="pointrange",color=pointRangeColor,linewidth=pointRangeLineWidth,size=pointRangeSize,position = position_dodge(width=dodgeWidth)) +
    scale_y_continuous(name="Divergence (%)",limits=ylimits) +
    scale_x_discrete(name="",labels=thisLabels) +
    scale_fill_manual(values = longColors)+
    scale_color_manual(values = longColors)+
    theme(legend.position="none") +
    labs(title="SNV Divergence")
  
  #Fig4R2D
  #CNA
  divergenceCNAaov=aov(CNADivergence ~ Cohort,data=aim4CNADataByPatient,na.action=na.omit)
  summary(divergenceCNAaov)
  TukeyHSD(divergenceCNAaov)
  tukeyDivergence=TukeyHSD(divergenceCNAaov)
  
  Ns=aim4CNADataByPatient[,.N,keyby=Cohort][,`:=`(label=paste(sep="\n",longShortNames[as.character(Cohort)],N))][]
  thisLabels=Ns[,label]
  names(thisLabels)=Ns[,Cohort]
  ylimits=c(min(aim4CNADataByPatient[,CNADivergence*100],na.rm=T),max(aim4CNADataByPatient[,CNADivergence*100],na.rm=T) * topYlimModFig41S)
  
  divergenceCNAPlot=ggplot(aim4CNADataByPatient,aes(x=Cohort,group=Cohort,y=CNADivergence*100,fill=Cohort))+
    geom_violin(aes(fill=Cohort),color=NA)+
    geom_quasirandom(alpha=1,size=quasirandomSize,shape=20,aes(color=Cohort),show.legend = F) +
    geom_signif(comparisons = list(c("0","2")),annotation = c(sprintf("%.1g",tukeyDivergence$Cohort["2-0","p adj"])), 
                step_increase = ggsignifStepIncrease,
                tip_length = ggsignifTipLength,
                textsize = ggsignifTextSize,
                vjust = ggsignifTextVjust)+
    stat_summary(fun.min=function(z) { quantile(z,0.25) },fun.max=function(z) { quantile(z,0.75) },fun=median,geom="pointrange",color=pointRangeColor,linewidth=pointRangeLineWidth,size=pointRangeSize,position = position_dodge(width=dodgeWidth)) +
    #stat_summary(geom = "text",fun.data=myMedianAnnotation,fun.args=list(yjust=0,label="ñ: ",ndigits=2),size=5,hjust=-0.3) +
    annotate("text",label=sprintf("Omnibus p = %.1g",summary(divergenceCNAaov)[[1]][["Pr(>F)"]][1]), y= Inf, x =Inf,vjust=1, hjust=1, size = ommnibusTextSize) +
    scale_y_continuous(name="% of non-overlapping\naltered genome",limits = ylimits) +
    scale_x_discrete(name="",labels=thisLabels) +
    scale_fill_manual(values = longColors)+
    scale_color_manual(values = longColors)+
    theme(legend.position="none") +
    labs(title="CNA Divergence")
  
  fig4=plot_grid(fig4A,fig4B,divergenceSNVPlot,divergenceCNAPlot,labels=c("A","B","C","D"),label_size = gridLabelFontSize)
  fig4
  
  ncol <- 2
  nrow <- 2
  save_plot(fig4,file=paste0(outDir,"fig4.pdf"),base_width = columnWidthInch, base_height = columnWidthInch*nrow/ncol/aspectRatio.plot)
  save_plot(fig4,file=paste0(outDir,"fig4.png"),base_width = columnWidthInch, base_height = columnWidthInch*nrow/ncol/aspectRatio.plot)
}

#Figure 5
#########
#Mod params
topYlimModFig5 <- 1.45

{
  ylimits=c(0,1*topYlimModFig5)
  
  dunnTest=dunn.test(x=aim4GLUTData$intScore,g=aim4GLUTData$Cohort,kw = T,method = "hs", altp = T)
  
  Ns=aim4GLUTData[,.(.N,Cohort=first(Cohort)),keyby=OriginalCohort][,`:=`(label=paste(sep="\n",longShortNames[as.character(Cohort)],N))][]
  thisLabels=Ns[,label]
  names(thisLabels)=Ns[,Cohort]
  
  glut1IntensityPlot=ggplot(aim4GLUTData,aes(x=Cohort,group=Cohort,y=intScore,fill=Cohort))+
    geom_violin(aes(fill=Cohort),color=NA)+
    geom_quasirandom(alpha=1,size=quasirandomSize,shape=20,aes(color=Cohort),show.legend = F) +
    geom_signif(comparisons = list(c("1","2"),
                                   c("0","2")                                 ),
                annotation = c(sprintf("%.1g",dunnTest$altP.adjusted[which(dunnTest$comparisons=="1 - 2")]),
                               sprintf("%.1g",dunnTest$altP.adjusted[which(dunnTest$comparisons=="0 - 2")])),
                step_increase = ggsignifStepIncrease,
                tip_length = ggsignifTipLength,
                textsize = ggsignifTextSize,
                vjust = ggsignifTextVjust)+
    stat_summary(fun.min=function(z) { quantile(z,0.25) },fun.max=function(z) { quantile(z,0.75) },fun=median,geom="pointrange",color=pointRangeColor,linewidth=pointRangeLineWidth,size=pointRangeSize,position = position_dodge(width=dodgeWidth)) +
    annotate("text",label=sprintf("Omnibus p = %.1g",kruskal.test(intScore~Cohort,data=aim4GLUTData)$p.value), y= Inf, x =Inf,vjust=1, hjust=1, size = ommnibusTextSize) +
    scale_y_continuous(name="Mean Intensity Score (MIS)", limits = ylimits) +
    scale_x_discrete(name="",labels=thisLabels) +
    scale_fill_manual(values = longColors)+
    scale_color_manual(values = longColors)+
    theme(legend.position="none")+
    labs(title="GLUT1")
  
  glut1IntensityPlot
  
  dunnTest=dunn.test(x=aim4ERPosERData$intScore,g=aim4ERPosERData$Cohort,kw = T,method = "hs", altp = T)
  
  Ns=aim4ERPosERData[,.(.N,Cohort=first(Cohort)),keyby=OriginalCohort][,`:=`(label=paste(sep="\n",longShortNames[as.character(Cohort)],N))][]
  thisLabels=Ns[,label]
  names(thisLabels)=Ns[,Cohort]
  
  erIntensityPlotERPos=ggplot(aim4ERPosERData,aes(x=Cohort,group=Cohort,y=intScore,fill=Cohort))+
    geom_violin(aes(fill=Cohort),color=NA)+
    geom_quasirandom(alpha=1,size=quasirandomSize,shape=20,aes(color=Cohort),show.legend = F) +
    geom_signif(comparisons = list(c("0","1"),c("0","2")),annotation = c(sprintf("%.1g",dunnTest$altP.adjusted[which(dunnTest$comparisons=="0 - 1")]),
                                                                         sprintf("%.1g",dunnTest$altP.adjusted[which(dunnTest$comparisons=="0 - 2")])), 
                step_increase = ggsignifStepIncrease,
                tip_length = ggsignifTipLength,
                textsize = ggsignifTextSize,
                vjust = ggsignifTextVjust)+  
    stat_summary(fun.min=function(z) { quantile(z,0.25) },fun.max=function(z) { quantile(z,0.75) },fun=median,geom="pointrange",color=pointRangeColor,linewidth=pointRangeLineWidth,size=pointRangeSize,position = position_dodge(width=dodgeWidth)) +
    annotate("text",label=sprintf("Omnibus p = %.1g",kruskal.test(intScore~Cohort,data=aim4ERPosERData)$p.value), y= Inf, x =Inf,vjust=1, hjust=1, size = ommnibusTextSize) +
    scale_y_continuous(name="Mean Intensity Score (MIS)", limits = ylimits) +
    scale_x_discrete(name="",labels=thisLabels) + 
    scale_fill_manual(values = longColors)+
    scale_color_manual(values = longColors)+
    theme(legend.position="none")+
    labs(title="ER")
  
  erIntensityPlotERPos
  
  fig5=plot_grid(glut1IntensityPlot,erIntensityPlotERPos,labels=c("A","B"),label_size = gridLabelFontSize)
  fig5
  
  ncol <- 2
  nrow <- 1
  save_plot(fig5,file=paste0(outDir,"fig5.pdf"),base_width = columnWidthInch, base_height = columnWidthInch*nrow/ncol/aspectRatio.plot)
  save_plot(fig5,file=paste0(outDir,"fig5.png"),base_width = columnWidthInch, base_height = columnWidthInch*nrow/ncol/aspectRatio.plot)
}

# ylimits=c(0,1*1.2)#WARNING HARDCODED with padding for geom_signif
# 
# dunnTest=dunn.test(x=aim4ERPosGLUTData$intScore,g=aim4ERPosGLUTData$Cohort,kw = T,method = "hs", altp = T)
# 
# Ns=aim4ERPosGLUTData[,.(.N,Cohort=first(Cohort)),keyby=OriginalCohort][,`:=`(label=paste(sep="\nn = ",OriginalCohort,N))][]
# thisLabels=Ns[,label]
# names(thisLabels)=Ns[,Cohort]
# 
# glut1IntensityPlotERPos=ggplot(aim4ERPosGLUTData,aes(x=Cohort,group=Cohort,y=intScore,fill=Cohort))+
#   geom_violin(aes(fill=Cohort),color=NA)+
#   geom_quasirandom(alpha=1,size=quasirandomSize,shape=20,aes(color=Cohort),show.legend = F) +
#   #geom_signif(comparisons = list(c("0","2")),annotation = c(sprintf("%.1g",dunnTest$altP.adjusted[which(dunnTest$comparisons=="0 - 2")])), step_increase = 0.075,tip_length = 0.025,)+
#   stat_summary(fun.min=function(z) { quantile(z,0.25) },fun.max=function(z) { quantile(z,0.75) },fun=median,geom="pointrange",color=pointRangeColor,linewidth=pointRangeLineWidth,size=pointRangeSize,position = position_dodge(width=dodgeWidth)) +
#   annotate("text",label=sprintf("Omnibus p = %.1g",kruskal.test(intScore~Cohort,data=aim4ERPosGLUTData)$p.value), y= Inf, x =Inf,vjust=1, hjust=1) +
#   scale_y_continuous(name="Mean Intensity Score (MIS)", limits = ylimits) +
#   scale_x_discrete(name="",labels=thisLabels) +
#   scale_fill_manual(values = longColors)+
#   scale_color_manual(values = longColors)+
#   theme(legend.position="none")+
#   labs(title="GLUT1")
# 
# glut1IntensityPlotERPos

#Supplementary Figure: IHC ER
#Mod params
topYlimModFigIHCER <- 1.1
{
  ylimits=c(0,1*topYlimModFigIHCER)
  
  dunnTest=dunn.test(x=aim4ERData$intScore,g=aim4ERData$Cohort,kw = T,method = "hs", altp = T)
  
  Ns=aim4ERData[,.(.N,Cohort=first(Cohort)),keyby=OriginalCohort][,`:=`(label=paste(sep="\n",OriginalCohort,N))][]
  thisLabels=Ns[,label]
  names(thisLabels)=Ns[,Cohort]
  
  erIntensityPlot=ggplot(aim4ERData,aes(x=Cohort,group=Cohort,y=intScore,fill=Cohort))+
    geom_violin(aes(fill=Cohort),color=NA)+
    geom_quasirandom(alpha=1,size=quasirandomSize,shape=20,aes(color=Cohort),show.legend = F) +
    geom_signif(comparisons = list(c("0","2")),annotation = c(sprintf("%.1g",dunnTest$altP.adjusted[which(dunnTest$comparisons=="0 - 2")])), 
                step_increase = ggsignifStepIncrease,
                tip_length = ggsignifTipLength,
                textsize = ggsignifTextSize,
                vjust = ggsignifTextVjust) +  
    stat_summary(fun.min=function(z) { quantile(z,0.25) },fun.max=function(z) { quantile(z,0.75) },fun=median,geom="pointrange",color=pointRangeColor,linewidth=pointRangeLineWidth,size=pointRangeSize) +
    annotate("text",label=sprintf("Omnibus p = %.1g",kruskal.test(intScore~Cohort,data=aim4ERData)$p.value), y= Inf, x =Inf,vjust=1, hjust=1, size=ommnibusTextSize) +
    scale_y_continuous(name="Mean Intensity Score (MIS)", limits = ylimits) +
    scale_x_discrete(name="",labels=thisLabels) +
    scale_fill_manual(values = longColors)+
    scale_color_manual(values = longColors)+
    theme(legend.position="none")+
    labs(title="ER: including ER-")
  
  erIntensityPlot
  
  ncol <- 1
  nrow <- 1
  save_plot(erIntensityPlot,file=paste0(outDir,"erIntensityPlot.pdf"),base_width = columnWidthInch, base_height = columnWidthInch*nrow/ncol/aspectRatio.plot)
  save_plot(erIntensityPlot,file=paste0(outDir,"erIntensityPlot.png"),base_width = columnWidthInch, base_height = columnWidthInch*nrow/ncol/aspectRatio.plot)
}

#Supplementary Figure: Longitudinal Phenotypic Divergence
{
  dunnTest=dunn.test(x=aim4GLUTData$EMDist,g=aim4GLUTData$Cohort,kw = T,method = "hs", altp = T)
  
  Ns=aim4GLUTData[,.(.N,Cohort=first(Cohort)),keyby=OriginalCohort][,`:=`(label=paste(sep="\n",OriginalCohort,N))][]
  thisLabels=Ns[,label]
  names(thisLabels)=Ns[,Cohort]
  
  glut1DivergencePlot=ggplot(aim4GLUTData,aes(x=Cohort,group=Cohort,y=EMDist,fill=Cohort))+
    geom_violin(aes(fill=Cohort),color=NA)+
    geom_quasirandom(alpha=1,size=quasirandomSize,shape=20,aes(color=Cohort),show.legend = F) +
    stat_summary(fun.min=function(z) { quantile(z,0.25) },fun.max=function(z) { quantile(z,0.75) },fun=median,geom="pointrange",color=pointRangeColor,linewidth=pointRangeLineWidth,size=pointRangeSize) +
    scale_y_continuous(name="Between-sample divergence (EMD)") +
    scale_x_discrete(name="",labels=thisLabels) +
    scale_fill_manual(values = longColors)+
    scale_color_manual(values = longColors)+
    theme(legend.position="none")+
    labs(title="GLUT1")
  
  glut1DivergencePlot
  
  dunnTest=dunn.test(x=aim4ERData$cdiScore,g=aim4ERData$Cohort,kw = T,method = "hs", altp = T)
  
  Ns=aim4ERData[,.(.N,Cohort=first(Cohort)),keyby=OriginalCohort][,`:=`(label=paste(sep="\n",OriginalCohort,N))][]
  thisLabels=Ns[,label]
  names(thisLabels)=Ns[,Cohort]
  
  erDivergencePlot=ggplot(aim4ERData,aes(x=Cohort,group=Cohort,y=cdiScore,fill=Cohort))+
    geom_violin(aes(fill=Cohort),color=NA)+
    geom_quasirandom(alpha=1,size=quasirandomSize,shape=20,aes(color=Cohort),show.legend = F) +
    stat_summary(fun.min=function(z) { quantile(z,0.25) },fun.max=function(z) { quantile(z,0.75) },fun=median,geom="pointrange",color=pointRangeColor,linewidth=pointRangeLineWidth,size=pointRangeSize,position = position_dodge(width=dodgeWidth)) +
    scale_y_continuous(name="Within-sample divergence (CDI)") +
    scale_x_discrete(name="",labels=thisLabels) +
    scale_fill_manual(values = longColors)+
    scale_color_manual(values = longColors)+
    theme(legend.position="none")+
    labs(title="ER")
  
  erDivergencePlot
  
  figLongitudinalPhenotypicDivergence=plot_grid(glut1DivergencePlot,erDivergencePlot,labels=c("A","B"),label_size = gridLabelFontSize)
  figLongitudinalPhenotypicDivergence
  
  ncol <- 2
  nrow <- 1
  save_plot(figLongitudinalPhenotypicDivergence,file=paste0(outDir,"longitudinalPhenotypicDivergencePlot.pdf"),base_width = pageWidthInch, base_height = pageWidthInch*nrow/ncol/aspectRatio.plot)
  save_plot(figLongitudinalPhenotypicDivergence,file=paste0(outDir,"longitudinalPhenotypicDivergencePlot.png"),base_width = pageWidthInch, base_height = pageWidthInch*nrow/ncol/aspectRatio.plot)
}

#Data prep for next steps
CNAvariables=c("CNADivergence","MeanAlteredP")
IHCvariables=c('MIS_GLUT1', 'EMD_GLUT1', 'MIS_ER', 'CDI_ER')
SNVvariables=c('SNVBurden', 'SNVDivergence')
humanReadableVariables=c(CNAvariables,IHCvariables,SNVvariables)
names(humanReadableVariables)=c("CNA divergence", "CNA burden", "GLUT1 intensity", "GLUT1 divergence", "ER intensity", "ER divergence","SNV burden","SNV divergence")

#Univariate cox tables
{
  coxTableRecurrence=rbindlist(lapply(humanReadableVariables,getCox,theData=allAim4DataWClinical,timeToEventVar = "TimeToEvent_mo",eventVar = "EventRecurrence"))
  coxTableNoninvasiveRecurrence=rbindlist(lapply(humanReadableVariables,getCox,theData=allAim4DataWClinical[EventProgression==F,],timeToEventVar = "TimeToEvent_mo",eventVar = "EventRecurrence"))
  #In principle, we are not planning to use logistic models for simplicity, but I have here the tools to bring them back if needed
  #logisticTableNoninvasiveRecurrence <- rbindlist(lapply(humanReadableVariables,getLogisticModel,theData=allAim4DataWClinical[EventProgression==F,][fullDataPatients,],eventVar = "EventRecurrence"))
  
  coxTableProgressionFromNonprogressors=rbindlist(lapply(humanReadableVariables,getCox,theData=allAim4DataWClinical,timeToEventVar = "TimeToEvent_mo",eventVar = "EventProgression"))
  coxTableProgression=rbindlist(lapply(humanReadableVariables,getCox,theData=allAim4DataWClinical[EventProgression==EventRecurrence,],timeToEventVar = "TimeToEvent_mo",eventVar = "EventProgression"))
  #logisticTableProgression <- rbindlist(lapply(humanReadableVariables,getLogisticModel,theData=allAim4DataWClinical[EventProgression==EventRecurrence,][fullDataPatients,],eventVar = "EventProgression"))
  
  coxTableProgressionFromNonprogressors[,`:=`(adj.p.val=(p.adjust(p.val,method = "holm")))]
  coxTableProgression[,`:=`(adj.p.val=(p.adjust(p.val,method = "holm")))]
  #logisticTableProgression[,`:=`(adj.p.val=(p.adjust(p.val,method = "holm")))]
  
  coxTableRecurrence[,`:=`(adj.p.val=(p.adjust(p.val,method = "holm")))]
  coxTableNoninvasiveRecurrence[,`:=`(adj.p.val=(p.adjust(p.val,method = "holm")))]
  #logisticTableNoninvasiveRecurrence[,`:=`(adj.p.val=(p.adjust(p.val,method = "holm")))]
  
  
  
  finalCoxTableProgression=coxTableProgression[,.(Variable=factor(Variable,levels=humanReadableVariables,labels=names(humanReadableVariables)),
                                                  p.val,
                                                  adj.p.val,
                                                  C, 
                                                  "HR (95% CI)"=sprintf("%0.2g (%0.2g - %0.2g)",HR,HRlowCI,HRupperCI))]
  
  finalCoxTableProgressionFromNonprogressors=coxTableProgressionFromNonprogressors[,.(Variable=factor(Variable,levels=humanReadableVariables,labels=names(humanReadableVariables)),
                                                       p.val,
                                                       adj.p.val,
                                                       C, 
                                                       "HR (95% CI)"=sprintf("%0.2g (%0.2g - %0.2g)",HR,HRlowCI,HRupperCI))]
  
  finalCoxTableRecurrence=coxTableRecurrence[,.(Variable=factor(Variable,levels=humanReadableVariables,labels=names(humanReadableVariables)),
                                                p.val,
                                                adj.p.val,
                                                C, 
                                                "HR (95% CI)"=sprintf("%0.2g (%0.2g - %0.2g)",HR,HRlowCI,HRupperCI))]
  
  finalCoxTableNoninvasiveRecurrence=coxTableNoninvasiveRecurrence[,.(Variable=factor(Variable,levels=humanReadableVariables,labels=names(humanReadableVariables)),
                                                p.val,
                                                adj.p.val,
                                                C, 
                                                "HR (95% CI)"=sprintf("%0.2g (%0.2g - %0.2g)",HR,HRlowCI,HRupperCI))]
  
  textTable(finalCoxTableRecurrence,outfile = paste0(outDir,"univariateCoxRecAll.tex"),alpha=1, 
             tableNumber = 7,
             supplementary = T,
             caption="Univariate Proportional Hazard Regressions of Time to Recurrence",
             footer="Adj.p.val: adjusted p.value using the Holm correction. C: concordance. HR (95\\%CI): hazard ratio and 95\\% confidence interval for the standard score (z-scores) of the variable of interest.")
  textTable(finalCoxTableProgressionFromNonprogressors,outfile = paste0(outDir,"univariateCoxProgFromNonprogAll.tex"),alpha = 1,
            tableNumber = 8,
            supplementary = T, 
            caption="Univariate Proportional Hazard Regressions of Time to Progression from Nonprogresors",
            footer="Adj.p.val: adjusted p.value using the Holm correction. C: concordance. HR (95\\%CI): hazard ratio and 95\\% confidence interval for the standard score (z-scores) of the variable of interest.")
  textTable(finalCoxTableNoninvasiveRecurrence,outfile = paste0(outDir,"univariateCoxNoninvRecAll.tex"),alpha=1,
            tableNumber = 9,
            supplementary = T, 
            caption="Univariate Proportional Hazard Regressions of Time to Noninvasive Recurrence",
            footer="Adj.p.val: adjusted p.value using the Holm correction. C: concordance. HR (95\\%CI): hazard ratio and 95\\% confidence interval for the standard score (z-scores) of the variable of interest.")
  textTable(finalCoxTableProgression,outfile = paste0(outDir,"univariateCoxProgAll.tex"),alpha = 1, 
            tableNumber = 10,
            supplementary = T, 
            caption="Univariate Proportional Hazard Regressions of Time to Progression",
            footer="Adj.p.val: adjusted p.value using the Holm correction. C: concordance. HR (95\\%CI): hazard ratio and 95\\% confidence interval for the standard score (z-scores) of the variable of interest.")

}

#Figure 6
#########
#KM plots for SNV burden
#Mod params
relHeightFig6 <- 0.2
aspectRatioFig6Mod <- 1.75

{
  # Previous version in which we take the threshold from a binomial regression. I think it is better to do everything with Cox to make it simpler (results are equivalent)
  # binomialSNVBurdenRecurrence=glm(Cohort!=0 ~ SNVBurden, data=allAim4DataWClinical,family=binomial)
  # binomialSNVBurdenProgression=glm(Cohort==2 ~ SNVBurden, data=allAim4DataWClinical,family=binomial)
  # 
  # progressionThreshold=getThreshold(binomialSNVBurdenProgression,"J")
  # recurrenceThreshold=getThreshold(binomialSNVBurdenRecurrence,"J")
  # 
  # allAim4DataWClinical[,`:=` (groupProgression=as.factor(ifelse(SNVBurden>=progressionThreshold,"high","low")))]
  # allAim4DataWClinical[,`:=` (groupRecurrence=as.factor(ifelse(SNVBurden>=recurrenceThreshold,"high","low")))]
  # 
  # recurrenceKMModel=surv_fit(Surv(TimeToEvent_mo,Cohort!=0) ~ groupRecurrence,data = allAim4DataWClinical[!is.na(groupRecurrence),])
  # abs(diff(surv_median(recurrenceKMModel)$median))
  # surv_pvalue(recurrenceKMModel)
  # 
  # progressionKMModel=surv_fit(Surv(TimeToEvent_mo,Cohort==2) ~ groupProgression,data = allAim4DataWClinical[!is.na(groupProgression),])
  # abs(diff(surv_median(progressionKMModel)$median))
  # surv_pvalue(progressionKMModel)
  #
  # recurrenceKMPlot=getKMPlot(recurrenceKMModel,threshold = recurrenceThreshold,labelConstant = "SNVs",ylab = "Recurrence-free Survival",colors = c('#8F770099','#0073C299'),title="Recurrence")
  # progressionKMPlot=getKMPlot(progressionKMModel,threshold = progressionThreshold,labelConstant = "SNVs",ylab = "Progression-free Survival",colors = c('#A7303099','#0073C299'),title="Progression")
  
  recurrenceKMPlot_Data <- allAim4DataWClinical[!is.na(SNVBurden),]
  recurrenceCoxModel <- coxph(Surv(TimeToEvent_mo,EventRecurrence) ~ SNVBurden, 
                              data = recurrenceKMPlot_Data)
  
  recurrenceKMPlot <- getKMPlotUnivariateCox(theseData = recurrenceKMPlot_Data,
                                             thesePatients = recurrenceKMPlot_Data$patient,
                                             coxModel = recurrenceCoxModel,
                                             timeToEventVar = "TimeToEvent_mo",
                                             eventVar = "EventRecurrence",
                                             theseColors = c('#8F770099','#0073C299'),
                                             title = "Recurrence",
                                             xname = "SNVs")
  
  recurrenceOnlyKMPlot_Data <- allAim4DataWClinical[!is.na(SNVBurden),][EventProgression==F,]
  recurrenceOnlyCoxModel <- coxph(Surv(TimeToEvent_mo,EventRecurrence) ~ SNVBurden, 
                                  data = recurrenceOnlyKMPlot_Data)
  
  recurrenceOnlyKMPlot <- getKMPlotUnivariateCox(theseData = recurrenceOnlyKMPlot_Data,
                                                 thesePatients = recurrenceOnlyKMPlot_Data$patient,
                                                 coxModel = recurrenceOnlyCoxModel,
                                                 timeToEventVar = "TimeToEvent_mo",
                                                 eventVar = "EventRecurrence",
                                                 theseColors = c('#8F770099','#0073C299'),
                                                 title = "Non-invasive Recurrence",
                                                 xname = "SNVs")
  
  progressionKMPlot_Data <- allAim4DataWClinical[!is.na(SNVBurden),]
  progressionCoxModel <- coxph(Surv(TimeToEvent_mo,EventProgression) ~ SNVBurden, 
                               data = progressionKMPlot_Data)
  
  progressionKMPlot <- getKMPlotUnivariateCox(theseData = progressionKMPlot_Data,
                                              thesePatients = progressionKMPlot_Data$patient,
                                              coxModel = progressionCoxModel,
                                              timeToEventVar = "TimeToEvent_mo",
                                              eventVar = "EventProgression",
                                              theseColors = c('#A7303099','#0073C299'),
                                              title = "Progression from non-progressors",
                                              xname = "SNVs")
  
  progressionOnlyKMPlot_Data <- allAim4DataWClinical[!is.na(SNVBurden),][EventRecurrence==EventProgression,]
  progressionOnlyCoxModel <- coxph(Surv(TimeToEvent_mo,EventProgression) ~ SNVBurden, 
                                   data = progressionOnlyKMPlot_Data)
  
  progressionOnlyKMPlot <- getKMPlotUnivariateCox(theseData = progressionOnlyKMPlot_Data,
                                                  thesePatients = progressionOnlyKMPlot_Data$patient,
                                                  coxModel = progressionOnlyCoxModel,
                                                  timeToEventVar = "TimeToEvent_mo",
                                                  eventVar = "EventProgression",
                                                  theseColors = c('#A7303099','#0073C299'),
                                                  title = "Progression",
                                                  xname = "SNVs")
  
  ncol <- 1
  nrow <- 2
  
  fig6S=plot_grid(recurrenceKMPlot$plot,recurrenceKMPlot$table,progressionKMPlot$plot,progressionKMPlot$table,labels=c("A","","B",""),rel_heights = c(1,relHeightFig6,1,relHeightFig6),align = 'v',label_size = gridLabelFontSize,nrow = 4) #cowplot to the rescue!
  save_plot(fig6S,file=paste0(outDir,"fig6S.pdf"),base_width = columnWidthInch, base_height = columnWidthInch*nrow/ncol/aspectRatio.plot/aspectRatioFig6Mod,device=cairo_pdf)
  save_plot(fig6S,file=paste0(outDir,"fig6S.png"),base_width = columnWidthInch, base_height = columnWidthInch*nrow/ncol/aspectRatio.plot/aspectRatioFig6Mod)
  
  fig6=plot_grid(recurrenceOnlyKMPlot$plot,recurrenceOnlyKMPlot$table,progressionOnlyKMPlot$plot,progressionOnlyKMPlot$table,labels=c("A","","B",""),rel_heights = c(1,relHeightFig6,1,relHeightFig6),align = 'v',label_size = gridLabelFontSize,nrow = 4) #cowplot to the rescue!
  save_plot(fig6,file=paste0(outDir,"fig6.pdf"),base_width = columnWidthInch, base_height = columnWidthInch*nrow/ncol/aspectRatio.plot/aspectRatioFig6Mod,device=cairo_pdf)
  save_plot(fig6,file=paste0(outDir,"fig6.png"),base_width = columnWidthInch, base_height = columnWidthInch*nrow/ncol/aspectRatio.plot/aspectRatioFig6Mod)
}

#Figure 7
#########
## Full clinical LASSO
relHeightFig7 <- 0.2
aspectRatioFig7Mod <- 1.75


#Data prep
fullDataPatients=allAim4Data[,.(patient,fullData=all(!is.na(.SD))),by=.I][fullData==T,patient]

imputedAllAim4DataOriginal=makeX(allAim4DataWClinical[,.SD,.SDcols=c(humanReadableVariables,clinicalVariablesForNow)],na.impute = T)
scaledImputedAllAim4DataOriginal=scale(imputedAllAim4DataOriginal)
scaledImputedAllAim4Data=data.table(scaledImputedAllAim4DataOriginal,
                                    allAim4DataWClinical[,.(EventRecurrence,EventProgression,TimeToEvent_mo,patient)])
imputedAllAim4Data=data.table(imputedAllAim4DataOriginal,
                                    allAim4DataWClinical[,.(EventRecurrence,EventProgression,TimeToEvent_mo,patient)])

NAAllAim4DataOriginal=makeX(allAim4DataWClinical[,.SD,.SDcols=c(humanReadableVariables,clinicalVariablesForNow)],na.impute = F)
scaledNAAllAim4DataOriginal=scale(NAAllAim4DataOriginal)
scaledNAAllAim4Data=data.table(scaledNAAllAim4DataOriginal,
                               allAim4DataWClinical[,.(EventRecurrence,EventProgression,TimeToEvent_mo,patient)])
NAAllAim4Data=data.table(NAAllAim4DataOriginal,
                               allAim4DataWClinical[,.(EventRecurrence,EventProgression,TimeToEvent_mo,patient)])

xColumnNames=colnames(scaledImputedAllAim4DataOriginal)

#Recurrence

#ATTENTION: along this section, the recurrenceOnly tag means the proportional hazards version of cohort 0 vs cohort1, while recurrence means cohort 0 vs cohort 1 + cohort 2

#We always use imputation in the covariates for covariate selection using LASSO
recurrenceFullDataLassoDevianceFreqs=getLassoFreq(seed=theSeed,
                                                  n=nLasso,
                                                  x=scaledImputedAllAim4Data[fullDataPatients,.SD,.SDcols=xColumnNames],
                                                  y=Surv(scaledImputedAllAim4Data[fullDataPatients,TimeToEvent_mo],scaledImputedAllAim4Data[fullDataPatients,EventRecurrence]),
                                                  family = "cox", 
                                                  standarize = F, 
                                                  type.measure = "deviance",
                                                  lambda="lambda.min")

recurrenceOnlyFullDataLassoDevianceFreqs=getLassoFreq(seed=theSeed,
                                                      n=nLasso,
                                                      x=scaledImputedAllAim4Data[fullDataPatients][EventProgression==F,.SD,.SDcols=xColumnNames],
                                                      y=Surv(scaledImputedAllAim4Data[fullDataPatients][EventProgression==F,TimeToEvent_mo],scaledImputedAllAim4Data[fullDataPatients][EventProgression==F,EventRecurrence]),
                                                      family = "cox", 
                                                      standarize = F, 
                                                      type.measure = "deviance",
                                                      lambda="lambda.min")

#Progression
progressionFullDataLassoDevianceFreqs=getLassoFreq(seed=theSeed,
                                                   n=nLasso,
                                                   x=scaledImputedAllAim4Data[fullDataPatients,.SD,.SDcols=xColumnNames],
                                                   y=Surv(scaledImputedAllAim4Data[fullDataPatients,TimeToEvent_mo],scaledImputedAllAim4Data[fullDataPatients,EventProgression]),
                                                   family = "cox", 
                                                   standarize = F, 
                                                   type.measure = "deviance",
                                                   lambda="lambda.min")

#ProgressionOnly 
progressionOnlyFullDataLassoDevianceFreqs=getLassoFreq(seed=theSeed,
                                                       n=nLasso,
                                                       x=scaledImputedAllAim4Data[fullDataPatients][EventRecurrence==EventProgression,.SD,.SDcols=xColumnNames],
                                                       y=Surv(scaledImputedAllAim4Data[fullDataPatients][EventRecurrence==EventProgression,TimeToEvent_mo],scaledImputedAllAim4Data[fullDataPatients][EventRecurrence==EventProgression,EventProgression]),
                                                       family = "cox", 
                                                       standarize = F, 
                                                       type.measure = "deviance",
                                                       lambda="lambda.min")

#ProgressionFromRecurrence (probably to remove)
progressionFromRecurrenceFullDataLassoDevianceFreqs=getLassoFreq(seed=theSeed,
                                                                 n=nLasso,
                                                                 x=scaledImputedAllAim4Data[fullDataPatients][EventRecurrence==T,.SD,.SDcols=xColumnNames],
                                                                 y=Surv(scaledImputedAllAim4Data[fullDataPatients][EventRecurrence==T,TimeToEvent_mo],scaledImputedAllAim4Data[fullDataPatients][EventRecurrence==T,EventProgression]),
                                                                 family = "cox",
                                                                 standarize = F,
                                                                 type.measure = "deviance",
                                                                 lambda="lambda.min")

#Recurrence subplots
{
  #Recurrence 
  figS7AB <- getForestAndKMCoxFromData(imputedData = scaledImputedAllAim4Data,
                                       originalData = allAim4Data, #We filter patients based on the non-imputed measured variables
                                       freqs = recurrenceFullDataLassoDevianceFreqs,
                                       minFreq = minLassoFreq,
                                       eventVar = "EventRecurrence",
                                       colors = c('#8F770099','#0073C299'),
                                       kmTitle = "Recurrence",
                                       forestTitle = "    Recurrence Proportional Hazards Model",
                                       termDict = humanReadableHotEncodedVariables,
                                       rel_height = relHeightFig7)
  
  #Noninvasive recurrence
  fig7AB <- getForestAndKMCoxFromData(imputedData = scaledImputedAllAim4Data,
                                      originalData = allAim4Data[EventProgression==F,],
                                      freqs = recurrenceOnlyFullDataLassoDevianceFreqs,
                                      minFreq = 0.65,
                                      eventVar = "EventRecurrence",
                                      colors = c('#8F770099','#0073C299'),
                                      kmTitle = "Non-invasive Recurrence",
                                      forestTitle = "    Non-invasive Recurrence Proportional Hazards Model",
                                      termDict = humanReadableHotEncodedVariables,
                                      rel_height = relHeightFig7,
                                      originalSDs = attr(scaledImputedAllAim4DataOriginal,"scaled:scale"))
  
  ## Estimated margin manual test
  ## recurrenceFig7Patients <- getFullDataPatients(allAim4Data,recurrenceFullDataLassoDevianceFreqs,minLassoFreq)
  ## summary(coxph(Surv(TimeToEvent_mo,EventRecurrence) ~ EstimatedMargin, data=scaledImputedAllAim4Data[recurrenceFig7Patients,]))
  ## summary(coxph(Surv(TimeToEvent_mo,EventRecurrence) ~ EstimatedMargin, data=scaledNAAllAim4Data[recurrenceFig7Patients,]))
  ## recurrenceOnlyFig7Patients <- getFullDataPatients(allAim4Data[EventProgression==F,],recurrenceOnlyFullDataLassoDevianceFreqs,minLassoFreq)
  ## summary(coxph(Surv(TimeToEvent_mo,EventRecurrence) ~ EstimatedMargin, data=scaledImputedAllAim4Data[recurrenceOnlyFig7Patients,]))
  ## summary(coxph(Surv(TimeToEvent_mo,EventRecurrence) ~ EstimatedMargin, data=scaledNAAllAim4Data[recurrenceOnlyFig7Patients,]))
  
  #Significant-only Alternatives
  #Recurrence 
  figS7AB_SigOnly <- getForestAndKMCoxFromData(imputedData = scaledImputedAllAim4Data,
                                               originalData = allAim4Data,
                                               freqs = recurrenceFullDataLassoDevianceFreqs,
                                               minFreq = minLassoFreq,
                                               eventVar = "EventRecurrence",
                                               colors = c('#8F770099','#0073C299'),
                                               kmTitle = "Recurrence",
                                               forestTitle = "    Recurrence Proportional Hazards Model",
                                               termDict = humanReadableHotEncodedVariables,
                                               rel_height = relHeightFig7,
                                               sigOnly = T)
  
  #Noninvasive recurrence
  fig7AB_SigOnly <- getForestAndKMCoxFromData(imputedData = scaledImputedAllAim4Data,
                                              originalData = allAim4Data[EventProgression==F,],
                                              freqs = recurrenceOnlyFullDataLassoDevianceFreqs,
                                              minFreq = minLassoFreq,
                                              eventVar = "EventRecurrence",
                                              colors = c('#8F770099','#0073C299'),
                                              kmTitle = "Non-invasive Recurrence",
                                              forestTitle = "    Non-invasive Recurrence Proportional Hazards Model",
                                              termDict = humanReadableHotEncodedVariables,
                                              rel_height = relHeightFig7,
                                              sigOnly = T)
}

#Progression subplots
{
  #Progression from Nonprogressors
  figS7CD <- getForestAndKMCoxFromData(imputedData = scaledImputedAllAim4Data,
                                      originalData = allAim4Data,
                                      freqs = progressionFullDataLassoDevianceFreqs,
                                      minFreq = minLassoFreq,
                                      eventVar = "EventProgression",
                                      colors = c('#A7303099','#0073C299'),
                                      kmTitle = "Progression from Nonprogressors",
                                      forestTitle = "    Progression from Nonprogressors Proportional Hazards Model",
                                      termDict = humanReadableHotEncodedVariables,
                                      rel_height = relHeightFig7)
  
  #ProgressionOnly
  fig7CD <- getForestAndKMCoxFromData(imputedData = scaledImputedAllAim4Data,
                                      originalData = allAim4Data[EventProgression==EventRecurrence,],
                                      freqs = progressionOnlyFullDataLassoDevianceFreqs,
                                      minFreq = minLassoFreq,
                                      eventVar = "EventProgression",
                                      colors = c('#A7303099','#0073C299'),
                                      kmTitle = "Progression",
                                      forestTitle = "    Progression Proportional Hazards Model",
                                      termDict = humanReadableHotEncodedVariables,
                                      rel_height = relHeightFig7,
                                      originalSDs = attr(scaledImputedAllAim4DataOriginal,"scaled:scale"))
  
  ## Estimated margin manual test
  ## progressionFig7Patients <- getFullDataPatients(allAim4Data,progressionFullDataLassoDevianceFreqs,minLassoFreq)
  ## summary(coxph(Surv(TimeToEvent_mo,EventProgression) ~ EstimatedMargin, data=scaledImputedAllAim4Data[progressionFig7Patients,]))
  ## summary(coxph(Surv(TimeToEvent_mo,EventProgression) ~ EstimatedMargin, data=scaledNAAllAim4Data[progressionFig7Patients,]))
  ## progressionOnlyFig7Patients <- getFullDataPatients(allAim4Data[EventProgression==EventRecurrence,],progressionOnlyFullDataLassoDevianceFreqs,minLassoFreq)
  ## summary(coxph(Surv(TimeToEvent_mo,EventProgression) ~ EstimatedMargin, data=scaledImputedAllAim4Data[progressionOnlyFig7Patients,]))
  ## summary(coxph(Surv(TimeToEvent_mo,EventProgression) ~ EstimatedMargin, data=scaledNAAllAim4Data[progressionOnlyFig7Patients,]))
  
  #ProgressionFromRecurrence (probably to remove)
  figSFromRecurrenceCD <- getForestAndKMCoxFromData(imputedData = scaledImputedAllAim4Data,
                                                    originalData = allAim4Data[EventRecurrence==T,],
                                                    freqs = progressionFromRecurrenceFullDataLassoDevianceFreqs,
                                                    minFreq = minLassoFreq,
                                                    eventVar = "EventProgression",
                                                    colors = c('#A7303099','#0073C299'),
                                                    kmTitle = "Progression from Recurrence",
                                                    forestTitle = "    Progression from Recurrence Proportional Hazards Model",
                                                    termDict = humanReadableHotEncodedVariables,
                                                    rel_height = relHeightFig7)
  #Significant-only Alternatives
  #Progression from Nonprogressors
  fig7CD_SigOnly <- getForestAndKMCoxFromData(imputedData = scaledImputedAllAim4Data,
                                              originalData = allAim4Data,
                                              freqs = progressionFullDataLassoDevianceFreqs,
                                              minFreq = minLassoFreq,
                                              eventVar = "EventProgression",
                                              colors = c('#A7303099','#0073C299'),
                                              kmTitle = "Progression from Nonprogressors",
                                              forestTitle = "    Progression from Nonprogressors Proportional Hazards Model",
                                              termDict = humanReadableHotEncodedVariables,
                                              rel_height = relHeightFig7,
                                              sigOnly = T)
  
  #ProgressionOnly
  fig7CD_SigOnly <- getForestAndKMCoxFromData(imputedData = scaledImputedAllAim4Data,
                                              originalData = allAim4Data[EventProgression==EventRecurrence,],
                                              freqs = progressionOnlyFullDataLassoDevianceFreqs,
                                              minFreq = minLassoFreq,
                                              eventVar = "EventProgression",
                                              colors = c('#A7303099','#0073C299'),
                                              kmTitle = "Progression",
                                              forestTitle = "    Progression Proportional Hazards Model",
                                              termDict = humanReadableHotEncodedVariables,
                                              rel_height = relHeightFig7,
                                              sigOnly = T,
                                              originalSDs = attr(scaledImputedAllAim4DataOriginal,"scaled:scale"))
  
  #ProgressionFromRecurrence (probably to remove)
  figSFromRecurrenceCD_SigOnly <- getForestAndKMCoxFromData(imputedData = scaledImputedAllAim4Data,
                                                            originalData = allAim4Data[EventRecurrence==T,],
                                                            freqs = progressionFromRecurrenceFullDataLassoDevianceFreqs,
                                                            minFreq = minLassoFreq,
                                                            eventVar = "EventProgression",
                                                            colors = c('#A7303099','#0073C299'),
                                                            kmTitle = "Progression from Recurrence",
                                                            forestTitle = "    Progression from Recurrence Proportional Hazards Model",
                                                            termDict = humanReadableHotEncodedVariables,
                                                            rel_height = relHeightFig7,
                                                            sigOnly = T)
  
}

#Plot composition and saving
{
  ncol <- 2
  nrow <- 2
  figS7 <- plot_grid(figS7AB[[1]],figS7AB[[2]],figS7CD[[1]],figS7CD[[2]],labels=c("A","B","C","D"),label_size = gridLabelFontSize)
  
  save_plot(figS7,file=paste0(outDir,"figS7.pdf"),base_width = pageWidthInch, base_height = pageWidthInch*nrow/ncol/aspectRatio.plot/aspectRatioFig7Mod,device=cairo_pdf)
  save_plot(figS7,file=paste0(outDir,"figS7.png"),base_width = pageWidthInch, base_height = pageWidthInch*nrow/ncol/aspectRatio.plot/aspectRatioFig7Mod)
  
  fig7 <- plot_grid(fig7AB[[1]],fig7AB[[2]],fig7CD[[1]],fig7CD[[2]],labels=c("A","B","C","D"),label_size = gridLabelFontSize)
  
  save_plot(fig7,file=paste0(outDir,"fig7.pdf"),base_width = pageWidthInch, base_height = pageWidthInch*nrow/ncol/aspectRatio.plot/aspectRatioFig7Mod,device=cairo_pdf)
  save_plot(fig7,file=paste0(outDir,"fig7.png"),base_width = pageWidthInch, base_height = pageWidthInch*nrow/ncol/aspectRatio.plot/aspectRatioFig7Mod)
  
  ncol <- 2
  nrow <- 1
  figSFromRecurrence <- plot_grid(figSFromRecurrenceCD[[1]],figSFromRecurrenceCD[[2]],labels=c("A","B"),label_size = gridLabelFontSize)
  save_plot(figSFromRecurrence,file=paste0(outDir,"figSFromRecurrence.pdf"),base_width = pageWidthInch, base_height = pageWidthInch*nrow/ncol/aspectRatio.plot/aspectRatioFig7Mod,device=cairo_pdf)
  save_plot(figSFromRecurrence,file=paste0(outDir,"figSFromRecurrence.png"),base_width = pageWidthInch, base_height = pageWidthInch*nrow/ncol/aspectRatio.plot/aspectRatioFig7Mod)
}

##Alternatives without imputation
#Recurrence subplots
{
  #Recurrence
  figS7AB_NoImputation <- getForestAndKMCoxFromData(imputedData = scaledNAAllAim4Data,
                                                    originalData = scaledNAAllAim4Data, #WARNING: When we do not want imputation, we need to use the clinical variables for patient selection too
                                                    freqs = recurrenceFullDataLassoDevianceFreqs,
                                                    minFreq = minLassoFreq,
                                                    eventVar = "EventRecurrence",
                                                    colors = c('#8F770099','#0073C299'),
                                                    kmTitle = "Recurrence",
                                                    forestTitle = "    Recurrence Proportional Hazards Model",
                                                    termDict = humanReadableHotEncodedVariables,
                                                    rel_height = relHeightFig7)
  
  #Noninvasive recurrence
  fig7AB_NoImputation <- getForestAndKMCoxFromData(imputedData = scaledNAAllAim4Data,
                                                   originalData = scaledNAAllAim4Data[EventProgression==F],
                                                   freqs = recurrenceOnlyFullDataLassoDevianceFreqs,
                                                   minFreq = minLassoFreq,
                                                   eventVar = "EventRecurrence",
                                                   colors = c('#8F770099','#0073C299'),
                                                   kmTitle = "Non-invasive Recurrence",
                                                   forestTitle = "    Non-invasive Recurrence Proportional Hazards Model",
                                                   termDict = humanReadableHotEncodedVariables,
                                                   rel_height = relHeightFig7)
}

#Progression subplots
{
  #Progression from nonprogressors
  figS7CD_NoImputation <- getForestAndKMCoxFromData(imputedData = scaledNAAllAim4Data,
                                                    originalData = scaledNAAllAim4Data,
                                                    freqs = progressionFullDataLassoDevianceFreqs,
                                                    minFreq = minLassoFreq,
                                                    eventVar = "EventProgression",
                                                    colors = c('#A7303099','#0073C299'),
                                                    kmTitle = "Progression from Nonprogressors",
                                                    forestTitle = "    Progression from Nonprogressors Proportional Hazards Model",
                                                    termDict = humanReadableHotEncodedVariables,
                                                    rel_height = relHeightFig7)
  #Progression only
  fig7CD_NoImputation <- getForestAndKMCoxFromData(imputedData = scaledNAAllAim4Data,
                                                   originalData = scaledNAAllAim4Data[EventRecurrence==EventProgression,],
                                                   freqs = progressionOnlyFullDataLassoDevianceFreqs,
                                                   minFreq = minLassoFreq,
                                                   eventVar = "EventProgression",
                                                   colors = c('#A7303099','#0073C299'),
                                                   kmTitle = "Progression",
                                                   forestTitle = "    Progression Proportional Hazards Model",
                                                   termDict = humanReadableHotEncodedVariables,
                                                   rel_height = relHeightFig7)
}

#Plot composition and saving
{
  ncol <- 2
  nrow <- 2
  
  figS7_NoImputation=plot_grid(figS7AB_NoImputation[[1]],figS7AB_NoImputation[[2]],figS7CD_NoImputation[[1]],figS7CD_NoImputation[[2]],labels=c("A","B","C","D"),label_size = gridLabelFontSize)
  
  save_plot(figS7_NoImputation,file=paste0(outDir,"figS7_NoImputation.pdf"),base_width = pageWidthInch, base_height = pageWidthInch*nrow/ncol/aspectRatio.plot/aspectRatioFig7Mod,device=cairo_pdf)
  save_plot(figS7_NoImputation,file=paste0(outDir,"figS7_NoImputation.png"),base_width = pageWidthInch, base_height = pageWidthInch*nrow/ncol/aspectRatio.plot/aspectRatioFig7Mod,base_asp = aspectRatio.plot)
  
  fig7_NoImputation=plot_grid(fig7AB_NoImputation[[1]],fig7AB_NoImputation[[2]],fig7CD_NoImputation[[1]],fig7CD_NoImputation[[2]],labels=c("A","B","C","D"),label_size = gridLabelFontSize)
  
  save_plot(fig7_NoImputation,file=paste0(outDir,"fig7_NoImputation.pdf"),base_width = pageWidthInch, base_height = pageWidthInch*nrow/ncol/aspectRatio.plot/aspectRatioFig7Mod,device=cairo_pdf)
  save_plot(fig7_NoImputation,file=paste0(outDir,"fig7_NoImputation.png"),base_width = pageWidthInch, base_height = pageWidthInch*nrow/ncol/aspectRatio.plot/aspectRatioFig7Mod)
  
}

#Figure 7 Sup with AdequateMargin
#################################
minLassoFreqAdequateMargin <- 0.1

#Using AdequateMargin instead of EstimatedMargin
clinicalVariablesThatMayBeConsideredAdequateMargin=c("AgeAtDCISDiagnosis", 
                                                     "MenopausalStatusCompleted", #recoded MenopausalStatus to use the Age at Diagnosis to reduce the missing data (threshold = 55 years, only for those without this information)
                                                     "RaceN",#recoded Race
                                                     "Height_cm", 
                                                     "Weight_kg",
                                                     "BMI",
                                                     "AxillaryDissection",
                                                     "NNodesExaminedN",# recoded NNodesExamined 
                                                     "DCISWasAtSurgicalMargin",
                                                     #"EstimatedMargin", #recoded "FinalMargin" + Margin_mm
                                                     "AdequateMargin", #recoded "FinalMargin" + Margin_mm with a 2mm threshold
                                                     "DCISSize_cm",
                                                     "NuclearGrade",
                                                     "NecrosisType",
                                                     #"Necrosis",
                                                     "MicrocalcType",
                                                     #"MicroCalcs",
                                                     #"DCISMicroCalcs",
                                                     #"BenignMicroCalcs",
                                                     "ER",
                                                     "PR",
                                                     "DCISTreatment",
                                                     #"SurgicalProcedure", #Recoded as DCISTreatment following Shelley's advice
                                                     #"RadiationDCIS", #Recoded as DCISTreatment following Shelley's advice
                                                     "HormonalTherapyDCIS",
                                                     #"TamoxifenDCIS",
                                                     #"HormonalTherapyType",
                                                     #To discard below this line, reason as a comment
                                                     "DCISDiagnosis_y", #Discarded because we expect it to be associated with the outcome, since the study finished at a given time
                                                     #"SurgicalProcedure_y", #chosen DCISDiagnosis_y, r=0.99 between the two
                                                     "HispanicOrLatino", #too few to take into account (7 only) 
                                                     "Laterality", #In principle not important
                                                     "Quadrant" #~41% of missing data
)

clinicalVariablesForNowAdequateMargin=clinicalVariablesThatMayBeConsideredAdequateMargin[1:(length(clinicalVariablesThatMayBeConsideredAdequateMargin)-4)]

#DataAdequateMargin prep
scaledImputedAllAim4DataAdequateMarginOriginal=scale(makeX(allAim4DataWClinical[,.SD,.SDcols=c(humanReadableVariables,clinicalVariablesForNowAdequateMargin)],na.impute = T))
scaledNAAllAim4DataAdequateMarginOriginal=scale(makeX(allAim4DataWClinical[,.SD,.SDcols=c(humanReadableVariables,clinicalVariablesForNowAdequateMargin)],na.impute = F))
xColumnNamesAdequateMargin=colnames(scaledImputedAllAim4DataAdequateMarginOriginal)
scaledImputedAllAim4DataAdequateMargin=data.table(scaledImputedAllAim4DataAdequateMarginOriginal,
                                                  allAim4DataWClinical[,.(EventRecurrence,EventProgression,TimeToEvent_mo,patient)])
scaledNAAllAim4DataAdequateMargin=data.table(scaledNAAllAim4DataAdequateMarginOriginal,
                                             allAim4DataWClinical[,.(EventRecurrence,EventProgression,TimeToEvent_mo,patient)])

#Recurrence

#ATTENTION: along this section, the recurrenceOnly tag means the proportional hazards version of cohort 0 vs cohort1, while recurrence means cohort 0 vs cohort 1 + cohort 2

#We always use imputation in the covariates for covariate selection using LASSO
recurrenceFullDataAdequateMarginLassoDevianceFreqs=getLassoFreq(seed=theSeed,
                                                                n=nLasso,
                                                                x=scaledImputedAllAim4DataAdequateMargin[fullDataPatients,.SD,.SDcols=xColumnNamesAdequateMargin],
                                                                y=Surv(scaledImputedAllAim4DataAdequateMargin[fullDataPatients,TimeToEvent_mo],scaledImputedAllAim4DataAdequateMargin[fullDataPatients,EventRecurrence]),
                                                                family = "cox", 
                                                                standarize = F, 
                                                                type.measure = "deviance",
                                                                lambda="lambda.min")

recurrenceOnlyFullDataAdequateMarginLassoDevianceFreqs=getLassoFreq(seed=theSeed,
                                                                    n=nLasso,
                                                                    x=scaledImputedAllAim4DataAdequateMargin[fullDataPatients][EventProgression==F,.SD,.SDcols=xColumnNamesAdequateMargin],
                                                                    y=Surv(scaledImputedAllAim4DataAdequateMargin[fullDataPatients][EventProgression==F,TimeToEvent_mo],scaledImputedAllAim4DataAdequateMargin[fullDataPatients][EventProgression==F,EventRecurrence]),
                                                                    family = "cox", 
                                                                    standarize = F, 
                                                                    type.measure = "deviance",
                                                                    lambda="lambda.min")

#Progression
progressionFullDataAdequateMarginLassoDevianceFreqs=getLassoFreq(seed=theSeed,
                                                                 n=nLasso,
                                                                 x=scaledImputedAllAim4DataAdequateMargin[fullDataPatients,.SD,.SDcols=xColumnNamesAdequateMargin],
                                                                 y=Surv(scaledImputedAllAim4DataAdequateMargin[fullDataPatients,TimeToEvent_mo],scaledImputedAllAim4DataAdequateMargin[fullDataPatients,EventProgression]),
                                                                 family = "cox", 
                                                                 standarize = F, 
                                                                 type.measure = "deviance",
                                                                 lambda="lambda.min")

#ProgressionOnly 
progressionOnlyFullDataAdequateMarginLassoDevianceFreqs=getLassoFreq(seed=theSeed,
                                                                     n=nLasso,
                                                                     x=scaledImputedAllAim4DataAdequateMargin[fullDataPatients][EventRecurrence==EventProgression,.SD,.SDcols=xColumnNamesAdequateMargin],
                                                                     y=Surv(scaledImputedAllAim4DataAdequateMargin[fullDataPatients][EventRecurrence==EventProgression,TimeToEvent_mo],scaledImputedAllAim4DataAdequateMargin[fullDataPatients][EventRecurrence==EventProgression,EventProgression]),
                                                                     family = "cox", 
                                                                     standarize = F, 
                                                                     type.measure = "deviance",
                                                                     lambda="lambda.min")

#ProgressionFromRecurrence (probably to remove)
progressionFromRecurrenceFullDataAdequateMarginLassoDevianceFreqs=getLassoFreq(seed=theSeed,
                                                                               n=nLasso,
                                                                               x=scaledImputedAllAim4DataAdequateMargin[fullDataPatients][EventRecurrence==T,.SD,.SDcols=xColumnNamesAdequateMargin],
                                                                               y=Surv(scaledImputedAllAim4DataAdequateMargin[fullDataPatients][EventRecurrence==T,TimeToEvent_mo],scaledImputedAllAim4DataAdequateMargin[fullDataPatients][EventRecurrence==T,EventProgression]),
                                                                               family = "cox",
                                                                               standarize = F,
                                                                               type.measure = "deviance",
                                                                               lambda="lambda.min")

#Recurrence subplots
{
  #Recurrence 
  figS7AB_AdequateMargin <- getForestAndKMCoxFromData(imputedData = scaledImputedAllAim4DataAdequateMargin,
                                                      originalData = allAim4Data, #We filter patients based on the non-imputed measured variables
                                                      freqs = recurrenceFullDataAdequateMarginLassoDevianceFreqs,
                                                      minFreq = minLassoFreqAdequateMargin,
                                                      eventVar = "EventRecurrence",
                                                      colors = c('#8F770099','#0073C299'),
                                                      kmTitle = "Recurrence",
                                                      forestTitle = "    Recurrence Proportional Hazards Model",
                                                      termDict = humanReadableHotEncodedVariables,
                                                      rel_height = relHeightFig7)
  
  #Noninvasive recurrence
  fig7AB_AdequateMargin <- getForestAndKMCoxFromData(imputedData = scaledImputedAllAim4DataAdequateMargin,
                                                     originalData = allAim4Data[EventProgression==F,],
                                                     freqs = recurrenceOnlyFullDataAdequateMarginLassoDevianceFreqs,
                                                     minFreq = minLassoFreqAdequateMargin,
                                                     eventVar = "EventRecurrence",
                                                     colors = c('#8F770099','#0073C299'),
                                                     kmTitle = "Non-invasive Recurrence",
                                                     forestTitle = "    Non-invasive Recurrence Proportional Hazards Model",
                                                     termDict = humanReadableHotEncodedVariables,
                                                     rel_height = relHeightFig7)
  
}

#Progression subplots
{
  
  #Progression from Nonprogressors
  figS7CD_AdequateMargin <- getForestAndKMCoxFromData(imputedData = scaledImputedAllAim4DataAdequateMargin,
                                                     originalData = allAim4Data,
                                                     freqs = progressionFullDataAdequateMarginLassoDevianceFreqs,
                                                     minFreq = minLassoFreqAdequateMargin,
                                                     eventVar = "EventProgression",
                                                     colors = c('#A7303099','#0073C299'),
                                                     kmTitle = "Progression from Nonprogressors",
                                                     forestTitle = "    Progression from Nonprogressors Proportional Hazards Model",
                                                     termDict = humanReadableHotEncodedVariables,
                                                     rel_height = relHeightFig7)
  
  #ProgressionOnly
  fig7CD_AdequateMargin <- getForestAndKMCoxFromData(imputedData = scaledImputedAllAim4DataAdequateMargin,
                                                     originalData = allAim4Data[EventProgression==EventRecurrence,],
                                                     freqs = progressionOnlyFullDataAdequateMarginLassoDevianceFreqs,
                                                     minFreq = minLassoFreqAdequateMargin,
                                                     eventVar = "EventProgression",
                                                     colors = c('#A7303099','#0073C299'),
                                                     kmTitle = "Progression",
                                                     forestTitle = "    Progression Proportional Hazards Model",
                                                     termDict = humanReadableHotEncodedVariables,
                                                     rel_height = relHeightFig7)
  
  #ProgressionFromRecurrence (probably to remove)
  figSFromRecurrenceCD_AdequateMargin <- getForestAndKMCoxFromData(imputedData = scaledImputedAllAim4DataAdequateMargin,
                                                                   originalData = allAim4Data[EventRecurrence==T,],
                                                                   freqs = progressionFromRecurrenceFullDataAdequateMarginLassoDevianceFreqs,
                                                                   minFreq = minLassoFreqAdequateMargin,
                                                                   eventVar = "EventProgression",
                                                                   colors = c('#A7303099','#0073C299'),
                                                                   kmTitle = "Progression from Recurrence",
                                                                   forestTitle = "    Progression from Recurrence Proportional Hazards Model",
                                                                   termDict = humanReadableHotEncodedVariables,
                                                                   rel_height = relHeightFig7)
}

#Plot composition and saving
{
  ncol <- 2
  nrow <- 2
  figS7_AdequateMargin <- plot_grid(figS7AB_AdequateMargin[[1]],figS7AB_AdequateMargin[[2]],figS7CD_AdequateMargin[[1]],figS7CD_AdequateMargin[[2]],labels=c("A","B","C","D"),label_size = gridLabelFontSize)
  
  save_plot(figS7_AdequateMargin,file=paste0(outDir,"figS7_AdequateMargin.pdf"),base_width = pageWidthInch, base_height = pageWidthInch*nrow/ncol/aspectRatio.plot/aspectRatioFig7Mod,device=cairo_pdf)
  save_plot(figS7_AdequateMargin,file=paste0(outDir,"figS7_AdequateMargin.png"),base_width = pageWidthInch, base_height = pageWidthInch*nrow/ncol/aspectRatio.plot/aspectRatioFig7Mod)
  
  fig7_AdequateMargin <- plot_grid(fig7AB_AdequateMargin[[1]],fig7AB_AdequateMargin[[2]],fig7CD_AdequateMargin[[1]],fig7CD_AdequateMargin[[2]],labels=c("A","B","C","D"),label_size = gridLabelFontSize)
  
  save_plot(fig7_AdequateMargin,file=paste0(outDir,"fig7_AdequateMargin.pdf"),base_width = pageWidthInch, base_height = pageWidthInch*nrow/ncol/aspectRatio.plot/aspectRatioFig7Mod,device=cairo_pdf)
  save_plot(fig7_AdequateMargin,file=paste0(outDir,"fig7_AdequateMargin.png"),base_width = pageWidthInch, base_height = pageWidthInch*nrow/ncol/aspectRatio.plot/aspectRatioFig7Mod)
  
  ncol <- 2
  nrow <- 1
  figSFromRecurrence_AdequateMargin <- plot_grid(figSFromRecurrenceCD_AdequateMargin[[1]],figSFromRecurrenceCD_AdequateMargin[[2]],labels=c("A","B"),label_size = gridLabelFontSize)
  save_plot(figSFromRecurrence_AdequateMargin,file=paste0(outDir,"figSFromRecurrence_AdequateMargin.pdf"),base_width = pageWidthInch, base_height = pageWidthInch*nrow/ncol/aspectRatio.plot/aspectRatioFig7Mod,device=cairo_pdf)
  save_plot(figSFromRecurrence_AdequateMargin,file=paste0(outDir,"figSFromRecurrence_AdequateMargin.png"),base_width = pageWidthInch, base_height = pageWidthInch*nrow/ncol/aspectRatio.plot/aspectRatioFig7Mod)
}

##Alternatives without imputation
#Recurrence subplots
{
  #Recurrence
  figS7AB_AdequateMargin_NoImputation <- getForestAndKMCoxFromData(imputedData = scaledNAAllAim4DataAdequateMargin,
                                                                   originalData = scaledNAAllAim4DataAdequateMargin, #WARNING: When we do not want imputation, we need to use the clinical variables for patient selection too
                                                                   freqs = recurrenceFullDataAdequateMarginLassoDevianceFreqs,
                                                                   minFreq = minLassoFreqAdequateMargin,
                                                                   eventVar = "EventRecurrence",
                                                                   colors = c('#8F770099','#0073C299'),
                                                                   kmTitle = "Recurrence",
                                                                   forestTitle = "    Recurrence Proportional Hazards Model",
                                                                   termDict = humanReadableHotEncodedVariables,
                                                                   rel_height = relHeightFig7)
  
  #Noninvasive recurrence
  fig7AB_AdequateMargin_NoImputation <- getForestAndKMCoxFromData(imputedData = scaledNAAllAim4DataAdequateMargin,
                                                                  originalData = scaledNAAllAim4DataAdequateMargin[EventProgression==F],
                                                                  freqs = recurrenceOnlyFullDataAdequateMarginLassoDevianceFreqs,
                                                                  minFreq = minLassoFreqAdequateMargin,
                                                                  eventVar = "EventRecurrence",
                                                                  colors = c('#8F770099','#0073C299'),
                                                                  kmTitle = "Non-invasive Recurrence",
                                                                  forestTitle = "    Non-invasive Recurrence Proportional Hazards Model",
                                                                  termDict = humanReadableHotEncodedVariables,
                                                                  rel_height = relHeightFig7)
}

#Progression subplots
{
  #Progression from nonprogressors
  figS7CD_AdequateMargin_NoImputation <- getForestAndKMCoxFromData(imputedData = scaledNAAllAim4DataAdequateMargin,
                                                                   originalData = scaledNAAllAim4DataAdequateMargin,
                                                                   freqs = progressionFullDataAdequateMarginLassoDevianceFreqs,
                                                                   minFreq = minLassoFreqAdequateMargin,
                                                                   eventVar = "EventProgression",
                                                                   colors = c('#A7303099','#0073C299'),
                                                                   kmTitle = "Progression from Nonprogressors",
                                                                   forestTitle = "    Progression from Nonprogressors Proportional Hazards Model",
                                                                   termDict = humanReadableHotEncodedVariables,
                                                                   rel_height = relHeightFig7)
  #Progression only
  fig7CD_AdequateMargin_NoImputation <- getForestAndKMCoxFromData(imputedData = scaledNAAllAim4DataAdequateMargin,
                                                                  originalData = scaledNAAllAim4DataAdequateMargin[EventRecurrence==EventProgression,],
                                                                  freqs = progressionOnlyFullDataAdequateMarginLassoDevianceFreqs,
                                                                  minFreq = minLassoFreqAdequateMargin,
                                                                  eventVar = "EventProgression",
                                                                  colors = c('#A7303099','#0073C299'),
                                                                  kmTitle = "Progression",
                                                                  forestTitle = "    Progression Proportional Hazards Model",
                                                                  termDict = humanReadableHotEncodedVariables,
                                                                  rel_height = relHeightFig7)
}

#Plot composition and saving
{
  ncol <- 2
  nrow <- 2
  
  figS7_AdequateMargin_NoImputation=plot_grid(figS7AB_AdequateMargin_NoImputation[[1]],figS7AB_AdequateMargin_NoImputation[[2]],figS7CD_AdequateMargin_NoImputation[[1]],figS7CD_AdequateMargin_NoImputation[[2]],labels=c("A","B","C","D"),label_size = gridLabelFontSize)
  
  save_plot(figS7_AdequateMargin_NoImputation,file=paste0(outDir,"figS7_AdequateMargin_NoImputation.pdf"),base_width = pageWidthInch, base_height = pageWidthInch*nrow/ncol/aspectRatio.plot/aspectRatioFig7Mod,device=cairo_pdf)
  save_plot(figS7_AdequateMargin_NoImputation,file=paste0(outDir,"figS7_AdequateMargin_NoImputation.png"),base_width = pageWidthInch, base_height = pageWidthInch*nrow/ncol/aspectRatio.plot/aspectRatioFig7Mod,base_asp = aspectRatio.plot)
  
  fig7_AdequateMargin_NoImputation=plot_grid(fig7AB_AdequateMargin_NoImputation[[1]],fig7AB_AdequateMargin_NoImputation[[2]],fig7CD_AdequateMargin_NoImputation[[1]],fig7CD_AdequateMargin_NoImputation[[2]],labels=c("A","B","C","D"),label_size = gridLabelFontSize)
  
  save_plot(fig7_AdequateMargin_NoImputation,file=paste0(outDir,"fig7_AdequateMargin_NoImputation.pdf"),base_width = pageWidthInch, base_height = pageWidthInch*nrow/ncol/aspectRatio.plot/aspectRatioFig7Mod,device=cairo_pdf)
  save_plot(fig7_AdequateMargin_NoImputation,file=paste0(outDir,"fig7_AdequateMargin_NoImputation.png"),base_width = pageWidthInch, base_height = pageWidthInch*nrow/ncol/aspectRatio.plot/aspectRatioFig7Mod)
  
}

#Other supplementary figures and paper stats
############################################

#Correlogram?
# cor(scaledImputedAllAim4Data[,-c("EventRecurrence","EventProgression","TimeToEvent_mo","patient")],use="pairwise.complete.obs") %>%
#   ggcorrplot(show.diag=F,type="lower",lab=T,lab_size=2,hc.order=T, insig = "blank") +
#   labs(title="Correlogram of (one-hot encoded) clinical variables")
# 
# cor(scaledNAAllAim4Data[,-c("EventRecurrence","EventProgression","TimeToEvent_mo","patient")],use="pairwise.complete.obs") %>%
#   ggcorrplot(show.diag=F,type="lower",lab=T,lab_size=2,hc.order=T, insig = "blank") +
#   labs(title="Correlogram of (one-hot encoded) clinical variables")

#Median follow-up time using the reverse KM method
reverseKM <- surv_fit(Surv(TimeToEvent_mo, noEvent) ~ 1, data = allAim4DataWClinical[,.(TimeToEvent_mo,noEvent=!(EventRecurrence | EventProgression))])
surv_median(reverseKM)
ggsurvplot(reverseKM)

#Correlation between CNA burden and divergence
cor.test(~CNADivergence + MeanAlteredP,data = allAim4DataWClinical,method="spearman")
cor.test(~CNADivergence + MeanAlteredP,data = allAim4DataWClinical[Cohort==0,],method="spearman")
cor.test(~CNADivergence + MeanAlteredP,data = allAim4DataWClinical[Cohort==1,],method="spearman")
cor.test(~CNADivergence + MeanAlteredP,data = allAim4DataWClinical[Cohort==2,],method="spearman")

