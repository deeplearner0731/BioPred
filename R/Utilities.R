#' Tutorial Data
#'
#' A dataset containing sample data for demonstrating the functionalities of the BioPred package.
#'
#' @format A data frame with the following columns:
#' \describe{
#'   \item{x1}{Numeric. A biomarker variable.}
#'   \item{x2}{Numeric. A biomarker variable.}
#'   \item{x3}{Numeric. A biomarker variable.}
#'   \item{x4}{Numeric. A biomarker variable.}
#'   \item{x5}{Numeric. A biomarker variable.}
#'   \item{x6}{Numeric. A biomarker variable.}
#'   \item{x7}{Numeric. A biomarker variable.}
#'   \item{x8}{Numeric. A biomarker variable.}
#'   \item{x9}{Numeric. A biomarker variable.}
#'   \item{x10}{Numeric. A biomarker variable.}
#'   \item{y.con}{Numeric. A continuous outcome variable.}
#'   \item{y.bin}{Binary. A binary outcome variable, where 0 represents one class and 1 represents another class.}
#'   \item{y.time}{Numeric. The time in months, used for survival analysis.}
#'   \item{y.event}{Binary. Event indicator variable, where 0 indicates censoring and 1 indicates the event of interest occurred.}
#'   \item{subgroup_label}{Binary. Ground truth of subgroup label. In real-world scenarios, this information is typically unavailable.}
#'   \item{treatment}{Binary. Treatment indicator variable, where 0 represents control and 1 represents treatment.}
#'   \item{treatment_categorical}{Factor. A categorical version of the treatment variable, with levels "Placebo" and "Treatment".}
#'   \item{risk_category}{Factor.}
#' }
#'
#' @details
#' This dataset is used to illustrate various functions within the BioPred package, including predictive modeling and subgroup analysis. The columns represent different types of data typically encountered in clinical studies.
#'
#' @usage data(tutorial_data)
#'
#' @keywords datasets
#'
#' @examples
#' data(tutorial_data)
#' head(tutorial_data)
"tutorial_data"











#### outlier detection based on boxplot ####

#' Outlier Detection Based on Boxplot
#'
#' This function identifies outliers in a numeric vector based on the interquartile range (IQR) method used in boxplots.
#'
#' @param x A numeric vector of values for outlier detection.
#' @return A logical vector indicating which elements of the input vector are outliers (TRUE if an outlier, FALSE otherwise).
#' @keywords internal
isout <- function(x) {
  iqr <- stats::IQR(x, na.rm = TRUE)
  up <- stats::quantile(x, probs = 3/4, na.rm = TRUE)
  lo <- stats::quantile(x, probs = 1/4, na.rm = TRUE)
  up <- up + 1.5 * iqr
  lo <- lo - 1.5 * iqr
  out <- x > up | x < lo
  out
}



#' Plot Predictive Biomarker Importance based on XGBoost-based Subgroup Model
#'
#' This function calculates and plots the importance of biomarkers in a trained XGBoostSub_con, XGBoostSub_bin or XGBoostSub_sur model.
#'
#' @param model The trained XGBoost-based model.
#' @return A barplot showing the biomarker importance.
#' @import xgboost
#' @export
predictive_biomarker_imp <- function(model) {
  # Calculate feature importance
  importance <- xgb.importance(model = model)

  # Plot feature importance
  graphics::barplot(importance$Gain,
          main = "Biomarker Importance",
          xlab = "Gain",
          ylab = "Biomarker",
          col = "skyblue",
          horiz = TRUE,
          names.arg = importance$Feature)

  # Return data frame containing feature importance
  importance_df <- data.frame(Biomarker = importance$Feature, Importance = importance$Gain)
  return(importance_df)
}




#' Get Subgroup Results
#'
#' This function predicts the treatment assignment for each patient based on a cutoff value.
#'
#' @param model The trained XGBoost-based subgroup model.
#' @param X_feature The data matrix containing patient features.
#' @param subgroup_label (Optional) The subgroup labels. In real-world data, this information is typically unknown and only available in simulated data. If provided, the prediction accuracy will also be returned.
#' @param cutoff The cutoff value for treatment assignment, defaulted to 0.5.
#' @return A data frame containing each subject and assigned treatment (1 for treatment, 0 for control). If subgroup labels are provided, it also returns the prediction accuracy of the subgroup labels.
#' @import xgboost
#' @export
get_subgroup_results <- function(model, X_feature, subgroup_label = NULL, cutoff = 0.5) {
  # Convert data to xgb.DMatrix
  if (is.null(subgroup_label)) {
    subgroup_label <- sample(c(0, 1), nrow(X_feature), replace = TRUE)
  }
  dtrain <- xgboost::xgb.DMatrix(data = as.matrix(X_feature), label = subgroup_label)
  # Predict treatment assignment for each patient
  pred <- stats::predict(model, dtrain)
  pred_prob <- 1.0/(1.0+exp(-pred))
  # Assign treatment based on the cutoff
  prediction <- ifelse(pred_prob > cutoff, 1, 0)

  acc <- NULL
  if (!is.null(subgroup_label)) {
    acc <- mean(subgroup_label == prediction)
  }
  # Create a data frame with patient ID and assigned treatment
  results <- data.frame(Patient_ID = seq_len(nrow(X_feature)),
                        Treatment_Assignment = prediction)

  if (!is.null(acc)) {
    return(list(assignment=results, acc = acc))
  }
  return(results)
}

#' CDF  Plot for a biomarker
#'
#' Cumulative Distribution Function (CDF) plot for a biomarker.
#'
#' @param xvar The biomarker name.
#' @param data The dataset.
#' @param y.int Increasement interval on the y.
#' @param xlim cdf plot range for xvar, when NULL, c(min(x), max(x)) will be used.
#' @param xvar.display Display name of the biomarker.
#' @param group A separate CDF line will be plotted for each group.
#' @return A ggplot object representing the CDF inverse plot.
#' @import ggplot2
#' @export
cdf_plot = function(xvar, data, y.int=5, xlim=NULL, xvar.display=xvar, group=NULL){
   num=NULL
   if(is.null(group)||is.na(group)){
    x = data[,xvar]
    x = x[!is.na(x)]

    if(is.null(xlim)) xlim=c(min(x),max(x))

    plot.data = NULL
    x.axis.vals = sort(unique(c(seq(from=xlim[1],to=xlim[2],length.out=1000), x)))
    for(x.uni.i in x.axis.vals){
      plot.data.i = data.frame(x=x.uni.i, num = mean(x<=x.uni.i)*100)
      plot.data = rbind(plot.data, plot.data.i)
    }


    fig = ggplot(data=plot.data, aes(x=x, y=num)) +
      geom_line(size=0.8) +
      labs(title=paste("Cumulative Distribution of", xvar.display) , x=xvar.display, y="% of Subjects <= Cutoff")+
      scale_y_continuous(breaks=seq(0, 100, y.int),limits = c(0,100))+
      scale_x_continuous(breaks=round(seq(from=xlim[1],to=xlim[2],length.out=11), digits=2))+
      theme_bw(base_size=15)+
      theme(text=element_text(size=15), axis.text=element_text(size=15),legend.text=element_text(size=15))

  }else{
    data = data[,c(xvar, group)]
    data = data[stats::complete.cases(data),]

    if(is.null(xlim)) xlim=c(min(data[,xvar]), max(data[,xvar]))

    plot.data = NULL
    x.axis.vals = sort(unique(c(seq(from=xlim[1],to=xlim[2],length.out=1000), data[,xvar])))

    for(grp.i in unique(data[,group])){
      x.i = data[data[,group]==grp.i, xvar]
      for(x.uni.i in x.axis.vals){
        plot.data.i = data.frame(x=x.uni.i, num = mean(x.i<=x.uni.i)*100, group=grp.i)
        plot.data = rbind(plot.data, plot.data.i)
      }
    }

    fig = ggplot2::ggplot(data=plot.data, aes(x=x, y=num, group=group, colour=group)) +
      geom_line(size=0.8) +
      labs(title=paste("Cumulative Distribution of", xvar.display) , x=xvar.display, y="% of Subjects <= Cutoff", colour = group)+
      scale_y_continuous(breaks=seq(0, 100, y.int),limits = c(0,100))+
      scale_x_continuous(breaks=round(seq(from=xlim[1],to=xlim[2],length.out=11), digits=2))+
      theme_bw(base_size=15)+
      theme(text=element_text(size=15), axis.text=element_text(size=15),legend.text=element_text(size=15))
  }
  fig
}

#' AUC ROC Table for Biomarkers Associated with Binary Outcomes
#'
#' Computes the area under the receiver operating characteristic (ROC) curve for Biomarkers Associated with Binary Outcomes,
#' and returns the results as a table.
#'
#' @param yvar Binary response variable name, where 0 represents controls and 1 represents cases.
#' @param xvars A vector of biomarker names.
#' @param dirs A vector of directions for the biomarkers. Options are "auto", ">", or "<".
#'   - "auto" (default): automatically determines in which group the median is higher and takes the direction accordingly.
#'   - ">": indicates that the biomarkers for the control group are higher than those for the case group (controls > t >= cases).
#'   - "<": indicates that the biomarkers for the control group are lower or equal to those for the case group (controls < t <= cases).
#' @param data The dataset containing the variables.
#' @param yvar.display Display name for the binary response variable.
#' @param xvars.display Display names for the biomarkers.
#' @return A table containing the AUC values for each biomarker.
#' @import pROC
#' @export
roc_bin = function(yvar, xvars, dirs, data, yvar.display=yvar, xvars.display=xvars){

  if(!all(data[,yvar]%in%c(0,1,NA)) | !all(c(0,1)%in%data[,yvar])) stop("Response needs to be 0/1.")
  if(length(xvars)!=length(dirs)) stop("xvars and dirs should be of the same length.")
  data[,yvar] = factor(data[,yvar], levels=c(0, 1))

  res.all = NULL
  for(i in 1:length(xvars)){
    xvar = xvars[i]
    dir = dirs[i]

    if(dir=="<") {
      dir=">"
    }else if(dir==">") {
      dir="<"
    }

    data.temp = data[,c(yvar,xvar)]
    data.temp = data.temp[stats::complete.cases(data.temp),]

    fml = stats::as.formula(paste(yvar,"~", xvar))
    res.i = roc(fml, data=data.temp, percent=F, na.rm=T, direction=dir, algorithm=1, quiet=T, smooth=F, auc=T, ci=F, plot=F)
    if(res.i$direction == ">") {
      res.i$direction = "<"
    }else if(res.i$direction == "<") {
      res.i$direction = ">"
    }
    res.i = data.frame(Response = yvar.display, Predictor = xvars.display[i], direction=paste("case",res.i$direction,"ctrl"),
                       AUC = as.numeric(res.i$auc))
    res.all = rbind(res.all, res.i)

  }

  res.all
}

#' ROC Plot Biomarkers Associated with Binary Outcomes
#'
#' Generates ROC plots for different biomarkers associated with binary outcomes.
#'
#' @param yvar Binary response variable name, where 0 represents controls and 1 represents cases.
#' @param xvars A vector of biomarker names.
#' @param dirs A vector of directions for the biomarkers. Options are "auto", ">", or "<".
#'   - "auto" (default): automatically determines in which group the median is higher and takes the direction accordingly.
#'   - ">" indicates that the biomarkers for the control group are higher than those for the case group (controls > t >= cases).
#'   - "<" indicates that the biomarkers for the control group are lower or equal to those for the case group (controls < t <= cases).
#' @param data The dataset containing the variables.
#' @param yvar.display Display name for the binary response variable.
#' @param xvars.display Display names for the biomarkers.
#' @return ROC plots for different biomarkers associated with binary outcomes.
#' @import pROC
#' @export
roc_bin_plot = function(yvar, xvars, dirs, data, yvar.display=yvar, xvars.display=xvars){

  if(!all(data[,yvar]%in%c(0,1,NA)) | !all(c(0,1)%in%data[,yvar])) stop("Response needs to be 0/1.")
  if(length(xvars)!=length(dirs)) stop("xvars and dirs should be of the same length.")
  data[,yvar] = factor(data[,yvar], levels=c(0, 1))

  for(i in 1:length(xvars)){
    xvar = xvars[i]
    dir = dirs[i]

    if(dir=="<") {
      dir=">"
    }else if(dir==">") {
      dir="<"
    }

    data.temp = data[,c(yvar,xvar)]
    data.temp = data.temp[stats::complete.cases(data.temp),]

    fml = stats::as.formula(paste(yvar,"~", xvar))
    res.i = roc(fml, data=data.temp, percent=F, na.rm=T, direction=dir, algorithm=1, quiet=T, smooth=F, auc=T, ci=F, plot=F)
    plot.roc(res.i, add=F, reuse.auc=T, axes=T, print.auc=T, grid=T, main=paste(yvar.display, "~", xvars.display[i]))
  }

}



#' Scatter Plot for a Biomarker Associated with Continuous Outcome
#'
#' Generates a scatter plot for exploring the relationship between a continuous response variable
#' and a biomarker variable.
#'
#' @param yvar Continuous response variable name.
#' @param xvar biomarker name.
#' @param data The dataset containing the variables.
#' @param ybreaks Breaks on the y-axis.
#' @param xbreaks Breaks on the x-axis.
#' @param yvar.display Display name for the response variable.
#' @param xvar.display Display name for the biomarker variable.
#' @return A list containing correlation coefficients, scatter plot, slope, and intercept.
#' @export
scat_cont_plot = function(yvar, xvar, data, ybreaks=NULL, xbreaks=NULL, yvar.display=yvar, xvar.display=xvar){
  data = data[,c(yvar,xvar)]
  data = data[stats::complete.cases(data),]
  names(data) = c("y","x")

  y=data[["y"]]
  x=data[["x"]]
  N = length(x)

  ## pearson ##
  cor.pears = round(stats::cor(x, y, use="complete.obs",method="pearson"),digits=2)
  cor.pears.p = round(stats::cor.test(x, y, alternative="two.sided", method="pearson")$p.value, digits=4)
  ## spearman ##
  cor.spear = round(stats::cor(x, y, use="complete.obs",method="spearman"),digits=2)
  cor.spear.p = round(stats::cor.test(x, y, alternative="two.sided", method="spearman")$p.value, digits=4)

  lm.fit = stats::lm(y~x)
  slope = summary(lm.fit)$coefficients["x","Estimate"]
  intercept = summary(lm.fit)$coefficients["(Intercept)","Estimate"]
  jitter.width = (max(x)-min(x))/80
  jitter.height = (max(y)-min(y))/80


  if(is.null(xbreaks)) xbreaks=round(seq(min(x), max(x), length.out = 10),digits=2)
  if(is.null(ybreaks)) ybreaks=round(seq(min(y), max(y), length.out = 10),digits=2)

  scat.plot = ggplot(data, aes(x=x, y=y)) +
    geom_point(size=3,shape=1, stroke=1.5,position=position_jitter(width=jitter.width,height=jitter.height))+
    theme_bw(base_size=15)+
    theme(text=element_text(size=15), axis.text=element_text(size=15),legend.text=element_text(size=15))+
    geom_abline(aes(slope=slope,intercept =intercept, color="Fitted Line"), size=0.8)+
    labs(color="") +
    xlab(xvar.display)+
    ylab(yvar.display)+
    scale_y_continuous(breaks=ybreaks)+
    scale_x_continuous(breaks=xbreaks)+
    ggtitle(paste(xvar.display," vs. ",yvar.display, " (N = ", N, ")\nCor.Pearson = ",cor.pears,
                  ", Cor.Spearman = ",cor.spear, "\nPearson.p = ",cor.pears.p,
                  ", Spearman.p = ",cor.spear.p, sep=""))+
    theme(plot.title = element_text(hjust = 0.5, size=15))


  list(Cor.Pearson = cor.pears, Cor.Spearman = cor.spear, fig=scat.plot, slope=slope, intercept=intercept)


}



#' GAM Plot
#'
#' Generates a generalized additive model (GAM) plot for exploring the relationship between a response
#' variable and a biomarker.
#'
#' @param yvar Response variable name.
#' @param censorvar Censoring variable name for survival analysis (0-censored, 1-event).
#' @param xvar Biomarker name.
#' @param xvars.adj Potential confounding variables to adjust for using linear terms.
#' @param sxvars.adj Potential confounding variables to adjust for using curve terms.
#' @param type "c" for continuous, "s" for survival, and "b" for binary response.
#' @param data The dataset containing the variables.
#' @param k Upper limit on the degrees of freedom associated with an s smooth.
#' @param pred.type "iterms" for trend of xvar, "response" for Y at the original scale.
#' @param link.scale Whether to show the plot in the scale of the link function.
#' @param title Title of the plot.
#' @param ybreaks Breaks on the y-axis.
#' @param xbreaks Breaks on the x-axis.
#' @param rugcol.var Variable name defining the color of the rug and points.
#' @param add.points Whether to add data points to the plot.
#' @param prt.sum Whether to print summary or not.
#' @param prt.chk Whether to print model diagnosis.
#' @param outlier.rm Whether to remove outliers based on 1.5IQR.
#' @param newdat User-supplied customized data for prediction and plotting.
#' @return A list containing p-table, s-table, GAM summary, GAM check, and the plot.
#' @import mgcv
#' @export
gam_plot = function(yvar, censorvar=NULL, xvar, xvars.adj=NULL, sxvars.adj=NULL, type, data, k,
                    pred.type="iterms", link.scale=TRUE, title="Trend Plot", ybreaks=NULL, xbreaks=NULL,
                    rugcol.var=NULL, add.points=FALSE, prt.sum=TRUE, prt.chk=FALSE, outlier.rm=FALSE, newdat=NULL){
  x=NULL
  yfit=NULL
  yfit.Lo=NULL
  yfit.Up=NULL
  rug=NULL
  rugcol=NULL
  y=NULL

  if(!pred.type%in%c("iterms","response")) stop("pred.type can only be iterms or response.")
  if(type=="s"){
    if(is.null(censorvar)) stop("Censoring variable missing!")
  }else if(type=="c"){
    type="c"
  }else if(type=="b"){
    if(!all(data[,yvar]%in%c(0,1,NA)) | !all(c(0,1)%in%data[,yvar])) stop("Response needs to be 0/1.")
  }else{
    stop("type can only be s, c or b.")
  }

  data = data[, c(yvar,censorvar,xvar, xvars.adj, sxvars.adj, rugcol.var)]
  data = data[stats::complete.cases(data[, c(yvar,censorvar,xvar, xvars.adj, sxvars.adj)]),]
  if(outlier.rm){
    data = data[!isout(data[,xvar]),]
  }


  if(is.null(xvars.adj)){
    par.fml = NULL
  }else{
    par.fml = paste(xvars.adj,collapse=" + ", sep="")
    par.fml = paste(par.fml,"+", sep="")
  }

  if(is.null(sxvars.adj)) {
    s.fml = NULL
  }else{
    s.fml = paste("s(",sxvars.adj,",k=",k,")",collapse="+",sep="")
    s.fml = paste(s.fml,"+",sep="")
  }


  s.fml = paste(s.fml,"s(",xvar,",k=",k,")",sep="")

  fml = stats::as.formula(paste(yvar,"~",par.fml,s.fml))

  if(type=="s"){
    res.gam = gam(fml, family=cox.ph(), data=data, weights=data[,censorvar], method="REML", select=F)
  }else if(type=="c"){
    res.gam = gam(fml, family=stats::gaussian(), data=data, method="REML", select=F)
  }else if(type=="b"){
    res.gam = gam(fml, family=stats::binomial(link = "logit"), data=data, method="REML", select=F)
  }

  ptable = summary(res.gam)$p.table
  stable = summary(res.gam)$s.table

  ptable.display = ptable
  if(!is.null(ptable.display)) ptable.display = round(ptable.display, digits=4)

  stable.display = stable
  if(!is.null(stable.display)) stable.display = round(stable.display, digits=4)

  if(prt.sum) {
    cat("GAM SUMMARY")
    print(summary(res.gam))
    cat("\n")
  }
  if(prt.chk) {
    cat("GAM CHECK")
    gam.check(res.gam)
    cat("\n")

  }

  if(is.null(newdat)){
    newdat = data.frame(x=seq(min(data[,xvar]), max(data[,xvar]), length.out=999))
    names(newdat)= xvar
    for(i in c(xvars.adj,sxvars.adj)){
      newdat[,i] = stats::median(data[,i])
    }
  }


  if(pred.type=="iterms"){
    x.term = paste("s(",xvar,")",sep="")
    dat.pred = stats::predict(object=res.gam, newdata=newdat, type="iterms",se.fit=T,
                       terms=x.term)

    plot.data = data.frame(yfit = dat.pred$fit[,x.term],
                           yfit.Up = dat.pred$fit[,x.term]+1.96*dat.pred$se.fit[,x.term],
                           yfit.Lo = dat.pred$fit[,x.term]-1.96*dat.pred$se.fit[,x.term], x=newdat[,xvar])

    if(type=="c"){
      ylab = "Mean Difference from Population Average"
      if(!link.scale) message("link.scale is ignored for type==c as link.scale can only be T for type==c.")
    }else if(type=="b"){
      ylab = "Log(OR) against Population Average"
      if(!link.scale){
        plot.data = data.frame(yfit = exp(plot.data$yfit), yfit.Up = exp(plot.data$yfit.Up),
                               yfit.Lo = exp(plot.data$yfit.Lo), x=plot.data[,"x"])
        ylab="OR against Population Average"
      }

    }else if(type=="s"){
      ylab = "Log(HR) against Population Average"
      if(!link.scale){
        plot.data = data.frame(yfit = exp(plot.data$yfit), yfit.Up = exp(plot.data$yfit.Up),
                               yfit.Lo = exp(plot.data$yfit.Lo), x=plot.data[,"x"])
        ylab="HR against Population Average"
      }

    }

  }else if(pred.type=="response"){

    if(type=="s") stop("time to event cannnot do pred.type=response.")


    dat.pred = stats::predict(object=res.gam, newdata=newdat, type="response",se.fit=T)

    plot.data = data.frame(yfit = as.vector(dat.pred$fit),
                           yfit.Up = as.vector(dat.pred$fit)+1.96*as.vector(dat.pred$se.fit),
                           yfit.Lo = as.vector(dat.pred$fit)-1.96*as.vector(dat.pred$se.fit),
                           x=newdat[,xvar])
    ylab = yvar
  }



  if(is.null(xbreaks)) xbreaks=round(seq(min(plot.data$x), max(plot.data$x), length.out = 10),digits=2)
  if(is.null(ybreaks)) ybreaks=round(seq(min(plot.data$yfit.Lo), max(plot.data$yfit.Up), length.out = 10),digits=2)

  if(is.null(rugcol.var)){
    data.rug = data.frame(rug=data[,xvar], rugcol="Individual Record", y=data[,yvar])
  }else{
    data.rug = data.frame(rug=data[,xvar], rugcol=data[,rugcol.var], y=data[,yvar])
  }

  rugcol.levels = as.character(sort(unique(data.rug$rugcol), na.last = T))
  rugcol.levels[is.na(rugcol.levels)] = "N/A"
  data.rug$rugcol = as.character(data.rug$rugcol)
  data.rug$rugcol[is.na(data.rug$rugcol)] = "N/A"
  data.rug$rugcol=factor(data.rug$rugcol, levels=rugcol.levels)

  if(pred.type=="response" & add.points){
    fig=ggplot(data=plot.data, aes(x=x, y=yfit)) +
      geom_ribbon(aes(ymin=yfit.Lo, ymax=yfit.Up),alpha=0.2) +
      geom_line(size=0.8)+
      geom_rug(aes(x=rug, colour=rugcol),data=data.rug, alpha=0.7, sides="b", inherit.aes=F)+
      geom_point(aes(x=rug, y=y, colour=rugcol),data=data.rug, size=2.5, inherit.aes=F)+
      labs(color="") +
      ylab(ylab)+
      xlab(xvar)+
      scale_y_continuous(breaks=ybreaks)+
      scale_x_continuous(breaks=xbreaks)+
      ggtitle(title)+
      theme_bw(base_size=15)+
      theme(text=element_text(size=15), axis.text=element_text(size=15),legend.text=element_text(size=15))


  }else{
    fig=ggplot(data=plot.data, aes(x=x, y=yfit)) +
      geom_ribbon(aes(ymin=yfit.Lo, ymax=yfit.Up),alpha=0.2) +
      geom_line(size=0.8)+
      geom_rug(aes(x=rug, colour=rugcol),data=data.rug, alpha=0.7, sides="b", inherit.aes=F)+
      labs(color="") +
      ylab(ylab)+
      xlab(xvar)+
      scale_y_continuous(breaks=ybreaks)+
      scale_x_continuous(breaks=xbreaks)+
      ggtitle(title)+
      theme_bw(base_size=15)+
      theme(text=element_text(size=15), axis.text=element_text(size=15),legend.text=element_text(size=15))
  }



  res = list(ptable=ptable, ptable.display=ptable.display, stable=stable, stable.display=stable.display,
             fig=fig, plot.data = plot.data)
  res
}





#' GAM Contrast Plot
#'
#' Computes and plots the contrasts between treatment and control group based on a GAM for exploring the relationship be-tween treatment benefit and biomarker.
#'
#' @param yvar Response variable name.
#' @param censorvar Censoring variable name (0-censored, 1-event). Required if type is "s" (survival).
#' @param xvar Biomarker name.
#' @param xvars.adj Potential confounding variables to adjust for using linear terms.
#' @param sxvars.adj Potential confounding variables to adjust for using curves.
#' @param trtvar Treatment variable that the contrast will build upon (treatment-control).
#' @param type Type of response variable. Options are "c" for continuous, "s" for survival, and "b" for binary response.
#' @param data The dataset containing the variables.
#' @param k Upper limit on the degrees of freedom associated with an s smooth.When this k is too large, program will report error saying
#   "Model has more coefficients than data."
#' @param title Title of the plot.
#' @param ybreaks Breaks on the y-axis.
#' @param xbreaks Breaks on the x-axis.
#' @param rugcol.var Variable name that defines the color of the rug.
#' @param link.scale Whether to show the plot (y-axis) in the scale of the link function (linear predictor).
#' @param prt.sum Whether to print summary or not.
#' @param prt.chk Whether to print model diagnosis.
#' @param outlier.rm Whether to remove outliers based on 1.5IQR.
#' @return A list containing the p-value table, summarized p-value table, s-value table, summarized s-value table, and the plot.
#' @import mgcv
#' @export
gam_ctr_plot = function(yvar, censorvar=NULL, xvar, xvars.adj=NULL, sxvars.adj=NULL, trtvar=NULL, type,
                        data, k, title="Group Contrast", ybreaks=NULL, xbreaks=NULL, rugcol.var=NULL,
                        link.scale=TRUE, prt.sum=TRUE, prt.chk=FALSE, outlier.rm=FALSE){
   x=NULL
   yfit=NULL
   yfit.Lo=NULL
   yfit.Up=NULL
   rug=NULL
   rugcol=NULL

  if(type=="s"){
    if(is.null(censorvar)) stop("Censoring variable missing!")
  }else if(type=="c"){
    type="c"
  }else if(type=="b"){
    if(!all(data[,yvar]%in%c(0,1,NA)) | !all(c(0,1)%in%data[,yvar])) stop("Response needs to be 0/1.")
  }else{
    stop("type can only be s, c or b.")
  }

  if(is.null(trtvar)) stop("trtvar is missing. Contrast is between treatments.")
  if(!is.numeric(data[,trtvar])) stop("Column of trtvar should be a numeric 0/1 vector.")
  if(!all(data[,trtvar]%in%c(0,1,NA)) || !all(c(0,1)%in%data[,trtvar])) stop("Column of trtvar needs to be 0/1.")


  data = data[, c(yvar,censorvar, xvar, xvars.adj, sxvars.adj, trtvar, rugcol.var)]
  data = data[stats::complete.cases(data[,c(yvar,censorvar, xvar, xvars.adj, sxvars.adj, trtvar)]),]

  if(outlier.rm){
    data = data[!isout(data[,xvar]),]
  }



  if(type!="s"){
    if("intercept"%in%names(data)) stop("Variable names should not include 'intercept'.")
    data$intercept = 1
    xvars.adj = c("intercept", xvars.adj)
  }


  if(is.null(xvars.adj)){
    par.fml = NULL
  }else{
    par.fml = paste(xvars.adj,collapse="+", sep="")
    par.fml = paste(par.fml,"+", sep="")
  }
  par.fml = paste(par.fml, trtvar, "+", sep="")


  if(is.null(sxvars.adj)) {
    s.fml = NULL
  }else{
    s.fml = paste("s(",sxvars.adj,",k=",k,")",collapse="+",sep="")
    s.fml = paste(s.fml,"+",sep="")
  }


  s.fml = paste(s.fml,"s(",xvar,",k=",k,")+s(",xvar,",k=",k, ",by=",trtvar,", pc=",
                mean(data[,xvar]),")",sep="")
  if(type!="s"){
    fml = stats::as.formula(paste(yvar,"~ -1+",par.fml,s.fml, sep=""))
  }else{
    fml = stats::as.formula(paste(yvar,"~",par.fml,s.fml, sep=""))
  }


  if(type=="s"){
    res.gam = gam(fml, family=cox.ph(), data=data, weights=data[,censorvar], method="REML", select=F)
  }else if(type=="c"){
    res.gam = gam(fml, family=stats::gaussian(), data=data, method="REML", select=F)
  }else if(type=="b"){
    res.gam = gam(fml, family=stats::binomial(link="logit"), data=data, method="REML", select=F)
  }

  ptable = summary(res.gam)$p.table
  stable = summary(res.gam)$s.table

  ptable.display = ptable
  if(!is.null(ptable.display)) ptable.display = round(ptable.display, digits=4)

  stable.display = stable
  if(!is.null(stable.display)) stable.display = round(stable.display, digits=4)

  if(prt.sum) {
    cat("GAM SUMMARY")
    print(summary(res.gam))
    cat("\n")
  }
  if(prt.chk) {
    cat("GAM CHECK")
    gam.check(res.gam)
    cat("\n")

  }

  newdat = data.frame(x=seq(min(data[,xvar]), max(data[,xvar]), length.out=1000), trt=1)
  names(newdat)= c(xvar, trtvar)
  for(i in xvars.adj){
    newdat[,i] = 0
  }
  for(i in sxvars.adj){
    newdat[,i] = stats::median(data[,i])
  }

  exclude.term = c(paste("s(",sxvars.adj,")",sep=""), paste("s(",xvar,")",sep=""))
  dat.pred = stats::predict(object=res.gam, newdata=newdat, type="link", se.fit=T, exclude=exclude.term)
  plot.data = data.frame(yfit = dat.pred$fit, yfit.Up = dat.pred$fit+1.96*dat.pred$se.fit,
                         yfit.Lo = dat.pred$fit-1.96*dat.pred$se.fit, x=newdat[,xvar])
  avg.line = ptable[trtvar,"Estimate"]


  if(type=="c"){
    ylab = "Mean Difference"
  }else if(type=="b"){
    ylab = "Log(OR)"
    if(!link.scale){
      plot.data = data.frame(yfit = exp(plot.data$yfit), yfit.Up = exp(plot.data$yfit.Up),
                             yfit.Lo = exp(plot.data$yfit.Lo), x=plot.data[,"x"])
      avg.line = exp(avg.line)
      ylab="OR"
    }
  }else if(type=="s"){
    ylab="Log(HR)"
    if(!link.scale){
      plot.data = data.frame(yfit = exp(plot.data$yfit), yfit.Up = exp(plot.data$yfit.Up),
                             yfit.Lo = exp(plot.data$yfit.Lo), x=plot.data[,"x"])
      avg.line = exp(avg.line)
      ylab="HR"
    }
  }

  if(is.null(xbreaks)) xbreaks=round(seq(min(plot.data$x), max(plot.data$x), length.out = 10),digits=2)
  if(is.null(ybreaks)) ybreaks=round(seq(min(plot.data$yfit.Lo), max(plot.data$yfit.Up), length.out = 10),digits=2)

  if(is.null(rugcol.var)){
    data.rug = data.frame(rug=data[,xvar], rugcol="Individual Record")
  }else{
    data.rug = data.frame(rug=data[,xvar], rugcol=data[,rugcol.var])
  }

  rugcol.levels = as.character(sort(unique(data.rug$rugcol), na.last = T))
  rugcol.levels[is.na(rugcol.levels)] = "N/A"
  data.rug$rugcol = as.character(data.rug$rugcol)
  data.rug$rugcol[is.na(data.rug$rugcol)] = "N/A"
  data.rug$rugcol=factor(data.rug$rugcol, levels=rugcol.levels)


  fig=ggplot(data=plot.data, aes(x=x, y=yfit)) +
    geom_ribbon(aes(ymin=yfit.Lo, ymax=yfit.Up),alpha=0.2) +
    geom_line(size=0.8)+
    geom_abline(aes(slope=0,intercept =avg.line, color="Mean"), linetype="dashed",size=0.9)+
    geom_rug(aes(x=rug, colour=rugcol),data=data.rug, alpha=0.7, sides="b", inherit.aes=F)+
    labs(color="") +
    ylab(ylab)+
    xlab(xvar)+
    scale_y_continuous(breaks=ybreaks)+
    scale_x_continuous(breaks=xbreaks)+
    ggtitle(title)+
    theme_bw(base_size=15)+
    theme(text=element_text(size=15), axis.text=element_text(size=15),legend.text=element_text(size=15))


  res = list(ptable=ptable, ptable.display=ptable.display, stable=stable, stable.display=stable.display, fig=fig)

  res


}

#' Fixed Cutoff Analysis for Individual Biomarker Associated with Binary Outcome Variables
#'
#' This function conducts fixed cutoff analysis for individual biomarker associated with binary outcome variables.
#'
#' @param yvar Binary response variable name. 0 represents controls and 1 represents cases.
#' @param xvar Biomarker name.
#' @param dir Cutoff direction for the desired subgroup. Options are ">", ">=", "<", or "<=".
#' @param cutoffs A vector of candidate cutoffs.
#' @param data The dataset containing the variables.
#' @param method Method for cutoff selection. Options are "Fisher", "Youden", "Conc.Prob", "Accuracy", or "Kappa".
#'        - "Fisher": Minimizes the Fisher test p-value.
#'        - "Youden": Maximizes the Youden index.
#'        - "Conc.Prob": Maximizes sensitivity * specificity.
#'        - "Accuracy": Maximizes accuracy.
#'        - "Kappa": Maximizes Kappa coefficient.
#' @param yvar.display Display name of the response variable.
#' @param xvar.display Display name of the predictor variable.
#' @param vert.x Whether to display the cutoff in a 90-degree angle when plotting (saves space).
#' @return A list containing statistical summaries, selected cutoff statistics, selected cutoff value, confusion matrix,
#'         and a ggplot object for visualization.
#' @import PropCIs
#' @export
fixcut_bin = function(yvar, xvar, dir, cutoffs, data, method="Fisher", yvar.display=yvar, xvar.display=xvar, vert.x=FALSE){
  Cutoff=NULL
  Value=NULL
  Measure=NULL
  Lower=NULL
  Upper=NULL
  data = data[,c(xvar, yvar)]
  data = data[stats::complete.cases(data),]
  row.names(data)=NULL
  if(!dir%in%c(">", "<", ">=", "<=")) stop("dir should be >, <, >=, or <=")
  if(!all(data[,yvar]%in%c(0,1,NA)) | !all(c(0,1)%in%data[,yvar])) stop("Response needs to be 0/1.")

  cutoffs = sort(unique(cutoffs))

  x=data[,xvar]
  y=data[,yvar]
  y=factor(y,levels=c(0,1))

  stat = NULL
  for(i in 1:length(cutoffs)){
    if(dir == "<"){
      pred = factor((x<cutoffs[i])*1,levels=c(0,1))
    }else if(dir == "<="){
      pred = factor((x<=cutoffs[i])*1,levels=c(0,1))
    }else if(dir == ">"){
      pred = factor((x>cutoffs[i])*1,levels=c(0,1))
    }else if(dir == ">="){
      pred = factor((x>=cutoffs[i])*1,levels=c(0,1))
    }

    cont = table(pred,y)
    TP = cont[2,2]
    FP = cont[2,1]
    TN = cont[1,1]
    FN = cont[1,2]
    sen = TP/(TP+FN)
    spe = TN/(TN+FP)
    ppv = TP/(TP+FP)
    npv = TN/(TN+FN)
    acc = (TP+TN)/(TP+FP+TN+FN)
    pe = mean(as.numeric(as.character(pred)))*mean(as.numeric(as.character(y)))+
      (1-mean(as.numeric(as.character(pred))))*(1-mean(as.numeric(as.character(y))))
    kap = (acc-pe)/(1-pe)


    sen.lo = scoreci(x=TP, n=TP+FN, conf.level = 0.95)$conf.int[1]
    sen.up = scoreci(x=TP, n=TP+FN, conf.level = 0.95)$conf.int[2]

    spe.lo = scoreci(x=TN, n=TN+FP, conf.level = 0.95)$conf.int[1]
    spe.up = scoreci(x=TN, n=TN+FP, conf.level = 0.95)$conf.int[2]

    ppv.lo = scoreci(x=TP, n=TP+FP, conf.level = 0.95)$conf.int[1]
    ppv.up = scoreci(x=TP, n=TP+FP, conf.level = 0.95)$conf.int[2]

    npv.lo = scoreci(x=TN, n=TN+FN, conf.level = 0.95)$conf.int[1]
    npv.up = scoreci(x=TN, n=TN+FN, conf.level = 0.95)$conf.int[2]

    acc.lo = scoreci(x=TP+TN, n=TP+FP+TN+FN, conf.level = 0.95)$conf.int[1]
    acc.up = scoreci(x=TP+TN, n=TP+FP+TN+FN, conf.level = 0.95)$conf.int[2]

    fish.p = stats::fisher.test(x=pred, y=y, or=1, alternative="two.sided")$p.value

    stat = rbind(stat, data.frame(Response = yvar.display, Predictor = xvar.display,
                                  Direction = dir,Cutoff = cutoffs[i], Sensitivity = sen, Specificity = spe,
                                  PPV = ppv, NPV = npv, Accuracy = acc, Kappa=kap, Fisher.pvalue = fish.p,
                                  Sensitivity.Lo = sen.lo, Sensitivity.Up = sen.up,
                                  Specificity.Lo = spe.lo, Specificity.Up = spe.up,
                                  PPV.Lo = ppv.lo, PPV.Up = ppv.up, NPV.Lo = npv.lo, NPV.Up = npv.up,
                                  Accuracy.Lo = acc.lo, Accuracy.Up = acc.up))

  }

  stat.display = stat
  temp = stat.display[,c("Sensitivity", "Specificity", "PPV", "NPV", "Accuracy", "Kappa", "Sensitivity.Lo", "Sensitivity.Up",
                         "Specificity.Lo", "Specificity.Up", "PPV.Lo", "PPV.Up", "NPV.Lo", "NPV.Up",
                         "Accuracy.Lo", "Accuracy.Up")]
  stat.display[,c("Sensitivity", "Specificity", "PPV", "NPV", "Accuracy", "Kappa","Sensitivity.Lo", "Sensitivity.Up",
                  "Specificity.Lo", "Specificity.Up", "PPV.Lo", "PPV.Up", "NPV.Lo", "NPV.Up",
                  "Accuracy.Lo", "Accuracy.Up")] = round(temp, digits=2)
  stat.display$Fisher.pvalue = round(stat.display$Fisher.pvalue, digits=4)

  stat.display$Sensitivity.CI = paste(stat.display$Sensitivity, " (", stat.display$Sensitivity.Lo, " , ", stat.display$Sensitivity.Up,")", sep="")
  stat.display$Specificity.CI = paste(stat.display$Specificity, " (", stat.display$Specificity.Lo, " , ", stat.display$Specificity.Up,")", sep="")
  stat.display$PPV.CI = paste(stat.display$PPV, " (", stat.display$PPV.Lo, " , ", stat.display$PPV.Up,")", sep="")
  stat.display$NPV.CI = paste(stat.display$NPV, " (", stat.display$NPV.Lo, " , ", stat.display$NPV.Up,")", sep="")
  stat.display$Accuracy.CI = paste(stat.display$Accuracy, " (", stat.display$Accuracy.Lo, " , ", stat.display$Accuracy.Up,")", sep="")
  stat.display = stat.display[,c("Response","Predictor","Direction", "Cutoff", "Sensitivity.CI",
                                 "Specificity.CI", "PPV.CI", "NPV.CI", "Accuracy.CI", "Kappa", "Fisher.pvalue")]

  if(method=="Fisher"){
    stat.sel = stat[which.min(stat$Fisher.pvalue), ]
    stat.sel.display = stat.display[which.min(stat$Fisher.pvalue),]
  }else if(method=="Youden"){
    stat.sel = stat[which.max(stat$Sensitivity+stat$Specificity), ]
    stat.sel.display = stat.display[which.max(stat$Sensitivity+stat$Specificity),]
  }else if(method=="Conc.Prob"){
    stat.sel = stat[which.max(stat$Sensitivity*stat$Specificity), ]
    stat.sel.display = stat.display[which.max(stat$Sensitivity*stat$Specificity),]
  }else if(method=="Accuracy"){
    stat.sel = stat[which.max(stat$Accuracy), ]
    stat.sel.display = stat.display[which.max(stat$Accuracy),]

  }else if(method=="Kappa"){
    stat.sel = stat[which.max(stat$Kappa), ]
    stat.sel.display = stat.display[which.max(stat$Kappa),]
  }


  row.names(stat.sel) = NULL
  row.names(stat.sel.display) = NULL
  cut.sel = stat.sel$Cutoff

  if(dir == "<"){
    pred = factor((x<cut.sel)*1,levels=c(0,1))
    dir.else = ">="
  }else if(dir == "<="){
    pred = factor((x<=cut.sel)*1,levels=c(0,1))
    dir.else = ">"
  }else if(dir == ">"){
    pred = factor((x>cut.sel)*1,levels=c(0,1))
    dir.else="<="
  }else if(dir == ">="){
    pred = factor((x>=cut.sel)*1,levels=c(0,1))
    dir.else="<"
  }

  cont = as.matrix(table(pred,y))
  colnames_cont = c(paste(yvar.display,"=", 0), paste(yvar.display,"=", 1))
  rownames_cont = c(paste(xvar.display,dir.else,cut.sel), paste(xvar.display,dir,cut.sel))
  temp.table = matrix(paste("(",round(cont/rowSums(cont), digits=3)*100, "%)",sep=""), nrow=dim(cont)[1])
  cont = matrix(paste(cont, temp.table), nrow=dim(cont)[1])
  rownames(cont) = rownames_cont
  colnames(cont) = colnames_cont

  #plot
  plot.sen = cbind(Measure="Sensitivity",stat[,c("Direction", "Cutoff", "Sensitivity", "Sensitivity.Lo", "Sensitivity.Up")])
  plot.spe = cbind(Measure="Specificity",stat[,c("Direction", "Cutoff", "Specificity", "Specificity.Lo", "Specificity.Up")])
  plot.ppv = cbind(Measure="PPV",stat[,c("Direction","Cutoff","PPV","PPV.Lo", "PPV.Up")])
  plot.npv = cbind(Measure="NPV",stat[,c("Direction","Cutoff","NPV","NPV.Lo", "NPV.Up")])
  plot.acc = cbind(Measure="Accuracy",stat[,c("Direction","Cutoff","Accuracy","Accuracy.Lo", "Accuracy.Up")])

  names(plot.sen) = c("Measure", "Direction", "Cutoff", "Value", "Lower", "Upper")
  names(plot.spe) = c("Measure", "Direction", "Cutoff", "Value", "Lower", "Upper")
  names(plot.ppv) = c("Measure", "Direction", "Cutoff", "Value", "Lower", "Upper")
  names(plot.npv) = c("Measure", "Direction", "Cutoff", "Value", "Lower", "Upper")
  names(plot.acc) = c("Measure", "Direction", "Cutoff", "Value", "Lower", "Upper")

  plot.data = rbind(plot.sen, plot.spe, plot.ppv, plot.npv, plot.acc)
  plot.data$Cutoff = factor(plot.data$Cutoff, levels=sort(unique(plot.data$Cutoff)))

  if(vert.x){
    fig = ggplot(plot.data, aes(x = Cutoff, y = Value, group=Measure, colour = Measure))+
      geom_errorbar(aes(ymin=Lower, ymax=Upper), width=0.3, size=0.8,position=position_dodge(0.2))+
      geom_line(size=0.8, position=position_dodge(0.2))+
      geom_point(size=1.8, position=position_dodge(0.2))+
      xlab(paste(xvar.display,dir,"Cutoff"))+
      ggtitle(paste(yvar.display, "~", xvar.display))+
      theme_bw(base_size=15)+
      theme(text=element_text(size=15), axis.text=element_text(size=15),legend.text=element_text(size=15),
            axis.text.x = element_text(angle = 90, hjust = 1, size=13))
  }else{
    fig = ggplot(plot.data, aes(x = Cutoff, y = Value, group=Measure, colour = Measure))+
      geom_errorbar(aes(ymin=Lower, ymax=Upper), width=0.3, size=0.8,position=position_dodge(0.2))+
      geom_line(size=0.8, position=position_dodge(0.2))+
      geom_point(size=1.8, position=position_dodge(0.2))+
      xlab(paste(xvar.display,dir,"Cutoff"))+
      ggtitle(paste(yvar.display, "~", xvar.display))+
      theme_bw(base_size=15)+
      theme(text=element_text(size=15), axis.text=element_text(size=15),legend.text=element_text(size=15))

  }

  list(stat=stat, stat.display=stat.display, stat.sel=stat.sel, stat.sel.display=stat.sel.display, cut.sel=cut.sel, confusion = cont, fig=fig)
}

#' Fixed Cutoff Analysis for Individual Biomarker Associated with Continuous Outcome
#'
#' This function conducts fixed cutoff analysis for individual biomarker associated with continuous outcome variables.
#'
#' @param yvar Continuous response variable name.
#' @param xvar Biomarker name.
#' @param dir Cutoff direction for the desired subgroup. Options are ">", ">=", "<", or "<=".
#' @param cutoffs A vector of candidate cutoffs.
#' @param data The dataset containing the variables.
#' @param method Method for cutoff selection. Currently only supports "t.test".
#'        - "t.test": Minimizes the t-test p-value.
#' @param yvar.display Display name of the response variable.
#' @param xvar.display Display name of the predictor variable.
#' @param vert.x Whether to display the cutoff in a 90-degree angle when plotting (saves space).
#' @return A list containing statistical summaries, selected cutoff statistics, selected cutoff value, group statistics,
#'         and a ggplot object for visualization.
#' @export
fixcut_con <- function(yvar, xvar, dir, cutoffs, data, method="t.test", yvar.display=yvar, xvar.display=xvar, vert.x=FALSE) {
  Cutoff <- NULL
  Mean <- NULL
  Group <- NULL
  Mean.Lo <- NULL
  Mean.Up <- NULL

  data <- data[, c(xvar, yvar)]
  data <- data[stats::complete.cases(data),]
  row.names(data) <- NULL
  if ("Group" %in% names(data)) stop("xvar and yvar should not have the name of Group.")
  if (!dir %in% c(">", "<", ">=", "<=")) stop("dir should be >, <, >=, or <=")
  cutoffs <- sort(unique(cutoffs))
  stat <- NULL
  for (i in seq_along(cutoffs)) {
    if (dir == ">") {
      data$Group <- (data[, xvar] > cutoffs[i]) * 1
    } else if (dir == ">=") {
      data$Group <- (data[, xvar] >= cutoffs[i]) * 1
    } else if (dir == "<") {
      data$Group <- (data[, xvar] < cutoffs[i]) * 1
    } else if (dir == "<=") {
      data$Group <- (data[, xvar] <= cutoffs[i]) * 1
    }

    # Check for empty groups
    if (nrow(data) == 0) {
      next
    }

    fml <- stats::as.formula(paste(yvar, "~ Group"))
    ttest <- try(stats::t.test(fml, data = data, alternative = "two.sided", mu = 0, var.equal = FALSE, conf.level = 0.95), silent = TRUE)

    if (!inherits(ttest, "try-error")) {
      pval.t <- ttest$p.value
      mean.diff.lo <- ttest$conf.int[1]
      mean.diff.up <- ttest$conf.int[2]
    } else {
      pval.t <- NA
      mean.diff.lo <- NA
      mean.diff.up <- NA
    }

    group0 <- data[data$Group == 0, ]
    group1 <- data[data$Group == 1, ]

    group0.N <- nrow(group0)
    group1.N <- nrow(group1)

    # Print group sizes for debugging
    #print(paste("Cutoff:", cutoffs[i], "Group0.N:", group0.N, "Group1.N:", group1.N))

    if (group0.N == 0 || group1.N == 0) {
      next
    }

    group0.t <- try(stats::t.test(group0[, yvar], alternative = "two.sided", mu = 0, conf.level = 0.95), silent = TRUE)
    group1.t <- try(stats::t.test(group1[, yvar], alternative = "two.sided", mu = 0, conf.level = 0.95), silent = TRUE)

    if (!inherits(group0.t, "try-error")) {
      group0.mean <- unname(group0.t$estimate)
      group0.lower <- group0.t$conf.int[1]
      group0.upper <- group0.t$conf.int[2]
      group0.median <- stats::median(group0[, yvar])
    } else {
      group0.mean <- mean(group0[, yvar], na.rm = TRUE)
      group0.lower <- NA
      group0.upper <- NA
      group0.median <- stats::median(group0[, yvar], na.rm = TRUE)
    }

    if (!inherits(group1.t, "try-error")) {
      group1.mean <- unname(group1.t$estimate)
      group1.lower <- group1.t$conf.int[1]
      group1.upper <- group1.t$conf.int[2]
      group1.median <- stats::median(group1[, yvar])
    } else {
      group1.mean <- mean(group1[, yvar], na.rm = TRUE)
      group1.lower <- NA
      group1.upper <- NA
      group1.median <- stats::median(group1[, yvar], na.rm = TRUE)
    }

    mean.diff <- group0.mean - group1.mean
    median.diff <- group0.median - group1.median

    stat <- rbind(stat, data.frame(Response = yvar.display, Predictor = xvar.display,
                                   Direction = dir, Cutoff = cutoffs[i], Mean.Diff = mean.diff,
                                   Mean.Diff.Lo = mean.diff.lo, Mean.Diff.Up = mean.diff.up,
                                   ttest.pval = pval.t, Median.Diff = median.diff,
                                   Group0.N = group0.N, Group0.Mean = group0.mean,
                                   Group0.Mean.Lo = group0.lower, Group0.Mean.Up = group0.upper,
                                   Group0.Median = group0.median,
                                   Group1.N = group1.N, Group1.Mean = group1.mean,
                                   Group1.Mean.Lo = group1.lower, Group1.Mean.Up = group1.upper,
                                   Group1.Median = group1.median))
  }

  if (is.null(stat) || nrow(stat) == 0) {
    stop("No valid data available after processing cutoffs.")
  }

  stat.display <- stat
  stat.display <- stat.display[, c("Response", "Predictor", "Direction", "Cutoff", "Group0.N", "Group1.N",
                                   "Median.Diff", "Mean.Diff", "Mean.Diff.Lo", "Mean.Diff.Up", "ttest.pval")]
  temp <- stat.display[, c("Median.Diff", "Mean.Diff", "Mean.Diff.Lo", "Mean.Diff.Up")]
  stat.display[, c("Median.Diff", "Mean.Diff", "Mean.Diff.Lo", "Mean.Diff.Up")] <- round(temp, digits = 2)
  stat.display$ttest.pval <- round(stat.display$ttest.pval, digits = 4)
  stat.display$Mean.Diff.CI <- paste(stat.display$Mean.Diff, " (", stat.display$Mean.Diff.Lo, " , ",
                                     stat.display$Mean.Diff.Up, ")", sep = "")
  stat.display <- stat.display[, c("Response", "Predictor", "Direction", "Cutoff", "Group0.N", "Group1.N",
                                   "Median.Diff", "Mean.Diff.CI", "ttest.pval")]

  if (method == "t.test") {
    stat.sel <- stat[which.min(stat$ttest.pval), ]
    stat.sel.display <- stat.display[which.min(stat$ttest.pval), ]
  }

  cut.sel <- stat.sel$Cutoff

  row.names(stat.sel) <- NULL
  row.names(stat.sel.display) <- NULL
  row.names(stat.display) <- NULL

  if (dir == "<") {
    dir.else <- ">="
  } else if (dir == "<=") {
    dir.else <- ">"
  } else if (dir == ">") {
    dir.else <- "<="
  } else if (dir == ">=") {
    dir.else <- "<"
  }

  group0.res <- stat.sel[, c("Response", "Group0.N", "Group0.Mean", "Group0.Mean.Lo", "Group0.Mean.Up", "Group0.Median")]
  group1.res <- stat.sel[, c("Response", "Group1.N", "Group1.Mean", "Group1.Mean.Lo", "Group1.Mean.Up", "Group1.Median")]

  colnames(group0.res) <- c("Response", "N", "Mean", "Mean.Lo", "Mean.Up", "Median")
  colnames(group1.res) <- c("Response", "N", "Mean", "Mean.Lo", "Mean.Up", "Median")
  group.res <- rbind(group0.res, group1.res)


  row.names(group.res) <- c(paste(xvar.display, dir.else, cut.sel), paste(xvar.display, dir, cut.sel))



  group.res.display <- group.res
  group.res.display[, c("Mean", "Mean.Lo", "Mean.Up", "Median")] <- round(group.res.display[, c("Mean", "Mean.Lo", "Mean.Up", "Median")], digits = 2)
  group.res.display$Mean.CI <- paste(group.res.display$Mean, " (", group.res.display$Mean.Lo, " , ", group.res.display$Mean.Up, ")")
  group.res.display <- group.res.display[, c("Response", "N", "Mean.CI", "Median")]

  ## plot ##
  data.plot0 <- stat[, c("Cutoff", "Group0.Mean", "Group0.Mean.Lo", "Group0.Mean.Up")]
  data.plot0$Group <- paste(dir.else, "Cutoff")
  row.names(data.plot0) <- NULL
  colnames(data.plot0) <- c("Cutoff", "Mean", "Mean.Lo", "Mean.Up", "Group")
  data.plot1 <- stat[, c("Cutoff", "Group1.Mean", "Group1.Mean.Lo", "Group1.Mean.Up")]
  data.plot1$Group <- paste(dir, "Cutoff")
  row.names(data.plot1) <- NULL
  colnames(data.plot1) <- c("Cutoff", "Mean", "Mean.Lo", "Mean.Up", "Group")

  data.plot <- rbind(data.plot0, data.plot1)

  if (nrow(data.plot) == 0) {
    stop("No data to plot.")
  }

  data.plot$Cutoff <- factor(data.plot$Cutoff, levels = sort(unique(data.plot$Cutoff)))

  fig <- ggplot(data.plot, aes(x = Cutoff, y = Mean, group = Group, colour = Group)) +
    geom_errorbar(aes(ymin = Mean.Lo, ymax = Mean.Up), width = 0.3, size = 0.8, position = position_dodge(0.2)) +
    geom_line(size = 0.8, position = position_dodge(0.2)) +
    geom_point(size = 1.8, position = position_dodge(0.2)) +
    xlab(paste("Cutoff of", xvar.display)) +
    ggtitle(paste(yvar.display, "~", xvar.display)) +
    theme_bw(base_size = 15) +
    theme(text = element_text(size = 15), axis.text = element_text(size = 15), legend.text = element_text(size = 15))

  if (vert.x) {
    fig <- fig + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 13))
  }

  list(stat = stat, stat.display = stat.display, stat.sel = stat.sel, stat.sel.display = stat.sel.display, cut.sel = cut.sel,
       group.res = group.res, group.res.display = group.res.display, fig = fig)
}



#' Fixed Cutoff Analysis for Individual Biomarker Associated with Survival Outcome
#'
#' This function conducts fixed cutoff analysis for Individual Biomarker Associated with survival outcome variables.
#'
#' @param yvar Survival response variable name.
#' @param censorvar Censoring variable. 0 indicates censored, 1 indicates an event.
#' @param xvar Biomarker name.
#' @param dir Cutoff direction for the desired subgroup.
#'        Options are ">", ">=", "<", or "<=".
#' @param cutoffs A vector of candidate cutoffs.
#' @param data The dataset containing the variables.
#' @param method Method for cutoff selection.
#'        Currently only supports "logrank".
#'        - "logrank": Minimizes the logrank test p-value.
#' @param yvar.display Display name of the response variable.
#' @param xvar.display Display name of the predictor variable.
#' @param vert.x Whether to display the cutoff in a 90-degree angle when plotting (saves space).
#' @return A list containing statistical summaries, selected cutoff statistics, selected cutoff value, group statistics,
#'         and a ggplot object for visualization.
#' @export
fixcut_sur = function(yvar, censorvar, xvar, dir, cutoffs, data, method="logrank", yvar.display=yvar, xvar.display=xvar, vert.x=FALSE){
  Cutoff=NULL
  HR.Lo=NULL
  HR.Up=NULL

  if(is.null(censorvar)) stop("Censoring variable missing!")
  data = data[,c(xvar, yvar, censorvar)]
  data = data[stats::complete.cases(data),]
  row.names(data)=NULL
  if("Group" %in% names(data)) stop("xvar, censorvar and yvar should not have the name of Group.")
  if(!dir%in%c(">", "<", ">=", "<=")) stop("dir should be >, <, >=, or <=")
  if(dir==">"){
    dir.else = "<="
  }else if(dir==">="){
    dir.else = "<"
  }else if(dir=="<"){
    dir.else = ">="
  }else if(dir=="<="){
    dir.else = ">"
  }


  cutoffs = sort(unique(cutoffs))

  stat = NULL
  for(i in 1:length(cutoffs)){
    if(dir==">"){
      data$Group = (data[,xvar]>cutoffs[i])*1
    }else if(dir==">="){
      data$Group = (data[,xvar]>=cutoffs[i])*1
    }else if(dir=="<"){
      data$Group = (data[,xvar]<cutoffs[i])*1
    }else if(dir=="<="){
      data$Group = (data[,xvar]<=cutoffs[i])*1
    }

    fml = stats::as.formula(paste("Surv(", yvar, ",", censorvar,") ~ Group"))
    res.cox = try(summary(coxph(fml, data=data)), silent=T)
    if(!inherits(res.cox,"try-error")){
      HR = res.cox$conf.int["Group", "exp(coef)"]
      HR.lower = res.cox$conf.int["Group", "lower .95"]
      HR.upper = res.cox$conf.int["Group", "upper .95"]
    }else{
      HR = as.numeric(NA)
      HR.lower = as.numeric(NA)
      HR.upper = as.numeric(NA)
    }

    data$Group = c(paste(xvar.display, dir.else, cutoffs[i]), paste(xvar.display, dir, cutoffs[i]))[data$Group+1]
    data$Group = factor(data$Group, levels=c(paste(xvar.display, dir.else, cutoffs[i]), paste(xvar.display, dir, cutoffs[i])))
    #logrank
    sdf = try(survdiff(fml, data=data, rho = 0),silent=T)
    if(!inherits(sdf,"try-error")){
      pval.logrank = 1 - stats::pchisq(sdf$chisq, length(sdf$n)-1)
    }else{
      pval.logrank = as.numeric(NA)
    }

    sfit = survfit(fml, data = data, conf.type="log-log")
    sfit.res = as.data.frame(summary(sfit, rmean="common")$table[,c("n.start","events", "rmean", "se(rmean)",
                                                                    "median","0.95LCL","0.95UCL")])
    sfit.res$rmean.lower = sfit.res[,"rmean"]-1.96*sfit.res[,"se(rmean)"]
    sfit.res$rmean.upper = sfit.res[,"rmean"]+1.96*sfit.res[,"se(rmean)"]
    sfit.res = sfit.res[, c("n.start","events", "rmean","rmean.lower", "rmean.upper",
                            "median","0.95LCL","0.95UCL")]
    names(sfit.res) = c("N", "Events", "Mean.Surv", "Mean.Surv.Lo", "Mean.Surv.Up", "Med.Surv", "Med.Surv.Lo", "Med.Surv.Up")

    sfit.group0.res = sfit.res[paste("Group=",paste(xvar.display, dir.else, cutoffs[i]),sep=""),]
    sfit.group1.res = sfit.res[paste("Group=",paste(xvar.display, dir, cutoffs[i]),sep=""),]
    names(sfit.group0.res) = paste("Group0.", names(sfit.res),sep="")
    names(sfit.group1.res) = paste("Group1.", names(sfit.res),sep="")

    stat.i = data.frame(Response = yvar.display, Predictor = xvar.display,
                        Direction = dir,Cutoff = cutoffs[i], HR =HR, HR.Lo=HR.lower, HR.Up=HR.upper, logrank.pvalue = pval.logrank)
    stat.i = cbind(stat.i, sfit.group0.res, sfit.group1.res)
    stat = rbind(stat, stat.i)

  }

  stat.display = stat
  temp.name = c("N", "Events", "Mean.Surv", "Mean.Surv.Lo", "Mean.Surv.Up", "Med.Surv", "Med.Surv.Lo", "Med.Surv.Up")
  temp = stat.display[,c("HR","HR.Lo","HR.Up",paste("Group0.", temp.name, sep=""), paste("Group1.", temp.name, sep=""))]
  stat.display[,c("HR","HR.Lo","HR.Up",paste("Group0.", temp.name, sep=""), paste("Group1.", temp.name, sep=""))] = round(temp, digits=2)
  stat.display$logrank.pvalue = round(stat.display$logrank.pvalue, digits=4)

  stat.display$HR.CI = paste(stat.display$HR, " (", stat.display$HR.Lo, " , ", stat.display$HR.Up, ")", sep="")
  stat.display$Group0.Mean.Surv.CI = paste(stat.display$Group0.Mean.Surv, " (", stat.display$Group0.Mean.Surv.Lo, " , ",
                                           stat.display$Group0.Mean.Surv.Up, ")", sep="")
  stat.display$Group0.Med.Surv.CI = paste(stat.display$Group0.Med.Surv, " (", stat.display$Group0.Med.Surv.Lo, " , ",
                                          stat.display$Group0.Med.Surv.Up, ")", sep="")
  stat.display$Group1.Mean.Surv.CI = paste(stat.display$Group1.Mean.Surv, " (", stat.display$Group1.Mean.Surv.Lo, " , ",
                                           stat.display$Group1.Mean.Surv.Up, ")", sep="")
  stat.display$Group1.Med.Surv.CI = paste(stat.display$Group1.Med.Surv, " (", stat.display$Group1.Med.Surv.Lo, " , ",
                                          stat.display$Group1.Med.Surv.Up, ")", sep="")
  stat.display = stat.display[,c("Response", "Predictor", "Direction", "Cutoff", "HR.CI",
                                 "Group0.N", "Group0.Events","Group0.Mean.Surv.CI","Group0.Med.Surv.CI",
                                 "Group1.N", "Group1.Events","Group1.Mean.Surv.CI","Group1.Med.Surv.CI",
                                 "logrank.pvalue")]

  if(method=="logrank"){
    stat.sel = stat[which.min(stat$logrank.pvalue), ]
    stat.sel.display = stat.display[which.min(stat$logrank.pval),]
  }

  row.names(stat.sel) = NULL
  row.names(stat.sel.display) = NULL
  cut.sel = stat.sel$Cutoff

  group0.res = stat.sel[,c("Response", paste("Group0.", temp.name, sep=""))]
  group1.res = stat.sel[,c("Response", paste("Group1.", temp.name, sep=""))]
  names(group0.res) = c("Response",temp.name)
  names(group1.res) = c("Response",temp.name)
  group.res = rbind(group0.res, group1.res)
  row.names(group.res) = c(paste(xvar.display, dir.else, cut.sel),paste(xvar.display, dir, cut.sel))


  group0.res.display = stat.sel.display[,c("Response","Group0.N","Group0.Events","Group0.Mean.Surv.CI","Group0.Med.Surv.CI")]
  group1.res.display = stat.sel.display[,c("Response","Group1.N","Group1.Events","Group1.Mean.Surv.CI","Group1.Med.Surv.CI")]
  names(group0.res.display) = c("Response","N","Events","Mean.Surv.CI","Med.Surv.CI")
  names(group1.res.display) = c("Response","N","Events","Mean.Surv.CI","Med.Surv.CI")
  group.res.display = rbind(group0.res.display, group1.res.display)
  row.names(group.res.display) = c(paste(xvar.display, dir.else, cut.sel),paste(xvar.display, dir, cut.sel))

  ## plot ##
  data.plot = stat[,c("Cutoff","HR","HR.Lo","HR.Up")]
  data.plot$Cutoff = factor(data.plot$Cutoff, levels=sort(unique(data.plot$Cutoff)))

  if(vert.x){
    fig = ggplot(data.plot, aes(x = Cutoff, y = HR))+
      geom_errorbar(aes(ymin=HR.Lo, ymax=HR.Up), width=0.3, size=0.8)+
      geom_line(size=0.8)+
      geom_point(size=1.8)+
      xlab(paste("Cutoff of", xvar.display))+
      ggtitle(paste(yvar.display, "~", xvar.display))+
      theme_bw(base_size=15)+
      theme(text=element_text(size=15), axis.text=element_text(size=15),legend.text=element_text(size=15),
            axis.text.x = element_text(angle = 90, hjust = 1, size=13))
  }else{
    fig = ggplot(data.plot, aes(x = Cutoff, y = HR))+
      geom_errorbar(aes(ymin=HR.Lo, ymax=HR.Up), width=0.3, size=0.8)+
      geom_line(size=0.8)+
      geom_point(size=1.8)+
      xlab(paste("Cutoff of", xvar.display))+
      ggtitle(paste(yvar.display, "~", xvar.display))+
      theme_bw(base_size=15)+
      theme(text=element_text(size=15), axis.text=element_text(size=15),legend.text=element_text(size=15))

  }

  list(stat=stat, stat.display=stat.display, stat.sel=stat.sel, stat.sel.display=stat.sel.display, cut.sel=cut.sel,
       group.res = group.res, group.res.display=group.res.display, fig=fig)
}


#' Cutoff Performance Evaluation
#'
#' This function evaluates the performance of a predictive model at a selected cutoff point.
#'
#' @param yvar Response variable name.
#' @param censorvar Censoring variable name (0-censored, 1-event).
#' @param xvar Biomarker name.
#' @param cutoff Selected cutoff value.
#' @param dir Direction for desired subgroup (">", ">=", "<", "<=").
#' @param xvars.adj Other covariates to adjust when evaluating the performance.
#' @param data Data frame containing the variables.
#' @param type Type of analysis: "c" for continuous, "s" for survival, and "b" for binary.
#' @param yvar.display Display name of response variable.
#' @param xvar.display Display name of biomarker variable.
#' @import survival
#' @import survminer
#' @return A list containing various performance metrics and optionally, plots.
#' @import survminer
#' @importFrom stats confint
#' @importFrom stats glm
#' @export
cut_perf = function(yvar, censorvar=NULL, xvar, cutoff, dir, xvars.adj=NULL, data, type, yvar.display=yvar, xvar.display=xvar){
  Group=NULL

  if(is.null(xvars.adj)){
    data = data[, c(yvar,censorvar,xvar)]
    data = data[stats::complete.cases(data),]
    if(!dir%in%c(">", "<", ">=", "<=")) stop("dir should be >, <, >=, or <=")
    if("Group" %in% names(data)) stop("Group cannot be one of the column names of data.")

    if(dir==">"){
      data$Group = (data[,xvar]>cutoff)*1
      dir.else = "<="
    }else if(dir==">="){
      data$Group = (data[,xvar]>=cutoff)*1
      dir.else = "<"
    }else if(dir=="<"){
      data$Group = (data[,xvar]<cutoff)*1
      dir.else = ">="
    }else if(dir=="<="){
      data$Group = (data[,xvar]<=cutoff)*1
      dir.else = ">"
    }

    if(type=="s"){
      if(is.null(censorvar)) stop("Censoring variable missing!")

      #coxph
      fml = stats::as.formula(paste("Surv(", yvar, ",", censorvar,") ~ Group"))
      res.cox = try(summary(coxph(fml, data=data)), silent=T)
      if(!inherits(res.cox,"try-error")){
        HR = res.cox$conf.int["Group", "exp(coef)"]
        HR.lower = res.cox$conf.int["Group", "lower .95"]
        HR.upper = res.cox$conf.int["Group", "upper .95"]
      }else{
        HR = as.numeric(NA)
        HR.lower = as.numeric(NA)
        HR.upper = as.numeric(NA)
      }


      data$Group = c(paste(xvar.display, dir.else, cutoff), paste(xvar.display, dir, cutoff))[data$Group+1]
      data$Group = factor(data$Group, levels=c(paste(xvar.display, dir.else, cutoff), paste(xvar.display, dir, cutoff)))
      #logrank
      sdf = try(survdiff(fml, data=data, rho = 0),silent=T)
      if(!inherits(sdf,"try-error")){
        pval.logrank = 1 - stats::pchisq(sdf$chisq, length(sdf$n)-1)
      }else{
        pval.logrank = as.numeric(NA)
      }

      comp.res = data.frame(HR=HR, HR.lower=HR.lower, HR.upper=HR.upper, logrank.pvalue = pval.logrank)
      comp.res.display = comp.res
      comp.res.display[,c("HR","HR.lower","HR.upper")] = round(comp.res.display[,c("HR","HR.lower","HR.upper")], digits=2)
      comp.res.display$logrank.pvalue = round(comp.res.display$logrank.pvalue, digits=4)

      comp.res.display$HR.CI = paste(comp.res.display$HR, " (", comp.res.display$HR.lower, " , ",
                                     comp.res.display$HR.upper, ")", sep="")
      comp.res.display = comp.res.display[,c("HR.CI", "logrank.pvalue")]
      header = data.frame(Response=yvar.display, Predictor=xvar.display, Direction=dir, Cutoff=cutoff)
      comp.res = cbind(header, comp.res)
      comp.res.display = cbind(header, comp.res.display)


      #median/mean survival time
      sfit = survfit(fml, data = data, conf.type="log-log")
      sfit.res = as.data.frame(summary(sfit, rmean="common")$table[,c("n.start","events", "rmean", "se(rmean)",
                                                                      "median","0.95LCL","0.95UCL")])
      sfit.res$rmean.lower = sfit.res[,"rmean"]-1.96*sfit.res[,"se(rmean)"]
      sfit.res$rmean.upper = sfit.res[,"rmean"]+1.96*sfit.res[,"se(rmean)"]
      sfit.res = sfit.res[, c("n.start","events", "rmean","rmean.lower", "rmean.upper",
                              "median","0.95LCL","0.95UCL")]
      names(sfit.res) = c("N", "Events", "Mean.Surv", "Mean.Surv.Lo", "Mean.Surv.Up", "Med.Surv", "Med.Surv.Lo", "Med.Surv.Up")


      sfit.res.display = sfit.res
      sfit.res.display = round(sfit.res.display, digits=2)
      sfit.res.display$Mean.Surv.CI = paste(sfit.res.display$Mean.Surv, " (", sfit.res.display$Mean.Surv.Lo, " , ",
                                            sfit.res.display$Mean.Surv.Up, ")", sep="")
      sfit.res.display$Med.Surv.CI = paste(sfit.res.display$Med.Surv, " (", sfit.res.display$Med.Surv.Lo, " , ",
                                           sfit.res.display$Med.Surv.Up, ")", sep="")
      sfit.res.display = sfit.res.display[,c("N", "Events", "Mean.Surv.CI", "Med.Surv.CI")]

      sfit.res = cbind(data.frame(Response=yvar.display), sfit.res)
      sfit.res.display = cbind(data.frame(Response=yvar.display), sfit.res.display)

      #plot
      fig = ggsurvplot(surv_fit(fml, data = data, conf.type="log-log"), data = data, risk.table = TRUE,
                       surv.median.line="none", pval=F, pval.method=F,  xlab=yvar.display, legend.title=xvar.display,
                       legend.labs = c(paste(dir.else, cutoff), paste(dir, cutoff)),
                       font.legend= 14, title=paste(yvar.display, "~", xvar.display))


      res = list(comp.res=comp.res, comp.res.display=comp.res.display, group.res=sfit.res, group.res.display=sfit.res.display, fig=fig)


    }else if(type=="c"){
      fml = stats::as.formula(paste(yvar,"~ Group"))
      wrstest = try(stats::wilcox.test(fml, data=data, alternative = "two.sided", mu=0, paired=FALSE), silent=T)
      ttest = try(stats::t.test(fml, data=data, alternative = "two.sided", mu=0, paired = FALSE, var.equal = FALSE,
                         conf.level = 0.95),silent=T)
      if(!inherits(ttest,"try-error")){
        pval.t = ttest$p.value
        mean.diff.lo = ttest$conf.int[1]
        mean.diff.up = ttest$conf.int[2]
      }else{
        pval.t = as.numeric(NA)
        mean.diff.lo = as.numeric(NA)
        mean.diff.up = as.numeric(NA)
      }

      if(!inherits(wrstest,"try-error")){
        pval.wrs = wrstest$p.value
      }else{
        pval.wrs = as.numeric(NA)
      }


      group0 = data[data$Group%in%0,]
      group1 = data[data$Group%in%1,]

      group0.N = dim(group0)[1]
      group1.N = dim(group1)[1]

      group0.t = try(stats::t.test(group0[,yvar], alternative = "two.sided", mu = 0, conf.level = 0.95), silent=T)
      group1.t = try(stats::t.test(group1[,yvar], alternative = "two.sided", mu = 0, conf.level = 0.95), silent=T)

      if(!inherits(group0.t,"try-error")){
        group0.mean = unname(group0.t$estimate)
        group0.lower = group0.t$conf.int[1]
        group0.upper = group0.t$conf.int[2]
        group0.median = stats::median(group0[,yvar])
      }else{
        group0.mean = mean(group0[,yvar])
        group0.lower = as.numeric(NA)
        group0.upper = as.numeric(NA)
        group0.median = stats::median(group0[,yvar])
      }
      group0.res = data.frame(Response=yvar.display, N=group0.N, Mean=group0.mean, Mean.Lo=group0.lower, Mean.Up=group0.upper, Median=group0.median)

      if(!inherits(group1.t,"try-error")){
        group1.mean = unname(group1.t$estimate)
        group1.lower = group1.t$conf.int[1]
        group1.upper = group1.t$conf.int[2]
        group1.median = stats::median(group1[,yvar])
      }else{
        group1.mean = mean(group1[,yvar])
        group1.lower = as.numeric(NA)
        group1.upper = as.numeric(NA)
        group1.median = stats::median(group1[,yvar])
      }
      group1.res = data.frame(Response=yvar.display, N=group1.N, Mean=group1.mean, Mean.Lo=group1.lower, Mean.Up=group1.upper, Median=group1.median)

      mean.diff = group0.mean - group1.mean
      median.diff = group0.median - group1.median

      comp.res = data.frame(Response=yvar.display, Predictor=xvar.display, Direction=dir, Cutoff=cutoff, Mean.Diff = mean.diff,
                            Mean.Diff.Lo = mean.diff.lo, Mean.Diff.Up = mean.diff.up, t.pvalue = pval.t, Median.Diff = median.diff,
                            WRS.pvalue = pval.wrs)
      comp.res.display = comp.res
      temp = comp.res.display[,c("Median.Diff", "Mean.Diff", "Mean.Diff.Lo", "Mean.Diff.Up")]
      comp.res.display[,c("Median.Diff", "Mean.Diff", "Mean.Diff.Lo", "Mean.Diff.Up")] = round(temp, digits = 2)
      comp.res.display[,c("t.pvalue","WRS.pvalue")] = round(comp.res.display[,c("t.pvalue","WRS.pvalue")], digits=4)
      comp.res.display$Mean.Diff.CI = paste(comp.res.display$Mean.Diff, " (", comp.res.display$Mean.Diff.Lo, " , ",
                                            comp.res.display$Mean.Diff.Up, ")", sep="")

      comp.res.display = comp.res.display[,c("Response", "Predictor", "Direction", "Cutoff",
                                             "Median.Diff", "Mean.Diff.CI","t.pvalue", "WRS.pvalue" )]

      group.res = rbind(group0.res, group1.res)
      rownames(group.res) = c(paste(xvar.display, dir.else, cutoff),paste(xvar.display, dir, cutoff))

      group.res.display = group.res
      group.res.display[,c("Mean", "Mean.Lo", "Mean.Up", "Median")] = round(group.res.display[,c("Mean", "Mean.Lo", "Mean.Up", "Median")], digits=2)
      group.res.display$Mean.CI = paste(group.res.display$Mean, " (", group.res.display$Mean.Lo, " , ", group.res.display$Mean.Up, ")")
      group.res.display = group.res.display[,c("Response","N", "Mean.CI", "Median")]

      #plot
      data.plot = data[,c(yvar,"Group")]
      data.plot$Group = c(paste(xvar.display, dir.else, cutoff), paste(xvar.display, dir, cutoff))[data.plot$Group+1]
      names(data.plot) = c("y", "Group")
      fig = ggplot(data.plot, aes(x=Group, y=y)) + geom_boxplot(outlier.shape=NA) +
        stat_summary(fun=mean, geom="point", shape=5, size=4) +
        geom_point(size=2.5, position=position_jitter(width=0.18,height=0))+
        ylab(yvar.display)+
        theme_bw(base_size=15)+
        theme(text=element_text(size=15), axis.text=element_text(size=15),legend.text=element_text(size=15))

      res = list(comp.res=comp.res, comp.res.display=comp.res.display, group.res=group.res, group.res.display=group.res.display, fig=fig)



    }else if(type=="b"){
      if(!all(data[,yvar]%in%c(0,1,NA)) | !all(c(0,1)%in%data[,yvar])) stop("Response needs to be 0/1.")

      x=data[,"Group"]
      y=data[,yvar]
      x=factor(x, levels=c(0,1))
      y=factor(y, levels=c(0,1))

      cont = table(x,y)
      TP = cont[2,2]
      FP = cont[2,1]
      TN = cont[1,1]
      FN = cont[1,2]
      sen = TP/(TP+FN)
      spe = TN/(TN+FP)
      ppv = TP/(TP+FP)
      npv = TN/(TN+FN)
      acc = (TP+TN)/(TP+FP+TN+FN)

      sen.lo = scoreci(x=TP, n=TP+FN, conf.level = 0.95)$conf.int[1]
      sen.up = scoreci(x=TP, n=TP+FN, conf.level = 0.95)$conf.int[2]

      spe.lo = scoreci(x=TN, n=TN+FP, conf.level = 0.95)$conf.int[1]
      spe.up = scoreci(x=TN, n=TN+FP, conf.level = 0.95)$conf.int[2]

      ppv.lo = scoreci(x=TP, n=TP+FP, conf.level = 0.95)$conf.int[1]
      ppv.up = scoreci(x=TP, n=TP+FP, conf.level = 0.95)$conf.int[2]

      npv.lo = scoreci(x=TN, n=TN+FN, conf.level = 0.95)$conf.int[1]
      npv.up = scoreci(x=TN, n=TN+FN, conf.level = 0.95)$conf.int[2]

      acc.lo = scoreci(x=TP+TN, n=TP+FP+TN+FN, conf.level = 0.95)$conf.int[1]
      acc.up = scoreci(x=TP+TN, n=TP+FP+TN+FN, conf.level = 0.95)$conf.int[2]

      fish.p = stats::fisher.test(x=x, y=y, or=1, alternative="two.sided")$p.value

      stat = data.frame(Response = yvar.display, Predictor = xvar.display,
                        Direction = dir,Cutoff = cutoff, N.Sel = TP+FP, N.unSel=TN+FN,
                        Sensitivity = sen, Specificity = spe,
                        PPV = ppv, NPV = npv, Accuracy = acc, Fisher.pvalue = fish.p,
                        Sensitivity.Lo = sen.lo, Sensitivity.Up = sen.up,
                        Specificity.Lo = spe.lo, Specificity.Up = spe.up,
                        PPV.Lo = ppv.lo, PPV.Up = ppv.up, NPV.Lo = npv.lo, NPV.Up = npv.up,
                        Accuracy.Lo = acc.lo, Accuracy.Up = acc.up)


      stat.display = stat
      temp = stat.display[,c("Sensitivity", "Specificity", "PPV", "NPV", "Accuracy", "Sensitivity.Lo", "Sensitivity.Up",
                             "Specificity.Lo", "Specificity.Up", "PPV.Lo", "PPV.Up", "NPV.Lo", "NPV.Up",
                             "Accuracy.Lo", "Accuracy.Up")]
      stat.display[,c("Sensitivity", "Specificity", "PPV", "NPV", "Accuracy", "Sensitivity.Lo", "Sensitivity.Up",
                      "Specificity.Lo", "Specificity.Up", "PPV.Lo", "PPV.Up", "NPV.Lo", "NPV.Up",
                      "Accuracy.Lo", "Accuracy.Up")] = round(temp, digits=2)
      stat.display$Fisher.pvalue = round(stat.display$Fisher.pvalue, digits=4)

      stat.display$Sensitivity.CI = paste(stat.display$Sensitivity, " (", stat.display$Sensitivity.Lo, " , ", stat.display$Sensitivity.Up,")", sep="")
      stat.display$Specificity.CI = paste(stat.display$Specificity, " (", stat.display$Specificity.Lo, " , ", stat.display$Specificity.Up,")", sep="")
      stat.display$PPV.CI = paste(stat.display$PPV, " (", stat.display$PPV.Lo, " , ", stat.display$PPV.Up,")", sep="")
      stat.display$NPV.CI = paste(stat.display$NPV, " (", stat.display$NPV.Lo, " , ", stat.display$NPV.Up,")", sep="")
      stat.display$Accuracy.CI = paste(stat.display$Accuracy, " (", stat.display$Accuracy.Lo, " , ", stat.display$Accuracy.Up,")", sep="")
      stat.display = stat.display[,c("Response","Predictor","Direction", "Cutoff", "Sensitivity.CI",
                                     "Specificity.CI", "PPV.CI", "NPV.CI", "Accuracy.CI", "Fisher.pvalue")]
      #rownames(stat) = NULL
      #rownames(stat.display) = NULL


      cont = as.matrix(cont)
      colnames_cont = c(paste(yvar.display,"=", 0), paste(yvar.display,"=", 1))
      rownames_cont = c(paste(xvar.display,dir.else,cutoff), paste(xvar.display,dir,cutoff))
      temp.table = matrix(paste("(",round(cont/rowSums(cont), digits=3)*100, "%)",sep=""), nrow=dim(cont)[1])
      cont = matrix(paste(cont, temp.table), nrow=dim(cont)[1])
      rownames(cont) = rownames_cont
      colnames(cont) = colnames_cont

      res = list(stat=stat, stat.display=stat.display, confusion=cont)


    }else{
      stop("type can only be c, b or s.")
    }

  }else{
    data = data[, c(yvar,censorvar,xvar, xvars.adj)]
    data = data[stats::complete.cases(data),]

    if(!dir%in%c(">", "<", ">=", "<=")) stop("dir should be >, <, >=, or <=")
    if("Group" %in% names(data)) stop("Group cannot be one of the column names of data.")

    if(dir==">"){
      data$Group = (data[,xvar]>cutoff)*1
      dir.else = "<="
    }else if(dir==">="){
      data$Group = (data[,xvar]>=cutoff)*1
      dir.else = "<"
    }else if(dir=="<"){
      data$Group = (data[,xvar]<cutoff)*1
      dir.else = ">="
    }else if(dir=="<="){
      data$Group = (data[,xvar]<=cutoff)*1
      dir.else = ">"
    }

    if(type=="s"){
      if(is.null(censorvar)) stop("Censoring variable missing!")

      #coxph
      xvars.adj.fml = paste(xvars.adj,collapse="+")
      fml = stats::as.formula(paste("Surv(", yvar, ",", censorvar,") ~ Group +", xvars.adj.fml))
      res.cox = try(summary(coxph(fml, data=data)), silent=T)
      if(!inherits(res.cox,"try-error")){
        HR = res.cox$conf.int["Group", "exp(coef)"]
        HR.lower = res.cox$conf.int["Group", "lower .95"]
        HR.upper = res.cox$conf.int["Group", "upper .95"]
        pval = res.cox$coefficients["Group", "Pr(>|z|)"]
      }else{
        HR = as.numeric(NA)
        HR.lower = as.numeric(NA)
        HR.upper = as.numeric(NA)
        pval = as.numeric(NA)
      }

      comp.res = data.frame(HR=HR, HR.lower=HR.lower, HR.upper=HR.upper, cox.pvalue = pval)
      comp.res.display = comp.res
      comp.res.display[,c("HR","HR.lower","HR.upper")] = round(comp.res.display[,c("HR","HR.lower","HR.upper")], digits=2)
      comp.res.display$cox.pvalue = round(comp.res.display$cox.pvalue, digits=4)
      comp.res.display$HR.CI = paste(comp.res.display$HR, " (", comp.res.display$HR.lower, " , ",
                                     comp.res.display$HR.upper, ")", sep="")
      comp.res.display = comp.res.display[,c("HR.CI", "cox.pvalue")]


      header = data.frame(Response=yvar.display, Predictor=xvar.display, Direction=dir, Cutoff=cutoff)
      comp.res = cbind(header, comp.res)
      comp.res.display = cbind(header, comp.res.display)

      res = list(comp.res=comp.res, comp.res.display=comp.res.display)


    }else if(type=="c"){
      xvars.adj.fml = paste(xvars.adj,collapse="+")
      fml = stats::as.formula(paste(yvar,"~ Group +",xvars.adj.fml))
      lm.fit = try(stats::lm(fml, data=data), silent=T)
      if(!inherits(lm.fit,"try-error")){
        mean.diff = summary(lm.fit)$coefficients["Group","Estimate"]
        mean.diff.lower = confint(lm.fit)["Group","2.5 %"]
        mean.diff.upper = confint(lm.fit)["Group","97.5 %"]
        pval = summary(lm.fit)$coefficients["Group","Pr(>|t|)"]
      }else{
        mean.diff = as.numeric(NA)
        mean.diff.lower = as.numeric(NA)
        mean.diff.upper = as.numeric(NA)
        pval = as.numeric(NA)
      }

      comp.res = data.frame(Mean.Diff = mean.diff, Mean.Diff.Lo = mean.diff.lower, Mean.Diff.Up=mean.diff.upper,
                            lm.pvalue = pval)
      comp.res.display = comp.res
      temp = comp.res.display[,c("Mean.Diff", "Mean.Diff.Lo", "Mean.Diff.Up")]
      comp.res.display[,c("Mean.Diff", "Mean.Diff.Lo", "Mean.Diff.Up")] = round(temp, digits = 2)
      comp.res.display$lm.pvalue = round(comp.res.display$lm.pvalue, digits=4)
      comp.res.display$Mean.Diff.CI = paste(comp.res.display$Mean.Diff, " (", comp.res.display$Mean.Diff.Lo, " , ",
                                            comp.res.display$Mean.Diff.Up, ")", sep="")
      comp.res.display = comp.res.display[,c("Mean.Diff.CI","lm.pvalue" )]


      header = data.frame(Response=yvar.display, Predictor=xvar.display, Direction=dir, Cutoff=cutoff)
      comp.res = cbind(header, comp.res)
      comp.res.display = cbind(header, comp.res.display)

      res = list(comp.res=comp.res, comp.res.display=comp.res.display)

    }else if(type=="b"){

      if(!all(data[,yvar]%in%c(0,1,NA)) | !all(c(0,1)%in%data[,yvar])) stop("Response needs to be 0/1.")

      xvars.adj.fml = paste(xvars.adj,collapse="+")
      fml = stats::as.formula(paste(yvar,"~ Group +",xvars.adj.fml))
      glm.fit = try(glm(fml, family=stats::binomial, data=data), silent=T)
      if(!inherits(glm.fit,"try-error")){
        log.OR = summary(glm.fit)$coefficients["Group", "Estimate"]
        log.OR.lower = confint(glm.fit)["Group","2.5 %"]
        log.OR.upper = confint(glm.fit)["Group","97.5 %"]
        pval = summary(glm.fit)$coefficients["Group", "Pr(>|z|)"]
      }else{
        log.OR = as.numeric(NA)
        log.OR.lower = as.numeric(NA)
        log.OR.upper = as.numeric(NA)
        pval = as.numeric(NA)

      }

      comp.res = data.frame(log.OR = log.OR, log.OR.Lo = log.OR.lower, log.OR.Up = log.OR.upper,
                            glm.pvalue = pval)
      comp.res.display = comp.res
      comp.res.display[,c("log.OR", "log.OR.Lo", "log.OR.Up")] = round(comp.res.display[,c("log.OR", "log.OR.Lo", "log.OR.Up")], digits=2)
      comp.res.display$glm.pvalue= round(comp.res.display$glm.pvalue, digits=4)
      comp.res.display$log.OR.CI = paste(comp.res.display$log.OR, " (", comp.res.display$log.OR.Lo, " , ",
                                         comp.res.display$log.OR.Up, ")", sep="")
      comp.res.display = comp.res.display[,c("log.OR.CI","glm.pvalue" )]

      header = data.frame(Response=yvar.display, Predictor=xvar.display, Direction=dir, Cutoff=cutoff)
      comp.res = cbind(header, comp.res)
      comp.res.display = cbind(header, comp.res.display)

      res = list(comp.res=comp.res, comp.res.display=comp.res.display)

    }else{
      stop("type can only be c, b or s.")
    }


  }

  res

}


#' Summarize Categorical Variables in Subgroup
#'
#' @description This function provides a summary of categorical variables in a dataset.
#'
#' @param yvar Name of the variable for summary.
#' @param yname A vector of ordered y values.
#' @param xvars Names of the variables for grouping.
#' @param xname.list A list (same order as xvars) of ordered x values for each xvar.
#' @param data The dataset.
#' @param yvar.display Display name for yvar.
#' @param xvars.display Display name for xvars.
#' @return A list containing the contingency table, frequency table, and percentage table.
#' @export
cat_summary = function(yvar, yname, xvars, xname.list, data, yvar.display=yvar, xvars.display=xvars){
  #yvar: name of the variable for summary
  #yname: a vector of ordered y values
  #xvars: names of the variables for grouping
  #xname.list: a list (same order as xvars) of ordered x values for each xvar
  #yvar.display: display name for yvar
  #xvars.display: displayname for xvars

  if(!is.list(xname.list)) stop("xname.list has to be a list.")
  data = data[,c(yvar, xvars)]
  data = data[stats::complete.cases(data),]
  yname = unique(yname)
  if(!all(yname%in%data[,yvar]) | !all(data[,yvar]%in%yname))
    warning("yname should be the values for the yvar variable.")
  data[,yvar] = factor(data[,yvar], levels=yname)

  if(length(xvars)!=length(xname.list)) stop("xvars and xname.list do not have the same length.")
  for(i in 1:length(xvars)){
    xname.list[[i]]=unique(xname.list[[i]])
    if(!all(xname.list[[i]]%in%data[,xvars[i]]) | !all(data[,xvars[i]]%in%xname.list[[i]]))
      warning("xname should be the values for the xvar variable.")
    data[,xvars[i]] = factor(data[,xvars[i]], levels=xname.list[[i]])
  }

  data = data[,c(xvars,yvar)]

  freq_table = table(data)
  perc_table = round(prop.table(freq_table,margin=1:length(xvars)), digits=3)*100
  freq_table = stats::ftable(freq_table)
  perc_table = stats::ftable(perc_table)

  perc_table.display = matrix(paste("(",perc_table,"%)", sep=""), nrow=nrow(perc_table))
  cont.display = matrix(paste(freq_table, perc_table.display), nrow=nrow(freq_table))
  cont.display = as.data.frame(cont.display, stringsAsFactors = F)
  names(cont.display) = paste(yvar.display, "=", yname)

  header = expand.grid(xname.list[length(xname.list):1], stringsAsFactors=F)[length(xname.list):1]
  names(header) = xvars.display

  cont.display = cbind(header, cont.display)
  list(cont.display=cont.display, freq_table=freq_table, perc_table=perc_table)

}



#' Subgroup Performance Evaluation for Predictive Cases
#'
#' @description This function evaluates the performance of subgroups based on different types of response variables in predictive cases.
#'
#' @param yvar Response variable name.
#' @param censorvar Censoring variable name (0-censored, 1-event).
#' @param grpvar Subgroup variable name.
#' @param grpname A vector of ordered subgroup names (values in the column of grpvar).
#' @param trtvar Treatment variable name.
#' @param trtname A vector of ordered treatment names (values in the column of trtvar).
#' @param xvars.adj Other covariates to adjust when evaluating the performance.
#' @param data The dataset.
#' @param type "c" for continuous; "s" for "survival", and "b" for binary.
#' @param yvar.display Display name of the response variable.
#' @param grpvar.display Display name of the group variable.
#' @param trtvar.display Display name of the treatment variable.
#' @return A list containing the comparison results, group results, and possibly a plot.
#' @importFrom car Anova
#' @export
subgrp_perf_pred = function(yvar, censorvar=NULL, grpvar, grpname, trtvar, trtname, xvars.adj=NULL, data,
                            type, yvar.display=yvar, grpvar.display=grpvar, trtvar.display=trtvar){
  Group=NULL
  y=NULL
  Trt=NULL
  grp=NULL
  Rate=NULL
  trt=NULL
  Rate.Lo=NULL
  Rate.Up=NULL
  data = data[, c(yvar,censorvar,grpvar,trtvar,xvars.adj)]
  data = data[stats::complete.cases(data),]

  grpname = unique(grpname)
  if(!all(grpname%in%data[,grpvar]) | !all(data[,grpvar]%in%grpname))
    stop("grpname should be the values for the grpvar variable.")

  trtname = unique(trtname)
  if(!all(trtname%in%data[,trtvar]) | !all(data[,trtvar]%in%trtname))
    stop("trtname should be the values for the trtvar variable.")

  data[,grpvar] = factor(data[,grpvar], levels=grpname)
  data[,trtvar] = factor(data[,trtvar], levels=trtname)

  unicomb1 = paste(rep(trtvar.display, length(trtname)*length(grpname)),"=",
                   rep(trtname, each=length(grpname)), ",",
                   rep(grpvar.display, length(trtname)*length(grpname)),"=",
                   rep(grpname, length(trtname)), sep="")
  datacomb1 = paste(trtvar.display,"=",data[,trtvar], ",",
                    grpvar.display,"=", data[,grpvar], sep="")
  unicomb1 = unicomb1[unicomb1%in%datacomb1]

  unicomb2 = paste(rep(trtname, each=length(grpname)), ",",
                   rep(grpname, length(trtname)), sep="")
  datacomb2 = paste(data[,trtvar], ",", data[,grpvar], sep="")
  unicomb2 = unicomb2[unicomb2%in%datacomb2]


  if(type=="s"){
    if(is.null(censorvar)) stop("Censoring variable missing!")

    #coxph
    if(is.null(xvars.adj)){
      xvars.adj.fml=NULL
    }else{
      xvars.adj.fml = paste(c("",xvars.adj),collapse="+")
    }

    fml = stats::as.formula(paste("Surv(", yvar, ",", censorvar,") ~ ", trtvar, "*",grpvar, xvars.adj.fml, sep=""))
    res.cox = try(Anova(coxph(fml, data=data), type=3, test.statistic="Wald"), silent=T)
    if(!inherits(res.cox,"try-error")){
      pval.pred = res.cox[paste(trtvar,":",grpvar,sep=""),"Pr(>Chisq)"]
      pval.prog = res.cox[grpvar, "Pr(>Chisq)"]

    }else{
      pval.pred = as.numeric(NA)
      pval.prog = as.numeric(NA)
    }

    comp.res = data.frame(pvalue.pred = pval.pred, pvalue.prog=pval.prog)
    comp.res.display = comp.res
    comp.res.display= round(comp.res.display, digits=4)

    header = data.frame(Response=yvar.display, Group.Variable=grpvar.display, Test="Cox ANOVA")
    comp.res = cbind(header, comp.res)
    comp.res.display = cbind(header, comp.res.display)


    #median/mean survival time
    fml = stats::as.formula(paste("Surv(", yvar, ",", censorvar,") ~ ", trtvar, "+",grpvar, sep=""))
    sfit = survfit(fml, data = data, conf.type="log-log")
    sfit.res = as.data.frame(summary(sfit, rmean="common")$table[,c("n.start","events", "rmean", "se(rmean)",
                                                                    "median","0.95LCL","0.95UCL")])
    sfit.res$rmean.lower = sfit.res[,"rmean"]-1.96*sfit.res[,"se(rmean)"]
    sfit.res$rmean.upper = sfit.res[,"rmean"]+1.96*sfit.res[,"se(rmean)"]
    sfit.res = sfit.res[, c("n.start","events", "rmean","rmean.lower", "rmean.upper",
                            "median","0.95LCL","0.95UCL")]
    names(sfit.res) = c("N", "Events", "Mean.Surv", "Mean.Surv.Lo", "Mean.Surv.Up", "Med.Surv", "Med.Surv.Lo", "Med.Surv.Up")
    row.names(sfit.res) = unicomb1

    sfit.res.display = sfit.res
    sfit.res.display = round(sfit.res.display, digits=2)
    sfit.res.display$Mean.Surv.CI = paste(sfit.res.display$Mean.Surv, " (", sfit.res.display$Mean.Surv.Lo, " , ",
                                          sfit.res.display$Mean.Surv.Up, ")", sep="")
    sfit.res.display$Med.Surv.CI = paste(sfit.res.display$Med.Surv, " (", sfit.res.display$Med.Surv.Lo, " , ",
                                         sfit.res.display$Med.Surv.Up, ")", sep="")
    sfit.res.display = sfit.res.display[,c("N", "Events", "Mean.Surv.CI", "Med.Surv.CI")]

    sfit.res = cbind(data.frame(Response=yvar.display), sfit.res)
    sfit.res.display = cbind(data.frame(Response=yvar.display), sfit.res.display)

    #plot
    fig = ggsurvplot(surv_fit(fml, data = data, conf.type="log-log"), data = data, risk.table = TRUE,
                     surv.median.line="none", pval=F, pval.method=F,  xlab=yvar.display,
                     legend.title=paste(trtvar.display,",", grpvar.display),
                     legend.labs = unicomb2, font.legend= 14,
                     title=paste("KM Curves for", yvar.display))


    res = list(comp.res=comp.res, comp.res.display=comp.res.display,
               group.res=sfit.res, group.res.display=sfit.res.display, fig=fig)

  }else if(type=="c"){
    if(is.null(xvars.adj)){
      xvars.adj.fml=NULL
    }else{
      xvars.adj.fml = paste(c("",xvars.adj),collapse="+")
    }

    fml = stats::as.formula(paste(yvar, " ~ ", trtvar, "*",grpvar, xvars.adj.fml, sep=""))
    lm.fit = try(Anova(stats::lm(fml, data=data), type=3), silent=T)
    if(!inherits(lm.fit,"try-error")){
      pval.pred = lm.fit[paste(trtvar,":",grpvar,sep=""),"Pr(>F)"]
      pval.prog = lm.fit[grpvar, "Pr(>F)"]

    }else{
      pval.pred = as.numeric(NA)
      pval.prog = as.numeric(NA)
    }

    comp.res = data.frame(pvalue.pred = pval.pred, pvalue.prog=pval.prog)
    comp.res.display = comp.res
    comp.res.display= round(comp.res.display, digits=4)

    header = data.frame(Response=yvar.display, Group.Variable=grpvar.display, Test="LM ANOVA")
    comp.res = cbind(header, comp.res)
    comp.res.display = cbind(header, comp.res.display)


    ## group res
    group.res = NULL
    for(unicomb2.i in unicomb2){
      datai = data[datacomb2%in%unicomb2.i,]
      datai.N = dim(datai)[1]
      datai.t = try(stats::t.test(datai[,yvar], alternative = "two.sided", mu = 0, conf.level = 0.95), silent=T)

      if(!inherits(datai.t,"try-error")){
        datai.mean = unname(datai.t$estimate)
        datai.lower = datai.t$conf.int[1]
        datai.upper = datai.t$conf.int[2]
        datai.median = stats::median(datai[,yvar])
      }else{
        datai.mean = mean(datai[,yvar])
        datai.lower = as.numeric(NA)
        datai.upper = as.numeric(NA)
        datai.median = stats::median(datai[,yvar])
      }

      datai.res = data.frame(Response=yvar.display, N=datai.N, Mean=datai.mean, Mean.Lo=datai.lower,
                             Mean.Up=datai.upper, Median=datai.median)
      row.names(datai.res) = paste(trtvar.display,",",grpvar.display," = ",unicomb2.i, sep="")

      group.res = rbind(group.res, datai.res)
    }

    group.res.display = group.res
    group.res.display[,c("Mean", "Mean.Lo", "Mean.Up", "Median")] = round(group.res.display[,c("Mean", "Mean.Lo", "Mean.Up", "Median")], digits=2)
    group.res.display$Mean.CI = paste(group.res.display$Mean, " (", group.res.display$Mean.Lo, " , ", group.res.display$Mean.Up, ")")
    group.res.display = group.res.display[,c("Response","N", "Mean.CI", "Median")]

    #plot
    data.plot = data[, c(yvar,grpvar,trtvar)]
    names(data.plot) = c("y", "Group", "Trt")

    fig = ggplot(data.plot, aes(x=Group, y=y, fill=Trt)) +
      geom_boxplot(outlier.shape=NA, outlier.colour=NA, position=position_dodge(0.8)) +
      stat_summary(fun=mean, geom="point", shape=5, size=3, position=position_dodge(0.8)) +
      geom_point(size=1.8, position=position_jitterdodge(jitter.width=0.13, dodge.width =0.8))+
      labs(x=grpvar.display, y=yvar.display, fill=trtvar.display)+
      theme_bw(base_size=15)+
      theme(text=element_text(size=15), axis.text=element_text(size=15),legend.text=element_text(size=15))


    res = list(comp.res=comp.res, comp.res.display=comp.res.display,
               group.res=group.res, group.res.display=group.res.display, fig=fig)



  }else if(type=="b"){
    if(!all(data[,yvar]%in%c(0,1,NA)) | !all(c(0,1)%in%data[,yvar])) stop("Response needs to be 0/1.")

    if(is.null(xvars.adj)){
      xvars.adj.fml=NULL
    }else{
      xvars.adj.fml = paste(c("",xvars.adj),collapse="+")
    }

    fml = stats::as.formula(paste(yvar, " ~ ", trtvar, "*",grpvar, xvars.adj.fml, sep=""))
    glm.fit = try(Anova(glm(fml, family=stats::binomial, data=data), type=3, test.statistic="Wald"), silent=T)

    if(!inherits(glm.fit,"try-error")){
      pval.pred = glm.fit[paste(trtvar,":",grpvar,sep=""),"Pr(>Chisq)"]
      pval.prog = glm.fit[grpvar, "Pr(>Chisq)"]

    }else{
      pval.pred = as.numeric(NA)
      pval.prog = as.numeric(NA)
    }

    comp.res = data.frame(pvalue.pred = pval.pred, pvalue.prog=pval.prog)
    comp.res.display = comp.res
    comp.res.display= round(comp.res.display, digits=4)

    header = data.frame(Response=yvar.display, Group.Variable=grpvar.display, Test="GLM ANOVA")
    comp.res = cbind(header, comp.res)
    comp.res.display = cbind(header, comp.res.display)

    #Group Res
    group.res = NULL
    for(unicomb2.i in unicomb2){
      datai = data[datacomb2%in%unicomb2.i,]
      datai.N = dim(datai)[1]
      datai.n = sum(datai[,yvar])
      datai.rr = mean(datai[,yvar])
      datai.rrlo = scoreci(x=sum(datai[,yvar]), n=datai.N, conf.level = 0.95)$conf.int[1]
      datai.rrup = scoreci(x=sum(datai[,yvar]), n=datai.N, conf.level = 0.95)$conf.int[2]

      datai.res = data.frame(Response=yvar.display, Treatment=datai[1,trtvar],
                             Group=datai[1,grpvar], n=datai.n, N=datai.N, Rate=datai.rr,
                             Rate.Lo=datai.rrlo, Rate.Up=datai.rrup)
      names(datai.res) = c("Response", trtvar.display, grpvar.display, "n", "N", "Rate",
                           "Rate.Lo", "Rate.Up")
      group.res = rbind(group.res, datai.res)
    }

    group.res.display = group.res
    group.res.display[,c("Rate", "Rate.Lo", "Rate.Up")] = round(group.res.display[,c("Rate", "Rate.Lo", "Rate.Up")], digits=3)
    group.res.display$Rate.CI = paste(group.res.display$Rate*100, "% (", group.res.display$Rate.Lo*100, "%",
                                      ", ",group.res.display$Rate.Up*100, "%)", sep="")
    group.res.display = group.res.display[,c("Response", trtvar.display, grpvar.display, "n", "N", "Rate.CI")]


    ##PLOT
    data.plot = group.res[,c(trtvar.display, grpvar.display, "Rate",  "Rate.Lo", "Rate.Up")]
    names(data.plot) = c("trt", "grp", "Rate",  "Rate.Lo", "Rate.Up")
    fig = ggplot(data.plot, aes(x=grp, y=Rate, fill=trt)) +
      geom_bar(position=position_dodge(0.9), stat="identity") +
      geom_errorbar(aes(ymin=Rate.Lo, ymax=Rate.Up),width=.2, position=position_dodge(.9))+
      labs(x=grpvar.display, y=yvar.display, fill=trtvar.display)+
      theme_bw(base_size=15)+
      theme(text=element_text(size=15), axis.text=element_text(size=15),legend.text=element_text(size=15))

    res = list(comp.res=comp.res, comp.res.display=comp.res.display,
               group.res=group.res, group.res.display=group.res.display, fig=fig)

  }else{
    stop("type can only be c, b or s.")
  }

}

#' Subgroup Performance Evaluation for Prognostic Cases
#'
#' This function evaluates subgroup performance based on different types of response variables.
#'
#' @param yvar The response variable name.
#' @param censorvar (Optional) The censoring variable name (0-censored, 1-event).
#' @param grpvar The subgroup variable name.
#' @param grpname A vector of ordered subgroup names (values in the column of grpvar).
#' @param xvars.adj (Optional) Other covariates to adjust when evaluating the performance.
#' @param data The dataset containing the variables.
#' @param type The type of response variable: "c" for continuous, "s" for survival, and "b" for binary.
#' @param yvar.display Display name of the response variable.
#' @param grpvar.display Display name of the group variable.
#'
#' @return A list containing subgroup performance results including logrank p-value, median and mean survival, Cox model p-value, ANOVA p-value, and more based on the specified response variable type.
#' @importFrom onewaytests welch.test
#' @importFrom car Anova
#' @export
subgrp_perf = function(yvar, censorvar=NULL, grpvar, grpname, xvars.adj=NULL, data, type, yvar.display=yvar, grpvar.display=grpvar){
  Group=NULL
  y=NULL

  if(is.null(xvars.adj)){
    data = data[, c(yvar,censorvar,grpvar)]
    data = data[stats::complete.cases(data),]
    grpname = unique(grpname)
    if(!all(grpname%in%data[,grpvar]) | !all(data[,grpvar]%in%grpname))
      stop("grpname should be the values for the grpvar variable.")

    data[,grpvar] = factor(data[,grpvar], levels=grpname)

    if(type=="s"){
      if(is.null(censorvar)) stop("Censoring variable missing!")

      #logrank
      fml = stats::as.formula(paste("Surv(", yvar, ",", censorvar,") ~ ",grpvar))
      sdf = try(survdiff(fml, data=data, rho = 0),silent=T)
      if(!inherits(sdf,"try-error")){
        pval.logrank = 1 - stats::pchisq(sdf$chisq, length(sdf$n)-1)
      }else{
        pval.logrank = as.numeric(NA)
      }

      comp.res = data.frame(logrank.pvalue = pval.logrank)
      comp.res.display = comp.res
      comp.res.display$logrank.pvalue = round(comp.res.display$logrank.pvalue, digits=4)

      header = data.frame(Response=yvar.display, Group.Variable=grpvar.display)
      comp.res = cbind(header, comp.res)
      comp.res.display = cbind(header, comp.res.display)


      #median/mean survival time
      sfit = survfit(fml, data = data, conf.type="log-log")
      sfit.res = as.data.frame(summary(sfit, rmean="common")$table[,c("n.start","events", "rmean", "se(rmean)",
                                                                      "median","0.95LCL","0.95UCL")])
      sfit.res$rmean.lower = sfit.res[,"rmean"]-1.96*sfit.res[,"se(rmean)"]
      sfit.res$rmean.upper = sfit.res[,"rmean"]+1.96*sfit.res[,"se(rmean)"]
      sfit.res = sfit.res[, c("n.start","events", "rmean","rmean.lower", "rmean.upper",
                              "median","0.95LCL","0.95UCL")]
      names(sfit.res) = c("N", "Events", "Mean.Surv", "Mean.Surv.Lo", "Mean.Surv.Up", "Med.Surv", "Med.Surv.Lo", "Med.Surv.Up")
      row.names(sfit.res) = paste(grpvar.display,"=",grpname, sep=" ")


      sfit.res.display = sfit.res
      sfit.res.display = round(sfit.res.display, digits=2)
      sfit.res.display$Mean.Surv.CI = paste(sfit.res.display$Mean.Surv, " (", sfit.res.display$Mean.Surv.Lo, " , ",
                                            sfit.res.display$Mean.Surv.Up, ")", sep="")
      sfit.res.display$Med.Surv.CI = paste(sfit.res.display$Med.Surv, " (", sfit.res.display$Med.Surv.Lo, " , ",
                                           sfit.res.display$Med.Surv.Up, ")", sep="")
      sfit.res.display = sfit.res.display[,c("N", "Events", "Mean.Surv.CI", "Med.Surv.CI")]

      sfit.res = cbind(data.frame(Response=yvar.display), sfit.res)
      sfit.res.display = cbind(data.frame(Response=yvar.display), sfit.res.display)

      #plot
      fig = ggsurvplot(surv_fit(fml, data = data, conf.type="log-log"), data = data, risk.table = TRUE,
                       surv.median.line="none", pval=F, pval.method=F,  xlab=yvar.display,
                       legend.title=grpvar.display, legend.labs = grpname,
                       font.legend= 14, title=paste(yvar.display, "~", grpvar.display))


      res = list(comp.res=comp.res, comp.res.display=comp.res.display, group.res=sfit.res, group.res.display=sfit.res.display, fig=fig)


    }else if(type =="c"){
      fml = stats::as.formula(paste(yvar, " ~ ",grpvar))
      welchtest = try(welch.test(fml, data=data, rate = 0, alpha = 0.05, na.rm = TRUE, verbose = F), silent=T)
      kwtest = try(stats::kruskal.test(fml, data=data), silent=T)

      if(!inherits(welchtest,"try-error")){
        pval.welch = welchtest$p.value
      }else{
        pval.welch = as.numeric(NA)
      }

      if(!inherits(kwtest,"try-error")){
        pval.kw = kwtest$p.value
      }else{
        pval.kw = as.numeric(NA)
      }

      comp.res = data.frame(Response=yvar.display, Group.Variable=grpvar.display,
                            Welch.pvalue=pval.welch, KW.pvalue = pval.kw)
      comp.res.display = comp.res
      comp.res.display[,c("Welch.pvalue","KW.pvalue")] = round(comp.res.display[,c("Welch.pvalue","KW.pvalue")], digits=4)


      group.res = NULL
      for(grp.i in grpname){
        datai = data[data[,grpvar]%in%grp.i,]
        datai.N = dim(datai)[1]
        datai.t = try(stats::t.test(datai[,yvar], alternative = "two.sided", mu = 0, conf.level = 0.95), silent=T)

        if(!inherits(datai.t,"try-error")){
          datai.mean = unname(datai.t$estimate)
          datai.lower = datai.t$conf.int[1]
          datai.upper = datai.t$conf.int[2]
          datai.median = stats::median(datai[,yvar])
        }else{
          datai.mean = mean(datai[,yvar])
          datai.lower = as.numeric(NA)
          datai.upper = as.numeric(NA)
          datai.median = stats::median(datai[,yvar])
        }

        datai.res = data.frame(Response=yvar.display, N=datai.N, Mean=datai.mean, Mean.Lo=datai.lower,
                               Mean.Up=datai.upper, Median=datai.median)
        row.names(datai.res) = paste(grpvar.display,"=",grp.i, sep=" ")

        group.res = rbind(group.res, datai.res)
      }

      group.res.display = group.res
      group.res.display[,c("Mean", "Mean.Lo", "Mean.Up", "Median")] = round(group.res.display[,c("Mean", "Mean.Lo", "Mean.Up", "Median")], digits=2)
      group.res.display$Mean.CI = paste(group.res.display$Mean, " (", group.res.display$Mean.Lo, " , ", group.res.display$Mean.Up, ")")
      group.res.display = group.res.display[,c("Response","N", "Mean.CI", "Median")]


      #plot
      data.plot = data[,c(yvar,grpvar)]
      names(data.plot) = c("y", "Group")

      fig = ggplot(data.plot, aes(x=Group, y=y)) + geom_boxplot(outlier.shape=NA) +
        stat_summary(fun=mean, geom="point", shape=5, size=4) +
        geom_point(size=2.5, position=position_jitter(width=0.18,height=0))+
        ylab(yvar.display)+xlab(grpvar.display)+
        theme_bw(base_size=15)+
        theme(text=element_text(size=15), axis.text=element_text(size=15),legend.text=element_text(size=15))

      res = list(comp.res=comp.res, comp.res.display=comp.res.display, group.res=group.res, group.res.display=group.res.display, fig=fig)


    }else if(type=="b"){
      if(!all(data[,yvar]%in%c(0,1,NA)) | !all(c(0,1)%in%data[,yvar])) stop("Response needs to be 0/1.")

      data[,yvar] = factor(data[,yvar], levels=c(0,1))
      cont = table(data[,grpvar], data[,yvar])
      colnames(cont) = paste(yvar.display, "=", colnames(cont))
      rownames(cont) = paste(grpvar.display, "=", rownames(cont))

      fish.p = stats::fisher.test(x=data[,grpvar], y=data[,yvar])$p.value

      comp.res =  data.frame(Response = yvar.display, Group.Variable=grpvar.display,
                             Fisher.pvalue = fish.p)
      comp.res.display = comp.res
      comp.res.display$Fisher.pvalue = round(comp.res.display$Fisher.pvalue, digits=4)

      temp.table = matrix(paste("(",round(cont/rowSums(cont), digits=3)*100, "%)",sep=""), nrow=dim(cont)[1])
      cont.display = matrix(paste(cont, temp.table), nrow=dim(cont)[1])
      rownames(cont.display) = rownames(cont)
      colnames(cont.display) = colnames(cont)

      res = list(comp.res=comp.res, comp.res.display=comp.res.display,
                 group.res=cont, group.res.display=cont.display)


    }else{
      stop("type can only be c, b or s.")
    }
  }else{
    data = data[, c(yvar,censorvar,grpvar,xvars.adj)]
    data = data[stats::complete.cases(data),]
    grpname = unique(grpname)
    if(!all(grpname%in%data[,grpvar]) | !all(data[,grpvar]%in%grpname))
      stop("grpname should be the values for the grpvar variable.")

    data[,grpvar] = factor(data[,grpvar], levels=grpname)

    if(type=="s"){
      if(is.null(censorvar)) stop("Censoring variable missing!")

      #coxph
      xvars.adj.fml = paste(xvars.adj,collapse="+")
      fml = stats::as.formula(paste("Surv(", yvar, ",", censorvar,") ~ ", grpvar, "+",xvars.adj.fml))
      res.cox = try(Anova(coxph(fml, data=data), type=3, test.statistic="Wald"), silent=T)
      if(!inherits(res.cox,"try-error")){
        pval = res.cox[grpvar,"Pr(>Chisq)"]

      }else{
        pval = as.numeric(NA)
      }

      comp.res = data.frame(cox.anova.pvalue = pval)
      comp.res.display = comp.res
      comp.res.display$cox.anova.pvalue = round(comp.res.display$cox.anova.pvalue, digits=4)

      header = data.frame(Response=yvar.display, Group.Variable=grpvar.display)
      comp.res = cbind(header, comp.res)
      comp.res.display = cbind(header, comp.res.display)

      res = list(comp.res=comp.res, comp.res.display=comp.res.display)


    }else if(type=="c"){
      xvars.adj.fml = paste(xvars.adj,collapse="+")
      fml = stats::as.formula(paste(yvar,"~", grpvar, "+",xvars.adj.fml))
      lm.fit = try(Anova(stats::lm(fml, data=data), type=3), silent=T)

      if(!inherits(lm.fit,"try-error")){
        pval = lm.fit[grpvar, "Pr(>F)"]
      }else{
        pval = as.numeric(NA)
      }

      comp.res = data.frame(Response=yvar.display, Group.Variable=grpvar.display,
                            lm.anova.pvalue = pval)
      comp.res.display = comp.res
      comp.res.display$lm.anova.pvalue = round(comp.res.display$lm.anova.pvalue, digits=4)

      res = list(comp.res=comp.res, comp.res.display=comp.res.display)


    }else if(type=="b"){
      if(!all(data[,yvar]%in%c(0,1,NA)) | !all(c(0,1)%in%data[,yvar])) stop("Response needs to be 0/1.")

      xvars.adj.fml = paste(xvars.adj,collapse="+")
      fml = stats::as.formula(paste(yvar,"~", grpvar, "+",xvars.adj.fml))
      glm.fit = try(Anova(glm(fml, family=stats::binomial, data=data), type=3, test.statistic="Wald"), silent=T)

      if(!inherits(glm.fit,"try-error")){
        pval = glm.fit[grpvar, "Pr(>Chisq)"]
      }else{
        pval = as.numeric(NA)
      }

      comp.res = data.frame(Response=yvar.display, Group.Variable=grpvar.display,
                            glm.anova.pvalue = pval)
      comp.res.display = comp.res
      comp.res.display$glm.anova.pvalue = round(comp.res.display$glm.anova.pvalue, digits=4)

      res = list(comp.res=comp.res, comp.res.display=comp.res.display)



    }else{
      stop("type can only be c, b or s.")
    }

  }

  res
}



