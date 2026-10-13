###############################################################################
# Manuscript: Integrating asymmetry, geography and timescale analyses to uncover fish population synchrony patterns in a tropical reservoir
# Author's: Emanuel R. Cararo
# Last update: 06/08/2026
# CODE 01: This code provide all function and packages necessary to run all following analyses. 
################################################################################

# Clean work space ----
rm(list = ls())

# Set seed to reproducibility ----
set.seed(263)

# Load required packages ----

library(tidyverse)
library(readxl)
library(wsyn)
library(ecodist)
library(lm.beta)
library(phytools)
library(ade4)
library(ape)
library(nlme)
library(ggplot2)
library(ggrepel)
library(sf)
library(ggeasy)
library(ggstatsplot)
library(patchwork)

# Packages used only through pkg::function() and therefore also required:
# pracma, ncf, psych, copula, matrixcalc, Matrix, mvtnorm, scales, stringr

# Create and input median species matrices ----

create_species_matrices = function(wide_data, local_col, year_col, month_col, first_species_col,impute_method = "none"){
  species_names = colnames(wide_data)[first_species_col:ncol(wide_data)]
  # 1. Standardize the timeline
  existing_time_ids = unique(paste(wide_data[[year_col]], sprintf("%02d", wide_data[[month_col]]), sep = "-"))
  start_point = min(existing_time_ids)
  end_point   = max(existing_time_ids)
  all_years  = min(wide_data[[year_col]]):max(wide_data[[year_col]])
  all_months = sort(unique(wide_data[[month_col]])) 
  full_grid  = expand.grid(month = all_months, year = all_years)
  full_grid$time_id = paste(full_grid$year, sprintf("%02d", full_grid$month), sep = "-")
  timeline = sort(full_grid$time_id[full_grid$time_id >= start_point & full_grid$time_id <= end_point])
  all_locals = sort(unique(wide_data[[local_col]]))
  all_matrices = list()
  for (species in species_names) {
    agg_df = aggregate(wide_data[[species]] ~ wide_data[[local_col]] +                 
                         paste(wide_data[[year_col]], sprintf("%02d", wide_data[[month_col]]), sep = "-"), FUN = sum)
    colnames(agg_df) = c("local", "time_id", "value")
    mat = matrix(NA, nrow = length(timeline), ncol = length(all_locals), 
                 dimnames = list(timeline, all_locals))
    for (i in 1:nrow(agg_df)) {
      t_id = as.character(agg_df$time_id[i])
      l_id = as.character(agg_df$local[i])
      if (t_id %in% timeline) {
        mat[t_id, l_id] = agg_df$value[i]
      }
    }
    if (impute_method != "none") {
      for (j in 1:ncol(mat)) {
        column_data = mat[, j]
        if (any(is.na(column_data)) && any(!is.na(column_data))) {
          val = if(impute_method == "median") median(column_data, na.rm = TRUE) else mean(column_data, na.rm = TRUE)
          column_data[is.na(column_data)] = val
        }
        mat[, j] = column_data
      }
    }
    all_matrices[[species]] = mat
  }
  return(all_matrices)
}

# Compute partial Spearman correlation Ghosh et al. 2020 Ecosphere ----

partialSpearman=function(u, v, bounds){
  if(!is.numeric(u) | !is.numeric(v)){
    stop("u and v must be numeric")}
  if(length(u)!=length(v)){
    stop("u and v have different length")}
  if(min(bounds)<0 | max(bounds)>1){
    stop("bounds must be between 0 and 1")}
  n=length(u)
  x=rank(u)/(n+1)
  y=rank(v)/(n+1)
  x_b=x[x+y > (2*min(bounds)) & x+y < (2*max(bounds))]
  y_b=y[x+y > (2*min(bounds)) & x+y < (2*max(bounds))]
  r=sum((x_b-mean(x))*(y_b-mean(y)))/((n-1)*sqrt(var(x)*var(y)))
  return(r)
}

# mSynch function from ncf package with Spearman correlation ----
mSynch2 = function (x, y = NULL, resamp = 1000, na.rm = FALSE,
                     quiet = FALSE, use.cor = "pairwise.complete.obs"){
  NAO <- FALSE
  if (any(!is.finite(x))) {
    if (na.rm) {
      warning("Missing values exist; Pairwise deletion will be used")
      NAO <- TRUE
    }
    else {
      stop("Missing values exist; use na.rm = TRUE for pairwise deletion")
    }
  }
  Sbar <- list(real = NA)
  n <- dim(x)[1]
  p <- dim(x)[2]
  x <- as.matrix(x) + 0
  if (!is.null(y)) {
    m <- dim(y)[1]
    y <- as.matrix(y) + 0
  }
  if (!is.null(y)) {
    synch <- cor(t(x), t(y), method="spearman", use = use.cor)
  }
  else {
    synch <- cor(t(x), method="spearman", use = use.cor)
  }
  if (is.null(y)) {
    triang <- lower.tri(synch)
    v <- synch[triang]
  }
  else {
    v <- synch
  }
  Sbar$real <- mean(v, na.rm = TRUE)
  Sbar$boot <- NULL
  if (resamp != 0) {
    Sbar$boot <- matrix(NA, nrow = resamp, ncol = 1)
    for (i in 1:resamp) {
      if (!quiet) {
        cat(i, " of ", resamp, "\n")
      }
      trekkx <- sample(1:n, replace = TRUE)
      if (!is.null(y)) {
        trekky <- sample(1:m, replace = TRUE)
        synchb <- synch[trekkx, trekky]
      }
      else {
        trekky <- trekkx
        synchb <- synch[trekkx, trekkx][diag(n)[trekkx,
                                                trekkx] != 1]
      }
      Sbar$boot[i, 1] <- mean(synchb, na.rm = TRUE)
    }
  }
  res <- Sbar
  class(res) <- "mSynch"
  res
}

# Rearranges data so the ranks match with simulations Walter et al. 2022 Ecology Letters  ----

alignranks=function(dat,sims)
{
  simsrk=apply(FUN=rank,MARGIN=c(2,3),X=sims)
  res=array(NA,dim(sims))
  for (counter in 1:dim(dat)[2])
  {
    res[,counter,]=dat[simsrk[,counter,],counter]
  }
  
  return(res)
}

# Create surrogate data-set Walter et al 2022 Ecology Letters ----

ncsurrog=function(m,numsurrog,plotcheckloc=NA)
{
  N=dim(m)[1]
  n=dim(m)[2]
  if (!all(is.finite(m)))
  {
    stop("Error in ncsurrog: NA, Inf, NaN, etc. not allowed in m")
  }
  scor=cor(m,method="spearman")
  ncov=copula::iRho(copula::normalCopula(0,2),scor)
  ev=eigen(ncov)$values
  if(!matrixcalc::is.positive.semi.definite(ncov)){
    warning("ncov is not +ve semi-definite, replacing with nearest +ve semi-definite matrix",call.=T,immediate.=T)
    
    ncovn=Matrix::nearPD(ncov,corr=TRUE, maxit=1000)
    ncovn=as.matrix(ncovn$mat)
    evn=eigen(ncovn)$values
    comp_ev=cbind(ev,evn)
    colnames(comp_ev)=c("evals_of_old_mat","evals_of_new_mat")
    cat("----This is comparing eigen values of the old and replacement matrices------","\n",comp_ev,"\n")
    ncov=ncovn  
  }
  
  surrogs=array(mvtnorm::rmvnorm(N*numsurrog,mean=rep(0,n),sigma=ncov),
                 dim=c(N,numsurrog,n))
  surrogs=aperm(surrogs,c(1,3,2)) 
  sm=apply(FUN=sort,X=m,MARGIN=2)
  surrogs=alignranks(sm,surrogs)
  if(!is.na(plotcheckloc))
  {
    scor_sur=array(NA,c(n,n,numsurrog))
    for (counter in 1:numsurrog)
    {
      scor_sur[,,counter]=cor(surrogs[,,counter],method="spearman")
    }
    scor_quant=apply(FUN=quantile,X=scor_sur,MARGIN=c(1,2),prob=c(0.005,0.995))
    scor_quant=aperm(scor_quant,c(2,3,1))
    scor=scor[lower.tri(scor)]
    scor_lq=scor_quant[,,1]
    scor_lq=scor_lq[lower.tri(scor_lq)]
    scor_uq=scor_quant[,,2]
    scor_uq=scor_uq[lower.tri(scor_uq)]
    xlimits=range(scor)
    ylimits=range(scor_lq,scor_uq)
    panwd=4
    panht=4
    xht=1
    ywd=1
    gap=.1
    totwd=ywd+panwd+gap
    totht=xht+panht+gap
    grDevices::pdf(file=paste0(plotcheckloc,".pdf"),width=totwd,height=totht)
    par(fig=c(ywd/totwd,
              (ywd+panwd)/totwd,
              (xht)/totht,
              (xht+panht)/totht),
        mai=c(0,0,0,0),mgp=c(3,0.75,0))
    plot(c(-1,1),c(-1,1),type="l",lty="dashed",
         xlim=xlimits,ylim=ylimits)
    for (counter in 1:length(scor))
    {
      lines(rep(scor[counter],2),c(scor_lq[counter],scor_uq[counter]),type="l",col="grey")
    }
    mtext("Data correlation",side=1,line=2)
    mtext("Surrogate correlation quantiles",side=2,line=2)
    grDevices::dev.off()
  }
  return(surrogs)
}

