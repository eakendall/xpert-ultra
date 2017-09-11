#This script is called by the primary script, defines several key functions
# (Reads in relevant parameters, creates the cohort to which xpert and ultra will be applied, etc)


# read in parameters
readparams <- function(tag, setting, cohortsize, date)
{
  if (tag==paste0("lowinc_",setting,"_",cohortsize, "_", date)) paramfile <- read.csv("lowinc_params.csv")[1:78,1:10] else
    if (tag==paste0("highinc_",setting,"_",cohortsize, "_", date)) paramfile <- read.csv("highinc_params.csv")[1:78,1:10] else
      paramfile <- read.csv("ultra_params.csv")[1:78,1:10]
  rownames(paramfile) <- paste0(paramfile$parameter,  paramfile$regimen, paramfile$resistance, paramfile$setting, paramfile$hiv)
  
  params <- paramfile[,c("paramtype","shape","estimate","low","high")]
  
  # modify parameters for sensitivity analyses, if specified by 'tag'
  if (tag==paste0("pessimisticrr_",setting,"_",cohortsize, "_", date)) 
    { params['deathrx_mdr_hiv',c(3:5)] <- params['deathrx_mdr_nohiv',c(3:5)] <- c(0.14,0.1,0.2)
      params['cureprob_mdr',c(3:5)] <- c(0.63, 0.54, 0.72)
  }
  
  if (tag==paste0("changebehavior_",setting,"_",cohortsize, "_", date)) 
  { 
    params['usens.trace.hiv.xneg',c(3:5)] <- params['usens.trace.hiv.xneg',c(3:5)]*0.9
    params['usens.trace.hiv.xpos',c(3:5)] <- params['usens.trace.hiv.xpos',c(3:5)]*0.9
    params['usens.trace.nohiv.xneg',c(3:5)] <- params['usens.trace.nohiv.xneg',c(3:5)]*0.95
    params['usens.trace.nohiv.xpos',c(3:5)] <- params['usens.trace.nohiv.xpos',c(3:5)]*0.95
    
    params['unspec.trace.nohist.xneg',c(3:5)] <- params['unspec.trace.nohist.xneg',c(3:5)]*0.5
    params['unspec.trace.nohist.xpos',c(3:5)] <- params['unspec.trace.nohist.xpos',c(3:5)]*0.5
    params['unspec.trace.hist.xneg',c(3:5)] <- params['unspec.trace.hist.xneg',c(3:5)]*0.5
    params['unspec.trace.hist.xpos',c(3:5)] <- params['unspec.trace.hist.xpos',c(3:5)]*0.5
  }
  
  if (tag==paste0("changebehavior2_",setting,"_",cohortsize, "_", date)) 
  { 
    params['usens.trace.hiv.xneg',c(3:5)] <- params['usens.trace.hiv.xneg',c(3:5)]*0.99
    params['usens.trace.hiv.xpos',c(3:5)] <- params['usens.trace.hiv.xpos',c(3:5)]*0.99
    params['usens.trace.nohiv.xneg',c(3:5)] <- params['usens.trace.nohiv.xneg',c(3:5)]*0.99
    params['usens.trace.nohiv.xpos',c(3:5)] <- params['usens.trace.nohiv.xpos',c(3:5)]*0.99
    
    params['unspec.trace.nohist.xneg',c(3:5)] <- params['unspec.trace.nohist.xneg',c(3:5)]*0.01
    params['unspec.trace.nohist.xpos',c(3:5)] <- params['unspec.trace.nohist.xpos',c(3:5)]*0.01
    params['unspec.trace.hist.xneg',c(3:5)] <- params['unspec.trace.hist.xneg',c(3:5)]*0.01
    params['unspec.trace.hist.xpos',c(3:5)] <- params['unspec.trace.hist.xpos',c(3:5)]*0.01
  }
  
  if (tag==paste0("chinaempiric_",setting,"_",cohortsize, "_", date)) 
  { 
    params['empiric_China',] <- params['empiric_India',]
  }
  return(params)
}


######################################
# create setting-specific cohorts of TB suspects #
######################################
# using point estimates for primary analysis #

make.cohort <- function(setting, params)
{
  cohortparams <- list()
  cohortparams$rerxfrac <- params[paste0("rerxfrac_",setting),"estimate"]
  cohortparams$rr_new <- params[paste0("rr_new_",setting),"estimate"]
  cohortparams$rr_rerx <- params[paste0("rr_rerx_",setting),"estimate"]
  cohortparams$hiv <- params[paste0("hiv_",setting),"estimate"]
  if (tag==paste0("stochastic_highHIV_",setting,"_",date)) cohortparams$hiv <- cohortparams$hiv*4
  cohortparams$negs_per_case <- params[paste0("negs_per_case_",setting),"estimate"]
  cohortparams$negwithhistory <- params[paste0("negwithhistory_",setting),"estimate"]
  
  cohort <- numeric(10)
  
  cohort <- with(cohortparams, c(rep(1/(negs_per_case+1),8), rep(negs_per_case/(negs_per_case+1),4)) * #tb or no 
                   c( 1-rerxfrac, 1-rerxfrac, rerxfrac, rerxfrac, 1-rerxfrac, 1-rerxfrac, rerxfrac, rerxfrac, 1-negwithhistory, 1-negwithhistory, negwithhistory, negwithhistory)  * # new or prev treated 
                   c(rep(c(hiv, 1-hiv),6)) * # hiv + or no
                   c(rep(c(1-rr_new, 1-rr_rerx, rr_new, rr_rerx), each=2), 1, 1, 1, 1) # mdr or no
  )
  names(cohort) <- c("new.ds.p", "new.ds.n", "ret.ds.p", "ret.ds.n", "new.rr.p", "new.rr.n", "ret.rr.p", "ret.rr.n", "new.neg.p", "new.neg.n", "ret.neg.p", "ret.neg.n")
  
  return(cohort)
}

###########################################################
## convert cohort to long format, w/ individual age, sex ##
###########################################################
make.longcohort <- function(cohortsize, setting, params)
{
  sized <- round(cohortsize*make.cohort(setting,params))
  sized[which.max(sized)] <- sized[which.max(sized)] + (cohortsize-sum(sized)) # if had rounding error add/subtract to the largest group (new negatives)
  
  # pull age and sex distribution of cases into matrix form:
  agesex <- read.csv("agesexdistrib.csv")[2:19,]
  rownames(agesex) <- agesex$lowbound_age
  pmat <- unlist(agesex[,paste0("p_age_",c("f_","m_"), setting)])
  
  # convert a matrix cohort into a list of individuals
  longcohort <- list()
  longcohort$ID <- 1:cohortsize
  longcohort$casetype <- rep(names(sized), sized)
  longcohort$TB <- rep(c(rep(T,8), rep(F,4)), sized) ## TB T/F
  longcohort$RR <- rep(c(rep(F,4), rep(T,4), rep(NA,4)), sized) ## RR T/F/NA
  longcohort$HIV <- rep(rep(c(T,F), times=6), times=sized) ## 
  longcohort$rerx <- rep(c(F,F,T,T,F,F,T,T,F,F,T,T), sized) ## rerx T/F
  psample <- sample(length(unlist(pmat)),size=cohortsize, replace=T, c(unlist(pmat)))
  longcohort$age <- floor(10*runif(cohortsize,0,1))/2 + 
    as.numeric(rownames(agesex))[psample - ifelse(psample>nrow(agesex),nrow(agesex),0) ] #age to nearest half-year (will use 6-mo time steps, and base mortality rate on floor(age))
  longcohort$sex <- ifelse(psample<=length(pmat)/2,1,2) #pull from age/sex distribution matrix. 1=f, 2=m
  longcohort$setting <- setting
  
  return(longcohort)
}


###############
## diagnosis ##
###############
assayresults <- function(tb, rr, hiv, hist, params)
{
  tbresult <- rrresult <- array(dim=c(length(tb),5))
   rrvect <- array(0, dim=c(length(tb),1))
  with(params, {
    # assign result with xpert, ultra trace, ultra no trace, and ultra conditional trace
    #if (tb) 
    xpert <- notrace <- trace <- conditional <- repeattrace <- logical(sum(tb))
    
    xpert <- rbinom(sum(tb),1,prob=c(xsens.hiv, xsens.nohiv)[2-hiv[tb]])
    notrace[xpert==1] <- rbinom(sum(xpert),1,prob=ifelse(hiv[tb], usens.notrace.hiv.xpos, usens.notrace.nohiv.xpos))
    notrace[xpert==0] <- rbinom(sum(!xpert),1,prob=ifelse(hiv[tb], usens.notrace.hiv.xneg, usens.notrace.nohiv.xneg))
    trace[xpert==1] <- rbinom(sum(xpert),1,prob=ifelse(hiv[tb], usens.trace.hiv.xpos, usens.trace.nohiv.xpos))
    trace[xpert==0] <- rbinom(sum(!xpert),1,prob=ifelse(hiv[tb], usens.trace.hiv.xneg, usens.trace.nohiv.xneg))
    conditional <- ifelse(hist[tb],notrace, trace)
    repeattrace[xpert==1] <- rbinom(sum(xpert),1,prob=ifelse(hiv[tb], usens.repeattrace.hiv.xpos, usens.repeattrace.nohiv.xpos))
    repeattrace[xpert==0] <- rbinom(sum(!xpert),1,prob=ifelse(hiv[tb], usens.repeattrace.hiv.xneg, usens.repeattrace.nohiv.xneg))
    tbresult[tb,] <- cbind(xpert, notrace, trace, conditional, repeattrace)

    #if (tb&rr) sensitivity for rr 
    rrvect[tb&rr&tbresult[,1]] <- rbinom(n = sum(tb&rr&tbresult[,1]), size = 1, prob = xrrsens)
    rrvect[tb&rr&!(tbresult[,1])&(tbresult[,2]|tbresult[,3])] <- rbinom(sum(tb&rr&!(tbresult[,1])&(tbresult[,2]|tbresult[,3])), 1, prob = urrsens)
    rrresult[tb&rr,] <- matrix(rrvect[tb&rr,], ncol=5, nrow=sum(tb&rr))*1*tbresult[tb&rr,]
    #if (tb& !rr) # now rr specificity also correlated
    rrresult[tb&!rr,] <- matrix(array(rbinom(sum(tb&!rr),1,(1-rrspec)), dim=c(sum(tb&!rr),1)), ncol=5, nrow=sum(tb&!rr)) *1*tbresult[tb&!rr,]
    
    #if (!tb) tb specificity
    xpert <- notrace <- trace <- conditional <- repeattrace <- logical(sum(!tb))
    
    xpert <- rbinom(sum(!tb),1,prob=c(xnspec.nohist, xnspec.hist)[hist[!tb]+1])
    notrace[xpert==1] <- rbinom(sum(xpert),1,prob=ifelse(hist[!tb], unspec.notrace.hist.xpos, unspec.notrace.nohist.xpos))
    notrace[xpert==0] <- rbinom(sum(!xpert),1,prob=ifelse(hist[!tb], unspec.notrace.hist.xneg, unspec.notrace.nohist.xneg))
    trace[xpert==1] <- rbinom(sum(xpert),1,prob=ifelse(hist[!tb], unspec.trace.hist.xpos, unspec.trace.nohist.xpos))
    trace[xpert==0] <- rbinom(sum(!xpert),1,prob=ifelse(hist[!tb], unspec.trace.hist.xneg, unspec.trace.nohist.xneg))
    conditional <- ifelse(hist[!tb],notrace, trace)
    repeattrace[xpert==1] <- rbinom(sum(xpert),1,prob=ifelse(hist[!tb], unspec.repeattrace.hist.xpos, unspec.repeattrace.nohist.xpos))
    repeattrace[xpert==0] <- rbinom(sum(!xpert),1,prob=ifelse(hist[!tb], unspec.repeattrace.hist.xneg, unspec.repeattrace.nohist.xneg))
    tbresult[!tb,] <- cbind(xpert, notrace, trace, conditional, repeattrace)
    
    rrresult[!tb,] <- matrix(array(rbinom(sum(!tb),1,(rr_falsenegs)), dim=c(sum(!tb),1)), ncol=5, nrow=sum(!tb)) *1*tbresult[!tb,]
    
  results <- array("noTB", dim=c(length(tb),5))
  results[tbresult==T] <- c("TBRS", "TBRR")[rrresult[tbresult==T]+1]
  return(list(tbresult, rrresult, results))
  
  })
}

##########################
## treatment and beyond ##
##########################

# Function to choose each patient's treatment, based on the result of the assay used
choosetreat <- function(assayresult, TB, empiric)
{
  
  e <- rbinom(TB, size = 1, prob = empiric)
  
  treat <-  matrix(c("treatDS","treatMDR", "unnecessaryRx","ongoing","truenegative"),nrow=5, ncol=length(TB))[rbind( 
    ((assayresult=="TBRS"|(assayresult=="noTB"&e))&TB),
    ((assayresult=="TBRR")&TB), 
    ((e|(assayresult=="TBRS"|assayresult=="TBRR"))&!TB), 
    ((assayresult=="noTB")&TB&!e), 
    ((assayresult=="noTB")&!TB&!e) )]
  
  return(treat)
}

# determine ultimate TB outcome for those who will have TB after the current diagnostic/tretment episode
ongoingoutcome <- function(ongcohort, vparams)
{
  # subtract off remaining TB death
  ongcohort$elapsedtime <- ongcohort$elapsedtime + 0.5
  ongcohort$currentstate[!ongcohort$currentRR] <- sample(c("ongoing","tbdeath"), size = sum(ongcohort$currentRR==F), replace = T, prob = c(1-vparams[paste0("tbdeath_untreated_rs_",setting)], vparams[paste0("tbdeath_untreated_rs_",setting)]))
  ongcohort$currentstate[ongcohort$currentRR] <- sample(c("ongoing","tbdeath"), size = sum(ongcohort$currentRR==T), replace = T, prob = c(1-vparams[paste0("tbdeath_untreated_rr")], vparams[paste0("tbdeath_untreated_rr")]))
  return(ongcohort)
}




#####################
#create mortality rate function, for individual age- sex- and setting-dependent mortality risks

morttable <- read.csv("agesexmort.csv")[1:22,]; rownames(morttable)<- morttable$lowbound_age

mort <- function(age, sex, hiv, setting, error=1, step=0.5)
{
  setting <- ifelse(hiv,"ZA",setting)
  age <- as.numeric(age); sex <- as.numeric(sex)
  rates <- error*as.numeric(morttable[cbind(as.character(pmin(120,age-age%%5)), paste0("mortrate_",ifelse(sex==1,"f_","m_"),setting))])*step
  rates[rates>1] <- 1
  return(rates)
}


##########################
# function for use in setting beta parameter distributions
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}