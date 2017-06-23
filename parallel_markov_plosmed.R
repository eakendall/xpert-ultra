setting <- commandArgs(trailingOnly=TRUE)[1] # ZA or India or China
fix <- commandArgs(trailingOnly=TRUE)[2] # NA, or any parameters to be assigned fixed values, e.g. "tbdeath_untreated_rs_ZA"
tag <- commandArgs(trailingOnly=TRUE)[3] # for distinguishing filenames, e.g. "samecfr", "lowprev", "overall", "pessimisticrr", "changebehavior", "chinaempiric", 'error'
if (tag=='error') taskid <- paste0('.', commandArgs(trailingOnly=TRUE)[4]) else taskid <- ''

print(paste0("tag is ",tag))

date <- "20170623"

require('triangle'); require('lhs'); require(parallel)

cohortsize <- 100000 
nsims <- 5000
varyassay <- T
cores <-detectCores(); savelocation <- "../scratch/ultradata/"
#   cores <- 1; savelocation <- ""
runname <- paste0(tag,"_",setting,"_",cohortsize,"_",date) 

source("Xpert Ultra cohorts plosmed.R")

readcohort <- data.frame(longcohort, stringsAsFactors = F)

vary <- (params$paramtype=="outcome"|params$paramtype=="outcome2"); if (varyassay) vary <- vary|params$paramtype=="assay"
varybeta <- vary&params$shape=="beta"
  
if(!is.na(fix)) vary[which(rownames(params)==fix)] <- FALSE
lhs <- randomLHS(nsims, sum(vary))

header <- c(rownames(params), 
            paste0(rep(c("assaydetected", "treated", "tbrxdeath" , "otherrxdeath", "rxcure", "rxfail", "rxfailtbdeath", "othertbdeath", "tbdeaths", "tbdeaths2", 
              "unnecessary", "unnecessary2", "assayfalsepos", "unnecessaryRR", "unnecessaryRR2", "missed", "missedRR", "cost", "ongoingtransmission", "ongoingtransmissionRR"), each=5), 
                   c(".x", ".notrace", ".u", ".condtrace", ".reptrace")))
write(header, file=paste0(savelocation,"markovoutput_",runname,".csv"), ncolumns=length(header), append = FALSE, sep=",")

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

ongoingoutcome <- function(ongcohort, vparams)
{
  # subtract off remaining TB death
  ongcohort$elapsedtime <- ongcohort$elapsedtime + 0.5
  ongcohort$currentstate[!ongcohort$currentRR] <- sample(c("ongoing","tbdeath"), size = sum(ongcohort$currentRR==F), replace = T, prob = c(1-vparams[paste0("tbdeath_untreated_rs_",setting)], vparams[paste0("tbdeath_untreated_rs_",setting)]))
  ongcohort$currentstate[ongcohort$currentRR] <- sample(c("ongoing","tbdeath"), size = sum(ongcohort$currentRR==T), replace = T, prob = c(1-vparams[paste0("tbdeath_untreated_rr")], vparams[paste0("tbdeath_untreated_rr")]))
  return(ongcohort)
}


# imperfectculture <- function(vparams)
# {
#   #Suppose culture is 98% sensitive, but some of the supposed added non-specificity of ultra (half of ultra-only positives in those without history) are culture-negative true cases.
#   vparams2 <- vparams
#   vparams2["xsens.hiv"] <- vparams["xsens.hiv"]*0.98
#   vparams2["xsens.nohiv"] <- vparams["xsens.nohiv"]*0.98
#   # and ultra sensitivity becomes (old ultra/(1- oldxpert) - 2% + (unspec/2)/2% )( 1 - new xpert).
#   vparams2["usens.trace.hiv"] <- min(1, ( (vparams["xsens.hiv"] + vparams["usens.trace.hiv"]*(1-vparams["xsens.hiv"]))*0.98 +
#                                             (vparams["unspec.notrace.nohist"]/2) -
#                                             vparams2["xsens.hiv"])/(1-vparams2["xsens.hiv"]))
#   vparams2["usens.trace.nohiv"] <- min(1, ( (vparams["xsens.nohiv"] + vparams["usens.trace.nohiv"]*(1-vparams["xsens.hiv"]))*0.98 + (vparams["unspec.notrace.nohist"]/2) - vparams2["xsens.nohiv"])/(1-vparams2["xsens.nohiv"]))
#   # And unspec becomes (unspec/2)/(1-unspec/2)
#   vparams2["unspec.trace.nohist"] <- vparams["unspec.trace.nohist"]/2/(1-vparams["unspec.trace.nohist"]/2)
# 
#   return(vparams2)
# }

params[,"alpha"] <- params[,"beta"] <- NA; 
params$alpha[params$shape=="beta"] <- estBetaParams(params$estimate[params$shape=="beta"], 
                                                    (((params$high-params$low)/4)^2)[params$shape=="beta"])$alpha
params$beta[params$shape=="beta"] <- estBetaParams(params$estimate[params$shape=="beta"], 
                                                   (((params$high-params$low)/4)^2)[params$shape=="beta"])$beta


onerun <- function(n, lhs, params, vary, readcohort, fix)
{
  ##################################
  ## sample parameters for each run ##
  vparams <- qtriangle(p = lhs[n,], a = params$low[vary], b=params$high[vary],c = params$estimate[vary])
  vparams[varybeta[vary]] <- qbeta(p = lhs[n,varybeta[vary]], shape1 = params$alpha[varybeta], shape2=params$beta[varybeta])
  names(vparams) <- rownames(params)[vary]

 # if (imperfcx) vparams <- imperfectculture(vparams)
  tempparams <- params
  tempparams$estimate[vary] <- vparams
  estimates <- as.list(tempparams$estimate); names(estimates) <- rownames(tempparams)
  if(!is.na(fix)) {vparams <- append(vparams, params[rownames(params)==fix,"estimate"]); if (sum(rownames(params)==fix)==1) names(vparams)[length(vparams)] <- fix}
  out <- character()
  
  #assay result
  readcohort[,10:14] <- assayresults(readcohort$TB, readcohort$RR, readcohort$HIV, readcohort$rerx, estimates)[[3]]
  colnames(readcohort)[10:14] <- paste0(c("xpert", "notrace", "trace", "conditional","repeattrace"),"result")
  
  #initial treatment decision. 
  readcohort[,15:19] <- NA
  colnames(readcohort)[15:19] <- paste0(c("xpert", "notrace", "trace", "conditional","repeattrace"),"treat")
  
  empiric <- vparams[paste0("empiric_",setting)]

readcohort$xperttreat <- choosetreat(readcohort$xpertresult, readcohort$TB, empiric)

readcohort$tracetreat[readcohort$traceresult!=readcohort$xpertresult] <- choosetreat(
    readcohort$traceresult[readcohort$traceresult!=readcohort$xpertresult], readcohort$TB[readcohort$traceresult!=readcohort$xpertresult], empiric)
  readcohort$tracetreat[readcohort$traceresult==readcohort$xpertresult] <- readcohort$xperttreat[readcohort$traceresult==readcohort$xpertresult]
  
readcohort$notracetreat[readcohort$notraceresult==readcohort$xpertresult] <- readcohort$xperttreat[readcohort$notraceresult==readcohort$xpertresult]
  readcohort$notracetreat[readcohort$notraceresult==readcohort$traceresult] <- readcohort$tracetreat[readcohort$notraceresult==readcohort$traceresult]
  if(sum((readcohort$notraceresult!=readcohort$traceresult)&(readcohort$notraceresult!=readcohort$xpertresult))>0)
    readcohort$notracetreat[(readcohort$notraceresult!=readcohort$traceresult)&(readcohort$notraceresult!=readcohort$xpertresult)] <- choosetreat(
    readcohort$notraceresult[(readcohort$notraceresult!=readcohort$traceresult)&(readcohort$notraceresult!=readcohort$xpertresult)], 
    readcohort$TB[(readcohort$notraceresult!=readcohort$traceresult)&(readcohort$notraceresult!=readcohort$xpertresult)], empiric)

readcohort$conditionaltreat[readcohort$rerx] <- readcohort$notracetreat[readcohort$rerx]; readcohort$conditionaltreat[!readcohort$rerx] <- readcohort$tracetreat[!readcohort$rerx]
  
readcohort$repeattracetreat[readcohort$repeattraceresult==readcohort$xpertresult] <- readcohort$xperttreat[readcohort$repeattraceresult==readcohort$xpertresult]
readcohort$repeattracetreat[readcohort$repeattraceresult==readcohort$traceresult] <- readcohort$tracetreat[readcohort$repeattraceresult==readcohort$traceresult]
if(sum((readcohort$repeattraceresult!=readcohort$traceresult)&(readcohort$repeattraceresult!=readcohort$xpertresult))>0)
  readcohort$repeattracetreat[(readcohort$repeattraceresult!=readcohort$traceresult)&(readcohort$repeattraceresult!=readcohort$xpertresult)] <- choosetreat(
    readcohort$repeattraceresult[(readcohort$repeattraceresult!=readcohort$traceresult)&(readcohort$repeattraceresult!=readcohort$xpertresult)], 
    readcohort$TB[(readcohort$repeattraceresult!=readcohort$traceresult)&(readcohort$repeattraceresult!=readcohort$xpertresult)], empiric)

readcohort$currentstate <- NA;
  readcohort$unnecessaryRx <- readcohort$unnecessaryMDR <- readcohort$missedTB <- readcohort$missedRR <- readcohort$empirictreatment <- readcohort$ongoingtransmission <- 
  readcohort$treated <- readcohort$cured <- readcohort$failed <- readcohort$otherrxdeath <- NA
  readcohort$currentRR<- readcohort$RR
  readcohort$DSRxs <- readcohort$MDRRxs <- readcohort$elapsedtime <- 0
  
  #for those with TB, determine outcomes if treated versus untreated, and then assign outcomes to assays accordingly
  DStreated <- MDRtreated <- untreated <- readcohort[readcohort$TB,] #ignoring true negatives for now, will assign their outcomes as separate step
  
  #assuming no empiricdeath for now:   if("deathrx_empiric" %in% names(vparams)) DStreated$currentstate[DStreated$result] <- sample(c("treatDS","tbdeath"),size = sum(cohort$empiricdiagnosis), replace=T, prob = c(1-vparams["deathrx_empiric"], vparams["deathrx_empiric"]))
  
  # if treatDS (hypothetical for whole case cohort for now)
  DStreated$elapsedtime <- 0.5; DStreated$DSRxs <- 1
  #DS treated as DS
  deathrx <- vparams[paste0("deathrx_ds_",setting,c("_hiv","_nohiv")[2-DStreated$HIV])]
  
  DStreated$currentstate[!DStreated$RR] <- apply(  cbind( 
    (deathrx), 
    (1-deathrx) * mort(DStreated$age + DStreated$elapsedtime, DStreated$sex, DStreated$HIV, setting, error=vparams["mortalityrate_error"]), 
    (1-deathrx) * (1-mort(DStreated$age + DStreated$elapsedtime, DStreated$sex, DStreated$HIV, setting, error=vparams["mortalityrate_error"]) ) * 
      vparams[paste0("cureprob_ds_rs_",setting)], 
    (1-deathrx) * (1-mort(DStreated$age + DStreated$elapsedtime, DStreated$sex, DStreated$HIV, setting, error=vparams["mortalityrate_error"]) ) * 
      (1-vparams[paste0("cureprob_ds_rs_",setting)]) )[!DStreated$RR,] , 
    MARGIN = 1, FUN = function(x) sample(c("tbrxdeath", "otherdeath", "cure", "ongoing"), size=1, prob = x))
  
  # with prob of acqres:
  DStreated$currentRR[DStreated$currentstate=="ongoing"&!DStreated$RR] <- sample(c(T,F), size = sum(DStreated$currentstate=="ongoing"&!DStreated$RR), replace = T, prob=c(vparams["acqres"], 1-vparams["acqres"]))
  
  #RR treated as DS # need to use ds mortality, otherwise being treated is worse than not treated
  DStreated$currentstate[DStreated$RR] <- apply(  cbind( 
    deathrx, 
    (1-deathrx) * mort(DStreated$age  + DStreated$elapsedtime, DStreated$sex, DStreated$HIV, setting, error=vparams["mortalityrate_error"]), 
    (1-deathrx) * (1-mort(DStreated$age + DStreated$elapsedtime, DStreated$sex, DStreated$HIV, setting, error=vparams["mortalityrate_error"])) * vparams[paste0("cureprob_ds_rr")], 
    (1-deathrx) * (1-mort(DStreated$age + DStreated$elapsedtime, DStreated$sex, DStreated$HIV, setting, error=vparams["mortalityrate_error"])) * (1-vparams[paste0("cureprob_ds_rr")]) )[DStreated$RR,] , 
    MARGIN = 1, FUN = function(x) sample(c("tbrxdeath", "otherdeath", "cure", "ongoing"), size=1, prob = x))
  DStreated$missedRR[DStreated$RR] <- TRUE
    
  DStreated[DStreated$currentstate=="ongoing","ongoingtransmission"] <- TRUE
  DStreated$treated <- TRUE
  DStreated$cured <- (DStreated$currentstate=="cure")
  DStreated$failed <- (DStreated$currentstate=="ongoing")
  DStreated$otherrxdeath <- (DStreated$currentstate=="otherdeath")
  
  
  # and get eventual TB death outcomes if ongoing after treatment
  DStreated[DStreated$currentstate=="ongoing",] <- ongoingoutcome(DStreated[DStreated$currentstate=="ongoing",], vparams)
  

  # if treat MDR
  MDRtreated$elapsedtime <- 0.5; MDRtreated$MDRRxs <- 1
  
  deathrx <- vparams[paste0("deathrx_mdr",c("_hiv","_nohiv")[2-MDRtreated$HIV])]
  #DS treated as RR   # assume same death and cure as RR on MDR for now
  MDRtreated$currentstate[!MDRtreated$RR] <- apply(  cbind( 
    deathrx, 
    (1-deathrx) * mort(MDRtreated$age + MDRtreated$elapsedtime, MDRtreated$sex, MDRtreated$HIV, setting, error=vparams["mortalityrate_error"]), 
    (1-deathrx) * (1-mort(MDRtreated$age + MDRtreated$elapsedtime, MDRtreated$sex, MDRtreated$HIV, setting, error=vparams["mortalityrate_error"])) * vparams["cureprob_mdr"],
    (1-deathrx) * (1-mort(MDRtreated$age + MDRtreated$elapsedtime, MDRtreated$sex, MDRtreated$HIV, setting, error=vparams["mortalityrate_error"])) * (1-vparams["cureprob_mdr"]) )[!MDRtreated$RR,] , 
    MARGIN = 1, FUN = function(x) sample(c("tbrxdeath", "otherdeath", "cure", "ongoing"), size=1, prob = x))
  MDRtreated$unnecessaryMDR[!MDRtreated$RR] <- TRUE
  
  #RR treated as RR 
  MDRtreated$currentstate[MDRtreated$RR] <- apply(  cbind( 
    deathrx, 
    (1-deathrx) * mort(MDRtreated$age + MDRtreated$elapsedtime, MDRtreated$sex, MDRtreated$HIV, setting, error=vparams["mortalityrate_error"]), 
    (1-deathrx) * (1-mort(MDRtreated$age + MDRtreated$elapsedtime, MDRtreated$sex, MDRtreated$HIV, setting, error=vparams["mortalityrate_error"])) * vparams["cureprob_mdr"], 
    (1-deathrx) * (1-mort(MDRtreated$age + MDRtreated$elapsedtime, MDRtreated$sex, MDRtreated$HIV, setting, error=vparams["mortalityrate_error"])) * (1-vparams["cureprob_mdr"]) )[MDRtreated$RR,] , 
    MARGIN = 1, FUN = function(x) sample(c("tbrxdeath", "otherdeath", "cure", "ongoing"), size=1, prob = x))
  
  MDRtreated[MDRtreated$currentstate=="ongoing","ongoingtransmission"] <- TRUE
  MDRtreated$treated <- TRUE
  MDRtreated$cured <- (MDRtreated$currentstate=="cure")
  MDRtreated$failed <- (MDRtreated$currentstate=="ongoing")
  MDRtreated$otherrxdeath <- (MDRtreated$currentstate=="otherdeath")

  # and set outcomes if ongoing after treatment, independent of those for DS treatment for now (could make the same if ongoig for both and same RR, but that's complicated)
  MDRtreated[MDRtreated$currentstate=="ongoing",] <- ongoingoutcome(MDRtreated[MDRtreated$currentstate=="ongoing",], vparams)
  
  # if not treated, apply ongoing function (which will also be used for those ongoing after unsuccessful treatment)
  untreated$currentstate <- "ongoing"
  untreated$ongoingtransmission <- TRUE
  untreated$missedTB <- TRUE
  untreated$missedRR[untreated$RR] <- TRUE
  untreated <- ongoingoutcome(untreated, vparams)
  untreated$treated <- FALSE
  untreated$cured <- NA
  untreated$failed <- NA
  untreated$otherrxdeath <- NA

  
  # now assign result depending on treatment, and results for negatives
  finalcohort <- list()
  for (assay in c("xpert", "notrace", "trace", "conditional", "repeattrace"))
  {
    finalcohort[[assay]] <- readcohort;
    finalcohort[[assay]][readcohort$TB,] [subset(readcohort,TB)[[paste0(assay,"treat")]]=="treatDS",] <- DStreated[readcohort[[paste0(assay,"treat")]]=="treatDS",]
    finalcohort[[assay]][readcohort$TB,] [subset(readcohort,TB)[[paste0(assay,"treat")]]=="treatMDR",] <- MDRtreated[readcohort[[paste0(assay,"treat")]]=="treatMDR",]
    finalcohort[[assay]][readcohort$TB,] [subset(readcohort,TB)[[paste0(assay,"treat")]]=="ongoing",] <- untreated[readcohort[[paste0(assay,"treat")]]=="ongoing",]
    
    finalcohort[[assay]][!readcohort$TB,"currentstate"] [subset(readcohort,!TB)[[paste0(assay,"treat")]]=="unnecessaryRx"] <- "unnecessaryRx"
    finalcohort[[assay]][!readcohort$TB,"unnecessaryRx"] <- finalcohort[[assay]][!readcohort$TB,"currentstate"]=="unnecessaryRx"
    finalcohort[[assay]][!readcohort$TB,"DSRxs"] [subset(readcohort,!TB)[[paste0(assay,"treat")]]=="unnecessaryRx"] <- 1
    finalcohort[[assay]][!readcohort$TB,"currentstate"] [subset(readcohort,!TB)[[paste0(assay,"treat")]]=="truenegative"] <- "truenegative"
    
    finalcohort[[assay]]["empirictreatment"] <- (readcohort[,paste0(assay,"result")]=="noTB")&(readcohort[,paste0(assay,"treat")]=="treatDS")
    
   }
  
  #  and track all outcomes
assaydetected <- treated <- tbrxdeath <- otherrxdeath <- rxcure <- rxfail <- rxfailtbdeath <- othertbdeath <- tbdeaths <- tbdeaths2 <- 
  unnecessary <- unnecessary2 <- assayfalsepos <- unnecessaryRR <- unnecessaryRR2 <- missed <- missedRR <- cost <- ongoingtransmission <- ongoingtransmissionRR <- numeric()

  for (assay in c("xpert", "notrace", "trace", "conditional", "repeattrace"))
  {
    assaydetected <- append(assaydetected, sum(readcohort$TB&finalcohort[[assay]][paste0(assay,"result")]!="noTB", na.rm=T))
    treated <- append(treated, sum(readcohort$TB & finalcohort[[assay]][paste0(assay,"treat")]!="ongoing", na.rm=T))
    tbrxdeath <-  append(tbrxdeath, sum(finalcohort[[assay]]$currentstate=="tbrxdeath") )
    otherrxdeath <- append(otherrxdeath, sum(readcohort$TB & finalcohort[[assay]]$otherrxdeath, na.rm=T) )
    rxcure <- append(rxcure, sum(readcohort$TB & finalcohort[[assay]]$cured, na.rm=T) )
    rxfail <- append(rxfail, sum(readcohort$TB & finalcohort[[assay]]$failed, na.rm=T) )
    rxfailtbdeath <- append(rxfailtbdeath, sum(readcohort$TB & finalcohort[[assay]]$failed & finalcohort[[assay]]$currentstate=="tbdeath", na.rm=T) )
    othertbdeath <- append(othertbdeath, sum(readcohort$TB & !finalcohort[[assay]]$treated  & finalcohort[[assay]]$currentstate=="tbdeath", na.rm=T) )
    
    tbdeaths <-  append(tbdeaths, 
                        sum(finalcohort[[assay]]$currentstate=="tbdeath")+ sum(finalcohort[[assay]]$currentstate=="tbrxdeath") )

    tbdeaths2 <-  append(tbdeaths2, 
                         sum(finalcohort[[assay]]$currentstate=="tbrxdeath") + 
                           sum(readcohort$TB & finalcohort[[assay]]$failed & finalcohort[[assay]]$currentstate=="tbdeath", na.rm=T)+
                           sum(readcohort$TB & !finalcohort[[assay]]$treated  & finalcohort[[assay]]$currentstate=="tbdeath", na.rm=T) )
    
    
    unnecessary <- append( unnecessary, sum(finalcohort[[assay]]$unnecessaryRx, na.rm = T))
    unnecessary2 <- append( unnecessary2, sum(finalcohort[[assay]][paste0(assay,"treat")]=="unnecessaryRx", na.rm = T))

    assayfalsepos <- append(assayfalsepos, sum(finalcohort[[assay]][paste0(assay,"result")]=="TBRS" & !readcohort$TB, na.rm = T))

    unnecessaryRR <- append(unnecessaryRR, sum(finalcohort[[assay]]$unnecessaryMDR, na.rm = T) )
    unnecessaryRR2 <- append( unnecessaryRR2, sum(finalcohort[[assay]][paste0(assay,"treat")]=="treatMDR" & (!readcohort$RR), na.rm = T))

    missed <- append(missed, sum(finalcohort[[assay]]$missedTB, na.rm = T) )
    
    missedRR <- append(missedRR, sum(finalcohort[[assay]]$missedRR, na.rm = T) )
    

    cost <- append(cost, vparams[paste0("rxcost_ds_",setting)]*(sum(finalcohort[[assay]]$DSRxs)) + 
                     vparams[paste0("rxcost_mdr_",setting)]*(sum(finalcohort[[assay]]$MDRRxs)) )
    
    ongoingtransmission <- append(ongoingtransmission, sum(finalcohort[[assay]]$ongoingtransmission, na.rm = T) )
    
    ongoingtransmissionRR <- append(ongoingtransmissionRR, sum(finalcohort[[assay]]$ongoingtransmission[finalcohort[[assay]]$currentRR==T], na.rm = T) )
    
  }
  
  out <- c(out, c(tempparams$estimate, 
                  assaydetected, treated, tbrxdeath, otherrxdeath, rxcure, rxfail, rxfailtbdeath, othertbdeath, tbdeaths, tbdeaths2, 
                    unnecessary, unnecessary2,assayfalsepos, unnecessaryRR, unnecessaryRR2, missed, missedRR, cost, ongoingtransmission, ongoingtransmissionRR))
  
  print(n)
  
  return(out)
  
}

pruns <- mclapply(X = 1:nsims, FUN = onerun, lhs = lhs, params=params, vary=vary, fix=fix, readcohort=readcohort, mc.cores=cores)

write(unlist(pruns), file=paste0(savelocation,"markovoutput_",runname, taskid,".csv"), ncolumns=length(header), append = TRUE, sep=",")
