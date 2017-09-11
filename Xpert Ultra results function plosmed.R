getmainresults <- function(setting, tag, date, chg=F, header=T)
{
  prob <- read.csv(paste0("markovoutput_",tag,"_",setting,"_1e+05_",date,".csv"), header = header)
  outputnames <- c("tbdeaths", "unnecessary", "assaydetected", "treated", "tbrxdeath","otherrxdeath","rxcure","rxfail","rxfailtbdeath","othertbdeath", "unnecessaryRR", "missed", "missedRR", "assayfalsepos")
  if (header==F) colnames(prob) <- colnames(read.csv(paste0("markovoutput_overall_",setting,"_1e+05_20170327.csv"), header = T))
  if (chg) {probx <- read.csv(paste0("markovoutput_overall_",setting,"_1e+05_20170327.csv"), header = T)} 
  if (chg==F) {probx <- prob}
  diffs <- (probx[,paste0(outputnames,".x")] - prob[,paste0(outputnames,".u")])
  names(diffs) <- outputnames; q <- c(0.025,0.1, 0.5, 0.9, 0.975)
  diffratio <- pmax((prob[,"unnecessary.u"] - probx[,"unnecessary.x"]),0)/(prob[,"tbdeaths.u"] - probx[,"tbdeaths.x"])
  diffratio[diffratio>0] <- -100000; mean(diffratio==-100000)
  
  return(list(
    'xpert' = 1000/100000*apply(probx[,paste0(c('tbdeaths','unnecessary'),".x")],2, quantile, q ,na.rm = T),
    'ultra' = 1000/100000*apply(prob[,paste0(c('tbdeaths','unnecessary'),".u")],2, quantile, q,na.rm = T), 
    'diffs' = (1000/100000*apply(diffs,2, quantile, q,na.rm = T))[,c('tbdeaths','unnecessary')], 
    'ratio' = quantile(diffratio,q,na.rm = T))) 
}


# for example: 
getmainresults("India", "highinc", "20170627")
