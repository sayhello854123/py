##
##Data pre-processing
Processing_1 <- function(clinical,CIBERSOFT){
  pFilter=0.05 
  cli <-data.table::fread(clinical,data.table = F)
  row.names(cli) <- cli[,1]
  cli$futime <- cli$OS.time/365
  cli$fustat <- cli$OS
  cli <- cli[,-c(1:4)]
  rt <- read.table(CIBERSOFT,sep = '\t',header = T,row.names = 1)
  data=rt[rt[,"P.value"]<pFilter,]
  data=data[,1:(ncol(rt)-4)]
  group=sapply(strsplit(rownames(data),"\\-"), "[", 4)
  table(group)
  group=sapply(strsplit(group,""), "[", 1)
  table(group)
  rt1=data[group==0,]
  sameSample=intersect(row.names(rt1), row.names(cli))
  rt=cbind(cli[sameSample,], data[sameSample,])
}

###univariate analysis
Processing_2 <- function(rt){
  colnames_sum <- colnames(rt)
  covariates <- colnames_sum[-which(colnames_sum %in% c("futime","fustat"))]
  univ_formulas <- sapply(covariates,
                          function(x) as.formula(paste("Surv(futime, fustat)~", x)))
  univ_models <- lapply(univ_formulas,
                        function(x) coxph(x, data = rt))
  univ_results <- lapply(univ_models,
                         function(x){ 
                           x <- summary(x)
                           #pvalue
                           p.value<-signif(x$wald["pvalue"], digits=3)
                           #HR
                           HR <-signif(x$coef[2], digits=3);
                           #95 % confidence interval
                           HR.confint.lower <- signif(x$conf.int[,"lower .95"],3)
                           HR.confint.upper <- signif(x$conf.int[,"upper .95"],3)
                           HR <- paste0(HR, " (", 
                                        HR.confint.lower, "-", 
                                        HR.confint.upper, ")")
                           res<-c(p.value,HR)
                           names(res)<-c("p.value","HR (95% CI for HR)")
                           return(res)
                         })
}
#survial_plot
survial_plot <- function(rt){
  rt$risk <- as.vector(ifelse(rt$riskScore>cutOp,"High", "Low"))
  kmfit <- survfit(Surv(futime, fustat) ~ risk, data = rt)
  train_survival_median <- ggsurvplot(kmfit,
                                      data=rt,
                                      pval = TRUE, 
                                      conf.int = F,
                                      xlab="Time(years)",
                                      legend.labs=c("High risk","Low risk" ),
                                      legend.title="Risk",
                                      title="Training set",
                                      font.main = c(15,"bold"),
                                      risk.table = TRUE, 
                                      risk.table.col = "strata", 
                                      linetype = "strata", 
                                      break.time.by = 1,
                                      # surv.median.line = "hv", 
                                      font.family = "Times",
                                      ggtheme = theme_bw(), 
                                      palette = c("#A73030FF", "#0073C2FF"))

}
survStatPlot=function(rt,outPdf=unll){
  rt=rt[order(rt$riskScore),]
  riskClass=rt[,"risk"]
  lowLength=length(riskClass[riskClass=="Low"])
  highLength=length(riskClass[riskClass=="High"])
  color=as.vector(rt$fustat)
  color[color==1]="red"
  color[color==0]="green"
  pdf(file=outPdf,width = 8,height = 9)
  plot(rt$futime,
       pch=19,
       xlab="Patients (increasing risk socre)",
       ylab="Survival time (years)",
       col=color)
  legend("topleft", c("Dead", "Alive"),bty="n",pch=19,col=c("red","green"),cex=1.2)
  abline(v=lowLength,lty=2)
  dev.off()
}

riskScorePlot <- function(rt,outPdf=unll){
  rt=rt[order(rt$riskScore),]
  riskClass=rt[,"risk"]
  lowLength=length(riskClass[riskClass=="Low"])
  highLength=length(riskClass[riskClass=="High"])
  line=rt[,"riskScore"]
  line[line>10]=10
  pdf(file=outPdf,width = 12,height = 5)
  plot(line,
       type="p",
       pch=20,
       xlab="Patients (increasing risk socre)",
       ylab="Risk score",
       col=c(rep("green",lowLength),
             rep("red",highLength)))
  trainMedianScore=median(rt$riskScore)
  abline(h=trainMedianScore,v=lowLength,lty=2)
  dev.off()
}
