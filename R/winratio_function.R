#' Sample Size Estimation for the Win Ratio 
#' 
#' A simulation-based sample size calculation approach, relying on the relationship between the win ratio and the rank distribution
#' @param n_arm_1 The sample size in arm 1 
#' @param n_arm_2 The sample size in arm 2 
#' @param alpha The two-sided type I error rate
#' @param p The probability of success in arm 1
#' @param n.iter The number of simulations 
#' @param seed The seed for random number generation
#' @param boot The number of bootstrap samples within each simulation 
#' @return A list with components including the median, 25th percentile, 75th percentile, mean, and standard deviation of Win Ratios across simulations
#' and the power of the proposed sample size under specified assumptions
#' @examples sim_example<-WinRatio_sampsize(n_arm_1=10,n_arm_2=10,alpha=0.05,p=0.6,n.iter=100,seed=1234,boot=100)
#' @export

WinRatio_sampsize <- function(n_arm_1, n_arm_2, alpha=0.05, p, n.iter=500, seed=1000, boot=500){
  
  n_total<-n_arm_1+n_arm_2
  odd<-ifelse(n_total%%2 == 1,1,0)
  
  if (odd==0){
    n_upper_rank<-n_total*0.5
    n_upper_rank_1<-n_upper_rank+1       
  }
  if (odd==1){ 
    n_upper_rank<-floor(n_total*0.5) + rbinom(1,1,0.5)
    n_upper_rank_1<-n_upper_rank+1          
  }
  
  no_cores <- detectCores() - 1
  cl <- makeCluster(no_cores)
  registerDoParallel(cl)
  
  n_hyper<-rep(NA,n.iter)
  win_ratio<-rep(NA,n.iter)
  
  mean_winratio_boot<-rep(NA,n.iter)
  sd_winratio_boot<-rep(NA,n.iter)
  lower_ci_normal<-rep(NA,n.iter)
  upper_ci_normal<-rep(NA,n.iter)
  power_normal<-rep(NA,n.iter)
  
  mean_log_winratio_boot<-rep(NA,n.iter)
  sd_log_winratio_boot<-rep(NA,n.iter)
  lower_log_ci_normal<-rep(NA,n.iter)
  upper_log_ci_normal<-rep(NA,n.iter)
  power_log_normal<-rep(NA,n.iter)
  
  lower_ci_percentile<-rep(NA,n.iter)
  upper_ci_percentile<-rep(NA,n.iter)
  power_percentile<-rep(NA,n.iter)
  
  lower_ci_percentile_pivot<-rep(NA,n.iter)
  upper_ci_percentile_pivot<-rep(NA,n.iter)
  power_pivot_percentile<-rep(NA,n.iter)
  
  lower_empirical<-rep(NA,n.iter)
  upper_empirical<-rep(NA,n.iter)
  power_empirical<-rep(NA,n.iter)
  
  for (j in 1:n.iter){
    set.seed(j+seed)
    
    if (p==0.5){
      #hypergeometric - under null
      n_hyper[j]<-rhyper(1,0.5*n_total,0.5*n_total,n_arm_1) 
    }
    if (p!=0.5){ 
      #alternative
      if (n_arm_1/n_arm_2 != 1){
        x<- -n_arm_1*0.5 + n_arm_1**2/n_total + (p/n_total)*2*n_arm_1*(n_total-n_arm_1) 
        p_avg<-x/n_arm_1
        odds<-oddsWNCHypergeo(p_avg*n_arm_1, 0.5*n_total, 0.5*n_total, n_arm_1, precision=0.1)
      }
      
      if (n_arm_1/n_arm_2 == 1){
        odds<-oddsWNCHypergeo(p*n_arm_1, 0.5*n_total, 0.5*n_total, n_arm_1, precision=0.1)
      }
      n_hyper[j]<-rWNCHypergeo(1,0.5*n_total,0.5*n_total,n_arm_1,odds) #non-central
    }
    
    ranks<-c(1:n_total) #all possible ranks - assuming 0.no ties
    treat_ranks<-c(sample(x=c(1:n_upper_rank),size=n_hyper[j],replace=F),sample(x=c(n_upper_rank_1:n_total),size=n_arm_1-n_hyper[j],replace=F))
    #control ranks are all those not included in treatment arm 
    control_ranks<-setdiff(ranks,treat_ranks)
    
    #estimate win ratio
    win_ratio[j]<-sum(sapply(treat_ranks,function(x){sum(x < control_ranks)}))/sum(sapply(treat_ranks,function(x){sum(x > control_ranks)}))
    
    #Bootstrap resample
    full<- c(treat_ranks,control_ranks)
    names(full)<-c(rep(1,n_arm_1),rep(0,n_arm_2))
    
    #estimate of win ratio for each bootstrap sample within one iteration
    winratio_est<-foreach(i=1:boot, .combine=c) %dopar% {
      set.seed(10000*j + i + seed)
      sample<-sample(full,replace=T,size=n_total)
      sample_treat<-sample[names(sample)==1]
      sample_control<-sample[names(sample)==0]
      #sample_treat<-sample(treat_ranks,replace=T,size=n_arm_1)
      #sample_control<-sample(control_ranks,replace=T,size=n_arm_2)
      #sum(sapply(sample_treat,function(x){sum(ifelse(x < sample_control,1,0))}))/sum(sapply(sample_treat,function(x){sum(ifelse(x > sample_control,1,0))}))
      sum(sapply(sample_treat,function(x){sum(x < sample_control)}))/sum(sapply(sample_treat,function(x){sum(x > sample_control)}))
    }
    
    #bootstrap CI based on percentiles (a/2, 1-a/2)
    lower_limit<-alpha/2
    upper_limit<-1-(alpha/2)
    quantiles<-quantile(winratio_est,c(lower_limit,upper_limit),na.rm=TRUE)
    lower_ci_percentile[j] <-quantiles[1]
    upper_ci_percentile[j] <-quantiles[2]
    power_percentile[j]<-ifelse(quantiles[1]>1|quantiles[2]<1,1,0) #reject if lower bound > 1 or upper bound < 1 (two-sided)
  
    #print(j)
  }
  stopImplicitCluster()
  

  Power<-length(power_percentile[power_percentile==1])/length(power_percentile)
 
  mylist<-list("WR_Median" = summary(win_ratio)[3], "WR_25th" = summary(win_ratio)[2], "WR_75th"=summary(win_ratio)[5], 
               "WR_Mean" = summary(win_ratio)[4], "WR_SD" = sd(win_ratio),"Power"=Power)
  return(mylist)
  WR<-multi_return()
}

