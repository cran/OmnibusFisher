#' @importFrom stats rnorm
#' @importFrom stats cov
#' @importFrom stats pchisq

Iteration_GR <- function(num, full_id, G_all, R_all, evectors1, evectors3, evals1, evals3, pskat1, pskat3, acc, method) { 
   p1<-p3<-vector()
   for (k in 1:num) {
    r_star = rnorm(dim(full_id)[1], 0, 1)
    r_star = cbind(full_id, r_star)
    G_all_star = merge(r_star, G_all)
    G_all_star = G_all_star[order(G_all_star$order),]
    R_all_star = merge(r_star, R_all)
    R_all_star = R_all_star[order(R_all_star$order),]
    r_star1 = t(evectors1) %*% G_all_star$r_star
    r_star3 = t(evectors3) %*% R_all_star$r_star
    if (length(evals1)==1) { Q1_star = r_star1 * evals1 * r_star1
    } else { Q1_star = t(r_star1) %*% diag(evals1) %*% r_star1 }
    if (length(evals3)==1) { Q3_star = r_star3 * evals3 * r_star3
    } else { Q3_star = t(r_star3) %*% diag(evals3) %*% r_star3 }
    if (method=="davies"){
      tmpout1_star<-davies(Q1_star, evals1, acc=acc)
      tmpout3_star<-davies(Q3_star, evals3, acc=acc)
      p1[k]<-tmpout1_star$Qq
      p3[k]<-tmpout3_star$Qq
    }
    else if (method=="kuonen"){
      p1[k]<-pchisqsum(Q1_star,rep(1,length(evals1)),evals1,lower.tail=F,method="sad")
      p3[k]<-pchisqsum(Q3_star,rep(1,length(evals3)),evals3,lower.tail=F,method="sad")
    }
   }

   p1[p1==0|p1<0] = 1e-4
   p1[p1>1] = 1
   p3[p3==0|p3<0] = 1e-4
   p3[p3>1] = 1

   p01_single = pskat1
   p03_single = pskat3
   T01 =  -2*(log(pskat1) + log(pskat3))
   E_T01 = 2*2
   V_T01 = 2*2*2 + 2*cov(-2*log(p1),-2*log(p3))
   v01 = 2*E_T01**2/V_T01
   c01 = v01/E_T01
   p01_dual = pchisq( c01*T01, df=v01, lower.tail=FALSE)
   min0 = min(c(p01_single,p03_single,p01_dual))
  
   p1_single = p1
   p3_single = p3
   T1 =  -2*(log(p1) + log(p3))
   E_T1 = 2*2
   V_T1 = 2*2*2 + 2*cov(-2*log(p1),-2*log(p3))
   v1 = 2*E_T1**2/V_T1
   c1 = v1/E_T1
   p1_dual = pchisq( c1*T1, df=v1, lower.tail=FALSE)
   min_pert = apply(cbind(p1_single,p3_single,p1_dual),1,min)
  
   pval_pert = sum(min_pert<min0)/(length(min_pert)+1)
   if (is.na(pval_pert)) pval_pert=1
   if (is.na(p01_dual)) p01_dual=1
   c(pval_pert, p01_dual, p01_single,p03_single)
}
