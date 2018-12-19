#' @importFrom stats rnorm
#' @importFrom stats cov
#' @importFrom stats pchisq

Iteration_MR <- function(num, full_id, M_all, R_all, evectors2, evectors3, evals2, evals3, pskat2, pskat3, acc, method) { 
   p2<-p3<-vector()
   for (k in 1:num) {
    r_star = rnorm(dim(full_id)[1], 0, 1)
    r_star = cbind(full_id, r_star)
    M_all_star = merge(r_star, M_all)
    M_all_star = M_all_star[order(M_all_star$order),]
    R_all_star = merge(r_star, R_all)
    R_all_star = R_all_star[order(R_all_star$order),]
    r_star2 = t(evectors2) %*% M_all_star$r_star
    r_star3 = t(evectors3) %*% R_all_star$r_star
    if (length(evals2)==1) { Q2_star = r_star2 * evals2 * r_star2
    } else { Q2_star = t(r_star2) %*% diag(evals2) %*% r_star2 }
    if (length(evals3)==1) { Q3_star = r_star3 * evals3 * r_star3
    } else { Q3_star = t(r_star3) %*% diag(evals3) %*% r_star3 }
    if (method=="davies"){
      tmpout2_star<-davies(Q2_star, evals2, acc=acc)
      tmpout3_star<-davies(Q3_star, evals3, acc=acc)
      p2[k]<-tmpout2_star$Qq
      p3[k]<-tmpout3_star$Qq
    }
    else if (method=="kuonen"){
      p2[k]<-pchisqsum(Q2_star,rep(1,length(evals2)),evals2,lower.tail=F,method="sad")
      p3[k]<-pchisqsum(Q3_star,rep(1,length(evals3)),evals3,lower.tail=F,method="sad")
    }
   }

   p2[p2==0|p2<0] = 1e-4
   p2[p2>1] = 1
   p3[p3==0|p3<0] = 1e-4
   p3[p3>1] = 1

   p02_single = pskat2
   p03_single = pskat3
   T01 =  -2*(log(pskat2) + log(pskat3))
   E_T01 = 2*2
   V_T01 = 2*2*2 + 2*cov(-2*log(p2),-2*log(p3))
   v01 = 2*E_T01**2/V_T01
   c01 = v01/E_T01
   p01_dual = pchisq( c01*T01, df=v01, lower.tail=FALSE)
   min0 = min(c(p02_single,p03_single,p01_dual))
  
   p2_single = p2
   p3_single = p3
   T1 =  -2*(log(p2) + log(p3))
   E_T1 = 2*2
   V_T1 = 2*2*2 + 2*cov(-2*log(p2),-2*log(p3))
   v1 = 2*E_T1**2/V_T1
   c1 = v1/E_T1
   p1_dual = pchisq( c1*T1, df=v1, lower.tail=FALSE)
   min_pert = apply(cbind(p2_single,p3_single,p1_dual),1,min)
  
   pval_pert = sum(min_pert<min0)/(length(min_pert)+1)
   if (is.na(pval_pert)) pval_pert=1
   if (is.na(p01_dual)) p01_dual=1
   c(pval_pert, p01_dual, p02_single,p03_single)
}
