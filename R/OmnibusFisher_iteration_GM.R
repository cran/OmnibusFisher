#' @importFrom stats rnorm
#' @importFrom stats cov
#' @importFrom stats pchisq

Iteration_GM <- function(num, full_id, G_all, M_all, evectors1, evectors2, evals1, evals2, pskat1, pskat2, acc, method) { 
   p1<-p2<-vector()
   for (k in 1:num) {
    r_star = rnorm(dim(full_id)[1], 0, 1)
    r_star = cbind(full_id, r_star)
    G_all_star = merge(r_star, G_all)
    G_all_star = G_all_star[order(G_all_star$order),]
    M_all_star = merge(r_star, M_all)
    M_all_star = M_all_star[order(M_all_star$order),]
    r_star1 = t(evectors1) %*% G_all_star$r_star
    r_star2 = t(evectors2) %*% M_all_star$r_star
    if (length(evals1)==1) { Q1_star = r_star1 * evals1 * r_star1
    } else { Q1_star = t(r_star1) %*% diag(evals1) %*% r_star1 }
    if (length(evals2)==1) { Q2_star = r_star2 * evals2 * r_star2
    } else { Q2_star = t(r_star2) %*% diag(evals2) %*% r_star2 }
    if (method=="davies"){
      tmpout1_star<-davies(Q1_star, evals1, acc=acc)
      tmpout2_star<-davies(Q2_star, evals2, acc=acc)
      p1[k]<-tmpout1_star$Qq
      p2[k]<-tmpout2_star$Qq
    }
    else if (method=="kuonen"){
      p1[k]<-pchisqsum(Q1_star,rep(1,length(evals1)),evals1,lower.tail=F,method="sad")
      p2[k]<-pchisqsum(Q2_star,rep(1,length(evals2)),evals2,lower.tail=F,method="sad")
    }
   }

   p1[p1==0|p1<0] = 1e-4
   p1[p1>1] = 1
   p2[p2==0|p2<0] = 1e-4
   p2[p2>1] = 1

   p01_single = pskat1
   p02_single = pskat2
   T01 =  -2*(log(pskat1) + log(pskat2))
   E_T01 = 2*2
   V_T01 = 2*2*2 + 2*cov(-2*log(p1),-2*log(p2))
   v01 = 2*E_T01**2/V_T01
   c01 = v01/E_T01
   p01_dual = pchisq( c01*T01, df=v01, lower.tail=FALSE)
   min0 = min(c(p01_single,p02_single,p01_dual))
  
   p1_single = p1
   p2_single = p2
   T1 =  -2*(log(p1) + log(p2))
   E_T1 = 2*2
   V_T1 = 2*2*2 + 2*cov(-2*log(p1),-2*log(p2))
   v1 = 2*E_T1**2/V_T1
   c1 = v1/E_T1
   p1_dual = pchisq( c1*T1, df=v1, lower.tail=FALSE)
   min_pert = apply(cbind(p1_single,p2_single,p1_dual),1,min)
  
   pval_pert = sum(min_pert<min0)/(length(min_pert)+1)
   if (is.na(pval_pert)) pval_pert=1
   if (is.na(p01_dual)) p01_dual=1
   c(pval_pert, p01_dual, p01_single,p02_single)
}
