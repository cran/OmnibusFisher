#' @importFrom stats rnorm
#' @importFrom stats cov
#' @importFrom stats pchisq

Iteration_GMR <- function(num, full_id, G_all, M_all, R_all, evectors1, evectors2, evectors3, evals1, evals2, evals3, pskat1, pskat2, pskat3, acc, method) { 
   p1<-p2<-p3<-vector()
   for (k in 1:num) {
    r_star = rnorm(dim(full_id)[1], 0, 1)
    r_star = cbind(full_id, r_star)
    G_all_star = merge(r_star, G_all)
    G_all_star = G_all_star[order(G_all_star$order),]
    M_all_star = merge(r_star, M_all)
    M_all_star = M_all_star[order(M_all_star$order),]
    R_all_star = merge(r_star, R_all)
    R_all_star = R_all_star[order(R_all_star$order),]
    r_star1 = t(evectors1) %*% G_all_star$r_star
    r_star2 = t(evectors2) %*% M_all_star$r_star
    r_star3 = t(evectors3) %*% R_all_star$r_star
    if (length(evals1)==1) { Q1_star = r_star1 * evals1 * r_star1
    } else { Q1_star = t(r_star1) %*% diag(evals1) %*% r_star1 }
    if (length(evals2)==1) { Q2_star = r_star2 * evals2 * r_star2
    } else { Q2_star = t(r_star2) %*% diag(evals2) %*% r_star2 }
    if (length(evals3)==1) { Q3_star = r_star3 * evals3 * r_star3
    } else { Q3_star = t(r_star3) %*% diag(evals3) %*% r_star3 }
    if (method=="davies"){
      tmpout1_star<-davies(Q1_star, evals1, acc=acc)
      tmpout2_star<-davies(Q2_star, evals2, acc=acc)
      tmpout3_star<-davies(Q3_star, evals3, acc=acc) 
      p1[k]<-tmpout1_star$Qq
      p2[k]<-tmpout2_star$Qq
      p3[k]<-tmpout3_star$Qq
    }
    else if (method=="kuonen"){
      p1[k]<-pchisqsum(Q1_star,rep(1,length(evals1)),evals1,lower.tail=F,method="sad")
      p2[k]<-pchisqsum(Q2_star,rep(1,length(evals2)),evals2,lower.tail=F,method="sad")
      p3[k]<-pchisqsum(Q3_star,rep(1,length(evals3)),evals3,lower.tail=F,method="sad")
    }
   }

   p1[p1==0|p1<0] = 1e-4
   p1[p1>1] = 1
   p2[p2==0|p2<0] = 1e-4
   p2[p2>1] = 1
   p3[p3==0|p3<0] = 1e-4
   p3[p3>1] = 1

   p01_single = pskat1
   p02_single = pskat2
   p03_single = pskat3
   T01 =  -2*(log(pskat1) + log(pskat2))
   T02 =  -2*(log(pskat1) + log(pskat3))
   T03 =  -2*(log(pskat2) + log(pskat3))
   T0tri = -2*(log(pskat1) + log(pskat2) + log(pskat3))
   E_T01 = E_T02 = E_T03 = 2*2
   E_T0tri = 2*3
   V_T01 = 2*2*2 + 2*cov(-2*log(p1),-2*log(p2))
   V_T02 = 2*2*2 + 2*cov(-2*log(p1),-2*log(p3))
   V_T03 = 2*2*2 + 2*cov(-2*log(p2),-2*log(p3))
   V_T0tri = 2*2*3 + 2*cov(-2*log(p1),-2*log(p2)) + 2*cov(-2*log(p1),-2*log(p3)) + 2*cov(-2*log(p2),-2*log(p3))
   v01 = 2*E_T01**2/V_T01
   v02 = 2*E_T02**2/V_T02
   v03 = 2*E_T03**2/V_T03
   v0tri = 2*E_T0tri**2/V_T0tri
   c01 = v01/E_T01
   c02 = v03/E_T02
   c03 = v02/E_T03
   c0tri = v0tri/E_T0tri
   p01_dual = pchisq( c01*T01, df=v01, lower.tail=FALSE)
   p02_dual = pchisq( c02*T02, df=v02, lower.tail=FALSE)
   p03_dual = pchisq( c03*T03, df=v03, lower.tail=FALSE)
   p0_tri = pchisq( c0tri*T0tri, df=v0tri, lower.tail=FALSE)
   min0 = min(c(p01_single,p02_single,p03_single,p01_dual,p02_dual,p03_dual,p0_tri))
  
   p1_single = p1
   p2_single = p2
   p3_single = p3
   T1 =  -2*(log(p1) + log(p2))
   T2 =  -2*(log(p1) + log(p3))
   T3 =  -2*(log(p2) + log(p3))
   Ttri = -2*(log(p1) + log(p2) + log(p3))
   E_T1 = E_T2 = E_T3 = 2*2
   E_Ttri = 2*3
   V_T1 = 2*2*2 + 2*cov(-2*log(p1),-2*log(p2))
   V_T2 = 2*2*2 + 2*cov(-2*log(p1),-2*log(p3))
   V_T3 = 2*2*2 + 2*cov(-2*log(p2),-2*log(p3))
   V_Ttri = 2*2*3 + 2*cov(-2*log(p1),-2*log(p2)) + 2*cov(-2*log(p1),-2*log(p3)) + 2*cov(-2*log(p2),-2*log(p3))
   v1 = 2*E_T1**2/V_T1
   v2 = 2*E_T2**2/V_T2
   v3 = 2*E_T3**2/V_T3
   vtri = 2*E_Ttri**2/V_Ttri
   c1 = v1/E_T1
   c2 = v3/E_T2
   c3 = v2/E_T3
   ctri = vtri/E_Ttri
   p1_dual = pchisq( c1*T1, df=v1, lower.tail=FALSE)
   p2_dual = pchisq( c2*T2, df=v2, lower.tail=FALSE)
   p3_dual = pchisq( c3*T3, df=v3, lower.tail=FALSE)
   p_tri = pchisq( ctri*Ttri, df=vtri, lower.tail=FALSE)
   min_pert = apply(cbind(p1_single,p2_single,p3_single,p1_dual,p2_dual,p3_dual,p_tri),1,min)
  
   pval_pert = sum(min_pert<min0)/(length(min_pert)+1)
   if (is.na(pval_pert)) pval_pert=1
   if (is.na(p0_tri)) p0_tri=1
   c(pval_pert, p0_tri, p01_single,p02_single,p03_single)
}
