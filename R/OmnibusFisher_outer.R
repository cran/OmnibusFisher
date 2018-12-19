#' A modified Fisher’s method (Omnibus-Fisher) to combine separate p-values of SNPs, RNA expressions and DNA methylations into an overall gene-level p-value
#'
#' * Each sample does not have to have all three types of omics data; each gene needs to have all three types of omics data mapped to it. \cr
#' * For example, 1,000 samples have SNPs mapped to 20,000 genes; 500 samples have methylated sites mapped to 18,000 genes; 300 samples have 16,000 expressed genes. \cr
#' * Then, all 1,000 samples (500 and 300 samples are included in the total 1,000 samples) contribute to the test. We are interested in the overlapped genes (e.g., 16,000 genes with SNPs, methylated sites and expression genes mapped to them)
#'
#' @name OmnibusFisher
#' @aliases OmnibusFisher
#' @param pheno A matrix of sample ID, trait (i.e., y) and covariates (class: data.frame).
#' @param full_id A vector of sample ID. This vector should include all IDs. In other words, samples with all 3 types of omics data, with 2 types and with 1 type should have their IDs in the vector (class: data.frame).
#' @param G A matrix of genotypes in a gene. The 1st column is sample ID. Each column is a SNP in the following columns (class: data.frame).
#' @param M A matrix of methylated sites in a gene. The 1st column is sample ID. Each column is a methylated site in the following columns (class: data.frame).
#' @param R A matrix of RNA expression probes in a gene (for microarray, one gene could have multiple probes mapped; for RNAseq, one gene always has one value). The 1st column is sample ID. For microarray, each column is a probe in the following columns; for RNAseq, the 2nd column is the expression value (only two columns) (class: data.frame).
#' @param exprs_G Regression model for SNPs under the null hypothesis (i.e., SNP effect is zero), y = cov1 + cov2 + ... + covp.
#' @param exprs_M Regression model for DNA methylation under the null hypothesis (i.e., methylation effect is zero), y = cov1 + cov2 + ... + covp.
#' @param exprs_R Regression model for RNA expression under the null hypothesis (i.e., RNA expression effect is zero), y = cov1 + cov2 + ... + covp.
#' @param type Either type="binary" or type="continuous" for binary or continuous traits.
#' @param method Method used to approximately calculate p-values: "kuonen" or "davies". Default "kuonen".
#' @param optimal Whether use optimal method to automatically search for the disease model (perturbation invovled).
#' @param perturb_iteration The number of perturbation iterations when using optimal method. For example, 1,000,000, then the lowest p-value can be obtained for perturbation method is 1/1,000,001.
#' @import CompQuadForm
#' @import stringr
#' @importFrom survey pchisqsum
#' @importFrom stats glm
#' @importFrom stats lm
#' @importFrom stats as.formula
#' @importFrom stats binomial
#' @include OmnibusFisher_iteration_GMR.R
#' @include OmnibusFisher_iteration_GM.R
#' @include OmnibusFisher_iteration_GR.R
#' @include OmnibusFisher_iteration_MR.R
#' @return 1. pval_GMR_pert:   the overall gene-level p-value automatically searching for the optimal disease model, when inputting three types of data. \cr
#'         2. pval_GMR_tri:   the overall gene-level p-value assuming all three types of data in the disease model, when inputting three types of data. \cr
#'         3. pval_GM_tri/pval_GR_tri/pval_MR_tri:   the overall gene-level p-value automatically searching for the optimal disease model, when inputting two types of data. \cr
#'         4. pval_GM_tri/pval_GR_tri/pval_MR_tri:   the overall gene-level p-value assuming two types of data in the disease model, when inputting two types of data. \cr
#'         5. pval_G/pval_M/pval_R:   the gene-level p-value for single type of data. \cr
#' @examples
#' ################
#' ### Examples ###
#' ################
#' data("example_data")
#' set.seed(123)
#' exprs_G = exprs_M = exprs_R = "aff ~ age + sex"
#' ### SNPs (G), DNA methylations (M) and RNA expressions (R) ###
#' results<-list()
#' for(i in 1:1){ #change to 1:3 for 3 genes
#'  results[[i]]<-OmnibusFisher(pheno=pheno, full_id=All_header, G=G[[i]], M=M[[i]],
#'  R=R[[i]], exprs_G=exprs_G, exprs_M=exprs_M, exprs_R=exprs_R, type="binary")
#'  # G[[1]] includes SNPs in gene1;
#'  # M[[1]] includes methylated sites in gene1;
#'  # R[[1]] includes gene expression probes in gene1 (or single gene1 expression value).
#' }
#' ### SNPs (G) and DNA methylations (M) ###
#' results<-list()
#' for(i in 1:1){ 
#'  results[[i]]<-OmnibusFisher(pheno=pheno, full_id=All_header, G=G[[i]], M=M[[i]],
#'  exprs_G=exprs_G, exprs_M=exprs_M, type="binary")
#' }
#' ### SNPs (G) and RNA expressions (R) ###
#' # results[[i]]<-OmnibusFisher(pheno=pheno, full_id=All_header, G=G[[i]], R=R[[i]],
#' # exprs_G=exprs_G, exprs_R=exprs_R, type="binary")
#' ### DNA methylations (M) and RNA expressions (R) ###
#' # results[[i]]<-OmnibusFisher(pheno=pheno, full_id=All_header, R=R[[i]], M=M[[i]],
#' # exprs_R=exprs_R, exprs_M=exprs_M, type="binary")
#' @export
OmnibusFisher <- function(pheno, full_id, G=NULL, M=NULL, R=NULL, exprs_G=NULL, exprs_M=NULL, exprs_R=NULL, type, optimal=FALSE, perturb_iteration=NULL, method="kuonen"){

  # SNP
  if(!is.null(G)){
   G_all = merge(G, pheno, by="ID")
   G_all$order = seq(1:dim(G_all)[1])
   exprs_G = str_replace_all(exprs_G, fixed(" "), "")
   cov_G = strsplit(exprs_G, "~")[[1]][2]
   cov_G = strsplit(cov_G, "\\+")
   X1 <- as.matrix(cbind(rep(1,dim(G_all)[1]), G_all[,cov_G[[1]]]))
   G_snp = G_all[,2:dim(G)[2]]
   G_snp = as.matrix(G_snp)
   W1 = diag(1, dim(G_snp)[2])
   y_G = strsplit(exprs_G, "~")[[1]][1]
   if(type=="binary"){
    model1 <- glm(as.formula(exprs_G), family=binomial(), data=G_all)
    beta1 <- summary(model1)$coefficients[,1]
    xbeta1 <- X1%*%beta1
    mu1 <- exp(xbeta1)/(1+exp(xbeta1))
    V1 <- diag(as.matrix(mu1*(1-mu1))[,1])
    res1 = as.matrix(G_all[,y_G] - mu1)
   }
   if(type=="continuous"){
    model1 <- lm(as.formula(exprs_G), data=G_all)
    beta1 <- summary(model1)$coefficients[,1]
    xbeta1 <- X1%*%beta1
    sigma1 <- (summary(model1)$sigma)
    V1 <- diag(sigma1**2, dim(X1)[1])
    res1 = as.matrix(G_all[,y_G] - xbeta1)
   }
   G_snp <- as.matrix(G_snp)
   Q1 = t(res1) %*% G_snp %*% W1 %*% t(G_snp) %*% res1
   P1 = V1 - V1 %*% X1 %*% solve(t(X1) %*% V1 %*% X1) %*% t(X1) %*% V1
   svdob1<-eigen(P1, symmetric=TRUE)
   P1_0.5<-svdob1$vectors%*%diag(sqrt(round(svdob1$values, 6)))%*%t(svdob1$vectors)
   eig1 = eigen(P1_0.5 %*% G_snp %*% W1 %*% t(G_snp) %*% P1_0.5, symmetric=T)
   evals1 = eig1$values[eig1$values>1e-6*eig1$values[1]]
   if (method=="davies"){
     power=3
     while(power<23){
       acc=1/(10^power)
       tmpout1<-davies(Q1, evals1, acc=acc, lim=1e+09)
       pskat1<-tmpout1$Qq
       power=power+1
       if (pskat1>0 & pskat1<1) {break}}
   }
   else if (method=="kuonen"){
     pskat1<-pchisqsum(Q1,rep(1,length(evals1)),evals1,lower.tail=F,method="sad")
     pskat1<-min(pskat1, 1)
     pskat1<-max(pskat1, 0)
   }
   evectors1 = eig1$vectors[,1:length(evals1)]
  }
  # Methylation
  if(!is.null(M)){
   M_all = merge(M, pheno, by="ID")
   M_all$order = seq(1:dim(M_all)[1])
   exprs_M = str_replace_all(exprs_M, fixed(" "), "")
   cov_M = strsplit(exprs_M, "~")[[1]][2]
   cov_M = strsplit(cov_M, "\\+")
   X2 <- as.matrix(cbind(rep(1,dim(M_all)[1]), M_all[,cov_M[[1]]]))
   M_loc = M_all[,2:dim(M)[2]]
   M_loc = as.matrix(M_loc)
   W2 = diag(1, dim(M_loc)[2])
   y_M = strsplit(exprs_M, "~")[[1]][1]
   if(type=="binary"){
    model2 <- glm(as.formula(exprs_M), family=binomial(), data=M_all)
    beta2 <- summary(model2)$coefficients[,1]
    xbeta2 <- X2%*%beta2
    mu2 <- exp(xbeta2)/(1+exp(xbeta2))
    V2 <- diag(as.matrix(mu2*(1-mu2))[,1])
    res2 = as.matrix(M_all[,y_M] - mu2)
   }
   if(type=="continuous"){
    model2 <- lm(as.formula(exprs_M), data=M_all)
    beta2 <- summary(model2)$coefficients[,1]
    xbeta2 <- X2%*%beta2
    sigma2 <- (summary(model2)$sigma)
    V2 <- diag(sigma2**2, dim(X2)[1])
    res2 = as.matrix(M_all[,y_M] - xbeta2)
   }
   M_loc <- as.matrix(M_loc)
   Q2 = t(res2) %*% M_loc %*% W2 %*% t(M_loc) %*% res2
   P2 = V2 - V2 %*% X2 %*% solve(t(X2) %*% V2 %*% X2) %*% t(X2) %*% V2
   svdob2<-eigen(P2, symmetric=TRUE)
   P2_0.5<-svdob2$vectors%*%diag(sqrt(round(svdob2$values, 6)))%*%t(svdob2$vectors)
   eig2 = eigen(P2_0.5 %*% M_loc %*% W2 %*% t(M_loc) %*% P2_0.5, symmetric=T)
   evals2 = eig2$values[eig2$values>1e-6*eig2$values[1]]
   if (method=="davies"){
     power=3
     while(power<23){
       acc=1/(10^power)
       tmpout2<-davies(Q2, evals2, acc=acc, lim=1e+09)
       pskat2<-tmpout2$Qq
       power=power+1
       if (pskat2>0 & pskat2<1) {break}}
   }
   else if (method=="kuonen"){
     pskat2<-pchisqsum(Q2,rep(1,length(evals2)),evals2,lower.tail=F,method="sad")
     pskat2<-min(pskat2, 1)
     pskat2<-max(pskat2, 0)
   }
   evectors2 = eig2$vectors[,1:length(evals2)]
  }
  # RNA
  if(!is.null(R)){
   R_all = merge(R, pheno, by="ID")
   R_all$order = seq(1:dim(R_all)[1])
   exprs_R = str_replace_all(exprs_R, fixed(" "), "")
   cov_R = strsplit(exprs_R, "~")[[1]][2]
   cov_R = strsplit(cov_R, "\\+")
   X3 <- as.matrix(cbind(rep(1,dim(R_all)[1]), R_all[,cov_R[[1]]]))
   R_snp = R_all[,2:dim(R)[2]]
   R_snp = as.matrix(R_snp)
   W3 = diag(1, dim(R_snp)[2])
   y_R = strsplit(exprs_R, "~")[[1]][1]
   if(type=="binary"){
    model3 <- glm(as.formula(exprs_R), family=binomial(), data=R_all)
    beta3 <- summary(model3)$coefficients[,1]
    xbeta3 <- X3%*%beta3
    mu3 <- exp(xbeta3)/(1+exp(xbeta3))
    V3 <- diag(as.matrix(mu3*(1-mu3))[,1])
    res3 = as.matrix(R_all[,y_R] - mu3)
   }
   if(type=="continuous"){
    model3 <- lm(as.formula(exprs_R), data=R_all)
    beta3 <- summary(model3)$coefficients[,1]
    xbeta3 <- X3%*%beta3
    sigma3 <- (summary(model3)$sigma)
    V3 <- diag(sigma3**2, dim(X3)[1])
    res3 = as.matrix(R_all[,y_R] - xbeta3)
   }
   R_snp <- as.matrix(R_snp)
   Q3 = t(res3) %*% R_snp %*% W3 %*% t(R_snp) %*% res3
   P3 = V3 - V3 %*% X3 %*% solve(t(X3) %*% V3 %*% X3) %*% t(X3) %*% V3
   svdob3<-eigen(P3, symmetric=TRUE)
   P3_0.5<-svdob3$vectors%*%diag(sqrt(round(svdob3$values, 6)))%*%t(svdob3$vectors)
   eig3 = eigen(P3_0.5 %*% R_snp %*% W3 %*% t(R_snp) %*% P3_0.5, symmetric=T)
   evals3 = eig3$values[eig3$values>1e-6*eig3$values[1]]
   if (method=="davies"){
     power=3
     while(power<23){
       acc=1/(10^power)
       tmpout3<-davies(Q3, evals3, acc=acc, lim=1e+09)
       pskat3<-tmpout3$Qq
       power=power+1
       if (pskat3>0 & pskat3<1) {break}}
   }
   else if (method=="kuonen"){
     pskat3<-pchisqsum(Q3,rep(1,length(evals3)),evals3,lower.tail=F,method="sad")
     pskat3<-min(pskat3, 1)
     pskat3<-max(pskat3, 0)
   }
   evectors3 = eig3$vectors[,1:length(evals3)]
  }

 if(!is.null(G) & !is.null(M) & !is.null(R)){
   It<-rep(9, 5)
   if (optimal & is.null(perturb_iteration)){
   tryCatch(It<-Iteration_GMR(num=200, full_id, G_all, M_all, R_all, evectors1, evectors2, evectors3, evals1, evals2, evals3, pskat1, pskat2, pskat3, acc=1e-4, method=method), error=function(e){cat("ERROR", "\n") })
     num2=round(1/min(It[-1]))*10
     if (num2>200 & num2<(1e6+1))  tryCatch(It<-Iteration_GMR(num=num2, full_id, G_all, M_all, R_all, evectors1, evectors2, evectors3, evals1, evals2, evals3, pskat1, pskat2, pskat3, acc=1e-4, method=method), error=function(e){cat("ERROR", "\n") })
     else if (num2>1e+6) {
       num2=1e+6
       tryCatch(It<-Iteration_GMR(num=num2, full_id, G_all, M_all, R_all, evectors1, evectors2, evectors3, evals1, evals2, evals3, pskat1, pskat2, pskat3, acc=1e-4, method=method), error=function(e){cat("ERROR", "\n") })
     }
   names(It)=c("pval_GMR_pert", "pval_GMR_tri", "pval_G", "pval_M", "pval_R")
   }
   else if (optimal & !is.null(perturb_iteration)){
     tryCatch(It<-Iteration_GMR(num=perturb_iteration, full_id, G_all, M_all, R_all, evectors1, evectors2, evectors3, evals1, evals2, evals3, pskat1, pskat2, pskat3, acc=1e-4, method=method), error=function(e){cat("ERROR", "\n") })
   names(It)=c("pval_GMR_pert", "pval_GMR_tri", "pval_G", "pval_M", "pval_R")
   }
   else if (!optimal){
   tryCatch(It<-Iteration_GMR(num=200, full_id, G_all, M_all, R_all, evectors1, evectors2, evectors3, evals1, evals2, evals3, pskat1, pskat2, pskat3, acc=1e-4, method=method), error=function(e){cat("ERROR", "\n") })
   It<-It[-1]
   names(It)=c("pval_GMR_tri", "pval_G", "pval_M", "pval_R")
   }
 }

 if(!is.null(G) & !is.null(M) & is.null(R)){
   It<-rep(9, 4)
   if (optimal & is.null(perturb_iteration)){
   tryCatch(It<-Iteration_GM(num=200, full_id, G_all, M_all, evectors1, evectors2, evals1, evals2, pskat1, pskat2, acc=1e-4, method=method), error=function(e){cat("ERROR", "\n") })
     num2=round(1/min(It[-1]))*10
     if (num2>200 & num2<(1e6+1))  tryCatch(It<-Iteration_GM(num=num2, full_id, G_all, M_all, evectors1, evectors2, evals1, evals2, pskat1, pskat2, acc=1e-4, method=method), error=function(e){cat("ERROR", "\n") })
     else if (num2>1e+6) {
       num2=1e+6
       tryCatch(It<-Iteration_GM(num=num2, full_id, G_all, M_all, evectors1, evectors2, evals1, evals2, pskat1, pskat2, acc=1e-4, method=method), error=function(e){cat("ERROR", "\n") })
     }
   names(It)=c("pval_GM_pert", "pval_GM_dual", "pval_G", "pval_M")
   }
   else if (optimal & !is.null(perturb_iteration)){
     tryCatch(It<-Iteration_GM(num=perturb_iteration, full_id, G_all, M_all, evectors1, evectors2, evals1, evals2, pskat1, pskat2, acc=1e-4, method=method), error=function(e){cat("ERROR", "\n") })
   names(It)=c("pval_GM_pert", "pval_GM_dual", "pval_G", "pval_M")
   }
   else if (!optimal){
   tryCatch(It<-Iteration_GM(num=200, full_id, G_all, M_all, evectors1, evectors2, evals1, evals2, pskat1, pskat2, acc=1e-4, method=method), error=function(e){cat("ERROR", "\n") })
   It<-It[-1]
   names(It)=c("pval_GM_dual", "pval_G", "pval_M")
   }
 }

 if(!is.null(G) & is.null(M) & !is.null(R)){
   It<-rep(9, 4)
   if (optimal & is.null(perturb_iteration)){
   tryCatch(It<-Iteration_GR(num=200, full_id, G_all, R_all, evectors1, evectors3, evals1, evals3, pskat1, pskat3, acc=1e-4, method=method), error=function(e){cat("ERROR", "\n") })
     num2=round(1/min(It[-1]))*10
     if (num2>200 & num2<(1e6+1))  tryCatch(It<-Iteration_GR(num=num2, full_id, G_all, R_all, evectors1, evectors3, evals1, evals3, pskat1, pskat3, acc=1e-4, method=method), error=function(e){cat("ERROR", "\n") })
     else if (num2>1e+6) {
       num2=1e+6
       tryCatch(It<-Iteration_GR(num=num2, full_id, G_all, R_all, evectors1, evectors3, evals1, evals3, pskat1, pskat3, acc=1e-4, method=method), error=function(e){cat("ERROR", "\n") })
     }
   names(It)=c("pval_GR_pert", "pval_GR_dual", "pval_G", "pval_R")
   }
   else if (optimal & !is.null(perturb_iteration)){
     tryCatch(It<-Iteration_GR(num=perturb_iteration, full_id, G_all, R_all, evectors1, evectors3, evals1, evals3, pskat1, pskat3, acc=1e-4, method=method), error=function(e){cat("ERROR", "\n") })
   names(It)=c("pval_GR_pert", "pval_GR_dual", "pval_G", "pval_R")
   }
   else if (!optimal){
   tryCatch(It<-Iteration_GR(num=200, full_id, G_all, R_all, evectors1, evectors3, evals1, evals3, pskat1, pskat3, acc=1e-4, method=method), error=function(e){cat("ERROR", "\n") })
   It<-It[-1]
   names(It)=c("pval_GR_dual", "pval_G", "pval_R")
   }
 }

 if(is.null(G) & !is.null(M) & !is.null(R)){
   It<-rep(9, 4)
   if (optimal & is.null(perturb_iteration)){
   tryCatch(It<-Iteration_MR(num=200, full_id, M_all, R_all, evectors2, evectors3, evals2, evals3, pskat2, pskat3, acc=1e-4, method=method), error=function(e){cat("ERROR", "\n") })
     num2=round(1/min(It[-1]))*10
     if (num2>200 & num2<(1e6+1))  tryCatch(It<-Iteration_MR(num=num2, full_id, M_all, R_all, evectors2, evectors3, evals2, evals3, pskat2, pskat3, acc=1e-4, method=method), error=function(e){cat("ERROR", "\n") })
     else if (num2>1e+6) {
       num2=1e+6
       tryCatch(It<-Iteration_MR(num=num2, full_id, M_all, R_all, evectors2, evectors3, evals2, evals3, pskat2, pskat3, acc=1e-4, method=method), error=function(e){cat("ERROR", "\n") })
     }
   names(It)=c("pval_MR_pert", "pval_MR_dual", "pval_M", "pval_R")
   }
   else if (optimal & !is.null(perturb_iteration)){
     tryCatch(It<-Iteration_MR(num=perturb_iteration, full_id, M_all, R_all, evectors2, evectors3, evals2, evals3, pskat2, pskat3, acc=1e-4, method=method), error=function(e){cat("ERROR", "\n") })
   names(It)=c("pval_MR_pert", "pval_MR_dual", "pval_M", "pval_R")
   }
   else if (!optimal){
   tryCatch(It<-Iteration_MR(num=200, full_id, M_all, R_all, evectors2, evectors3, evals2, evals3, pskat2, pskat3, acc=1e-4, method=method), error=function(e){cat("ERROR", "\n") })
   It<-It[-1]
   names(It)=c("pval_MR_dual", "pval_M", "pval_R")
   }
 }

return(It)
}
