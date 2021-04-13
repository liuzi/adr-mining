library(Rlda)

pres_matrix = read.csv(
    "/data/liu/mimic3/CLAMP_NER/single_drug_analysis/FEATURE/PRE_PROCESS/pres_rxnorm_matrix.csv",
    row.names=1)
print("pres_rxnorm matrix:\n")
print(dim(pres_matrix))
diag_matrix = read.csv(
    "/data/liu/mimic3/CLAMP_NER/single_drug_analysis/FEATURE/PRE_PROCESS/diag_matrix.csv",
    row.names=1)
print("diag matrix:\n")
print(dim(diag_matrix))


res2file <- function(res, n_gibbs, n_community, gamma=0.01,alpha0=0.01,alpha1=0.01,filename = 'pres'){
    # prefix = "/data/liu/LDA/lda_R_result"
    prefix = "/data/liu/mimic3/LDA_MODEL/JOINT_LDA"
    args = paste(c("ngib",n_gibbs,"_ncomp",n_community,"_gama",gamma,"_alpha",alpha0), collapse = '')
    file_path = file.path(prefix,args)
    dir.create(file_path, showWarnings = FALSE)

    theta_df  = getTheta.rlda(res)
    phi_df = getPhi.rlda(res)
    write.csv(theta_df, file = file.path(file_path, paste(filename,"_theta.csv",sep = '')),row.names=FALSE)
    write.csv(phi_df, file = file.path(file_path, paste(filename,"_phi.csv",sep = '')),row.names=FALSE)

    pdf(file = file.path(file_path, paste(filename,"_rlda.pdf",sep = '')))
    plot(res)
    dev.off()

}

runlda <- function(data,filename,n_gibbs, n_community,gamma=0.01,alpha0=0.01,alpha1=0.01){
    # Set seed
    set.seed(2021)
    # Hyperparameters for each prior distribution
#     gamma <- 0.01
#     alpha0 <- 0.01
#     alpha1 <- 0.01
    # Execute the LDA for the Bernoulli entry
    res <- rlda.bernoulli(data = data, n_community = n_community,
                          alpha0 = alpha0, alpha1 = alpha1, gamma = gamma,
                          n_gibbs = n_gibbs,ll_prior = TRUE, display_progress = TRUE)
   
    res2file(res, n_gibbs, n_community, gamma, alpha0, alpha1,filename)
   
    return(res)
}

gamma=0.01
alpha0=0.01
alpha1=0.01

res = runlda(pres_matrix,"pres_rxnorm_Full",700,10,gamma,alpha0,alpha1)
res = runlda(diag_matrix,"diag_Full",700,10,gamma,alpha0,alpha1)