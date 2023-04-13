opt_simulation_3region <- function(k) {
  cat(paste0("Begin sim:", k, "\n"))

  Y <- simulated_data[[k]]
  L1 <- L
  L2 <- L
  L3 <- L
  M <- M
  
  X_Region1 <- Y$X_Region1
  X_Region2 <- Y$X_Region2
  X_Region3 <- Y$X_Region3
  
  # Correlation of Averages from Chao's code
  # Vanilla
  X_R1_mat <- matrix(X_Region1, nrow=M)
  X_R2_mat <- matrix(X_Region2, nrow=M)
  X_R3_mat <- matrix(X_Region3, nrow=M)
  
  X_R1_avg <- apply(X_R1_mat, 1, mean)
  X_R2_avg <- apply(X_R2_mat, 1, mean)
  X_R3_avg <- apply(X_R3_mat, 1, mean)
  
  results_pearson <- data.frame(cor_vanilla_12 = cor(X_R1_avg, X_R2_avg),
                                cor_vanilla_13 = cor(X_R1_avg, X_R3_avg),
                                cor_vanilla_23 = cor(X_R2_avg, X_R3_avg))
  
  # Bspline
  X_R1_df <- data.frame(signal=X_Region1, time=rep(1:M, times=L1))
  X_R2_df <- data.frame(signal=X_Region2, time=rep(1:M, times=L2))
  X_R3_df <- data.frame(signal=X_Region3, time=rep(1:M, times=L2))
  
  fitR1 <- lm(signal ~  -1 + bs(time, df=K, intercept=T), data=X_R1_df)
  fitR2 <- lm(signal ~  -1 + bs(time, df=K, intercept=T), data=X_R2_df)
  fitR3 <- lm(signal ~  -1 + bs(time, df=K, intercept=T), data=X_R3_df)
  
  mean_plus_eta_R1 <- predict(fitR1, data.frame(time=1:M))
  mean_plus_eta_R2 <- predict(fitR2, data.frame(time=1:M))
  mean_plus_eta_R3 <- predict(fitR3, data.frame(time=1:M))
  
  results_pearson$cor_bspline_12 <-  cor(mean_plus_eta_R1, mean_plus_eta_R2)
  results_pearson$cor_bspline_13 <-  cor(mean_plus_eta_R1, mean_plus_eta_R3)
  results_pearson$cor_bspline_23 <-  cor(mean_plus_eta_R2, mean_plus_eta_R3)
  
  
  # Intra-regional model
  ## Initilize fixed effect parameters
  G_region <- bs(rep(1:M, L1), df=K, intercept=T)
  ## Optimization
  ### Region 1
  parametersInit_region1 <-c(0, 0, 0, coef(fitR1))
  result_region1 <- opt_intra(parametersInit_region1,
                                              X_Region1, G_region,
                                              dist_sqrd_mat_region1, time_sqrd_mat,
                                              L1, M, kernel_type)
  
  ### Region 1 results
  region1_para <- result_region1$theta
  nu_1 <- result_region1$nu
  eta_1 <- c(G_region[1:M,] %*% nu_1)
  var_eta_1 <- var(eta_1)
  ############
  
  ### Region 2
  parametersInit_region2 <- c(0, 0, 0, coef(fitR2))
  result_region2 <- opt_intra(parametersInit_region2,
                                              X_Region2, G_region,
                                              dist_sqrd_mat_region2, time_sqrd_mat,
                                              L2, M, kernel_type)
  
  #Region 2 results
  region2_para <- result_region2$theta
  nu_2 = result_region2$nu
  eta_2 <- c(G_region[1:M,] %*% nu_2)
  var_eta_2 <- var(eta_2)
  ############
  
  ### Region 3
  parametersInit_region3 <- c(0, 0, 0, coef(fitR3))
  result_region3 <- opt_intra(parametersInit_region3,
                                              X_Region3, G_region,
                                              dist_sqrd_mat_region3, time_sqrd_mat,
                                              L3, M, kernel_type)
  
  #Region 2 results
  region3_para <- result_region3$theta
  nu_3 = result_region3$nu
  eta_3 <- c(G_region[1:M,] %*% nu_3)
  var_eta_3 <- var(eta_3)
  ############
  
  # Cor Intra
  results_pearson$cor_intra_12 <-  cor(eta_1, eta_2)
  results_pearson$cor_intra_13 <-  cor(eta_1, eta_3)
  results_pearson$cor_intra_23 <-  cor(eta_2, eta_3)
  
  
  # Inter-regional model
  Z <- matrix(c(rep(c(1,0), each = (L1*M)),
                rep(c(0,1), each = (L2*M))), (L1+L2)*M, 2)
  
  ## Pair 1-2
  gamma_vec_12 <- c(region1_para, region2_para)
  parametersInit_pair_12 <- c(0, 0, 0, -2, mean(X_R1_avg), mean(X_R2_avg))

  result_pair_12 <- opt_inter(parametersInit_pair_12,
                                      matrix(c(X_Region1, X_Region2), ncol=1), 
                                      Z,
                                      dist_sqrd_mat_region1,
                                      dist_sqrd_mat_region2,
                                      time_sqrd_mat, gamma_vec_12,
                                      kernel_type)
  
  result_pair_12_theta <- result_pair_12$theta
  cat(paste0("REML rho_12: ", result_pair_12_theta[1]))
  cat("\n")
  
  ## Pair 1-3
  gamma_vec_13 <- c(region1_para, region3_para)
  parametersInit_pair_13 <- c(0, 0, 0, -2, mean(X_R1_avg), mean(X_R3_avg))
  result_pair_13 <- opt_inter(parametersInit_pair_13,
                                      matrix(c(X_Region1, X_Region3), ncol=1), 
                                      Z,
                                      dist_sqrd_mat_region1,
                                      dist_sqrd_mat_region3,
                                      time_sqrd_mat, gamma_vec_13,
                                      kernel_type)
  
  result_pair_13_theta <- result_pair_13$theta
  cat(paste0("REML rho_13: ", result_pair_13_theta[1]))
  cat("\n")
  
  ## Pair 1-3
  gamma_vec_23 <- c(region2_para, region3_para)
  parametersInit_pair_23 <- c(0, 0, 0, -2, mean(X_R2_avg), mean(X_R3_avg))
  result_pair_23 <- opt_inter(parametersInit_pair_23,
                                      matrix(c(X_Region2, X_Region3), ncol=1), 
                                      Z,
                                      dist_sqrd_mat_region2,
                                      dist_sqrd_mat_region3,
                                      time_sqrd_mat, gamma_vec_23,
                                      kernel_type)
  
  result_pair_23_theta <- result_pair_23$theta
  cat(paste0("REML rho_23: ", result_pair_23_theta[1]))
  cat("\n")
  
  
  
  rho_inter <- data.frame(rho_12=result_pair_12_theta[1], 
                          rho_13=result_pair_13_theta[1], 
                          rho_23=result_pair_23_theta[1])
  asymp_var_inter <- data.frame(var_12=result_pair_12$asymptotic_var, 
                                var_13=result_pair_13$asymptotic_var, 
                                var_23=result_pair_23$asymptotic_var)
  inter_para <- rbind(result_pair_12_theta[2:4],
                      result_pair_13_theta[2:4],
                      result_pair_23_theta[2:4])
  
  colnames(inter_para) <- c("tauEta", "kEta", "nugget")

  cat(paste0("Done sim:", k, "\n"))

  list(rho_pearson=results_pearson, 
              rho_inter=rho_inter,
              asymp_var_inter=asymp_var_inter,
              result_pair_12=result_pair_12,
              result_pair_13=result_pair_13,
              result_pair_23=result_pair_23)
}
##########################################################################################