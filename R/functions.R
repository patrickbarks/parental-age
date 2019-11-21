



#### Gam models to smooth transition rates ------------------------------------
GamFnSurv <- function(data, plot = FALSE) {
  # data:  flat df with columns for individual survival (surv) and age (age)
  # plot:  if TRUE use high-resolution sequence of ages (for nicer plotting)
  mod <- gam(surv ~ s(age, k = 4), family = 'binomial', data = data)
  if (plot == TRUE) {
    age <- seq(1, o, length.out = 100)
  } else {
    age <- 1:o
  }
  pred <- predict(mod, newdata = list(age = age), type = 'response')
  return(data.frame(i = age, P_ij = pred))
}

GamFnFert <- function(data, plot = FALSE) {
  # data:  flat df with columns for individual identity (ident.focal), age-specific
  #   fecundity (fecund), and age (age)
  # plot:  if TRUE use high-resolution sequence of ages (for nicer plotting)
  mod <- gamm(fecund ~ s(age), family = 'poisson', data = data,
              random = list(ident.focal = ~1))
  if (plot == TRUE) {
    age <- seq(1, o, length.out = 100)
  } else {
    age <- 1:o
  }
  pred <- predict(mod$gam, newdata = list(age = age), type = 'response')
  return(data.frame(i = age, F_ij = pred))
}




#### Split reproductive age classes into equally-spaced parental age classes --
map_ij <- function(F_i, s) {
  # F_i:  age-specific fecundity vector
  # s:  number of parental age classes
  rep_stage <- F_i > 0
  n_prerep <- length(which(rep_stage == FALSE))
  
  j <- c(
    rep(1L, n_prerep),
    as.integer(as.factor(cut(seq_along(which(rep_stage)), s)))
  )
  
  return(j)
}




#### Construct A_par -----------------------------------------------------------
# vec-permutation; from Caswell (2012) Theor Ecol 5:403–417, Equation 15
vec_perm <- function(o, s) {
  # o:  number of age classes (omega)
  # s:  number of parental age classes
  K <- matrix(0, nrow = o*s, ncol = o*s)
  
  for(i in 1:s) {
    for(j in 1:o) {
      E <- matrix(0, nrow = s, ncol = o)
      E[i,j] <- 1
      K <- K + (E %x% t(E))
    }
  }
  
  return(K)
}


# construct matrix A_par
build_Apar <- function(Tr, q_i, stasis = FALSE, checks = TRUE) {
  # Tr:   flat df with transition rates (P_ij and F_ij) by age (i) and parental age (j)
  # q_i:  vector giving mapping between age classes (indices) and parental age classes (values)
  # stasis: if TRUE, include stasis loop at maximum age class omega
  # checks: if TRUE, check whether A_par ergodic and irreducible
  o <- length(unique(Tr$i))
  s <- length(unique(Tr$j))
  
  P_mat <- arrange(Tr, j, i)$P_ij %>% matrix(ncol = s)
  F_mat <- arrange(Tr, j, i)$F_ij %>% matrix(ncol = s)
  
  if (stasis == FALSE) P_mat[o,] <- 0
  
  # create stage-classified projection matrices for each age class
  F_a <- list()
  P_a <- list()
  
  P_block <- matrix(0, nrow = o*s, ncol = o*s)
  F_block <- matrix(0, nrow = o*s, ncol = o*s)
  
  for(i in 1:o) {
    P_a[[i]] <- diag(P_mat[i,])
    
    F_a[[i]] <- matrix(0, nrow = s, ncol = s)
    F_a[[i]][q_i[i],] <- F_mat[i,]
    
    E <- matrix(0, nrow = o, ncol = o)
    E[i,i] <- 1
    
    P_block <- P_block + (E %x% P_a[[i]])
    F_block <- F_block + (E %x% F_a[[i]])
  }
  
  # create age-transition matrices
  D_P <- matrix(0, nrow = o, ncol = o)
  D_P[2:o, 1:(o-1)] <- diag(o-1)
  if (stasis == TRUE) D_P[o,o] <- 1
  D_P_block <- diag(s) %x% D_P
  
  D_F <- matrix(0, nrow = o, ncol = o)
  D_F[1,] <- rep(1, o)
  D_F_block <- diag(s) %x% D_F
  
  # create vec-permutation  matrix
  K <- vec_perm(o, s)
  
  # convert to sparse
  K <- as(K, "sparseMatrix")
  D_P_block <- as(D_P_block, "sparseMatrix")
  D_F_block <- as(D_F_block, "sparseMatrix")
  P_block <- as(P_block, "sparseMatrix")
  F_block <- as(F_block, "sparseMatrix")
  
  # projection matrix structured by both age and parental age
  A_par <- (t(K) %*% D_P_block %*% K %*% P_block) +
    (t(K) %*% D_F_block %*% K %*% F_block)
  
  # check whether A_par ergodic and irreducible
  if (checks == TRUE) {
    A_par_check <- as.matrix(A_par)
    check1 <- popdemo::isErgodic(A_par_check)
    check2 <- popdemo::isIrreducible(A_par_check)
    if (check1 == FALSE) warning('A_par is non-ergodic')
    if (check2 == FALSE) warning('A_par is reducible')
  }
  
  return(list(o = o, s = s, q_i = q_i, P_mat = P_mat, F_mat = F_mat,
              D_P_block = D_P_block, D_F_block = D_F_block,
              P_block = P_block, F_block = F_block,
              K = K, A_par = A_par))
}


# construct matrix A_par (light version, used in optimization to set lambda) ---
build_Apar_lite <- function(Tr, q_i, stasis = FALSE) {
  # Tr:   flat df with transition rates (P_ij and F_ij) by age (i) and parental age (j)
  # q_i:  vector giving mapping between age classes (indices) and parental age classes (values)
  # stasis: if TRUE, include stasis loop at maximum age class omega
  o <- length(unique(Tr$i))
  s <- length(unique(Tr$j))
  
  P_mat <- arrange(Tr, j, i)$P_ij %>% matrix(ncol = s)
  F_mat <- arrange(Tr, j, i)$F_ij %>% matrix(ncol = s)
  
  if (stasis == FALSE) P_mat[o,] <- 0
  
  # create stage-classified projection matrices for each age class
  F_a <- list()
  P_a <- list()
  
  P_block <- matrix(0, nrow = o*s, ncol = o*s)
  F_block <- matrix(0, nrow = o*s, ncol = o*s)
  
  for(i in 1:o) {
    P_a[[i]] <- diag(P_mat[i,])
    
    F_a[[i]] <- matrix(0, nrow = s, ncol = s)
    F_a[[i]][q_i[i],] <- F_mat[i,]
    
    E <- matrix(0, nrow = o, ncol = o)
    E[i,i] <- 1
    
    P_block <- P_block + (E %x% P_a[[i]])
    F_block <- F_block + (E %x% F_a[[i]])
  }
  
  # create age-transition matrices
  D_P <- matrix(0, nrow = o, ncol = o)
  D_P[2:o, 1:(o-1)] <- diag(o-1)
  if (stasis == TRUE) D_P[o,o] <- 1
  D_P_block <- diag(s) %x% D_P
  
  D_F <- matrix(0, nrow = o, ncol = o)
  D_F[1,] <- rep(1, o)
  D_F_block <- diag(s) %x% D_F
  
  # create vec-permutation  matrix
  K <- vec_perm(o, s)
  
  # convert to sparse
  K <- as(K, "sparseMatrix")
  D_P_block <- as(D_P_block, "sparseMatrix")
  D_F_block <- as(D_F_block, "sparseMatrix")
  P_block <- as(P_block, "sparseMatrix")
  F_block <- as(F_block, "sparseMatrix")
  
  # projection matrix structured by both age and parental age
  A_par <- (t(K) %*% D_P_block %*% K %*% P_block) +
    (t(K) %*% D_F_block %*% K %*% F_block)

  return(as.matrix(A_par))
}




#### Analyze A_par -------------------------------------------------------------

# extract age-by-parental-age sensitivities from A_par
extract_sens <- function(A_par_l, w, v, stasis) {
  # A_par_l:  list object returned from function build_Apar
  # w:  stable distribution of A_par (scaled right dominant eigenvector)
  # v:  reproductive value distribution of A_par (scaled left dominant eigenvector)
  # stasis:  if TRUE, include stasis loop at maximum age class omega
  o <- A_par_l$o
  s <- A_par_l$s
  
  v <- v / sum(v * w)
  w <- w / sum(v * w)
  
  I_sw <- as(diag(1, nrow = o*s, ncol = o*s), "sparseMatrix")
  vec_I_w <- as(c(diag(rep(1, o))), "sparseMatrix")
  I_s <- as(diag(rep(1, s)), "sparseMatrix")
  I_s2 <- as(diag(rep(1, s^2)), "sparseMatrix")
  
  w_t <- as(t(w), "sparseMatrix")
  v_t <- as(t(v), "sparseMatrix")
  
  term1 <- w_t %x% v_t
  term2_surv <- I_sw %x% (t(A_par_l$K) %*% A_par_l$D_P_block %*% A_par_l$K)
  term2_fert <- I_sw %x% (t(A_par_l$K) %*% A_par_l$D_F_block %*% A_par_l$K)
  term4 <- vec_I_w %x% I_s2
  
  sens_surv <- matrix(0, nrow = o, ncol = s)
  sens_fert <- matrix(0, nrow = o, ncol = s)
  
  for(i in 1:o) {
    E_ii <- Matrix(0, nrow = o, ncol = o, sparse = T)
    E_ii[i,i] <- 1
    
    term3 <- E_ii %x% A_par_l$K %x% I_s
    
    vec_sens_surv <- term1 %*% term2_surv %*% term3 %*% term4
    vec_sens_fert <- term1 %*% term2_fert %*% term3 %*% term4
    
    sens_surv[i,] <- diag(matrix(vec_sens_surv, nrow = s))
    sens_fert[i,] <- matrix(vec_sens_fert, nrow = s)[A_par_l$q_i[i],]
  }
  
  # sensitivty to mortality hazard
  sens_hazr <- sens_surv * -A_par_l$P_mat

  # remove sensitivities to survival for final age class, if stasis == FALSE
  if (stasis == FALSE) {
    sens_surv[o,] <- NA_real_
    sens_hazr[o,] <- NA_real_
  }
  
  return(list(sens_surv = sens_surv, sens_hazr = sens_hazr, sens_fert = sens_fert))
}


# analysis of A_par and A_ref
analyze_mod <- function(A_par_l, stasis = FALSE, checks = TRUE) {
  # A_par_l:  list object returned from function build_Apar
  # stasis:  if TRUE, include stasis loop at maximum age class omega
  # checks:  if TRUE, check whether lambda_par == lambda_ref and w_i_par == w_i_ref
  o <- A_par_l$o
  s <- A_par_l$s
  
  ## summary analysis of A_par
  eig_par <- eigen(A_par_l$A_par)
  lmax_par <- which.max(Re(eig_par$values))
  lamda_par <- Re(eig_par$values[lmax_par])
  W_par <- eig_par$vectors
  w_par <- abs(Re(W_par[,lmax_par]))
  w_par <- w_par / sum(w_par)
  v_par <- abs(Re(eigen(t(A_par_l$A_par))$vectors[,lmax_par]))
  stable_mat_par <- matrix(w_par, byrow = TRUE, ncol = s)
  
  ## age-by-parental-age sensitivities from A_par
  sens_par <- extract_sens(A_par_l, w = w_par, v = v_par, stasis = stasis)
  sens_surv_par <- sens_par$sens_surv
  sens_hazr_par <- sens_par$sens_hazr
  sens_fert_par <- sens_par$sens_fert
  
  ### age-specific sensitivities from A_par
  sens_surv_age_par <- rowSums(sens_surv_par)
  sens_hazr_age_par <- rowSums(sens_hazr_par)
  sens_fert_age_par <- rowSums(sens_fert_par)
  
  ## stable distributions from A_par
  stable_age_par <- rowSums(stable_mat_par)
  stable_stage_par <- colSums(stable_mat_par)
  stable_stage_par_rel <- stable_mat_par / stable_age_par
  
  ## expected age-specific survival and fecundity at stable distribution
  P_ref <- rowSums(stable_stage_par_rel * A_par_l$P_mat)
  F_ref <- rowSums(stable_stage_par_rel * A_par_l$F_mat)
  
  ## take an unweighted average (response to Assoc. Ed.)
  # P_ref <- rowMeans(A_par_l$P_mat)
  # F_ref <- rowMeans(A_par_l$F_mat)
  
  ## test: weight ref tr by relative repro value rather than stable stage
  v_mat_par <- matrix(v_par, byrow = TRUE, ncol = s)
  v_age_par <- rowSums(stable_stage_par_rel * v_mat_par)
  
  ## create A_ref (same as A_par but without parental age effect)
  A_ref <- matrix(0, nrow = o, ncol = o)
  A_ref[1,] <- F_ref
  A_ref[2:o, 1:(o-1)] <- diag(P_ref[1:(o-1)])
  if (stasis == TRUE) A_ref[o,o] <- P_ref[o]
  
  ## summary analysis of A_ref
  eig_ref <- eigen(A_ref)
  lmax_ref <- which.max(Re(eig_ref$values))
  lamda_ref <- Re(eig_ref$values[lmax_ref])
  W_ref <- eig_ref$vectors
  w_ref <- abs(Re(W_ref[,lmax_ref]))
  w_ref <- w_ref / sum(w_ref)
  v_ref <- abs(Re(eigen(t(A_ref))$vectors[,lmax_ref]))
  sens_ref <- v_ref %o% w_ref / sum(v_ref * w_ref)
  
  ## age-specific sensitivities from A_ref
  sens_fert_age_ref <- sens_ref[1,]
  sens_surv_age_ref <- sens_ref[cbind(c(2:o, o), 1:o)]
  sens_hazr_age_ref <- sens_surv_age_ref * -P_ref
  
  # remove P and sens_surv for maximum age class, if stasis == FALSE
  if (stasis == FALSE) {
    P_ref[o] <- NA_real_
    sens_surv_age_par[o] <- NA_real_
    sens_surv_age_ref[o] <- NA_real_
    sens_hazr_age_par[o] <- NA_real_
    sens_hazr_age_ref[o] <- NA_real_
  }
  
  ## check whether equilibrium properties of A_par and A_ref are equivalent
  if (checks == TRUE) {
    if (abs(lamda_par - lamda_ref) > 1e-8) warning('lambda_par != lambda_ref')
    if (any(abs(stable_age_par - w_ref) > 1e-8)) warning('w_i_par != w_i_ref')
  }
  
  ### Gather and summarize model outputs
  w_par_vec <- c(matrix(w_par, byrow = TRUE, ncol = s))
  v_par_vec <- c(matrix(v_par, byrow = TRUE, ncol = s))
  sens_surv_par_vec <- c(sens_surv_par)
  sens_hazr_par_vec <- c(sens_hazr_par)
  sens_fert_par_vec <- c(sens_fert_par)
  
  df_par <- expand.grid(i = 1:o, j = 1:s) %>% 
    mutate(w = w_par_vec,
           v = v_par_vec/v_age_par[1],
           sens_surv = sens_surv_par_vec,
           sens_hazr = sens_hazr_par_vec,
           sens_fert = sens_fert_par_vec)
  
  df_age <- data.frame(
    i = 1:o,
    surv_ref = P_ref,
    fert_ref = F_ref,
    w_ref = w_ref,
    w_par = stable_age_par,
    v_ref = v_ref/v_ref[1],
    v_par = v_age_par/v_age_par[1],
    sens_surv_ref = sens_surv_age_ref,
    sens_hazr_ref = sens_hazr_age_ref,
    sens_fert_ref = sens_fert_age_ref,
    sens_surv_par = sens_surv_age_par,
    sens_hazr_par = sens_hazr_age_par,
    sens_fert_par = sens_fert_age_par
  )
  
  return(list(df_age = df_age, df_par = df_par, A_ref = A_ref, P_ref = P_ref, F_ref = F_ref))
}




#### Simulate parental age effects ---------------------------------------------

# logit function
logit <- function(p) {
  # p: probability value to transform to log-odds scale 
  return(log(p / (1 - p)))
}

# inverse logit function
logit_inv <- function(a) {
  # a: log-odds value to transform to probability scale
  return(exp(a) / (exp(a) + 1))
}

# main analysis
SimulateParEffect <- function(P_i, F_i, s, P_range = 0.50, F_range = 0.15) {
  # P_i:     vector of baseline age-specific survival
  # F_i:     vector of baseline age-specific fecundity
  # s:       number of parental age classes
  # P_range: range (+/-) for simulated parental age effect on survival
  # F_range: range (+/-) for simulated parental age effect on fecundity
  o <- length(P_i)
  
  P_mat_raw <- matrix(rep(P_i, s), nrow = o)
  F_mat_raw <- matrix(rep(F_i, s), nrow = o)
  
  P_eta_mat <- matrix(rep(seq(P_range, -P_range, length.out = s), o), nrow = o, byrow = TRUE)
  P_logit <- logit(P_mat_raw)
  P_out <- logit_inv(P_eta_mat + P_logit)
  
  F_gamma_mat <- matrix(rep(seq(1+F_range, 1-F_range, length.out = s), o), nrow = o, byrow = TRUE)
  F_out <- F_gamma_mat * F_mat_raw
  
  out_df <- P_out %>% 
    as_tibble() %>% 
    setNames(1:s) %>% 
    mutate(i = 1:o) %>% 
    gather(j, P_ij, -i) %>% 
    mutate(F_ij = c(F_out)) %>% 
    mutate(j = as.integer(j))
  
  return(out_df)
}


# Scenario C (offspring quality increases with parental age)
SimulateParEffectRev <- function(P_i, F_i, s, P_range = 0.50, F_range = 0.15) {
  # P_i:     vector of baseline age-specific survival
  # F_i:     vector of baseline age-specific fecundity
  # s:       number of parental age classes
  # P_range: range (+/-) for simulated parental age effect on survival
  # F_range: range (+/-) for simulated parental age effect on fecundity
  o <- length(P_i)
  
  P_mat_raw <- matrix(rep(P_i, s), nrow = o)
  F_mat_raw <- matrix(rep(F_i, s), nrow = o)
  
  P_eta_mat <- matrix(rep(seq(-P_range, P_range, length.out = s), o), nrow = o, byrow = TRUE)
  P_logit <- logit(P_mat_raw)
  P_out <- logit_inv(P_eta_mat + P_logit)
  
  F_gamma_mat <- matrix(rep(seq(1-F_range, 1+F_range, length.out = s), o), nrow = o, byrow = TRUE)
  F_out <- F_gamma_mat * F_mat_raw
  
  out_df <- P_out %>% 
    as_tibble() %>% 
    setNames(1:s) %>% 
    mutate(i = 1:o) %>% 
    gather(j, P_ij, -i) %>% 
    mutate(F_ij = c(F_out)) %>% 
    mutate(j = as.integer(j))
  
  return(out_df)
}



#### Adjust fecundity transition rates to achieve given lambda -----------------
opt_fn <- function(par, s, P_i, F_i, P_range, F_range, lam, stasis = FALSE) {
  F_i_adj <- par * F_i
  Tr_sim <- SimulateParEffect(P_i, F_i_adj, s = s, P_range = P_range, F_range = F_range)
  q_i <- map_ij(F_i_adj, s = s)
  A_par <- build_Apar_lite(Tr_sim, q_i, stasis = stasis)
  return(abs(lam - popbio::lambda(A_par)))
}

opt_fn_rev <- function(par, s, P_i, F_i, P_range, F_range, lam, stasis = FALSE) {
  F_i_adj <- par * F_i
  Tr_sim <- SimulateParEffectRev(P_i, F_i_adj, s = s, P_range = P_range, F_range = F_range)
  q_i <- map_ij(F_i_adj, s = s)
  A_par <- build_Apar_lite(Tr_sim, q_i, stasis = stasis)
  return(abs(lam - popbio::lambda(A_par)))
}




#### Plotting functions --------------------------------------------------------
LabelFn <- function(x) {
  # x:  label
  return(formatC(x, format = 'fg'))
}

