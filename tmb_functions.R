


#' @title Build Stack for TMB MBG model
#' @description Organize Data and Parameter Stacks to fit a TMB MBG model
#' @author Roy Burstein
#'
#' @param d prepped model frame with no NAs
#' @param yl a vector of years for analysis (i.e. c(2001,2002,2003))
#' @param zl a vector of zcol for analysis (i.e. ages c(1,2,3,4,5). Must be integers starting with 1)
#' @param fes a string of model fixed effects of form 'var1 + var2 + var3' or as string vector
#'  corresponding to column names in d
#' @param indic indicator name, corresponding to the response variable in d
#' @param country_re TRUE/FALSE include country re.
#' @param nid_re TRUE/FALSE include re on survey. Need nid column in the data frame. defaults to use_nid_res from config. 
#' @param exclude_cs character vector of covariates to exclude from centrescaling
#' @param mesh an inla mesh object. NOTE! this gets remade inside this function... Roy said there was a good reason
#' @param s2mesh Logical. Should the mesh be created on the surface of
#'   a sphere? If TRUE, s2params is used to specify mesh parameters
#'   instead of max_edge and mesh_offset
#' @param s2params string of 3 element numeric vector in R
#'   notation. e.g. "c(25, 500, 1000)". The entries describe the
#'   minimum triangle edge length allowed, hos far to extend the mesh
#'   beyond the 'simple' boundary, and the maximum allowed triangle
#'   edge length, respectively. Units are in kilometers. Used only if
#'   s2mesh=TRUE.
#' @param cov_constraints named int vector. integer vector indexed by covariate names.
#'
#' @return returns a named list with Data and Parameters to be passed into fit_mbg_tmb()
#' 
build_mbg_data_stack_tmb <- function(d          = df,                 
                                     fes        = all_fixed_effects,  
                                     indic      = indicator, 
                                     mesh       = mesh_s,
                                     cov_constraints = covariate_constraint_vectorize(config)
                                     ){   
  
  # require some libraries
  require(INLA)
  require(data.table)
  
  # ensure d is a dt
  d <- setDT(d)

  # For a single time period and age, all observations are in the same grouping
  d$group <- 1
  
  # coordinates at data points. these are passed to TMB in long,lat so
  # we keep them and create another set for 3d coords
  coords   <- cbind(d$longitude,d$latitude)
  
  # if we have  mesh on s2, first convert coords to spherical to project to mesh
  data.locs <- coords ## long, lat
  if(mesh$manifold == "S2"){
    ## then the mesh is on the sphere and we need to use 3d coords
    data.locs <- lonlat3D(data.locs[, 1], data.locs[, 2])
  }
  
  # make a projection matrix from data to st mesh                                    
  A.proj <- inla.spde.make.A(mesh  = mesh,
                             loc   = data.locs,
                             group = d$group)
  
  # Build SPDE object using INLA (must pass mesh$idx$loc when supplying Boundary) to get matrices
  spde <- inla.spde2.matern(mesh, alpha = 2, constr = FALSE) 
                                       
  # make a clean design matrix. make sure all fes appear in d
  fes  <- unlist(strsplit(fes, ' \\+ ')) # i dont think this is safe and we might want to consider changing to formulas
  if(!all(fes %in% names(d)))  
    stop('Check your fes argument, not all covariate names appear in d.')
  if(length(fes)!=0) {
    X_xp <- as.matrix(cbind(int=1, d[,c(fes),with=FALSE]))
  } else {
    X_xp <- as.matrix(cbind(int=rep(1,nrow(d))))
  }

  # Center-scale data. imports seegMBG
  cs_df <- getCentreScale(X_xp, exclude ='int')
  X_xp  <- centreScale(X_xp, df = cs_df)

  # get data range in case we want to clamp for prediction later
  clamper <- data.table(apply(X_xp,2,range))
    
  # set GP RE vector, or matrix depending on if we have a z dimension
  epsilon_s <- rep(0, length=mesh$n)

  # run a quick regression to get starting values for the fixed effects 
  # This can speed up model fitting if iterations are slow.
  message('LM for starting fe vals')
  y <- (d[[indic]]+.0001)/d$N
  y[y<=0] <- 0.001
  y[y>=1] <- 0.999
  fe_start <- round( unname( lm(qlogis(y) ~ -1 + X_xp)$coefficients ), 4)
  message(sprintf('starting values for fixed effects: %s',paste0(fe_start,collapse=', ')))

  # check if user wants to scale gaussian variance by N, if not set them all to one in the gaussian rows
  n_i <- d$N

  use_poisson <- as.integer(tolower(indicator_family)=='poisson')
  
  # Construct a list of all Data necessary to TMB to fit
  Data <- list(
    num_i          = nrow(d),            # Total number of observations
    num_s          = mesh$n,             # Number of vertices in SPDE mesh
    y_i            = d[[indic]],         # Number of observed events in the cluster
    n_i            = d$N,                # Number of observed exposures in the cluster
    w_i            = d$weight,           # Data weight for each row
    X_ij           = X_xp,               # Covariate design matrix
    M0             = spde$param.inla$M0, # SPDE sparse matrix
    M1             = spde$param.inla$M1, # SPDE sparse matrix
    M2             = spde$param.inla$M2, # SPDE sparse matrix
    Aproj          = A.proj,             # mesh to prediction point projection matrix
    options = list(
      use_priors = 1,    # option1==1 use priors 
      adreport_off = 1,  # option2==1 ADREPORT off
      use_poisson = use_poisson # Option 3==1: Use a Poisson likelihood
    ),
    fconstraints = tmb_cov_constraint(colnames(X_xp), cov_constraints)
  )
  
  # Set staring values for parameters
  Parameters <- list(
    alpha_j          = fe_start,  # FE parameters alphas
    logtau           = -0.5,      # Matern/AR tau
    logkappa         = -0.5,	    # Matern Range
    epsilon_s        = epsilon_s  # Random Effects: GP locations
  )
  
  
  # put bounds on parameters (Note, this is especially important for rhos)
  L  <- c(rep(-10,ncol(X_xp)),-10,-10)
  U  <- c(rep( 10,ncol(X_xp)), 10, 10)
  pn <- c(rep('alpha_j',ncol(X_xp)),'logtau','logkappa')
  names(L) <- names(U) <- pn
  
  # return the list
  return(list(Data         = Data,
              Parameters   = Parameters,
              clamper      = clamper,
              cs_df        = cs_df,
              coords       = coords,
              mesh         = mesh,
              L            = L,
              U            = U))

}





#' @title Fit a TMB MBG model
#' @description Fit a TMB MBG model and pull jointPrecision report
#' @author Roy Burstein, adapted by Nat Henry for KIT TB project
#' 
#' @param lbdcorerepo core repo location
#' @param cpp_template name of cpp template file within ./<lbdcorerepo>/mbg_central/
#' @param tmb_input_stack object that results from build_mbg_data_stack_tmb() or build_mbg_data_stack(...,tmb=TRUE)
#' @param ADmap_list map parameter for ignoring parameters
#' @param control_list pass control list to nlminb()
#' @param optimizer which software to use for optimization (optim, or nlminb)
#'
#' @return list of tmb model objects: ADfun, opt, sdrep
#'
#' @useDynLib mbg_tmb_model
#'
fit_mbg_tmb <- function(lbdcorerepo     = core_repo,
                        cpp_template    = 'mbg_tmb_model',
                        tmb_input_stack = input_data,
                        ADmap_list      = NULL,
                        control_list    = NULL,
                        optimizer       = 'nlminb',
                        newton_steps    = 0,
                        sparse_ordering  = TRUE
                        ){
  library(optimx)
  message('WARNING: This TMB implementation does not support sum to one constraints yet.')
  
  # make sure TMB is loaded
  require(TMB)
  
  # compile the cpp file and dynload it
  message('compiling template')
  TMB::compile(sprintf('%s/mbg_central/%s.cpp', lbdcorerepo, cpp_template))
  dyn.load( dynlib(sprintf('%s/mbg_central/%s', lbdcorerepo, cpp_template)) )
  
  # deal with parallelization
  threads <- system('echo $OMP_NUM_THREADS', intern = TRUE)
  if(threads != '') {
    message(sprintf('Detected %s threads in OMP environmental variable.',threads))
    openmp(as.numeric(threads))
  } else {
    message('Did not detect environmental OMP variable, defaulting to 4 cores. \n
             You can set this using OMP_NUM_TREADS or when launching singularity image.')
    openmp(4)
  }
  
  if(is.null(ADmap_list)) ADmap_list <- list()

  # REs
  randompars <- c("epsilon_s")

  # set Data flag for Kaspers normalization fix
  tmb_input_stack$Data$flag <- 1 

  # print
  message(paste0('ADMAP_LIST: ',  paste0(names(ADmap_list),collapse=', ')))
  message(paste0('Random Pars: ', paste0(randompars,collapse=', ')))
  
  # make the AD object
  message('Making AD object')
  obj <- MakeADFun(data       = tmb_input_stack$Data, 
                   parameters = tmb_input_stack$Parameters,  
                   map        = ADmap_list, 
                   random     = randompars, 
                   hessian    = TRUE, 
                   DLL        = cpp_template)
  
  # Run optimizer
  message('Running MLE')
  # FIT THE MODEL
  # Try using multiple optimization algorithms
  converged <- FALSE
  opt_methods <- c('BFGS','nlminb')
  try_opt_num <- 1
  while( (converged==FALSE) & (try_opt_num <= length(opt_methods)) ){
    this_opt_method <- opt_methods[try_opt_num]
    opt0 <- optimx(
      par     = obj$par,
      fn      = function(x) as.numeric(obj$fn(x)),
      gr      = obj$gr,
      method  = this_opt_method,
      itnmax  = 500,
      control = control_list
    )
    if (opt0$convcode == 0){
      converged <- TRUE
    } else {
      # In case the model did not converge, go on to the next optimization method
      try_opt_num <- try_opt_num + 1
    } 
  }

  # Check to make sure that the distribution converged
  if (opt0$convcode != 0) stop("Model did not converge")

  message(sprintf("Running %i Newton steps...",newton_steps))
  for(i in seq_len(newton_steps)) {
    g = as.numeric( obj$gr(opt0$par) )
    h = optimHess(opt0$par, fn=obj$fn, gr=obj$gr)
    opt0$par = opt0$par - solve(h, g)
    opt0$objective = obj$fn(opt0$par)
  }
  
  # run sdreport to get joint precision of all parameters
  for(i in 1:20) message('** Getting Joint Precision **')
  SD0 <- TMB::sdreport(obj, getJointPrecision = TRUE, bias.correct = TRUE)

  # return
  return( list ( ADfun   = obj,
                 opt     = opt0,
                 sdrep   = SD0,
                 fenames = colnames(tmb_input_stack$Data$X_ij)) )
  
}





#' Take multivariate normal draws given a mean vector and precision matrix
#'
#' @param mu vector of parameter means
#' @param prec joint precision matrix
#' @param n.sims number of draws
#'
#' @return length(mu) by n.sims matrix of parameter draws
#' 
rmvnorm_prec <- function(mu, prec, n.sims) {
  z <- matrix(rnorm(length(mu) * n.sims), ncol=n.sims)
  L <- Cholesky(prec, super = TRUE)
  z <- solve(L, z, system = "Lt") ## z = Lt^-1 %*% z
  z <- solve(L, z, system = "Pt") ## z = Pt    %*% z
  z <- as.matrix(z)
  return(mu + z)
}









#' @title Predict MBG from a TMB Model
#'
#' @description Project out to a full sample space defined by the sr argument
#' 
#' @author Roy Burstein, adapted by Nat Henry for KIT TB project
#' 
#' @param samples Number of draws to take
#' @param seed Seed to set for RNG
#' @param model_fit_object object output from function fit_mbg_tmb()
#' @param tmb_input_data input_data output of build_mbg_data_stack_tmb
#' @param fes a string of model fixed effects of form 'var1 + var2 + var3' or as string vector
#'  corresponding to column names in d
#' @param sr simple raster object
#' @param yl a vector of years for analysis (i.e. c(2001,2002,2003))
#' @param zl a vector of z for analysis (i.e. c(1,2,3,4,5)). Config item z_list. No z-dimension it should just be zero
#' @param covs_list named list of covariate bricks, each should be length(yl) long
#'
#' @return a cell_preds object
#'
predict_mbg_tmb <- function(samples,
                            seed             = NULL,
                            tmb_input_stack  = input_data,
                            model_fit_object = model_fit,
                            fes              = all_fixed_effects,
                            sr               = simple_raster,
                            yl               = year_list,
                            zl               = z_list,
                            transform        = 'inverse-logit',
                            covs_list        = cov_list,
                            clamp_covs       = FALSE,
                            cov_constraints = covariate_constraint_vectorize(config)) {

  # libs
  require(raster)
  require(sp)
  
  #
  sdrep     <- model_fit_object$sdrep
  
  # pull a few useful things from the input data stack
  cs_transform         <- tmb_input_stack$cs_df
  mesh                 <- tmb_input_stack$mesh
  cntry_re_map         <- tmb_input_stack$cntry_re_map
  if(clamp_covs == TRUE) {
    clamper            <- tmb_input_stack$clamper
  } else {
    clamper            <- NULL
  }

  # set seed if it is requested
  if(!is.null(seed)) set.seed(seed)
  
  # vector of means
  mu    <- c(sdrep$par.fixed,sdrep$par.random)
  
  # simulate draws
  draws <- rmvnorm_prec(mu = mu , prec = sdrep$jointPrecision, n.sims = samples)
  
  ## separate out the draws
  parnames      <- c(names(sdrep$par.fixed), names(sdrep$par.random))
  epsilon_draws <- draws[parnames=='epsilon_s',]
  alpha_draws   <- draws[parnames=='alpha_j',]

  # remove non-predicted FEs from draws so they are discareded from here out
  alpha_draws   <- alpha_draws[which(!model_fit_object$fenames %in% non_pred_fes),]
  model_fit_object$fenames <- model_fit_object$fenames[!model_fit_object$fenames %in% non_pred_fes]
  # seperate out Z FE draws
  FE_z_draws    <- NULL
  if(length(zl) > 1){
    FE_z_draws    <- alpha_draws[ grepl('FE_z_level__',model_fit_object$fenames),] # separate out z-level fixed effects from other FEs
    alpha_draws   <- alpha_draws[!grepl('FE_z_level__',model_fit_object$fenames),] # remove any z-level fixed effects from other FEs
  }

  
  if(length(zl) > 1)
    if(dim(FE_z_draws)[1] != (length(zl)-1) )
      stop('Incorrect number of fixed effects for levels of z in the z_list')
  
  # names of fes
  tmb_const <- tmb_cov_constraint(model_fit_object$fenames, cov_constraints)
  fes       <- unlist(strsplit(fes, ' \\+ '))
  
  # mask covariates to simple raster
  for(l in 1:length(covs_list)) {
    covs_list[[l]]  <- crop(covs_list[[l]], extent(sr))
    covs_list[[l]]  <- setExtent(covs_list[[l]], sr)
    covs_list[[l]]  <- mask(covs_list[[l]], sr)
  }
  
  # keep only covariates used in model, typically stacking
  covs_list <- covs_list[names(covs_list) %in% fes]
  
  
  # get coordinates of full projection space
  f_orig <- data.table(cbind(xyFromCell(sr, seegSDM:::notMissingIdx(sr)), gaul_code=as.vector(sr[seegSDM:::notMissingIdx(sr)]))) # extract gaul
  f_orig$t <- f_orig$z <- 1 # set initial time and Z
  f_orig[,tmpord:=1:.N]
  
  # add time periods and z periods as needed
  grp <- setDT(expand.grid(1:length(yl), 1:length(zl)))
  setnames(grp,c('Var1','Var2'),c('t','z'))
  grp[,group := 1:.N]
  fullsamplespace <- data.table()
  for(g in 1:max(grp$group)){
    tmp <- f_orig
    tmp[,z  := grp$z[grp$group==g]]
    tmp[,t  := grp$t[grp$group==g]]
    tmp[,gp := g]
    fullsamplespace <- rbind(fullsamplespace,tmp)
  }
  fullsamplespace[,idx := 1:.N]

  # pull out covariates in format we expect them
  # a list of length periods with a brick of named covariates inside
  new_cl <- list()
  if(length(covs_list)==0){
    message('No covariates detected, predicting using intercept only.')
  } else {
    message(sprintf('%i covariates detected.',length(covs_list)))
    for(p in 1:length(yl)){
      new_cl[[p]] <- list()
      for(n in names(covs_list)){
        if(dim(covs_list[[n]])[3]==1) { # synoptic mean covariates
          new_cl[[p]][[n]] <- covs_list[[n]]
        } else if (dim(covs_list[[n]])[3]==length(yl)) { # time varying covariates
          new_cl[[p]][[n]] <- covs_list[[n]][[p]]
        } else { # error if there is some other weird non-conforming year thing
          stop(sprintf('Covariate %s is a brick with %i layers, while year_list has %i years',
                       n,dim(covs_list[[n]])[3],length(yl)))
        }
      }
      new_cl[[p]] <- brick(new_cl[[p]])
    }
  }
  
  # get surface locs to project on to
  pcoords        <- cbind(x=fullsamplespace[t==1,x], y=fullsamplespace[t==1,y]) ## used for cov raster extract

  ## setup coords for GP projection. convert coords to spherical if
  ## using spherical modeling mesh. used if you made a mesh on s2
  if(mesh_s$manifold == "S2"){
    gp_coords <- lonlat3D(pcoords[, 1], pcoords[, 2])
  } else {
    gp_coords <- pcoords
  }

  ## define grouping across periods
  groups_periods <- fullsamplespace[t==1,gp]
                  
  # extract cell values  from covariates, deal with timevarying covariates here
  cov_vals <- list()
  for(z in 1:length(zl)){
    cov_vals[[z]] <- list()
    for(p in 1:length(yl)){
      if(length(fes)>0) {
        
        # raster extract and keep only fes
        cov_vals[[z]][[p]] <- raster::extract(new_cl[[p]], pcoords[1:nrow(f_orig),])
        cov_vals[[z]][[p]] <- cov_vals[[z]][[p]][,colnames(cov_vals[[z]][[p]]) %in% c(fes)]
        # If there is only a single covariate, convert from vector to matrix
        if( (length(covs_list)==1) & !('matrix' %in% class(cov_vals[[z]][[p]]))){
          cov_vals[[z]][[p]] <- matrix(cov_vals[[z]][[p]], ncol=1)
        }
        colnames(cov_vals[[z]][[p]]) <- fes

        # transform raw covariate values (center scaled) if needed (i.e. if cs_tranform is not 1 0 for that variable)
        cov_vals[[z]][[p]] <- centreScale(cov_vals[[z]][[p]],cs_transform) 
        
        # clamp covariates if clamper is not null
        if(!is.null(clamper)){
          message('** Clamping **')
          for(fe in fes){
            tmpvec <- cov_vals[[z]][[p]][,colnames(cov_vals[[z]][[p]])==fe]
            old_max <- max(tmpvec, na.rm=T)
            old_min <- min(tmpvec, na.rm=T)
            mn <- as.numeric(clamper[,fe,with=FALSE][1])
            mx <- as.numeric(clamper[,fe,with=FALSE][2])
            tmpvec[tmpvec<mn] <- mn
            tmpvec[tmpvec>mx] <- mx
            message(sprintf(
              "Old range was [%s, %s], now is [%s, %s]",
              old_min, old_max, min(tmpvec,na.rm=T), max(tmpvec,na.rm=T)
            ))
            cov_vals[[z]][[p]][,colnames(cov_vals[[z]][[p]])==fe] <- tmpvec
          }
        }
        # add an intercept
        cov_vals[[z]][[p]] <- cbind(int = 1, cov_vals[[z]][[p]])
      } else {
        # if no covariates just do intercept only
        cov_vals[[z]][[p]] <- cbind(int = rep(1,nrow(f_orig)))
      }
    }
  }
  
  ## use inla helper functions to project the spatial effect.
  A.pred <- inla.spde.make.A(
    mesh  = mesh,
    loc   = gp_coords,
    group = groups_periods
  )
  
  ### values of GP ST surface at each cell (long by nperiods)
  # if we have multiple zs then do this by z since its possible to throw a SuiteSparse 'Problem too large' error here. 
  if(length(zl) > 1){
    cell_s <- list()
    for(zz in zl){
      cell_s[[zz]] <- as.matrix(A.pred[(which(fullsamplespace$z==zz)),] %*% epsilon_draws)
    }
    cell_s <- do.call('rbind',cell_s)
  } else{
    cell_s <- as.matrix(A.pred %*% epsilon_draws)
  }
  # Replicate across all years
  cell_s <- do.call(rbind, replicate(length(yl), cell_s, simplify=FALSE))

  # covariate values by alpha draws
  l_vals <- list()
  for(z in 1:length(zl)){
    l_vals[[z]] <- list()
    for(p in 1:length(yl))  
      l_vals[[z]][[p]] <- cov_vals[[z]][[p]] %*% apply_constraints(tmb_const, rbind(alpha_draws,FE_z_draws))
  }
  cell_l <- do.call("rbind",unlist(l_vals, recursive = FALSE))
  
  # add together linear and st components
  pred_tmb <- cell_l + cell_s
  
  # transform
  if(transform=='inverse-logit') { 
    pred_tmb <- plogis(as.matrix(pred_tmb))
  } else {
    pred_tmb <- eval(parse(text=sprintf('%s(as.matrix(pred_tmb))',transform)))
  }
  
  # if there is more than one z, then return a list of length zl cell_preds
  if(length(zl) > 1){
    pred_tmb_list <- list()
    chunklength <- dim(pred_tmb)[1]/length(zl) 
    for(z in 1:length(zl))
      pred_tmb_list[[z]] <- pred_tmb[((z-1)*chunklength+1):(chunklength*z),1:samples]
    pred_tmb <- pred_tmb_list
  }
  
  # return the predicted cell_pred object
  return(pred_tmb)

}




#' @title Get fitted model parameters from TMB object
#' @description Make a nice table of fitted model parameters from a geostat TMB object
#' @author Roy Burstein
#'
#' @param model_fit fitted TMB object that comes out of fit_mbg_tmb()
#' @param exp_fixed_effects Boolean, should the fixed effects be exponentiated (if model was fit with logit link). Defaults to TRUE
#' @param transform_hyperparams Boolean, should hyperparmeters be transformed from fitted to more natural space. Defaults to TRUE
#' @param draws Integer, number of draws to use. Defaults to 1000
#' @param calculate_range Boolean, should we calculate and report range using kappa in draws space. Defaults to TRUE
#' @param calculate_nom_var  Boolean, should we calculate and report nominal variance using kappa and tau in draws space. Defaults to TRUE
#'
#' @return formatted data.table object
fitted_param_table_tmb <- function(model_fit,
                                   exp_fixed_effects     = TRUE,
                                   transform_hyperparams = TRUE,
                                   draws                 = 1000,
                                   calculate_range       = TRUE,
                                   calculate_nom_var     = TRUE,
                                   cov_constraints = covariate_constraint_vectorize(config)) {
  
  # get draws of parameter values
  mu <- model_fit$sdrep$par.fixed
  pm <- model_fit$sdrep$jointPrecision[1:length(mu),1:length(mu)]
  draws <- rmvnorm_prec(mu,pm,draws)

  # deal with given names of parameters
  fn <- names(model_fit$sdrep$par.fixed)
  fn[fn == 'alpha_j'] <- model_fit$fenames

  # Apply constraint transformations
  tmb_const <- tmb_cov_constraint(model_fit$fenames, cov_constraints)
  draws[names(model_fit$sdrep$par.fixed) == "alpha_j",] <-
    apply_constraints(tmb_const, draws[names(model_fit$sdrep$par.fixed) == "alpha_j",])

  # Transform fixed effects
  if(exp_fixed_effects == TRUE){
    draws[which(fn %in% model_fit$fenames),] <- exp(draws[which(fn %in% model_fit$fenames),])
  }

  # Get the range parameter
  if(calculate_range == TRUE){
    # see equation 6.16 (pg 196) of Blangiardo and Cameletti Book 
    ranger <- sqrt(8) / exp(draws[which(fn == 'logkappa'),])
    draws  <- rbind(draws, ranger)
    fn     <- c(fn, 'range')
  }
  
  # Get the nominal variance parameter
  if(calculate_nom_var == TRUE){
    # see equation 6.17 (pg 196) of Blangiardo and Cameletti Book 
    nominal_variance <- 1 / (4 * pi * (exp(draws[which(fn == 'logkappa'),]))^2 * (exp(draws[which(fn == 'logtau'),]))^2)
    draws <- rbind(draws, nominal_variance)
    fn    <- c(fn, 'nominal_variance')
  }
  
  # transform hyperparmeters
  if(transform_hyperparams == TRUE){
    draws[which(fn == 'logtau'),]             <- exp(draws[which(fn == 'logtau'),])
    draws[which(fn == 'logkappa'),]           <- exp(draws[which(fn == 'logkappa'),])
    draws[which(fn == 'log_nugget_sigma'),]   <- exp(draws[which(fn == 'log_nugget_sigma'),])
    draws[which(fn == 'log_cre_sigma'),]      <- exp(draws[which(fn == 'log_cre_sigma'),])
    draws[which(fn == 'log_nidre_sigma'),]    <- exp(draws[which(fn == 'log_nidre_sigma'),])
    
    fn[fn == 'logtau']           <- 'tau'
    fn[fn == 'logkappa']         <- 'kappa'
    fn[fn == 'log_nugget_sigma'] <- 'nugget_SD'
    fn[fn == 'log_cre_sigma']    <- 'country_RE_SD'
    fn[fn == 'log_nidre_sigma']  <- 'NID_RE_SD'
    
    draws[which(fn == 'zrho'),] <- (exp( draws[which(fn == 'zrho'),] ) - 1) / (exp( draws[which(fn == 'zrho'),] ) + 1)
    draws[which(fn == 'trho'),] <- (exp( draws[which(fn == 'trho'),] ) - 1) / (exp( draws[which(fn == 'trho'),] ) + 1)
    fn[fn == 'zrho'] <- 'age_rho'
    fn[fn == 'trho'] <- 'year_rho'
    
  }
  
  # summarize draws and clean up data table
  su <- data.table(t(apply(draws,1,quantile,c(0.025,0.500,0.975))))
  su[, fn := fn]
  colnames(su) <- c('lower','median','upper','param_name')
  su <- su[,c('param_name','median','lower','upper'),with=FALSE]
  
  # return the final table
  return(su)
}

#' @title Read INLA prior for TMB
#'
#' @description Read in a prior specification that is suited for INLA and make 
#'   it TMB readable.
#'
#' @param prior_string character, character vec of length 1 specifying priors
#'
#' @return List specifying a TMB prior, containing three elements:
#'   - logNormal: Is the prior lognormal (1) or not (0)?
#'   - par1: The first shape parameter. In the lognormal case, the mean
#'   - par2: The second shape parameter. In the lognormal case, the variance
#'
read_inla_prior <- function(prior_string){
  prior_list <- eval(parse(text=prior_string[1]))
  if(!(prior_list$prior %in% c("normal", "loggamma"))){
    stop("TMB implementation only supports normal and loggamma priors.")
  }
  return(list(
    logNormal = ifelse(prior_list$prior == "normal", 1, 0),
    par1 = prior_list$param[1],
    par2 = prior_list$param[2]
  ))
}

#' @title Modify constraints to fit TMB formatting
#' 
#' @author Neal Marquez
#' 
#' @description Optionally add constraints to fixed effects in the TMB 
#'   optimization model. In the TMB model, a constraint label of '0' indicates
#'   that the fixed effect is unconstrained, a label of '1' constrains the 
#'   above zero, and and a label of '-1' constrains the variable below zero.
#'
#' @param fes character, fixed effects in model
#' @param constraints named int vector output from the 
#'   `covariate_constraint_vectorize()` function. Specifies how to constrain 
#'   each fixed effect.
#' @param zl int, additional constraints to pad on which will be 0
#' 
tmb_cov_constraint <- function(
  model_fes,
  constraints = covariate_constraint_vectorize(config)
  ){
  tmb_const <- sapply(unlist(strsplit(model_fes, " \\+ ")), function(m){
    if(m %in% names(constraints)){
      x <- unname(constraints[m])
    }
    else{
      x <- 0
    }
    x
  })

  return(tmb_const)
}

#' @title Apply constraints to fitted value draws
#' 
#' @author Neal Marquez
#' 
#' @description Given draws of fixed effect values generated from a fitted TMB
#'   model and a list of constraints applied to fixed effects in that model, 
#'   apply transformations on constrained fixed effects to reproduce how they
#'   were incorporated in the model (constraining them above or below zero). If
#'   a fixed effect had a constraint label of 0 (unconstrained, the default),
#'   the untransformed draws will be returned.
#' 
#' @param tmb_const int vector, tmb prepped constraints
#' @param alpha_draws matrix, fitted draws of coefficients
#' 
#' @return matrix of transformed beta coefficients for fixed effect draws
#' 
apply_constraints <- function(tmb_const, FE_draws){
  Nfe <- nrow(FE_draws)
  
  FEdraws_const <- sapply(1:Nfe, function(j){
    X <- FE_draws[j,]
    if(tmb_const[j] == 1){
      X <- exp(X)
    }
    else if(tmb_const[j] == -1){
      X <- -exp(X)
    }
    
    return(X)
  })
  
  return(t(FEdraws_const))
}

