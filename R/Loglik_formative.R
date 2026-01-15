Loglik_formative <- function(K , nD, mapping,nL,mapping2, paraOpt,  paraFixe , posfix , paras_k ,
                             sequence , type_int , ind_seq_i, MCnr , nmes ,
                             m_is , Mod_MatrixY , Mod_MatrixYprim , df,
                             x , z , q , nb_paraD ,
                             x0 , z0 , q0 , cholesky ,
                             data_surv , data_surv_intY , nYsurv , basehaz , knots_surv, 
                             np_surv , survival , assoc , truncation, 
                             nE, Xsurv1 , Xsurv2 ,
                             if_link , zitr, ide,
                             tau , tau_is, 
                             modA_mat, DeltaT, ii){
  
  # fix model.matrix for formative
  mappingLP2LP1 <- pmin(table(mapping, mapping2), 1)
  mappingLP2LP1_vec <- apply(mappingLP2LP1,2,function(x)which(x!=0))
  x <-model_matrix_nL(x,repeats_lp = mappingLP2LP1_vec)
  z <-model_matrix_nL(z,repeats_lp = mappingLP2LP1_vec)
  q <- if(all(q)==0)rep(0,nL)else rep(q,times= as.numeric(table(mappingLP2LP1_vec)))
  nb_paraD <- nb_paraD
  x0 <-model_matrix_nL(x0,repeats_lp = mappingLP2LP1_vec)
  z0 <-model_matrix_nL(z0,repeats_lp = mappingLP2LP1_vec)
  q0 <- if(all(q0)==0)rep(0,nL)else rep(q0,times= as.numeric(table(mappingLP2LP1_vec)))
  
  # LogLik
  res <- Loglik(
    K = K,
    nD = nL,
    mapping =  mapping2,
    paraOpt = paraOpt,
    paraFixe = paraFixe,
    posfix = posfix,
    paras_k = paras_k,
    sequence = sequence,
    type_int = type_int,
    ind_seq_i = ind_seq_i,
    MCnr = MCnr,
    nmes = nmes,
    m_is = m_is,
    Mod_MatrixY = Mod_MatrixY,
    Mod_MatrixYprim = Mod_MatrixYprim,
    df = df,
    x = x,
    z = z,
    q = q,
    nb_paraD = nb_paraD,
    x0 = x0,
    z0 = z0,
    q0 = q0,
    cholesky = cholesky,
    data_surv = data_surv,
    data_surv_intY = data_surv_intY,
    nYsurv = nYsurv,
    basehaz = basehaz,
    knots_surv = knots_surv,
    np_surv = np_surv,
    survival = survival,
    assoc =  assoc,
    truncation = truncation,
    nE = nE,
    Xsurv1 = Xsurv1,
    Xsurv2 = Xsurv2,
    if_link = if_link,
    zitr = zitr,
    ide = ide,
    tau = tau,
    tau_is = tau_is,
    modA_mat = modA_mat,
    DeltaT,
    ii = ii
  )
  
  return(res)
  
}