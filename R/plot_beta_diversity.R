#' Analysis and plot of Beta diversity (Turnover)
#'
#'

plot_beta_diversity <- function(){
  my_list <- list()
  fungo <- phyloseq::subset_samples(identps_bakt, Trtment != '1champ')
  bako <- phyloseq::subset_samples(identps_fungi, Trtment != '1champ')
  #---- Rarefy to even depth and Compute shannon diversity -----------------------
  set.seed(100)
  fungorare  <-  phyloseq::filter_taxa(fungo, function(x) sum(x) > 1, TRUE)
  fungorare <- phyloseq::rarefy_even_depth(fungorare, sample.size = 1500)
  fungocommun <- phyloseq::filter_taxa(fungo, function(x) sum(x >= 1) > 1, TRUE)
  fungocommun <- phyloseq::rarefy_even_depth(fungocommun, sample.size = 1500)
  bakorare <- phyloseq::filter_taxa(bako, function(x) sum(x) > 1, TRUE)
  bakorare <- phyloseq::rarefy_even_depth(bakorare, sample.size = 1500)
  bakocommun <- phyloseq::filter_taxa(bako, function(x) sum(x >= 1) > 1, TRUE)
  bakocommun <- phyloseq::rarefy_even_depth(bakocommun, sample.size = 1500)
  bc <- phyloseq::sample_data(bakocommun)$ID
  fungocommunf <- fungocommun
  oldDF <- as(phyloseq::sample_data(fungocommunf), "data.frame")
  newDF <- subset(oldDF, ID %in% bc)
  phyloseq::sample_data(fungocommunf) <- phyloseq::sample_data(newDF)
  fungoraref <- fungorare
  oldDF <- as(phyloseq::sample_data(fungoraref), "data.frame")
  newDF <- subset(oldDF, ID %in% phyloseq::sample_data(bakorare)$ID)
  phyloseq::sample_data(fungoraref) <- phyloseq::sample_data(newDF)

  #---- Distance matrix for fungi and bacteria and nematode
  prep_bc_dfs <- function(ps_df){
    otu_ident_f <- phyloseq::otu_table(ps_df)
    tmpdf <- data.frame(plot=phyloseq::sample_data(ps_df)$Plot,
                        trt=phyloseq::sample_data(ps_df)$Trtment)
    tmpdf$trt <- stringr::str_replace(tmpdf$trt,'^1(.*?)$', '\\1')
    tmpdf$trt <- stringr::str_replace(tmpdf$trt,'^2n$', '12n')
    #table(tmpdf$trt)
    rownames(otu_ident_f) <-toupper(paste0(tmpdf$plot, tmpdf$trt))
    sotu <- otu_ident_f[sort(rownames(otu_ident_f)),]
    bray_ident <- vegan::vegdist(as(sotu, 'matrix'), method='jaccard', binary = F)
    a <- a[!rownames(a) %in% 'APIST',]
    a <- a[!rownames(a) %in% 'BPIRE',]
    sa <- a[sort(rownames(sotu)),]
    bc_arbre <- vegan::vegdist(as.matrix(sa[,3:ncol(sa)]))
    # functional dissimilarity
    data_cwm <- FD::dbFD(as.data.frame(scale(x)),a[,c(-1, -2)], calc.FRic= T, calc.FDiv= T)$CWM
    data_cwm <- data_cwm[row.names(sotu),]
    eucl_func_arbre <- dist(data_cwm)

    # phylogenetic distance
    arbre_uni <- picante::unifrac(a[row.names(sotu),c(-1,-2)],tree_ar)

    # soil distance matrix
    env_df <- dplyr::select(phyloseq::sample_data(ps_df), pH, C, N,no3_med, nh4_med)
    env_df <- purrr::map_df(env_df, scale)
    env_df <- as.data.frame(env_df)
    rownames(env_df) <-toupper(paste0(tmpdf$plot, tmpdf$trt))
    env_df <- as.matrix(env_df[sort(rownames(sotu)),])
    eucl_env <- dist(env_df)

    return(list(bc_microbe = bray_ident, arbre = list(
                bc_arbre = bc_arbre,
                unifrac_phylo_arbre = arbre_uni,
                eucl_func_arbre =eucl_func_arbre,
                eucl_env = eucl_env)))
  }

  prep_bc_dfs_nem <- function(){
    sgenus <- nem_matrix_genus[sort(rownames(nem_matrix_genus)),]
    tokeep <- expand.grid(c('A','B','C','D'), c('ACRU','ACSA','BEAL','BEPA','PIGL','PIST','THOC','2NR2A','2N5','2NR5','2N8','4NR7','4N8'))
    tokeep <- paste0(tokeep[,1], tokeep[,2])
    sgenus <- na.omit(sgenus[tokeep,])
    sgenus <- sgenus[,!colSums(sgenus) == 0]
    sgenus <- sgenus[sort(rownames(sgenus)),]
    bray_ident <- vegan::vegdist(sgenus, method='bray')
    data(a)
    sa <- a[rownames(sgenus),]
    san <- dplyr::select(sa,ABBA:ncol(sa))
    san <- na.omit(san[,which(colSums(san, na.rm =T)>0)])
    bc_arbre <- vegan::vegdist(san)

    # bacteria dissbray_ident_f
    tmpdf <- data.frame(plot=phyloseq::sample_data(bakocommun)$Plot,
                        trt=phyloseq::sample_data(bakocommun)$Trtment)
    tmpdf$trt <- stringr::str_replace(tmpdf$trt,'^1(.*?)$', '\\1')
    tmpdf$trt <- stringr::str_replace(tmpdf$trt,'^2n$', '12n')
    phyloseq::sample_data(bakocommun)$trt <- toupper(paste0(tmpdf$plot,tmpdf$trt))

    bakonem <- bakocommun
    oldDF <- as(phyloseq::sample_data(bakonem), "data.frame")
    newDF <- subset(oldDF, trt %in% rownames(sgenus))
    phyloseq::sample_data(bakonem) <- phyloseq::sample_data(newDF)

    otu_ident_b <- phyloseq::otu_table(bakonem)
    rownames(otu_ident_b) <- phyloseq::sample_data(bakonem)$trt
    sotu <- otu_ident_b[sort(rownames(otu_ident_b)),]
    bray_ident_b <- vegan::vegdist(sotu, method='bray')
    sgenus_b <- nem_matrix_genus[sort(rownames(otu_ident_b)),]
    bray_ident_nb <- vegan::vegdist(sgenus_b, method='bray')

    # fungi diss
    tmpdf <- data.frame(plot=phyloseq::sample_data(fungocommun)$Plot,
                        trt=phyloseq::sample_data(fungocommun)$Trtment)
    tmpdf$trt <- stringr::str_replace(tmpdf$trt,'^1(.*?)$', '\\1')
    tmpdf$trt <- stringr::str_replace(tmpdf$trt,'^2n$', '12n')
    phyloseq::sample_data(fungocommun)$trt <- toupper(paste0(tmpdf$plot,tmpdf$trt))

    fungonem <- fungocommun
    oldDF <- as(sample_data(fungonem), "data.frame")
    newDF <- subset(oldDF, trt %in% rownames(sgenus))
    sample_data(fungonem) <- sample_data(newDF)

    otu_ident_f <- phyloseq::otu_table(fungonem)
    rownames(otu_ident_f) <- phyloseq::sample_data(fungonem)$trt
    sotu <- otu_ident_f[sort(rownames(otu_ident_f)),]
    bray_ident_f <- vegan::vegdist(sotu, method='bray')
    sgenus_f <- nem_matrix_genus[sort(rownames(otu_ident_f)),]
    bray_ident_nf <- vegan::vegdist(sgenus_f, method='bray')

    # functional dissimilarity
    data_cwm <- FD::dbFD(as.data.frame(scale(x)),a[,c(-1, -2)], calc.FRic= T, calc.FDiv= T)$CWM
    data_cwm <- data_cwm[row.names(sgenus),]
    data_cwm <- na.omit(data_cwm)
    eucl_func_arbre <- dist(data_cwm)

    # phylogenetic distance
    arbre_uni <- picante::unifrac(san,tree_ar)

    # soil distance matrix
    env_df <- dplyr::select(phyloseq::sample_data(fungocommun), pH, C, N, no3_med, nh4_med)
    tmpdf <- data.frame(plot=phyloseq::sample_data(fungocommun)$Plot,
                        trt=phyloseq::sample_data(fungocommun)$Trtment)
    tmpdf$trt <- stringr::str_replace(tmpdf$trt,'^1(.*?)$', '\\1')
    tmpdf$trt <- stringr::str_replace(tmpdf$trt,'^2n$', '12n')
    env_df <- purrr::map_df(env_df, scale)
    env_df <- as.data.frame(env_df)
    rownames(env_df) <-toupper(paste0(tmpdf$plot, tmpdf$trt))
    env_df <- as.matrix(env_df[sort(rownames(sgenus)),])
    eucl_env <- dist(env_df)

    return(list(bc_microbe_n = bray_ident,
                bc_microbe_n_b = bray_ident_nb,
                bc_microbe_n_f = bray_ident_nf,
                bc_microbe_b = bray_ident_b,
                bc_microbe_f = bray_ident_f,
                arbre = list(
                  bc_arbre = bc_arbre,
                  unifrac_phylo_arbre = arbre_uni,
                  eucl_func_arbre =eucl_func_arbre,
                  eucl_env = eucl_env)))
  }
  bc_df_b <- prep_bc_dfs(bakorare)
  bc_df_f <- prep_bc_dfs(fungorare)
  bc_df_n <- prep_bc_dfs_nem()
# my_list for bacteria
  my_list <- list()
  counter <- 1
  for(i in bc_df_b$arbre){
    my_list$b$tmp <- list(mantel_test_p = NULL)
    my_list$b$tmp$mantel_test_p <- ape::mantel.test(as.matrix(bc_df_b$bc_microbe), as.matrix(i))$p
    tmpdf <- data.frame(fung = c(as.matrix(bc_df_b$bc_microbe)),arbre = c( as.matrix(i)))
    tmpdf <- dplyr::filter(tmpdf, fung != 0, arbre != 0)
    mantel_reg <- phytools::multi.mantel(bc_df_b$bc_microbe, i, nperm = 100)
    my_list$b$tmp$plot <- ggplot2::ggplot(tmpdf, ggplot2::aes(x=arbre, y=fung)) +
      ggplot2::geom_point(alpha=0.1) +
      ggplot2::geom_abline(intercept = mantel_reg$coefficients[1], slope = mantel_reg$coefficients[2])
    my_list$b$tmp$mantel_reg_rsq <- mantel_reg$r.squared
    my_list$b$tmp$mantel_reg_coef <- mantel_reg$coefficients
    my_list$b$tmp$mantel_reg_probf <- mantel_reg$probF
    names(my_list$b)[counter]<-names(bc_df_b$arbre[counter])
    counter <- 1 + counter
  }
# my_list for fungi
  counter <- 1
  bc_df_n$bc_microbe_n
  for(i in bc_df_f$arbre){
    my_list$f$tmp <- list(mantel_test_p = NULL)
    my_list$f$tmp$mantel_test_p <- ape::mantel.test(as.matrix(bc_df_f$bc_microbe), as.matrix(i))$p
    tmpdf <- data.frame(fung = c(as.matrix(bc_df_f$bc_microbe)),arbre = c( as.matrix(i)))
    tmpdf <- dplyr::filter(tmpdf, fung != 0, arbre != 0)
    mantel_reg <- phytools::multi.mantel(bc_df_f$bc_microbe, i, nperm = 100)
    my_list$f$tmp$plot <- ggplot2::ggplot(tmpdf, ggplot2::aes(x=arbre, y=fung)) +
      ggplot2::geom_point(alpha=0.1) +
      ggplot2::geom_abline(intercept = mantel_reg$coefficients[1], slope = mantel_reg$coefficients[2])
    my_list$f$tmp$mantel_reg_rsq <- mantel_reg$r.squared
    my_list$f$tmp$mantel_reg_coef <- mantel_reg$coefficients
    my_list$f$tmp$mantel_reg_probf <- mantel_reg$probF
    names(my_list$f)[counter]<-names(bc_df_f$arbre[counter])
    counter <- 1 + counter
  }
# my_list for nematodes
  counter <- 1
  for(i in bc_df_n$arbre){
    my_list$n$tmp <- list(mantel_test_p = NULL)
    my_list$n$tmp$mantel_test_p <- ape::mantel.test(as.matrix(bc_df_n$bc_microbe_n), as.matrix(i))$p
    tmpdf <- data.frame(fung = c(as.matrix(bc_df_n$bc_microbe_n)),arbre = c( as.matrix(i)))
    tmpdf <- dplyr::filter(tmpdf, fung != 0, arbre != 0)
    mantel_reg <- phytools::multi.mantel(bc_df_n$bc_microbe_n, i, nperm = 100)
    my_list$n$tmp$plot <- ggplot2::ggplot(tmpdf, ggplot2::aes(x=arbre, y=fung)) +
      ggplot2::geom_point(alpha=0.1) +
      ggplot2::geom_abline(intercept = mantel_reg$coefficients[1], slope = mantel_reg$coefficients[2])
    my_list$n$tmp$mantel_reg_rsq <- mantel_reg$r.squared
    my_list$n$tmp$mantel_reg_coef <- mantel_reg$coefficients
    my_list$n$tmp$mantel_reg_probf <- mantel_reg$probF
    names(my_list$n)[counter]<-names(bc_df_f$arbre[counter])
    counter <- 1 + counter
  }
  #---- Correlation (rho) between nematode and bacteria and fungi --------------
  my_list$nbf$mantel_b <-  vegan::mantel(as.matrix(bc_df_n$bc_microbe_n_b), as.matrix(bc_df_n$bc_microbe_b))
  my_list$nbf$mantel_f <- vegan::mantel(as.matrix(bc_df_n$bc_microbe_n_f), as.matrix(bc_df_n$bc_microbe_f))
  tmpdf <- data.frame(nematode = c(c(as.matrix(bc_df_n$bc_microbe_n_b)), c(as.matrix(bc_df_n$bc_microbe_n_f))),
                      dissimilarity = c(c( as.matrix(bc_df_n$bc_microbe_b),c(as.matrix(bc_df_n$bc_microbe_f)))),
                      group = c(rep('bacteria', length(c(as.matrix(bc_df_n$bc_microbe_n_b)))),
                                rep('fungi', length(c(as.matrix(bc_df_n$bc_microbe_n_f))))))
  tmpdf <- dplyr::filter(tmpdf, nematode != 0, dissimilarity != 0)
  my_list$nbf$plot <- ggplot2::ggplot(tmpdf, ggplot2::aes(y=dissimilarity, x=nematode)) +
    ggplot2::geom_point(alpha=0.5)+
    facet_grid(.~group)

  #---- Correlation (rho) between fungi and bacteria ---------------------------
  bc_df_bako <- prep_bc_dfs(bakorare)
  bc_df_fung <- prep_bc_dfs(fungoraref)
  my_list$bf$mantel <- vegan::mantel(as.matrix(bc_df_bako$bc_microbe), as.matrix(bc_df_fung$bc_microbe))
  tmpdf <- data.frame(fungi = c(as.matrix(bc_df_fung$bc_microbe)),bacteria = c( as.matrix(bc_df_bako$bc_microbe)))
  tmpdf <- dplyr::filter(tmpdf, fungi != 0, bacteria != 0)
  my_list$bf$plot <- ggplot2::ggplot(tmpdf, ggplot2::aes(x=bacteria, y=fungi)) +
    ggplot2::geom_point(alpha=0.1)
  return(my_list)
}

#
# envdf_only_soil <- envdf.f[,c('id','C','N','pH','fun_div','no3_med','nh4_med',
#                               'percent_gymnosperm','mixity')]
# envdf_only_soil$mixity <- as.integer(envdf_only_soil$mixity)
# envdf.fna <- na.omit(envdf_only_soil)
#
#
# datatr_bakt <-decostand(sotu_prop_ident_ps_bakt,"hellinger")
# datatr_fung <-decostand(sotu_prop_ident_ps_fung,"hellinger")
#
# mod <- varpart(datatr_bakt[envdf.fna$id,], sa[envdf.fna$id,3:ncol(sa)],envdf.fna[,-1])
# mod
#
# mod <- varpart(datatr_fung[envdf.fna$id,], sa[envdf.fna$id,3:ncol(sa)],envdf.fna[,-1])
# mod
#
# my_list_test <- plot_beta_diversity()
