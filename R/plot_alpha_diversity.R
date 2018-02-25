#' Load environment
#'

plot_alpha_diversity <- function() {
  my_list <- list()
  fungo <- phyloseq::subset_samples(identps_bakt, Trtment != '1champ' & Trtment != '12n')
  bako <- phyloseq::subset_samples(identps_fungi, Trtment != '1champ' & Trtment != '12n')
#---- Rarefy to even depth and Compute shannon diversity -----------------------
  set.seed(100)
  fungorare = phyloseq::filter_taxa(fungo, function(x) sum(x) > 1, TRUE)
  fungorare = phyloseq::rarefy_even_depth(fungorare, sample.size = 1500)
  fungocommun = phyloseq::filter_taxa(fungo, function(x) sum(x >= 1) > 1, TRUE)
  fungocommun = phyloseq::rarefy_even_depth(fungocommun, sample.size = 1500)
  bakorare = phyloseq::filter_taxa(bako, function(x) sum(x) > 1, TRUE)
  bakorare = phyloseq::rarefy_even_depth(bakorare, sample.size = 1500)
  bakocommun = phyloseq::filter_taxa(bako, function(x) sum(x >= 1) > 1, TRUE)
  bakocommun = phyloseq::rarefy_even_depth(bakocommun, sample.size = 1500)

  shannon_baktrare <- vegan::diversity(phyloseq::otu_table(bakorare),'shannon')
  shannon_baktcommun <- vegan::diversity(phyloseq::otu_table(bakocommun),'shannon')
  shannon_fungrare <- vegan::diversity(phyloseq::otu_table(fungorare),'shannon')
  shannon_fungcommun <- vegan::diversity(phyloseq::otu_table(fungocommun),'shannon')

#----Relation arbre diversité phylogénétique - fungi/bakt  phylogénétique -----
  # Bakorare
  sha_bakorare <- data.frame(row.names = phyloseq::sample_data(bakorare)$id,
                             shannon_baktrare,
                             Plot = phyloseq::sample_data(bakorare)$Plot)
  pd_test_bakorare <- na.omit(cbind(sha_bakorare,
                                    pd_arbre[match(rownames(sha_bakorare), rownames(pd_arbre)),]))
  # Bakocommun
  sha_bakocommun <- data.frame(row.names = phyloseq::sample_data(bakocommun)$id,shannon_baktcommun,
                               Plot =phyloseq::sample_data(bakocommun)$Plot)
  pd_test_bakocommun <- na.omit(cbind(sha_bakocommun,
                                      pd_arbre[match(rownames(sha_bakocommun), rownames(pd_arbre)),]))
  # Fungorare
  sha_fungrare <- data.frame(row.names = phyloseq::sample_data(fungorare)$id,shannon_fungrare,
                             Plot = phyloseq::sample_data(fungorare)$Plot)
  pd_test_fungrare <- na.omit(cbind(sha_fungrare,
                                    pd_arbre[match(rownames(sha_fungrare), rownames(pd_arbre)),]))
  # Fungocommun
  sha_fungcommun <- data.frame(row.names = phyloseq::sample_data(fungocommun)$id,shannon_fungcommun,
                               Plot = phyloseq::sample_data(fungocommun)$Plot)
  pd_test_fungcommun <- na.omit(cbind(sha_fungcommun,
                                      pd_arbre[match(rownames(sha_fungcommun), rownames(pd_arbre)),]))

  # test linear relationship
  my_list$summary_lm_pd_br <- summary(nlme::lme(exp(shannon_baktrare) ~ PD, random = ~1| Plot,data=pd_test_bakorare))
  my_list$summary_lm_pd_bc <- summary(nlme::lme(data=pd_test_bakocommun, exp(shannon_baktcommun) ~ PD, random = ~1| Plot))
  my_list$summary_lm_pd_fr <- summary(nlme::lme(data=pd_test_fungrare, exp(shannon_fungrare) ~ PD, random = ~1| Plot))
  my_list$summary_lm_pd_fc <- summary(nlme::lme(data=pd_test_fungcommun, exp(shannon_fungcommun) ~ PD, random = ~1| Plot))
  # Nothing significant will plot tree richness and say that phylogenetic
  # did not produce different results

#----Relation arbre diversité functionel - fungi/bakt  phylogénétique ------
  data_fdis <- as.data.frame(FD::dbFD(as.data.frame(scale(x)),a[,c(-1, -2)], calc.FRic= T, calc.FDiv= T)$FDis)
  names(data_fdis) <- "fdis"

  # Bacteria rare
  sha_bakorare <- data.frame(row.names = phyloseq::sample_data(bakorare)$id,
                             shannon_baktrare,
                             Plot = phyloseq::sample_data(bakorare)$Plot)
  fdis_test_bakorare <- na.omit(cbind(sha_bakorare,
                                      fdis=data_fdis[match(rownames(sha_bakorare), rownames(data_fdis)),]))
  # Bacteria commun
  sha_bakocommun <- data.frame(row.names = sample_data(bakocommun)$id,
                               shannon_baktcommun,
                               Plot = sample_data(bakocommun)$Plot)
  fdis_test_bakocommun <- na.omit(cbind(sha_bakocommun,
                                        fdis=data_fdis[match(rownames(sha_bakocommun), rownames(data_fdis)),]))
  # Fungi rare
  sha_fungrare <- data.frame(row.names = sample_data(fungorare)$id,
                             shannon_fungrare,
                             Plot = sample_data(fungorare)$Plot)
  fdis_test_fungrare <- na.omit(cbind(sha_fungrare,
                                      fdis=data_fdis[match(rownames(sha_fungrare), rownames(data_fdis)),]))
  # Fungi commun
  sha_fungcommun <- data.frame(row.names = sample_data(fungocommun)$id,
                               shannon_fungcommun,
                               Plot = sample_data(fungocommun)$Plot)
  fdis_test_fungcommun <- na.omit(cbind(sha_fungcommun,
                                      fdis=data_fdis[match(rownames(sha_fungcommun), rownames(data_fdis)),]))

  # test linear relationship
  my_list$summary_lm_fdis_br <- summary(nlme::lme(data=fdis_test_bakorare, exp(shannon_baktrare) ~ fdis, random = ~1| Plot))
  my_list$summary_lm_fdis_bc <- summary(nlme::lme(data=fdis_test_bakocommun, exp(shannon_baktcommun) ~ fdis, random = ~1| Plot))
  my_list$summary_lm_fdis_fr <- summary(nlme::lme(data=fdis_test_fungrare, exp(shannon_fungrare) ~ fdis, random = ~1| Plot))
  my_list$summary_lm_fdis_fc <- summary(nlme::lme(data=fdis_test_fungcommun, exp(shannon_fungcommun) ~ fdis, random = ~1| Plot))

#----Relation arbre diversité spécifique - fungi/bakt ------
  # Bacteria rare
  sha_bakorare <- data.frame(row.names = phyloseq::sample_data(bakorare)$id,
                             shannon_baktrare,
                             Plot = phyloseq::sample_data(bakorare)$Plot,
                             sr = sample_data(bakorare)$DIV)
  # Bacteria commun
  sha_bakocommun <- data.frame(row.names = phyloseq::sample_data(bakocommun)$id,
                             shannon_baktcommun,
                             Plot = phyloseq::sample_data(bakocommun)$Plot,
                             sr = sample_data(bakocommun)$DIV)
  # Fungi rare
  sha_fungrare <- data.frame(row.names = phyloseq::sample_data(fungorare)$id,
                             shannon_fungrare,
                             Plot = phyloseq::sample_data(fungorare)$Plot,
                             sr = sample_data(fungorare)$DIV)
  # Fungi commun
  sha_fungcommun <- data.frame(row.names = phyloseq::sample_data(fungocommun)$id,
                               shannon_fungcommun,
                               Plot = phyloseq::sample_data(fungocommun)$Plot,
                               sr = sample_data(fungocommun)$DIV)
  # test linear relationship
  my_list$summary_lm_sr_br <- summary(nlme::lme(data=sha_bakorare, exp(shannon_baktrare) ~ sr, random = ~1| Plot))
  my_list$summary_lm_sr_bc <- summary(nlme::lme(data=sha_bakocommun, exp(shannon_baktcommun) ~ sr, random = ~1| Plot))
  my_list$summary_lm_sr_fr <- summary(nlme::lme(data=sha_fungrare, exp(shannon_fungrare) ~ sr, random = ~1| Plot))
  my_list$summary_lm_sr_fc <- summary(nlme::lme(data=sha_fungcommun, exp(shannon_fungcommun) ~ sr, random = ~1| Plot))
  MuMIn::r.squaredGLMM(nlme::lme(data=sha_bakorare, exp(shannon_baktrare) ~ sr, random = ~1| Plot))
  MuMIn::r.squaredGLMM(nlme::lme(data=sha_fungrare, exp(shannon_fungrare) ~ sr, random = ~1| Plot))

  #Nematodes
  sgenus <- nem_matrix_genus[sort(rownames(nem_matrix_genus)),]
  tokeep <- expand.grid(c('A','B','C','D'), c('ACRU','ACSA','BEAL','BEPA','PIGL','PIST','THOC','2NR2A','2N5','2NR5','2N8','4NR7','4N8'))
  tokeep <- paste0(tokeep[,1], tokeep[,2])
  sgenus <- na.omit(sgenus[tokeep,])
  df_nem <- data.frame(id = rownames(sgenus),
                       shannon = vegan::diversity(sgenus, index = 'shannon'),
                       stringsAsFactors = F)
  df_nem2 <- df_nem %>%
    mutate(Plot = substr(id, 1,1),
           sr = as.integer(stringr::str_replace_all(substr(id, 2,2), '[A-Z]', '1')))
  summary(nlme::lme(data=df_nem2, exp(shannon) ~ sr, random = ~1| Plot))
  MuMIn::r.squaredGLMM(nlme::lme(data=df_nem2, exp(shannon) ~ sr, random = ~1| Plot))
#---- Plot species richness ----------------------------------------------------
  alldf <- data.frame(treespdiv=c(df_nem2$sr,
                                  sample_data(fungorare)$DIV,
                                  sample_data(bakorare)$DIV),
                      shannon=c(df_nem2$shannon,
                                shannon_fungrare,
                                shannon_baktrare),
                      group=c(rep('nematode',nrow(df_nem2)),
                              rep('fungi',nrow(sample_data(fungorare))),
                              rep('bacteria',nrow(sample_data(bakorare)))),
                      Plot = c(df_nem2$Plot,
                               sample_data(fungorare)$Plot,
                               sample_data(bakorare)$Plot))
  alldf <- subset(alldf, treespdiv != 0)
  pt <- ggplot2::ggplot(dplyr::filter(alldf, treespdiv <5), ggplot2::aes(x=treespdiv, y=exp(shannon), color=group))+
    ggplot2::geom_jitter(alpha=0.6,width=0.1, height = 0)+
    ggplot2::stat_smooth(method='lm') +
    ggplot2::theme(legend.position="none")+
    ggplot2::xlab('Tree species \nrichness') +
    ggplot2::ylab('Shannon diversity of bacteria, fungi and nematode')

#---- Plot functional diversity ----------------------------------------------------
  alldf <- data.frame(treefundiv=c(fdis_test_fungrare$fdis,
                                  fdis_test_fungrare$fdis,
                                  fdis_test_bakorare$fdis),
                      shannon=c(sample_data(fungorare)$Hnem,
                                fdis_test_fungrare$shannon_fungrare,
                                fdis_test_bakorare$shannon_baktrare),
                      group=c(rep('nematode',nrow(sample_data(fungorare))),
                              rep('fungi',nrow(sample_data(fungorare))),
                              rep('bacteria',nrow(fdis_test_bakorare))),
                      Plot = c(rep(sample_data(fungorare)$Plot,2),
                               sample_data(bakorare)$Plot))
  # change all zeroes for ones.
  alldf <- dplyr::mutate(alldf, treefundiv = replace( treefundiv, treefundiv == 0, 1))
  ptfun <- ggplot2::ggplot(alldf, ggplot2::aes(x=treefundiv, y=exp(shannon), color=group))+
    #ggplot2::geom_jitter(alpha=0.2,width=0.1)+
    ggplot2::geom_point()+
    ggplot2::stat_smooth(method='lm') +
    ggplot2::theme(legend.position="none")+
    ggplot2::xlab('Tree functional \ndiversity') +
    ggplot2::ylab('Shannon diversity of bacteria, fungi and nematode')

  my_list$summary_lm_fdis1_fr <- summary(nlme::lme(exp(shannon) ~ treefundiv, random = ~1|Plot, alldf[alldf$group == 'fungi',]))
  my_list$summary_lm_fdis1_br <- summary(nlme::lme(exp(shannon) ~ treefundiv, random = ~1|Plot, alldf[alldf$group == 'bacteria',]))
  my_list$summary_lm_fdis1_n <- summary(nlme::lme(exp(shannon) ~ treefundiv, random = ~1|Plot, na.omit(alldf[alldf$group == 'nematode',])))
#----- Plot nematode against bako fung -----------------------------------------
  nfb <- data.frame(nematode=c(sample_data(fungorare)$Hnem,
                               sample_data(bakorare)$Hnem),
                    shannon=c(shannon_fungrare,shannon_baktrare),
                    group=c(rep('fungi',nrow(sample_data(fungorare))),
                            rep('bacteria',nrow(sample_data(bakorare)))),
                    Plot = c(sample_data(fungorare)$Plot,
                             sample_data(bakorare)$Plot))
  pn <- ggplot2::ggplot(nfb, ggplot2::aes(x=exp(nematode), y=exp(shannon), color=group))+
    ggplot2::geom_point(alpha = 0.6)+
    ggplot2::stat_smooth(method='lm') +
    ggplot2::theme(legend.position="none")+
    ggplot2::xlab('Shannon diversity of \nnematode genus') +
    ggplot2::ylab('Shannon diversity of bacteria and fungi')+
    ggplot2::scale_colour_discrete(drop=TRUE,
                                   limits = levels(alldf$group))

  nfb <- na.omit(nfb)
  cor(x=nfb$nematode, y=nfb$shannon)
  cor.test(x=nfb$nematode[nfb$group == 'fungi'], y=nfb$shannon[nfb$group == 'fungi'])
  cor.test(x=nfb$nematode[nfb$group == 'bacteria'], y=nfb$shannon[nfb$group == 'bacteria'])

  my_list$summary_lm_nem_fr <- summary(nlme::lme(exp(shannon) ~ nematode, random = ~1|Plot, nfb[nfb$group == 'fungi',]))
  my_list$summary_lm_nem_br <- summary(nlme::lme(exp(shannon) ~ nematode, random = ~1|Plot, nfb[nfb$group == 'bacteria',]))

  ## Put block-treatment rowname
  ##------------------------------------------------------------------------------
  bakt_fungi <- data.frame(shannon_baktrare =shannon_baktrare,
                           shannon_fungrare =shannon_fungrare[names(shannon_baktrare)])
  bakt_fungi$group <- 'fungi'
  bakt_fungi$group <- as.factor(bakt_fungi$group)


  (cortest <- cor.test(x=exp(shannon_baktcommun), y=exp(shannon_fungcommun[names(shannon_baktcommun)])))
  my_list$cortest_rare_bf <- (cortest <- cor.test(x=exp(shannon_baktrare), y=exp(shannon_fungrare[names(shannon_baktrare)])))

  pbf <- ggplot2::ggplot(bakt_fungi, ggplot2::aes(x=exp(shannon_baktrare), y=exp(shannon_fungrare), color = group))+
    ggplot2::geom_point(alpha=0.6)+
    ggplot2::stat_smooth(method='lm') +
    ggplot2::xlab('Shannon diversity \nof bacteria') +
    ggplot2::ylab('Shannon diversity of fungi') +
    ggplot2::annotate('text',x=50,y=1,label=paste('r = ',
                                                 round(cortest$estimate,digits=2),
                                                 ' ± ',
                                                 round(cortest$estimate-cortest$conf.int[1],digits=3),
                                                 '***',
                                                 sep='')
    ) +
    ggplot2::scale_colour_discrete(drop=TRUE,
                                   limits = levels(alldf$group))
  return(c(my_list, list('bfn_fundiv' = ptfun,'bfn_div' = pt,'bf_n' = pn,'b_f' = pbf)))
}
