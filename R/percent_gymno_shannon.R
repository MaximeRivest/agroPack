#' analyse and plot shannon diversity as a function of gymnosperm concentration

percent_gymno_shannon <- function(){
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

  my_list <- list()
  bkcf <- as(sample_data(fungorare), Class='data.frame')
  bkcf$shannon_bc <- shannon_fungrare

  #bkc <- as(sample_data(fungocommun), Class='data.frame')
  #bkc$shannon_bc <- shannon_fungcommun

  bkcb <- as(sample_data(bakorare), Class='data.frame')
  bkcb$shannon_bc <- shannon_baktrare

  #bkc <- as(sample_data(bakocommun), Class='data.frameexpression("Nematode's " * alpha * "-diversity")')
  #bkc$shannon_bc <- shannon_baktcommun

  #bkcf$percent_gymnosperm <- as.factor(bkcf$percent_gymnosperm)
  #bkcb$percent_gymnosperm <- as.factor(bkcb$percent_gymnosperm)
  #lme1 <- nlme::lme(exp(shannon_bc) ~ percent_gymnosperm, random  = ~1 | Plot, data=dplyr::filter(bkcf, shannon_bc < 5.5))
  lme1 <- nlme::lme(exp(shannon_bc) ~ percent_gymnosperm, random  = ~1 | Plot, data=bkcf)
  #anova(lme1)
  my_list$summary_lme1f <- summary(lme1)
  my_list$rsquaref <- MuMIn::r.squaredGLMM(lme1)
  lme1 <- nlme::lme(exp(shannon_bc) ~ percent_gymnosperm, random  = ~1 | Plot, data=bkcb)
  #anova(lme1)
  my_list$summary_lme1b <- summary(lme1)
  my_list$rsquareb <- MuMIn::r.squaredGLMM(lme1)
  #comp <- multcomp::glht(lme1, linfct=multcomp::mcp(percent_gymnosperm = "Tukey"))
  #print(summary(comp))
  bkcn <- bkcf[!is.na(bkcf$Hnem),]
  lme1 <- nlme::lme(exp(Hnem) ~ percent_gymnosperm, random = ~1 | Plot, data = bkcn)
  my_list$summary_lme1n <- summary(lme1)
  my_list$rsquaren <- MuMIn::r.squaredGLMM(lme1)

  bkcf$group <- 'fungi'
  bkcb$group <- 'bacteria'
  bkcn$group <- 'nematodes'

  bkc <- rbind(bkcb, bkcf, bkcn)
  bkc$group <- as.factor(bkc$group)
  my_list$plot <- ggplot2::ggplot(bkc, ggplot2::aes(x=percent_gymnosperm, y = exp(shannon_bc), color = group)) +
    ggplot2::geom_jitter(alpha=0.6, width = 5, height = 0)+
    ggplot2::stat_smooth(method='lm')
  return(my_list)
}
