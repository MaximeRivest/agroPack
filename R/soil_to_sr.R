#' analyse and plot shannon diversity as a function of gymnosperm concentration

soil_to_sr <- function(){
  bako <- phyloseq::subset_samples(identps_bakt, Trtment != '1champ' & Trtment != '12n')
  fungo <- phyloseq::subset_samples(identps_fungi, Trtment != '1champ' & Trtment != '12n')
  #---- Rarefy to even depth and Compute shannon diversity -----------------------
  set.seed(100)
  fungorare = phyloseq::filter_taxa(fungo, function(x) sum(x) > 1, TRUE)
  fungorare = phyloseq::rarefy_even_depth(fungorare, sample.size = 1500)
  fungocommun = phyloseq::filter_taxa(fungo, function(x) sum(x >= 1) > 1, TRUE)
  fungocommun = phyloseq::rarefy_even_depth(fungocommun, sample.size = 1500)
  bakorare = phyloseq::filter_taxa(bako, function(x) sum(x) > 1, TRUE)
  bakorare = phyloseq::rarefy_even_depth(bakorare, sample.size = 7000)
  bakocommun = phyloseq::filter_taxa(bako, function(x) sum(x >= 1) > 1, TRUE)
  bakocommun = phyloseq::rarefy_even_depth(bakocommun, sample.size = 7000)

  shannon_baktrare <- vegan::diversity(phyloseq::otu_table(bakorare),'shannon')
  shannon_baktcommun <- vegan::diversity(phyloseq::otu_table(bakocommun),'shannon')


  shannon_fungrare <- vegan::diversity(phyloseq::otu_table(fungorare),'shannon')
  shannon_fungcommun <- vegan::diversity(phyloseq::otu_table(fungocommun),'shannon')

  my_list <- list()
  bkcf <- as(sample_data(fungorare), Class='data.frame')
  bkcf$shannon_bc <- shannon_fungrare

  bkcb <- as(sample_data(bakorare), Class='data.frame')
  bkcb$shannon_bc <- shannon_baktrare

  bkcf_no_NA <- filter(bkcf, !is.na(percent_gymnosperm) & !is.na(C) & !is.na(pH) & !is.na(no3_med)) %>%
    select(shannon_bc, percent_gymnosperm, C, N, pH, no3_med, nh4_med, Plot, DIV, fun_div)
  lme1 <- nlme::lme(exp(shannon_bc) ~ C + pH + no3_med, random  = ~1 | Plot,data=bkcf_no_NA)
  my_list$summary_lme1f <- summary(lme1)
  my_list$rsquaref <- MuMIn::r.squaredGLMM(lme1)

  bkcb_no_NA <- filter(bkcb, !is.na(percent_gymnosperm) & !is.na(C) & !is.na(pH) & !is.na(no3_med)) %>%
    select(shannon_bc, C, N, pH, no3_med, nh4_med, Plot, DIV, fun_div)
  lme1 <- nlme::lme(exp(shannon_bc) ~  C + pH + no3_med, random  = ~1 | Plot,data=bkcb_no_NA)
  my_list$summary_lme1b <- summary(lme1)
  my_list$rsquareb <- MuMIn::r.squaredGLMM(lme1)

   return(my_list)
}
