#' Load environment
#'

plot_alpha_diversity <- function() {
  my_list <- list()
  bako <- baktps
  fungo <- fungips
  #---- Rarefy to even depth and Compute shannon diversity -----------------------
  set.seed(100)
  fungorare = phyloseq::filter_taxa(fungo, function(x) sum(x) > 1, TRUE)
  fungorare = phyloseq::rarefy_even_depth(fungorare, sample.size = 2200)
  fungocommun = phyloseq::filter_taxa(fungo, function(x) sum(x >= 1) > 1, TRUE)
  fungocommun = phyloseq::rarefy_even_depth(fungocommun, sample.size = 2200)
  bakorare = phyloseq::filter_taxa(bako, function(x) sum(x) > 1, TRUE)
  bakorare = phyloseq::rarefy_even_depth(bakorare, sample.size = 9500)
  bakocommun = phyloseq::filter_taxa(bako, function(x) sum(x >= 1) > 1, TRUE)
  bakocommun = phyloseq::rarefy_even_depth(bakocommun, sample.size = 9500)

  shannon_baktrare <- vegan::diversity(phyloseq::otu_table(bakorare),'shannon')
  shannon_baktcommun <- vegan::diversity(phyloseq::otu_table(bakocommun),'shannon')
  shannon_fungrare <- vegan::diversity(phyloseq::otu_table(fungorare),'shannon')
  shannon_fungcommun <- vegan::diversity(phyloseq::otu_table(fungocommun),'shannon')

  #----Relation arbre diversité spécifique - fungi/bakt ------
  # Bacteria rare
  sha_bakorare <- data.frame(row.names = phyloseq::sample_data(bakorare)$sample_id,
                             shannon_baktrare,
                             Plot = phyloseq::sample_data(bakorare)$Plot,
                             sr = sample_data(bakorare)$Trtment)
  # Bacteria commun
  sha_bakocommun <- data.frame(row.names = phyloseq::sample_data(bakocommun)$sample_id,
                               shannon_baktcommun,
                               Plot = phyloseq::sample_data(bakocommun)$Plot,
                               sr = sample_data(bakocommun)$Trtment)
  # Fungi rare
  sha_fungrare <- data.frame(row.names = phyloseq::sample_data(fungorare)$sample_id,
                             shannon_fungrare,
                             Plot = phyloseq::sample_data(fungorare)$Plot,
                             sr = sample_data(fungorare)$Trtment)
  # Fungi commun
  sha_fungcommun <- data.frame(row.names = phyloseq::sample_data(fungocommun)$sample_id,
                               shannon_fungcommun,
                               Plot = phyloseq::sample_data(fungocommun)$Plot,
                               sr = sample_data(fungocommun)$Trtment)
  # test linear relationship
  my_list$summary_lm_sr_br <- summary(nlme::lme(data=sha_bakorare, exp(shannon_baktrare) ~ sr, random = ~1| Plot))
  my_list$summary_lm_sr_bc <- summary(nlme::lme(data=sha_bakocommun, exp(shannon_baktcommun) ~ sr, random = ~1| Plot))
  my_list$summary_lm_sr_fr <- summary(nlme::lme(data=sha_fungrare, exp(shannon_fungrare) ~ sr, random = ~1| Plot))
  my_list$summary_lm_sr_fc <- summary(nlme::lme(data=sha_fungcommun, exp(shannon_fungcommun) ~ sr, random = ~1| Plot))
  MuMIn::r.squaredGLMM(nlme::lme(data=sha_bakorare, exp(shannon_baktrare) ~ sr, random = ~1| Plot))
  MuMIn::r.squaredGLMM(nlme::lme(data=sha_fungrare, exp(shannon_fungrare) ~ sr, random = ~1| Plot))

  my_list$my_plot_br <- ggplot2::ggplot(sha_bakorare,
                             ggplot2::aes(x = sr,
                                          y = exp(shannon_baktrare)
                                          #,color = as.factor(Treatment_meter)
                             ))+
    ggplot2::geom_boxplot(alpha = 1)+
    ggplot2::xlab("Distance from windbreak (m)") +
    ggplot2::ylab(expression(paste(e^{H*minute}," Bacteria")))+
    theme_Publication() +
    ggplot2::stat_summary(fun.y=mean,
                          colour="black",
                          geom="point",
                          shape=1,
                          size=2,
                          show.legend = FALSE) +
    ggplot2::theme(legend.position = "none")

  my_list$my_plot_br <- ggplot2::ggplot(sha_fungrare,
                                        ggplot2::aes(x = sr,
                                                     y = exp(shannon_fungrare)
                                                     #,color = as.factor(Treatment_meter)
                                        ))+
    ggplot2::geom_boxplot(alpha = 1)+
    ggplot2::xlab("Distance from windbreak (m)") +
    ggplot2::ylab(expression(paste(e^{H*minute}," Fungi")))+
    theme_Publication() +
    ggplot2::stat_summary(fun.y=mean,
                          colour="black",
                          geom="point",
                          shape=1,
                          size=2,
                          show.legend = FALSE) +
    ggplot2::theme(legend.position = "none")
}
