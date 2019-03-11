Using the agroPack package to generate the results
================

-   [Introduction](#introduction)
-   [loading the packages](#loading-the-packages)
-   [Give path to my document](#give-path-to-my-document)
-   [From sequencing data to richness](#from-sequencing-data-to-richness)
-   [From earthworm species occurence table to diversity metrics](#from-earthworm-species-occurence-table-to-diversity-metrics)
-   [Preparation of the list of variable](#preparation-of-the-list-of-variable)
-   [Statistical analysis](#statistical-analysis)
-   [Plots](#plots)

Introduction
------------

agroPack contains the data and the functions used to analyse the relationship between soil biodiversity and windbreaks in arable lands of southern Quebec

<a href="#top">Back to top</a>

loading the packages
--------------------

Throughout the code I use agroPack. Which can be installed with `devtools::install_github("maximerivest/agroPack")`. I also use `dplyr` and `ggplot2`, both available on `CRAN`.

``` r
library(agroPack)
library(ggplot2)
library(dplyr)
```

<a href="#top">Back to top</a>

Give path to my document
------------------------

I use two computers with exactly the some content in "Documents" but one is hosted on a different hardrive. Depending where I run this script I comment in/out either of the paths.

``` r
#cur_dir <- "/media/bigdrive/document_FedoraFall2018/"
cur_dir <- "~/Documents/"
```

<a href="#top">Back to top</a>

From sequencing data to richness
--------------------------------

I inserted data into the package. I tend to load raw data into the package and do the transformation here. So the following lines of codes transforms the OTU data.table into species diversity metric. I start from a `phyloseq` data structure, use phyloseq for some filtering and transformation and then I extract the otu\_table and match rownames. Finally, this data structure can be feed into a vegan function.

``` r
bako <- baktps
fungo <- fungips
#---- Rarefy to even depth and Compute shannon diversity -----------------------
set.seed(100)
fungo = phyloseq::filter_taxa(fungo, function(x) sum(x) > 1, TRUE)
fungo = phyloseq::rarefy_even_depth(fungo, sample.size = 2200)
```

    ## You set `rngseed` to FALSE. Make sure you've set & recorded
    ##  the random seed of your session for reproducibility.
    ## See `?set.seed`

    ## ...

    ## 5 samples removedbecause they contained fewer reads than `sample.size`.

    ## Up to first five removed samples are:

    ## 36288289290291   

    ## ...

    ## 61OTUs were removed because they are no longer 
    ## present in any sample after random subsampling

    ## ...

``` r
bako = phyloseq::filter_taxa(bako, function(x) sum(x) > 1, TRUE)
bako = phyloseq::rarefy_even_depth(bako, sample.size = 9500)
```

    ## You set `rngseed` to FALSE. Make sure you've set & recorded
    ##  the random seed of your session for reproducibility.
    ## See `?set.seed`
    ## 
    ## ...

    ## 4 samples removedbecause they contained fewer reads than `sample.size`.

    ## Up to first five removed samples are:

    ## 289290291348 

    ## ...

    ## 15OTUs were removed because they are no longer 
    ## present in any sample after random subsampling

    ## ...

``` r
shannon_bakt <- vegan::diversity(phyloseq::otu_table(bako),'shannon')
shannon_fung <- vegan::diversity(phyloseq::otu_table(fungo),'shannon')

# Bacteria
Sample_ID <- paste0(stringr::str_extract(phyloseq::sample_data(bako)$Site, "[1-2]"), "-",
                    toupper(stringr::str_extract(phyloseq::sample_data(bako)$Site, "[a-z]")),
                    phyloseq::sample_data(bako)$Plot, "-",
                    phyloseq::sample_data(bako)$Trtment)

sha_bako <- data.frame(row.names = phyloseq::sample_data(bako)$sample_id,
                       Sample_ID,
                       shannon_bakt)
# Fungi
Sample_ID <- paste0(stringr::str_extract(phyloseq::sample_data(fungo)$Site, "[1-2]"), "-",
                    toupper(stringr::str_extract(phyloseq::sample_data(fungo)$Site, "[a-z]")),
                    phyloseq::sample_data(fungo)$Plot, "-",
                    phyloseq::sample_data(fungo)$Trtment)
sha_fung <- data.frame(row.names = phyloseq::sample_data(fungo)$sample_id,
                       Sample_ID,
                       shannon_fung)
dplyr::left_join(main_data, sha_bako) %>%
  dplyr::left_join(sha_fung) -> main_data
```

    ## Joining, by = "Sample_ID"

    ## Warning: Column `Sample_ID` joining factors with different levels, coercing
    ## to character vector

    ## Joining, by = "Sample_ID"

    ## Warning: Column `Sample_ID` joining character vector and factor, coercing
    ## into character vector

<a href="#top">Back to top</a>

From earthworm species occurence table to diversity metrics
-----------------------------------------------------------

Same type of things as for bacteria and fungi but with a bit less filtering and without using `phyloseq` because earthworms data are based on actuall individual sampling rather than sequencing data.

``` r
# Earthworms
names(earthworm_community)[1] <- "Sample_ID"
mat_ea <- as.matrix(earthworm_community[,2:9])
row.names(mat_ea) <- earthworm_community[,1]
shannon_worm <- vegan::diversity(mat_ea,'shannon')
sha_worm <- data.frame(row.names = earthworm_community$Sample_ID,
                       Sample_ID = earthworm_community$Sample_ID,
                       shannon_worm)
main_data <- dplyr::left_join(main_data, sha_worm)
```

    ## Joining, by = "Sample_ID"

    ## Warning: Column `Sample_ID` joining character vector and factor, coercing
    ## into character vector

Preparation of the list of variable
-----------------------------------

Here, I will be preparing a list with all the data and parameters that I use for the statistical analysis. I am making very similar plots and statistics as I related several dependent variables to one independent variable, namely "distance from the windbreaks".

``` r
ls_dependent_variables_to_edge_distance <- list(
  water_content = list(varname = "Gravimetric humidity",
                       ylab = "Gravimetric humidity",
                       dataset = list(to_plot = main_data$Gravimetric_humidity,
                                      for_stat = main_data$Gravimetric_humidity)
  ),
  whc = list(varname = "Water holding capacity",
             ylab = "Water holding capacity",
             dataset = list(to_plot = main_data$Whc,
                            for_stat = main_data$Whc)
  ),
  total_c = list(varname = "Organic C (g/kg)",
                 ylab = expression(paste("Organic C (g ",kg^{-1},")")),
                 dataset = list(to_plot = main_data$Ctot,
                                for_stat = main_data$Ctot)
  ),
  total_n = list(varname = "Total nitrogen (g/kg)",
                 ylab = expression(paste("Total nitrogen (g ",kg^{-1},")")),
                 dataset = list(to_plot = main_data$Ntot,
                                for_stat = main_data$Ntot)
  ),
  c_n_ratio = list(varname = "C:N",
                   ylab = "C:N",
                 dataset = list(to_plot = main_data$C_N,
                                for_stat = main_data$C_N)
  ),
  pH = list(varname = "pH",
            ylab = "pH",
            dataset = list(to_plot = main_data$pH,
                           for_stat = main_data$pH)
  ),
  nem = list(varname = "log(Nematode density (no./kg))",
             ylab = expression(paste("Nematode density (no. ",kg^{-1},")")),
             dataset = list(to_plot = main_data$Nematode_count * 20,
                            for_stat = log((main_data$Nematode_count * 20)+1))
  ),
  ea_count = list(varname = "log(Earthworm density (no./m2))",
                  ylab = expression(paste("Earthworm density (no. ",m^{2},")")),
                  dataset = list(to_plot = main_data$Earthworm_count/0.09,
                                 for_stat = log((main_data$Earthworm_count/0.09)+1))
  ),
  ea_div = list(varname = "log(Earthworm species richness)",
                ylab = "Earthworm species richness",
                dataset = list(to_plot = main_data$Earthworm_div,
                               for_stat = log(main_data$Earthworm_div + 1))
  ),
  ea_biomass = list(varname = "log(Earthworm biomass (g/m2))",
                    ylab = expression(paste("Earthworm biomass g", m^{2}, ")")),
                    dataset = list(to_plot = main_data$Earthworm_biomass/0.09,
                                   for_stat = log((main_data$Earthworm_biomass/0.09)+1))
  ),
  bakt_div = list(varname = "Bacteria Shannon exponential",
                  ylab = expression(paste(e^{H*minute}," Bacteria")),
                  dataset = list(to_plot = main_data$shannon_bakt,
                                 for_stat = main_data$shannon_bakt)
  ),
  fung_div = list(varname = "Fungi Shannon exponential",
                  ylab = expression(paste(e^{H*minute}," fungi")),
                  dataset = list(to_plot = main_data$shannon_fung,
                                 for_stat = main_data$shannon_fung)
  )
)
```

<a href="#top">Back to top</a>

Statistical analysis
--------------------

Here, I use the function that I created to make mix effect linear regression, mixed effect Anova, and multiple comparaison of mix effect linear regression (with factors). Within this code chunk, I also make tables of results that I export to csv.

``` r
for(i in 1:length(ls_dependent_variables_to_edge_distance)){
  ls_dependent_variables_to_edge_distance[[i]]$stats <- 
    stats_for_soil_properties(ls_dependent_variables_to_edge_distance[[i]]$dataset$for_stat)
}

table_anova_mixed_effect <- function(my_list, varname = ""){
  res <- summary(aov(y ~ Treatment + Error(Field_block/Treatment), data=my_list$df))[[2]][[1]]
  class(res) <- "data.frame"
  res <- res[1,]
  res$variable_name <- varname
  return(res)
}

list_table_anova <- list(1:length(ls_dependent_variables_to_edge_distance))
for(i in 1:length(ls_dependent_variables_to_edge_distance)){
  list_table_anova[[i]] <- table_anova_mixed_effect(ls_dependent_variables_to_edge_distance[[i]]$stats, 
                                                  deparse(ls_dependent_variables_to_edge_distance[[i]]$varname))
}

all_anova_mixed <- bind_rows(list_table_anova)

all_anova_mixed %>%
  mutate(
    "Benjamini-Hochberg corrected p-value" = format(`Pr(>F)`/dplyr::percent_rank(`Pr(>F)`), scientific=F),
    "p-value" = format(`Pr(>F)`, scientific=F)
  ) -> all_anova_mixed
write.csv(all_anova_mixed, paste0(cur_dir, "/Msc/manuscript_nem_vdt_david_jan29/data/anova_distance_effect.csv"))


list_table_multcomp <- list(1:length(ls_dependent_variables_to_edge_distance))
for(i in 1:length(ls_dependent_variables_to_edge_distance)){
  list_table_multcomp[[i]] <- table_of_multiple_comparaison(ls_dependent_variables_to_edge_distance[[i]]$stats, 
                                                    deparse(ls_dependent_variables_to_edge_distance[[i]]$varname))
}

all_comps <- bind_rows(list_table_multcomp)
```

    ## Warning in bind_rows_(x, .id): Unequal factor levels: coercing to character

    ## Warning in bind_rows_(x, .id): binding character and factor vector,
    ## coercing into character vector

    ## Warning in bind_rows_(x, .id): binding character and factor vector,
    ## coercing into character vector

    ## Warning in bind_rows_(x, .id): binding character and factor vector,
    ## coercing into character vector

    ## Warning in bind_rows_(x, .id): binding character and factor vector,
    ## coercing into character vector

    ## Warning in bind_rows_(x, .id): binding character and factor vector,
    ## coercing into character vector

    ## Warning in bind_rows_(x, .id): binding character and factor vector,
    ## coercing into character vector

    ## Warning in bind_rows_(x, .id): binding character and factor vector,
    ## coercing into character vector

    ## Warning in bind_rows_(x, .id): binding character and factor vector,
    ## coercing into character vector

    ## Warning in bind_rows_(x, .id): binding character and factor vector,
    ## coercing into character vector

    ## Warning in bind_rows_(x, .id): binding character and factor vector,
    ## coercing into character vector

    ## Warning in bind_rows_(x, .id): binding character and factor vector,
    ## coercing into character vector

    ## Warning in bind_rows_(x, .id): binding character and factor vector,
    ## coercing into character vector

``` r
all_comps %>%
  mutate(
    "Benjamini-Hochberg corrected p-value" = format(`p-value`/dplyr::percent_rank(`p-value`), scientific=F),
    "p-value" = format(`p-value`, scientific=F)
  ) -> all_comps

write.csv(all_comps, paste0(cur_dir, "/Msc/manuscript_nem_vdt_david_jan29/data/multiple_comparaisons.csv"))
```

<a href="#top">Back to top</a>

Plots
-----

Here, I use the function that I created to make boxplots of the dependent variable as a function of windbreak distance. Within this code chunk, I also save the plots as png.

``` r
for(i in 1:length(ls_dependent_variables_to_edge_distance)){
  ls_dependent_variables_to_edge_distance[[i]]$boxplot_dist <- 
    plot_for_soil_properties(ls_dependent_variables_to_edge_distance[[i]]$stats, 
                             my_ylab = ls_dependent_variables_to_edge_distance[[i]]$ylab)
}

for(i in 1:length(ls_dependent_variables_to_edge_distance)){
  png(paste0(cur_dir,
             'Msc/manuscript_nem_vdt_david_jan29/figure/',
             stringr::str_replace(ls_dependent_variables_to_edge_distance[[i]]$varname,"/", " "),
             '.png'),
      width = 8.4, height = 8.4, units = 'cm', res = 300)
  print(ls_dependent_variables_to_edge_distance[[i]]$boxplot_dist)
  dev.off()
  
  setEPS()
  postscript(paste0(cur_dir,
                    'Msc/manuscript_nem_vdt_david_jan29/figure/',
                    stringr::str_replace(ls_dependent_variables_to_edge_distance[[i]]$varname,"/", " "),
                    '.png'), width = 3.30, height = 3.30)
  print(ls_dependent_variables_to_edge_distance[[i]]$boxplot_dist)
  dev.off()
}


# # Combine plots
# png(paste0(cur_dir,'Msc/manuscript_nem_vdt_david_jan29/figure/soil_combined.png'), width = 8.4, height = 8.4, units = 'cm', res = 300)
# gridExtra::grid.arrange(p_water_content,
#                         p_whc,
#                         p_c,
#                         p_n,
#                         p_c_n,
#                         p_ph, ncol = 3)
# dev.off()
# 
# setEPS()
# postscript(paste0(cur_dir,'Msc/manuscript_nem_vdt_david_jan29/figure/soil_combined.eps'), width = 3.3, height = 3.3)
# gridExtra::grid.arrange(p_water_content,
#                         p_whc,
#                         p_c,
#                         p_n,
#                         p_c_n,
#                         p_ph, ncol = 3)
# dev.off()
```

<a href="#top">Back to top</a>
