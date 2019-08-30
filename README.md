IBD肠道微生物整合代谢网络
================

``` r
library(curatedMetagenomicData)
library(tidyverse)

library(readxl)
library(clusterProfiler)
library(igraph)
library(mmnet)
library(furrr)
library(patchwork)
```

``` r
source_files <- list.files("R", "R$", full.names = TRUE)
for (file in source_files) source(file)
```

数据预处理
----------

数据源自Nielsen等人的文章[1]，用R包[curatedMetagenomicData](https://github.com/waldronlab/curatedMetagenomicData)获取功能注释数据。我们选择spanish样本进行分析，去除不一致样本，包括肥胖（BMI &gt;= 30)、同一个样本中第二次测序的数据、使用发酵乳产品 。

``` r
# 从nielsen文章中下载样本元数据nielsen_sample.xls
samples <- read_xls("data/nielsen_sample.xls")
spanish_sample <- filter(samples, Nationality == "spanish") %>% 
  arrange(`Sample ID`)
# 78 samples were sequenced twice, delete the sequence samples sequenced in the second time
individual <- unique(spanish_sample$`Individual ID`)
indx <- match(spanish_sample$`Individual ID`, individual)
selected_sample <- slice(spanish_sample, which(!duplicated(indx)))
selected_sample <- mutate(
  selected_sample, 
  `Sample ID` = str_replace_all(`Sample ID`, "\\.", "_")
)
# All samples did't consumed a defined fermented milk product containing the 
# previously sequenced Bifidobacterium animalis subsp. lactis CNCM I-2494
# all(selected_sample$`Known consumers of a defined fermented milk product (DFMP)` == "NA")
selected_sample$`Known consumers of a defined fermented milk product (DFMP)` <- NULL
# remove the sample with no bmi, and obese (BMI >= 30)
selected_sample <- filter(selected_sample, BMI != "NA" & BMI < 30)
```

用[curatedMetagenomicData](https://github.com/waldronlab/curatedMetagenomicData)下载功能注释数据

``` r
# set ExperimentHub cache dir
setExperimentHubOption("CACHE", "~/ExperimentHub/")
nile_ds <- curatedMetagenomicData(
  "NielsenHB_2014.genefamilies_relab.stool",
  counts = TRUE,
  dryrun = FALSE
)
```

    ## Working on NielsenHB_2014.genefamilies_relab.stool

    ## snapshotDate(): 2019-04-29

    ## see ?curatedMetagenomicData and browseVignettes('curatedMetagenomicData') for documentation

    ## downloading 0 resources

    ## loading from cache 
    ##     'EH1287 : 1287'

``` r
nile_ds <- nile_ds[[1]]
# subset the functional profilling according to the sample id
sample_ids <- selected_sample$`Sample ID`
sample_indx <- match(sampleNames(nile_ds), sample_ids)
selected_ds <- nile_ds[, !is.na(sample_indx)]
features <- featureNames(selected_ds)
gene_family_abundance <- exprs(selected_ds)
sample_meta <- pData(selected_ds)
```

[curatedMetagenomicData](https://github.com/waldronlab/curatedMetagenomicData)默认使用[humann2](https://bitbucket.org/biobakery/humann2/src/default/humann2/)参考Uniref90进行功能注释，我们使用humann2将其映射至KEGG orthology。

``` r
humann2_regroup_in <- as.data.frame(
  gene_family_abundance, 
  row.names = row.names(gene_family_abundance)
) %>% 
  rownames_to_column(var = "# Gene Family")

# save as humann2_regoup.py input
# write.table(humann2_regroup_in, "output/data/humann2_regroup_in.tsv", 
#   sep = "\t", quote = FALSE, row.names = FALSE
# )
# regoup using humann2
# humann2_regroup_table --input humann2_regroup_in.tsv \
# --output humann2_regoup_out.tsv --groups uniref90_ko

humann2_regroup_out <- read_tsv("output/data/humann2_regroup_out.tsv") %>% 
  slice(-1)
```

    ## Parsed with column specification:
    ## cols(
    ##   .default = col_double(),
    ##   `# Gene Family` = col_character()
    ## )

    ## See spec(...) for full column specifications.

``` r
index_remove <- select(humann2_regroup_out, -`# Gene Family`) %>% 
  rowSums(.) > 0
humann2_regroup_out <- slice(humann2_regroup_out, which(index_remove)) %>% 
  column_to_rownames("# Gene Family")
```

KO丰度的归一化，为了保证结果的可靠性，丰度在2以下的KO设为0。

``` r
# remove low confidence (abundance) kos， set to 0 while lower than 2, low prevelance
kos_table <- humann2_regroup_out
kos_table[kos_table < 2] = 0
kos_table <- kos_table[rowSums(kos_table) > 0, ]
kos_norm <- sweep(kos_table, 2, sample_meta$number_reads * 100, "/")
state <- sample_meta$disease
```

KO相对丰度在127个样本之间是一致的（平均相关系数为0.83, Spearman相关检验），健康样本间的KO相对丰度的一致性（平均相关系数0.86，Spearman相关检验）高于IBD样本间的一致性（平均相关系数为0.81，Spearman相关检验）

``` r
# ko correaltion across samples
load("data/sysdata.rda") # reference metabolic network based on KEGG metabolic pathway from mmnet
refnet <- RefDbcache$network
refnode <- V(refnet)$name
ann_kos_table <- kos_norm[row.names(kos_norm) %in% refnode, ]
sample_cor <- cor(ann_kos_table, method = "spearman")
ibd_cor <- sample_cor[state == "IBD", state == "IBD"]
ibd_cor[lower.tri(ibd_cor)] %>% mean()
```

    ## [1] 0.8116212

``` r
health_cor <- sample_cor[state != "IBD", state != "IBD"]
health_cor[lower.tri(health_cor)] %>% mean()
```

    ## [1] 0.8547604

``` r
inter_cor <- sample_cor[state != "IBD", state == "IBD"]
mean(inter_cor)
```

    ## [1] 0.8284058

整合代谢网络
------------

我们参考[mmnet](https://github.com/yiluheihei/mmnet)构建整合代谢网络。节点表示KO，如果a KO参与的代谢反应的产物可作为 b KO 参与代谢反应的底物，那么存在由a到b的连边。

``` r
# load("data/sysdata.rda") # reference metabolic network based on KEGG metabolic pathway from mmnet
# refnet <- RefDbcache$network
# refnode <- V(refnet)$name
subnodes <- intersect(refnode, row.names(kos_norm))

# construct integrated metabolic network
biom.data <- make_biom(kos_norm, observation_metadata = row.names(kos_norm))
biom.data$type <- "enzymatic genes abundance"
ssns <- constructSSN(biom.data)
diff_ssn <- diff_net(ssns, sample.state = state)
```

根据比值比OR预测IBD相关KO
-------------------------

OR&gt;2表示KO在IBD中显著富集，OR&lt;0.5表示KO在IBD中显著减少。共发现366个IBD相关KO，其中254个高表达，112个低表达

``` r
or <- vertex_attr(diff_ssn, "OR")
or_p <- vertex_attr(diff_ssn, "OR_p")
enrich_indx <- which(or > 2 & or_p < 0.05) 
deplete_indx <-which(or < 0.5 & or_p < 0.05) 
enrich_ko <- V(diff_ssn)$name[enrich_indx]
deplete_ko <- V(diff_ssn)$name[deplete_indx]
```

差异KO进行KEGG代谢通路富集分析

``` r
enrich_pathway <- enrichKEGG(gene = enrich_ko, organism = 'ko', pvalueCutoff = 0.05, qvalueCutoff = 0.05, 
  pAdjustMethod = "fdr")
deplete_pathway <- enrichKEGG(gene = deplete_ko, organism = "ko", pvalueCutoff = 0.05, qvalueCutoff = 0.05, pAdjustMethod = "fdr")
# remove low confidence pathway
enrich_pathway@result <- dplyr::filter(enrich_pathway@result, p.adjust < 0.05)
enrich_pathway@result <- dplyr::filter(enrich_pathway@result, Description != "Carbon metabolism") # 高表达和低表达都富集了，应该在低表达的KO中富集
enrich_pathway@result <- dplyr::arrange(enrich_pathway@result, p.adjust)
deplete_pathway@result <- dplyr::filter(deplete_pathway@result, p.adjust < 0.05)
deplete_pathway@result <- dplyr::filter(deplete_pathway@result, Description != "Glycolysis / Gluconeogenesis")

p_pathway_enrich <- enrichplot::dotplot(enrich_pathway, showCategory = 10) + 
  enrichplot::dotplot(deplete_pathway, showCategory = 10) +
  plot_layout(ncol = 1)
```

    ## wrong orderBy parameter; set to default `orderBy = "x"`
    ## wrong orderBy parameter; set to default `orderBy = "x"`

``` r
p_pathway_enrich
```

<img src="README_files/figure-markdown_github/pathway-enrich-1.png" width=".8\linewidth" />

``` r
# ggsave(p_pathway_enrich, filename = "output/pathway_enrich.pdf", width = 8, height = 9)
```

IBD相关KO与网络拓扑属性的关联
-----------------------------

介数中心度表示网络中心性，KO差异值与中心度负相关，IBD相关的KO中心度比其他KO的中心度低。

``` r
bc <- betweenness(diff_ssn, normalized = TRUE)
diff_score <- abs(log2(or))
cor(diff_score, bc, method = "spearman") # -0.1629
```

    ## [1] -0.1629375

``` r
wilcox.test(diff_score, bc) # 2.2e-16
```

    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  diff_score and bc
    ## W = 1418468, p-value < 2.2e-16
    ## alternative hypothesis: true location shift is not equal to 0

``` r
diff_bc <- bc[c(enrich_indx, deplete_indx)]
other_bc <- bc[-c(enrich_indx, deplete_indx)]
wilcox.test(diff_bc, other_bc, alternative = "less") 
```

    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  diff_bc and other_bc
    ## W = 141324, p-value = 1.142e-07
    ## alternative hypothesis: true location shift is less than 0

``` r
topo_tibble <- generate_topos_tibble(diff_ssn)

se <- function(x) sd(x)/sqrt(length(x))
summary_topo_tibble <- group_by(topo_tibble, method, state) %>% 
  summarise(
    m_value = mean(value),
    sd_value = sd(value),
    se = se(value)
  )

p_centrality <- ggplot(dplyr::filter(summary_topo_tibble, method == "Betweenness centrality"), 
  aes(state, m_value, fill = state)) +
  geom_col(width = 0.6) +
  geom_errorbar(aes(ymin = m_value - se, ymax = m_value + se), width = 0.2) +
  labs(x=NULL, y = "Beteenness Centrality") + 
  ggsci::scale_fill_npg() +
  theme_bw() +
  theme(legend.position = "none")
p_centrality
```

![](README_files/figure-markdown_github/p-centrality-1.png)

``` r
#ggsave(p_centrality, "output/p_centrality.pdf", width = 4, height = 3)
```

与介数中心度类似，IBD相关KO的度显著小于其KO的度；相反，IBD相关KO的聚类系数显著高于其他KO的聚类系数

``` r
cc <- transitivity(diff_ssn, type ="weighted", vids = V(diff_ssn), isolates=c("zero"))
diff_cc <- cc[c(enrich_indx, deplete_indx)]
other_cc <- cc[-c(enrich_indx, deplete_indx)]
wilcox.test(diff_cc, other_cc, alternative = "greater") # < 0.0003
```

    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  diff_cc and other_cc
    ## W = 193226, p-value = 0.0004353
    ## alternative hypothesis: true location shift is greater than 0

``` r
dg <- degree(diff_ssn, mode = "all", normalized = TRUE)
diff_dg <- dg[c(enrich_indx, deplete_indx)]
other_dg <- dg[-c(enrich_indx, deplete_indx)]
wilcox.test(diff_dg, other_dg, alternative = "less")
```

    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  diff_dg and other_dg
    ## W = 157296, p-value = 0.005419
    ## alternative hypothesis: true location shift is less than 0

``` r
p_other_topo <- ggplot(summary_topo_tibble, aes(state, m_value, fill = state)) + 
  geom_col(width = 0.6) +
  geom_errorbar(aes(ymin = m_value - se, ymax = m_value + se), width = 0.2) +
  facet_wrap(~method, scales = "free_y") +
  labs(x=NULL, y = "Centrality") + 
  ggsci::scale_fill_npg() +
  theme_bw() +
  theme(legend.position = "none")
p_other_topo
```

![](README_files/figure-markdown_github/p-other-topo-1.png)

``` r
# ggsave(p_other_topo, "output/p_other_topo.pdf", width = 7, height = 5)
```

IBD和healthy网络比较
--------------------

把所有样本分成IBD和healthy两组，分别构建整合代谢网络，比较其拓扑变化。

``` r
ibd_indx <- state == "IBD"
ibd_table <- kos_norm[, state == "IBD"]
ibd_table <- ibd_table[rowSums(ibd_table) > 0, ]
health_table <- kos_norm[, state != "IBD"]
health_table <- health_table[rowSums(health_table) > 0, ]

ibd_kos <- intersect(refnode, row.names(ibd_table))
health_kos <- intersect(refnode, row.names(health_table))
health_net <- induced_subgraph(refnet, health_kos) 
ibd_net <- induced_subgraph(refnet, ibd_kos) 
```

IBD减少网络的模块性和密度，可能是由于IBD引起肠道微生物多样性和功能减少引起的。随机打乱样本标签（IBD或者healthy）1000次，得到1000对网络作为零分布用于计算网络模块和密度差异是否显著。我们发现真正的healthy-network比零分布中大多数health-network的模块性强（87.5%），而真正的IBD-network比零分布中大多数IBD-network模块性弱（76.6%），同时healthy-network与IBD-network之间的模块化程度差异显著。

``` r
ibd_modularity <- cluster_leading_eigen(as.undirected(ibd_net)) %>% 
  modularity()
health_modularity <- cluster_leading_eigen(as.undirected(health_net)) %>% 
  modularity()
diff_modularity <- health_modularity - ibd_modularity

# modularity shuffled 1000 times
shuffled_modularity <- map(1:1000, ~ shuffle_modularity(kos_norm, state, refnet) %>% unlist()) %>% 
  do.call(rbind, .) %>% 
  as.data.frame()

shuffled_diff_modularity <- shuffled_modularity$health_modularity - shuffled_modularity$disease_modularity
sum(health_modularity > shuffled_modularity$health_modularity) 
```

    ## [1] 756

``` r
sum(ibd_modularity < shuffled_modularity$disease_modularity) 
```

    ## [1] 948

``` r
sum(diff_modularity > shuffled_diff_modularity)
```

    ## [1] 915

``` r
## density
health_density <- edge_density(health_net)
ibd_density <- edge_density(ibd_net)
shuffled_density <- purrr::rerun(1000, shuffle_density(kos_norm, state, refnet) %>% unlist()) %>% 
  do.call(rbind, .) %>% 
  as.data.frame()
diff_density <- health_density - ibd_density
shuffled_diff_density <- shuffled_density$health_density - shuffled_density$disease_density
sum(health_density > shuffled_density$health_density)
```

    ## [1] 888

``` r
sum(ibd_density < shuffled_density$disease_density)
```

    ## [1] 740

``` r
sum(diff_density > shuffled_diff_density)
```

    ## [1] 873

用移除网络中的节点模拟外界环境对代谢网络的扰动，分别按三种顺序移除节点：介数中心度、度和随机。通过比较发现，healthy-network比IBD-network更稳定，无论用何种方法移除节点，需要注意的是，随着移除节点数的增多，按介数中心度（约50%）和随机移除节点（约70%），IBD-network比healthy-network更稳健。

``` r
nc_health_degree <- nc_attack(as.undirected(health_net),method = "degree")
nc_ibd_degree <- nc_attack(as.undirected(ibd_net),method = "degree")
nc_health_bc <- nc_attack(as.undirected(health_net), method = "betweenness")
nc_ibd_bc <- nc_attack(as.undirected(ibd_net), method = "betweenness")
nc_health_random <- nc_attack(as.undirected(health_net), method = "random")
nc_ibd_random <- nc_attack(as.undirected(ibd_net), method = "random")

nc_tibble <- bind_rows(
  make_nc_tibble(unlist(nc_health_degree), "health", "degree"),
  make_nc_tibble(unlist(nc_ibd_degree), "IBD", "degree"),
  make_nc_tibble(unlist(nc_health_bc), "health", "betweenness"),
  make_nc_tibble(unlist(nc_ibd_bc), "IBD", "betweenness")
  , make_nc_tibble(unlist(nc_health_random), "health", "random"),
   make_nc_tibble(unlist(nc_ibd_random), "IBD", "random")
)
p_nc <- ggplot(nc_tibble, aes(x, y, color = state)) +
  geom_line() +
  theme_bw() + 
  ggsci::scale_color_npg() +
  theme(legend.position = c(0.07, 0.9), legend.title = element_blank()) +
  labs(x = "Percentage of removed nodes", y = "Natural connectivity") +
  facet_wrap(~method, scales = "free_y")

p_nc
```

![](README_files/figure-markdown_github/natural-connectivity-1.png)

``` r
# ggsave(p_nc, filename = "output/p_nc.pdf", width = 7, height = 4)
```

其他
----

### modularity的p值

除了随机打乱样本label之外，我们还可以随机打乱网络的边（保持每个节点的出入度）1000次，然后裁剪到1000对网络作为零分布

``` r
shuffled_edge_modularity <- map(1:100,
  ~ shuffle_edge_modularity2(diff_ssn, health_kos, ibd_kos) %>% unlist()) %>%
  do.call(rbind, .) %>%
  as.data.frame()
shuffled_edge_modularity <- map(1:100,
  ~ shuffle_edge_modularity(health_net, ibd_net) %>% unlist()) %>%
  do.call(rbind, .) %>%
  as.data.frame()
shuffled_diff_edge_modularity <- shuffled_edge_modularity$health_modularity -
  shuffled_edge_modularity$disease_modularity
sum(diff_modularity > shuffled_diff_edge_modularity)
sum(health_modularity > shuffled_edge_modularity$health_modularity)
sum(ibd_modularity > shuffled_edge_modularity$disease_modularity)
```

### 其他拓扑属性的关联和比较

如度分布，seed set，参考[RevEcoR](https://github.com/yiluheihei/RevEcoR)

``` r
health_net_dd <- degree_distribution(health_net, mode = "all")
ibd_net_dd <- degree_distribution(ibd_net, cumulative = FALSE)
n_health_dd <- length(health_net_dd)
n_ibd_dd <- length(ibd_net_dd)
dd_tibble <- tibble(
  x = c(0:(n_health_dd - 1), 0:(n_ibd_dd - 1)),
  y = c(health_net_dd, ibd_net_dd),
  state = rep(c("health", "IBD"), c(n_health_dd, n_ibd_dd))
)
```

### 另外我们可以对拓扑属性（模块度和密度等）做稀疏分析

比如我们发现IBD模块度较小，我们可以不断增大随机抽取样本数，计算这些拓扑属性，看它们是否趋于稳定（收敛），因为我们的例子中样本数较小（127:IBD 73,健康54），增大样本数目的稀疏曲线不够好。可以不断增加随机抽取序列数进行稀疏分析（因为我们直接利用别人的功能注释数据进行分析，就没有进行这一步尝试）。

``` r
rarefed_net <- raref_sample(kos_norm, state, refnet)
# 计算时间超长，不知道为何
rarefed_modu <- rerun(100,  map(rarefed_net, ~ map(.x, ~ as.undirected(.x) %>% 
    cluster_leading_eigen(options = list(maxiter = 1000000000, ncv = 8)) %>% 
    modularity())
))
rarefed_health_modu <- map(rarefed_modu, ~ unlist(.x[[1]])) %>% do.call(cbind, .) %>% rowMeans()
rarefed_ibd_modu <- map(rarefed_modu, ~ unlist(.x[[2]])) %>% do.call(cbind, .) %>% rowMeans()
rarefed_tibble <- tibble(
  x = c(1:length(rarefed_health_modu), 1:length(rarefed_ibd_modu)),
  modularity = c(rarefed_health_modu, rarefed_ibd_modu),
  state = rep(c("healthy", "IBD"), c(length(rarefed_health_omdu), length(rarefed_ibd_omdu)))
)
ggplot(rarefed_tibble, aes(x, modularity, color = state)) + geom_line()
# raref_n <- lengths(rarefed_modularity)
# raref_vcount <- rerun(10, map_depth(rarefed_net, 2, vcount))
# raref_modularity_tibble <- tibble(
#   x = c(1:raref_n[1], 1:raref_n[2]),
#   modularity = unlist(rarefed_modularity),
#   state = rep(c("healthy", "IBD"), raref_n)
# )
# ggplot(raref_modularity_tibble, aes(x, modularity, color = state)) + 
#   geom_line()
```

sessionInfo
-----------

``` r
sessionInfo()
```

    ## R version 3.6.0 (2019-04-26)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Ubuntu 16.04.6 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/openblas-base/libblas.so.3
    ## LAPACK: /usr/lib/libopenblasp-r0.2.18.so
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] parallel  stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] patchwork_0.0.1               furrr_0.1.0                  
    ##  [3] future_1.14.0                 mmnet_1.13.0                 
    ##  [5] biom_0.4.0                    igraph_1.2.4.1               
    ##  [7] clusterProfiler_3.12.0        readxl_1.3.1                 
    ##  [9] forcats_0.4.0                 stringr_1.4.0                
    ## [11] purrr_0.3.2                   readr_1.3.1                  
    ## [13] tidyr_0.8.3                   tibble_2.1.3                 
    ## [15] ggplot2_3.2.1                 tidyverse_1.2.1              
    ## [17] curatedMetagenomicData_1.14.1 ExperimentHub_1.10.0         
    ## [19] dplyr_0.8.3                   Biobase_2.44.0               
    ## [21] AnnotationHub_2.16.0          BiocFileCache_1.8.0          
    ## [23] dbplyr_1.4.2                  BiocGenerics_0.30.0          
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] backports_1.1.4               fastmatch_1.1-0              
    ##   [3] plyr_1.8.4                    lazyeval_0.2.2               
    ##   [5] splines_3.6.0                 BiocParallel_1.18.0          
    ##   [7] listenv_0.7.0                 urltools_1.7.3               
    ##   [9] digest_0.6.20                 htmltools_0.3.6              
    ##  [11] GOSemSim_2.10.0               viridis_0.5.1                
    ##  [13] GO.db_3.8.2                   magrittr_1.5                 
    ##  [15] memoise_1.1.0                 globals_0.12.4               
    ##  [17] modelr_0.1.4                  enrichplot_1.4.0             
    ##  [19] prettyunits_1.0.2             colorspace_1.4-1             
    ##  [21] blob_1.2.0                    rvest_0.3.4                  
    ##  [23] rappdirs_0.3.1                ggrepel_0.8.1                
    ##  [25] haven_2.1.1                   xfun_0.8                     
    ##  [27] crayon_1.3.4                  RCurl_1.95-4.12              
    ##  [29] jsonlite_1.6                  zeallot_0.1.0                
    ##  [31] glue_1.3.1                    polyclip_1.10-0              
    ##  [33] gtable_0.3.0                  UpSetR_1.4.0                 
    ##  [35] Rhdf5lib_1.6.0                scales_1.0.0                 
    ##  [37] DOSE_3.10.2                   DBI_1.0.0                    
    ##  [39] Rcpp_1.0.2                    viridisLite_0.3.0            
    ##  [41] xtable_1.8-4                  progress_1.2.2               
    ##  [43] gridGraphics_0.4-1            bit_1.1-14                   
    ##  [45] europepmc_0.3                 stats4_3.6.0                 
    ##  [47] httr_1.4.0                    fgsea_1.10.0                 
    ##  [49] RColorBrewer_1.1-2            modeltools_0.2-22            
    ##  [51] pkgconfig_2.0.2               XML_3.98-1.20                
    ##  [53] flexmix_2.3-15                farver_1.1.0                 
    ##  [55] nnet_7.3-12                   RJSONIO_1.3-1.2              
    ##  [57] ggplotify_0.0.4               tidyselect_0.2.5             
    ##  [59] labeling_0.3                  rlang_0.4.0                  
    ##  [61] reshape2_1.4.3                later_0.8.0                  
    ##  [63] AnnotationDbi_1.46.0          munsell_0.5.0                
    ##  [65] cellranger_1.1.0              tools_3.6.0                  
    ##  [67] cli_1.1.0                     generics_0.0.2               
    ##  [69] RSQLite_2.1.2                 broom_0.5.2                  
    ##  [71] ggridges_0.5.1                evaluate_0.14                
    ##  [73] yaml_2.2.0                    knitr_1.23                   
    ##  [75] bit64_0.9-7                   ggraph_1.0.2                 
    ##  [77] nlme_3.1-139                  mime_0.7                     
    ##  [79] DO.db_2.9                     xml2_1.2.0                   
    ##  [81] compiler_3.6.0                rstudioapi_0.10              
    ##  [83] curl_4.0                      interactiveDisplayBase_1.22.0
    ##  [85] tweenr_1.0.1                  stringi_1.4.3                
    ##  [87] lattice_0.20-38               Matrix_1.2-17                
    ##  [89] ggsci_2.9                     vctrs_0.2.0                  
    ##  [91] pillar_1.4.2                  BiocManager_1.30.4           
    ##  [93] triebeard_0.3.0               data.table_1.12.2            
    ##  [95] cowplot_1.0.0                 bitops_1.0-6                 
    ##  [97] httpuv_1.5.1                  qvalue_2.16.0                
    ##  [99] R6_2.4.0                      promises_1.0.1               
    ## [101] gridExtra_2.3                 IRanges_2.18.1               
    ## [103] codetools_0.2-16              MASS_7.3-51.1                
    ## [105] assertthat_0.2.1              rhdf5_2.28.0                 
    ## [107] withr_2.1.2                   S4Vectors_0.22.0             
    ## [109] hms_0.5.0                     grid_3.6.0                   
    ## [111] rmarkdown_1.14                rvcheck_0.1.3                
    ## [113] ggforce_0.2.2                 shiny_1.3.2                  
    ## [115] lubridate_1.7.4

[1] H B, Almeida M, Juncker A S, et al. Identification and assembly of genomes and genetic elements in complex metagenomic samples without using reference genomes. Nature Biotechnology, 2014, 32(8): 822–828.
