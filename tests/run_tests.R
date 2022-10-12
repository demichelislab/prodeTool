# library(devtools)
# install()
# library(prodeTool)
#
# # SAMPLE RUN ===================================================================
#
# N_GENES = 10000
# N_SAMPLES = 10
#
# gr <- data.frame(
#     n1 =  sample(paste0("YY", 1:(N_GENES-1)), N_GENES*100, replace=T),
#     n2 =  sample(paste0("YY", 1:N_GENES), N_GENES*100, replace=T)
# )
#
# mis <- paste0("YY", N_GENES)
#
# .getAdjMatr2 <- function(etab, gns){
#
#     gns2n <- 1:length(gns)
#     names(gns2n) <- gns
#
#     etab <- as.matrix(etab[which(etab[,1]%in%gns & etab[,2]%in%gns),])
#
#     mm <- matrix(0, length(gns), length(gns))
#     colnames(mm) <- rownames(mm) <- gns
#     mm[etab] <- 1
#     return(mm)
# }
#
# .getAdjMatr3 <- function(etab, gns){
#
#     gns2n <- 1:length(gns)
#     names(gns2n) <- gns
#
#     idxm <- cbind(
#         c(as.numeric(gns2n[etab[,1]]), as.numeric(gns2n[etab[,2]]), 1:length(gns)),
#         c(as.numeric(gns2n[etab[,2]]), as.numeric(gns2n[etab[,1]]), length(gns):1)
#     )
#
#     idxm <- idxm[complete.cases(idxm), ]
#
#     mm <- Matrix(0, length(gns), length(gns))
#     colnames(mm) <- rownames(mm) <- gns
#     mm[idxm] <- 1
#     return(mm)
# }
#
# etab <- gr
#
# mm1 <- .getAdjMatr3(etab, paste0("YY", 1:(N_GENES)))
# mm2 <- .getAdjMatr2(etab, paste0("YY", 1:(N_GENES)))
#
#
#
# adj <- .getAdjMatr(gr, paste0("YY", 1:(N_GENES-1)))
#
# dim(adj)
#
# adj2 <- Matrix(adj)
#
# ds <-
#     matrix(
#         rnorm(N_GENES*N_SAMPLES),
#         nrow=N_GENES
#     )
#
# rownames(ds) <- paste0("YY", 1:N_GENES)
# colnames(ds) <- paste0("S", 1:ncol(ds))
#
# dm <- data.frame(
#     a = rep(c(1, 0), each=N_SAMPLES/2),
#     b = sample(c("a", "b"), N_SAMPLES, replace=T)
# )
# rownames(dm) <- paste0("S", 1:ncol(ds))
#
# covs <- "b"
# cond <- "a"
#
# prodeInput <- prepareProdeObject( # i can remove cores and filterCtrl
#     scores_matrix  = ds,
#     samples_info   = dm,
#     condition      = cond,
#     covariates     = covs,
#     edge_table     = gr,
#     cores          = 2,
#     filterCtrl     = T
# )
#
# proout <- runProde( # This has to be changed so that it takes object+formula+Padj.method+filterThreshold
#     scores_matrix  = ds,
#     columns_data   = dm,
#     condition      = cond,
#     covariates     = covs,
#     edge_table     = gr,
#     cores          = 2,
#     filterCtrl     = T
# )
#
#
# # Test this prode version on real data -----------------------------------------
#
# roo <- c("~/shares/CIBIO-Storage/CO/SPICE/thomas/proDe_v2/")
#
# source(paste0(roo, "./code/utils/_imports_.R"))
# source(paste0(roo, "./code/utils/_utils_depmap_functions.R"))
# source(paste0(roo, "./code/utils/_utils_prode_functions.R"))
# source(paste0(roo, "./code/utils/_utils_ppin_functions.R"))
#
# # Load depmap data -------------------------------------------------------------
#
# dep_data <- readRDS(paste0(roo,"./db/depmap_data_processed/depmap22Q1_complete_all_models_expression.rds"))
#
# dttmp <- fread(paste0(roo, "./db/ppin_data/edge_list_combined_400_genes.txt"), data.table=F)
#
# mm <- .getAdjMatr3(etab = dttmp, gns=colnames(dep_data[[2]]))
#
# # Select 9p21.3 loss models ----------------------------------------------------
#
# getLineage <- function(depmap_data, lineage){
#
#     if (!lineage%in%depmap_data[[1]]$lineage){
#         stop("Lineage not found")
#     }
#
#     to_keep <- which(depmap_data[["meta_data"]]$lineage==lineage)
#
#     depmap_data[[1]] <- depmap_data[[1]][to_keep,]
#     depmap_data[[2]] <- depmap_data[[2]][to_keep,]
#     depmap_data[[3]] <- depmap_data[[3]][to_keep,]
#
#     return(depmap_data)
# }
#
# dep_data.tmp <- dep_data#<- getLineage(dep_data, "urinary_tract")
#
# tmp_xp <- dep_data.tmp[[3]][,c("CDKN2A","CDKN2B","MTAP")]
#
# wt_mods  <- rownames(tmp_xp)[
#     tmp_xp$MTAP >= 2 &
#         tmp_xp$CDKN2A >= 2 &
#         tmp_xp$CDKN2B >= 2 #&
#     # tmp_xp$MTAP <= 6 &
#     # tmp_xp$CDKN2A <= 8 &
#     # tmp_xp$CDKN2B <= 4
# ] #rownames(tmp_cn)[apply(tmp_cn, 1, function(x) all(x==2)) & apply(tmp_sn, 1, function(x) all(x==0))]
#
# mut_mods <- rownames(tmp_xp)[
#     tmp_xp$MTAP <= 1 &
#         tmp_xp$CDKN2A <= 1 &
#         tmp_xp$CDKN2B <= 1
# ]#rownames(tmp_cn)[apply(tmp_cn, 1, function(x) all(x==0))]
#
# # Run ProDe --------------------------------------------------------------------
#
# design_matr <- dep_data.tmp[[1]][c(wt_mods, mut_mods),]
#
# design_matr$group <- c(
#     rep(0, length(wt_mods)),
#     rep(1, length(mut_mods))
# )
#
# score_matr <- t(dep_data[[2]][c(wt_mods, mut_mods),])
#
# covs <- c("lineage")
#
# if (length(table(design_matr$msi)) == 2){
#     covs <- c(covs, "msi")
# }
#
# if (length(table(design_matr$culture_type)) == 2){
#     covs <- c(covs, "culture_type")
# }
#
# # subset graph and beta to same gns ............................................
# # This part can be included in the object creation before running prode.
#
# all_gns <- unique(unlist(dttmp))
# common <- intersect(all_gns, rownames(score_matr))
# score_matr <- score_matr[common,]
# dttmp <- dttmp[which(dttmp$gene1%in%common & dttmp$gene2%in%common),]
# score_matr[which(is.na(score_matr))] <- 0
#
# # ..............................................................................
#
# out <- prodeTool::runProde(
#     scores_matrix  = score_matr,
#     columns_data   = design_matr,
#     condition      = "group",
#     covariates     = covs,
#     edge_table     = dttmp,
#     cores          = 1,
#     filterCtrl     = T
# )
#
# dd <- fread(
#     paste0(
#         "~/shares/CIBIO-Storage/CO/SPICE/thomas/proDe_v2/",
#         "data/06_validation_inhouse_isogenic_lines/prode_rra_9p21.txt"), data.table=F)
#
# rownames(dd) <- dd$gene
# head(dd)
# head(out)
#
# cc <- intersect(dd$gene, out$gene)
#
# {
#     plot(
#         -log10(out[cc,"fdr"]),
#         -log10(dd[cc,"fdr"])
#     )
#     abline(a = 0, b=1, col="red")
# }
#
# {
#     plot(
#         -log10(out[cc,"p.value"]),
#         -log10(dd[cc,"pval_rra"])
#     )
#     abline(a = 0, b=1, col="red")
# }
#
# plot(
#     x=dd[cc,"u"],
#     y=out[cc,"u"]
# )
#
# plot(
#     x=dd[cc,"score_rra"],
#     y=out[cc,"rra_score"]
# )
#
# plot(
#     x=dd[cc,"beta"],
#     y=out[cc,"Estimate"]
# )
#
#
#
# adj <- fread(  "~/shares/CIBIO-Storage/CO/SPICE/thomas/proDe_v2/data/06_validation_inhouse_isogenic_lines/adj_matrix_9p21.txt")
# adj <- as.matrix(adj)
# rownames(adj) <- colnames(adj)
#
# head(adj[1:10,1:10])
# diag(adj) <- 1
# head(adj[,1:10])
#
# adj <- adj[-which(rownames(adj) == "GFY"), -which(rownames(adj) == "GFY")]
# #adj <- adj[rownames(adj_mat), colnames(adj_mat)]
#
# all(rownames(adj) == rownames(adj_mat))
# all(adj == adj_mat)
#
# nn <- rowSums(adj[cc,cc])
# mm <- rowSums(adj_mat[cc,cc])
#
# plot(mm, nn)
#
#
# # Check out computed pbs method ------------------------------------------------
#
# all(rownames(adj) == dd$gene)
#
# bet1 <- dd[-which(dd$gene=="GFY"),]
# bet2 <- out[dd$gene[-which(dd$gene=="GFY")],]
#
# plot(
#     bet1$beta_scaled,
#     bet2$t.value
# )
# abline(a=0, b=1, col="red")
#
# plot(rank(dd$beta_scaled)/nrow(dd), dd$u)
# abline(a=0, b=1, col="red")
#
# plot(rank(bet1$beta_scaled)/nrow(bet1), rank(bet2$t.value)/nrow(bet2))
# plot(bet1$u, bet2$u)
# abline(a=0, b=1, col="red")
#
# oo <- apply(adj, 1, function(nn){
#
#     re <- sort(bet$u[which(nn!=0)])
#
#     #rnks <- sort(bet$pval[which(nn!=0)])
#
#     ps <- pbeta(re, 1:length(re), length(re) - 1:length(re) + 1)
#
#     min(ps, na.rm=T)
#
# })
#
# hist(oo)
#
# {
#     plot(oo, dd[rownames(adj),"score_rra"])
#     abline(a=0, b=1, col="red")
# }
#
#
# dim(adj)
# dim(out)
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
