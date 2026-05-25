#!/usr/bin/env Rscript
library(DropletUtils)
library(Seurat)

args <- commandArgs(trailingOnly = TRUE)
sce_file <- args[1]          # 输入 SCE 文件路径
raw_seurat_file <- args[2]   # 输入原始 Seurat 文件路径
out_filtered_seurat <- args[3]  # 输出过滤后的 Seurat 文件路径
out_plot <- args[4]           # 输出 knee plot图片路径

# 读取 SCE 对象
sce <- readRDS(sce_file)

# 计算总 UMI
my.counts <- counts(sce)
total_umi <- colSums(my.counts)
total_umi_sorted <- sort(total_umi, decreasing = TRUE)

# top10%阈值
n_top10 <- ceiling(length(total_umi_sorted) * 0.10)
umi_threshold_top10 <- total_umi_sorted[n_top10]

# barcodeRanks
set.seed(0)
br.out <- barcodeRanks(my.counts, fit.bounds = c(umi_threshold_top10, max(total_umi)))

# emptyDrops 检测空载液滴
set.seed(100)
e.out <- emptyDrops(my.counts, retain = metadata(br.out)$knee, lower = umi_threshold_top10)

# 提取 knee 与 top10%之间的细胞
knee_val <- metadata(br.out)$knee
inflection_val <- metadata(br.out)$inflection
top10_val <- umi_threshold_top10
umi_per_bc <- colSums(my.counts)
mid_bc <- names(umi_per_bc)[umi_per_bc < knee_val & umi_per_bc >= top10_val]
mid_res <- e.out[mid_bc, ]
mid_res$UMI <- umi_per_bc[mid_bc]
mid_res <- mid_res[order(mid_res$UMI, decreasing = TRUE), ]

# 找到显著比例小于50%的窗口
window_size <- 20
threshold <- 0.5
last_umi <- NA
for (i in seq(window_size, nrow(mid_res))) {
  window_res <- mid_res[(i - window_size + 1):i, ]
  sig_ratio <- mean(window_res$FDR <= 0.01, na.rm = TRUE)
  if (!is.na(sig_ratio) && !is.nan(sig_ratio) && sig_ratio < threshold) {
    last_umi <- window_res$UMI[window_size]
    break
  }
}
if (is.na(last_umi)) {
  last_umi <- top10_val
}

# re knee
knee_val_raw <- knee_val
repeat {
  if (abs(last_umi - knee_val) > 5) {
    break
  }
  set.seed(100)
  br.out2 <- barcodeRanks(my.counts, fit.bounds = c(knee_val, max(total_umi)))
  knee_val2 <- metadata(br.out2)$knee
  e.out2 <- emptyDrops(my.counts, retain = metadata(br.out2)$knee, lower = knee_val)
  
  mid_bc2 <- names(umi_per_bc)[umi_per_bc < knee_val2 & umi_per_bc > knee_val-1]
  mid_res2 <- e.out2[mid_bc2, ]
  mid_res2$UMI <- umi_per_bc[mid_bc2]
  mid_res2 <- mid_res2[order(mid_res2$UMI, decreasing = TRUE), ]
  
  last_umi_new <- NA
  for (i in seq(window_size, nrow(mid_res2))) {
    window_res <- mid_res2[(i - window_size + 1):i, ]
    sig_ratio <- mean(window_res$FDR <= 0.01, na.rm = TRUE)
    if (!is.na(sig_ratio) && !is.nan(sig_ratio) && sig_ratio < threshold) {
      last_umi_new <- window_res$UMI[window_size]
      break
    }
  }
  if (is.na(last_umi_new)) {
    last_umi_new <- knee_val
  }
  knee_val <- knee_val2
  last_umi <- last_umi_new
}

# 标记最终保留的细胞
e.out$IsCell_byUMI <- e.out$Total >= last_umi
is.cell <- rownames(e.out)[e.out$IsCell_byUMI == TRUE]

# 绘制 knee plot
png(out_plot, width = 4000, height = 3200, res = 600)
col.vec <- ifelse(e.out$IsCell_byUMI, "red", "gray")
plot(br.out$rank, br.out$total, log="xy",
     xlab="Rank", ylab="Total", main="Knee plot",
     col=col.vec, pch=16, cex=0.3)
abline(h = knee_val_raw, col = "dodgerblue", lty = 2)
if (exists("knee_val2") && !is.null(knee_val2)) {
     abline(h = knee_val2, col = "blue4", lty = 2)
}
abline(h = metadata(br.out)$inflection, col = "forestgreen", lty = 2)
abline(h = umi_threshold_top10, col = "purple", lty = 2)
abline(h = last_umi, col = "red3", lty = 2)
text(x=max(br.out$rank)*0.003, y=knee_val_raw, 
     labels=paste0("knee=", round(knee_val_raw,2)), 
     col="dodgerblue", pos=3, cex=1.2)
if (exists("knee_val2") && !is.null(knee_val2)) {
    text(x = max(br.out$rank) * 0.003, y = knee_val2, 
         labels = paste0("knee last=", round(knee_val2, 2)), 
         col = "blue4", pos = 3, cex = 1.2)
}
text(x=max(br.out$rank)*0.003, y=inflection_val, 
     labels=paste0("inflection=", round(inflection_val,2)), 
     col="forestgreen", pos=3, cex=1.2)
text(x=max(br.out$rank)*0.003, y=umi_threshold_top10, 
     labels=paste0("top10%=", round(umi_threshold_top10,2)), 
     col="purple", pos=3, cex=1.2)
text(x=max(br.out$rank)*0.003, y=last_umi, 
     labels=paste0("last_umi=", round(last_umi,2)), 
     col="red3", pos=3, cex=1.2)
legend_labels <- c("knee", "inflection", "top10%", "last_umi")
legend_cols   <- c("dodgerblue", "forestgreen", "purple", "red3")
if (exists("knee_val2") && !is.null(knee_val2)) {
    legend_labels <- c(legend_labels, "knee last")
    legend_cols   <- c(legend_cols, "blue4")
}
legend("bottomleft", lty = 2, col = legend_cols, legend = legend_labels)
dev.off()

# 读取原始 Seurat 对象并筛选细胞
raw_seurat_obj <- readRDS(raw_seurat_file)
seurat_obj_filtered <- subset(raw_seurat_obj, cells = is.cell)

# 保存过滤后的 Seurat 对象
saveRDS(seurat_obj_filtered, out_filtered_seurat)
  