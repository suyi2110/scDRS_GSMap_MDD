#!/usr/bin/env Rscript
library(data.table)

# 1. 获取文件列表并排除汇总表
files <- list.files(pattern="\\.tsv$")
exclude_files <- c("summary_Nsample.tsv", "summary_Nsample_new.tsv")
files <- files[!files %in% exclude_files]

res_list <- vector("list", length(files))

for (i in seq_along(files)) {
  f <- files[i]
  out_name <- sub("\\.tsv$", ".snpp", f)
  
  # 增量跳过逻辑（如果想重跑 Vitiligo，可以临时删掉 Vitiligo.snpp）
  if (file.exists(out_name)) {
    message(sprintf("跳过已存在: %s", f))
    next 
  }

  message(sprintf("正在尝试处理: %s", f))
  
  # 使用 fill=TRUE 增加容错性，并读取首行
  dt_sample <- fread(f, nrows = 5, fill = TRUE)
  header <- names(dt_sample)
  
  # 清理列名：去除前后的空白字符并转大写
  header_clean <- trimws(toupper(header))
  
  # 匹配列索引
  rsid_idx <- which(header_clean %in% c("RSID", "SNP", "SNPID", "MARKER"))[1]
  p_idx    <- which(header_clean %in% c("P", "P-VALUE", "PVALUE", "P.VALUE"))[1]
  n_idx    <- which(header_clean %in% c("N", "NSAMPLE", "WEIGHT", "SAMPLE_SIZE"))[1]

  # 如果没找到列，打印诊断信息
  if (is.na(rsid_idx) || is.na(p_idx)) {
    message(sprintf("!!! 匹配失败: %s", f))
    message(sprintf("   检测到的列名为: [%s]", paste(header, collapse="|")))
    next
  }

  rsid_col <- header[rsid_idx]
  p_col    <- header[p_idx]
  
  # 2. 仅读入必要列，大幅减少内存和时间
  cols_to_read <- c(rsid_col, p_col)
  if (!is.na(n_idx)) cols_to_read <- c(cols_to_read, header[n_idx])
  
  dt <- fread(f, select = cols_to_read)
  
  # 3. 输出 .snpp (去掉 NA)
  out_snpp <- na.omit(dt[, .(get(rsid_col), get(p_col))])
  setnames(out_snpp, c("RSID", "P"))
  
  fwrite(out_snpp, file = out_name, sep = "\t", col.names = FALSE)

  # 4. 计算 N 众数
  mode_n <- NA
  if (!is.na(n_idx)) {
    n_col_real <- header[n_idx]
    mode_n <- dt[, .N, by = n_col_real][order(-N)][1, get(n_col_real)]
  }
  
  res_list[[i]] <- data.table(file = f, N_mode = mode_n)
}

# 汇总结果
summary_dt <- rbindlist(res_list, fill = TRUE)
if (nrow(summary_dt) > 0) {
    fwrite(summary_dt, file = "summary_Nsample_new.tsv", sep = "\t")
    message("--- 处理结束 ---")
} else {
    message("--- 没有新文件被处理 ---")
}
