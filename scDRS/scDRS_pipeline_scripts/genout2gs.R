library(data.table)
library(dplyr)
# 读取ref文件
ref <- fread("./NCBI37.3.gene.loc",header = F)
# 获取工作目录下所有以genes.out结尾的文件
file_names <- list.files(pattern = "*.genes.out")
# 创建一个空列表来存储数据
all_data <- list()
# 循环读取所有文件
for (file in file_names) {
  # 读取每个文件
  alt <- fread(file)
  # 使用left_join根据ref中的V1与alt中的GENE列进行匹配，并将V6列的基因名称加入到alt中
  # 然后使用mutate来调整列
  alt_transformed <- alt %>%
    left_join(ref, by = c("GENE" = "V1")) %>%
    mutate(GENE = V6) %>%
    select(-V6)
  # 清除NA并保存到列表中
  alt_transformed_clean <- na.omit(alt_transformed)
  # 选择GENE和P列
  da <- alt_transformed_clean %>%
    select(GENE, P)
  # 重命名P列为文件名称中的关键部分
  object_name <- gsub("^trait_","",gsub(".genes.out", "", file))
  da <- da %>%
    rename(!!object_name := P)
  # 将数据存储到all_data列表中
  all_data[[object_name]] <- da
}
# 初始化result_df为第一个文件的GENE列
result_df <- data.frame(GENE = all_data[[1]]$GENE, stringsAsFactors = FALSE)
# 合并所有数据
for (object_name in names(all_data)) {
  current_data <- all_data[[object_name]]
  # 使用GENE列进行左连接，确保所有基因都包含在结果中
  result_df <- merge(result_df, current_data, by = "GENE", all.x = TRUE)
  # 重命名合并列（前面已经做过了）
  # names(result_df)[names(result_df) == "P"] <- object_name
}
# 将Pvalue 缺失值替换为1 
result_df[is.na(result_df)] <- 1
# 输出最终结果
fwrite(as.data.table(result_df), file = "geneset.pval", sep = "\t", quote = FALSE)
# 查看结果
print(head(result_df))


