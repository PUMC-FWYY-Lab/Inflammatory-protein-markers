
#---------------------  导入芬兰数据作为发现队列循环每个蛋白--------------------------------------

library(data.table)
AF_duplicate_total <- readRDS("~/房颤与炎症性疾病/Data/AF_duplicate.RDS")
# 初始化一个空的列表来存储所有的exposure_to_AF数据框
exposure_to_AF_list <- list()

# 获取目录下所有的.RDS文件
rds_files <- list.files(path = "~/房颤与炎症性疾病/Data/59个cis-pQTLs", pattern = "_exp_dat.RDS", full.names = TRUE)

# 获取总基因数
total_genes <- length(rds_files)

# 遍历每个.RDS文件
for (i in 1:length(rds_files)) {
  # 读取RDS文件
  tryCatch({
    exposure_dat <- readRDS(rds_files[i])
    
    # 提取基因名
    gene_name <- gsub("_exp_dat.RDS", "", basename(rds_files[i]))
    
    # 打印进度
    cat("Processing gene", i, "out of", total_genes, "\n")
    
    # 生成harmonise_data
    AF_duplicate_subset <- AF_duplicate_total[AF_duplicate_total$SNP %in% exposure_dat$SNP, ]
    
    # 检查AF_duplicate_subset$SNP是否为空
    if (length(AF_duplicate_subset$SNP) == 0) {
      # 如果为空，创建一个空的数据框
      exposure_to_AF <- data.frame(matrix(NA, ncol = 13, nrow = 1))
    } else {
      AF_duplicate_harmonise_data <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = AF_duplicate_subset)
      
      # 进行MR分析
      data <- mr_heterogeneity(AF_duplicate_harmonise_data)
      
      # 根据条件选择MR方法
      if (any(data$method == "Inverse variance weighted" & data$Q_pval < 0.05)) {
        exposure_to_AF <- mr(AF_duplicate_harmonise_data, method_list = c('mr_ivw_mre'))
      } else {
        exposure_to_AF <- mr(AF_duplicate_harmonise_data)
      }
      
      # 生成odds ratios
      exposure_to_AF <- generate_odds_ratios(exposure_to_AF)
      
      # 去掉前两列
      exposure_to_AF <- exposure_to_AF[, -c(1, 2)]
    }
    
    # 添加基因名和其他信息
    exposure_to_AF$exposure <- gene_name
    exposure_to_AF$outcome <- "AF"
    
    if (exists("AF_duplicate_subset")) {
      exposure_to_AF$exposure_SNP_number <- length(exposure_dat$SNP)
      exposure_to_AF$exposure_SNP <- list(exposure_dat$SNP)
      exposure_to_AF$outcome_SNP_number <- length(AF_duplicate_subset$SNP)
      exposure_to_AF$outcome_SNP <- list(AF_duplicate_subset$SNP)
    } else {
      exposure_to_AF$exposure_SNP_number <- NA
      exposure_to_AF$exposure_SNP <- NA
      exposure_to_AF$outcome_SNP_number <- NA
      exposure_to_AF$outcome_SNP <- NA
    }
    
    # 将生成的exposure_to_AF添加到列表中
    exposure_to_AF_list <- c(exposure_to_AF_list, list(exposure_to_AF))
  }, error = function(e) {
    cat("Error occurred for gene:", gene_name, "\n", "Error message:", conditionMessage(e), "\n")
    next
  })
}

# 合并所有的exposure_to_AF数据框为最终的表格
final_table <- rbindlist(exposure_to_AF_list, fill = TRUE)

final_table_subset_0.05<-final_table[final_table$pval<0.05,]
final_table_subset_0.05_59<-final_table[final_table$pval<0.05/59,]

saveRDS(final_table, file = "~/房颤与炎症性疾病/Data/芬兰发现队列循环蛋白最终MR结果.RDS")
saveRDS(final_table_subset_0.05, file = "~/房颤与炎症性疾病/Data/芬兰发现队列循环蛋白最终MR结果0.05.RDS")
saveRDS(final_table_subset_0.05_59, file = "~/房颤与炎症性疾病/Data/芬兰发现队列循环蛋白最终MR结果0.05_59.RDS")

#---------------------  导入GCST006414数据作为验证队列循环每个蛋白--------------------------------------

AF_dat_total <- readRDS("~/房颤与炎症性疾病/Data/AF_dat.RDS")
# 初始化一个空的列表来存储所有的exposure_to_AF数据框
exposure_to_AF_list <- list()

# 获取目录下所有的.RDS文件
rds_files <- list.files(path = "~/房颤与炎症性疾病/Data/59个cis-pQTLs", pattern = "_exp_dat.RDS", full.names = TRUE)

# 获取总基因数
total_genes <- length(rds_files)

# 遍历每个.RDS文件
for (i in 1:length(rds_files)) {
  # 读取RDS文件
  tryCatch({
    exposure_dat <- readRDS(rds_files[i])
    
    # 提取基因名
    gene_name <- gsub("_exp_dat.RDS", "", basename(rds_files[i]))
    
    # 打印进度
    cat("Processing gene", i, "out of", total_genes, "\n")
    
    # 生成harmonise_data
    AF_duplicate_subset <- AF_dat_total[AF_dat_total$SNP %in% exposure_dat$SNP, ]
    
    # 检查AF_duplicate_subset$SNP是否为空
    if (length(AF_duplicate_subset$SNP) == 0) {
      # 如果为空，创建一个空的数据框
      exposure_to_AF <- data.frame(matrix(NA, ncol = 13, nrow = 1))
    } else {
      AF_duplicate_harmonise_data <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = AF_duplicate_subset)
      
      # 进行MR分析
      data <- mr_heterogeneity(AF_duplicate_harmonise_data)
      
      # 根据条件选择MR方法
      if (any(data$method == "Inverse variance weighted" & data$Q_pval < 0.05)) {
        exposure_to_AF <- mr(AF_duplicate_harmonise_data, method_list = c('mr_ivw_mre'))
      } else {
        exposure_to_AF <- mr(AF_duplicate_harmonise_data)
      }
      
      # 生成odds ratios
      exposure_to_AF <- generate_odds_ratios(exposure_to_AF)
      
      # 去掉前两列
      exposure_to_AF <- exposure_to_AF[, -c(1, 2)]
    }
    
    # 添加基因名和其他信息
    exposure_to_AF$exposure <- gene_name
    exposure_to_AF$outcome <- "AF"
    
    if (exists("AF_duplicate_subset")) {
      exposure_to_AF$exposure_SNP_number <- length(exposure_dat$SNP)
      exposure_to_AF$exposure_SNP <- list(exposure_dat$SNP)
      exposure_to_AF$outcome_SNP_number <- length(AF_duplicate_subset$SNP)
      exposure_to_AF$outcome_SNP <- list(AF_duplicate_subset$SNP)
    } else {
      exposure_to_AF$exposure_SNP_number <- NA
      exposure_to_AF$exposure_SNP <- NA
      exposure_to_AF$outcome_SNP_number <- NA
      exposure_to_AF$outcome_SNP <- NA
    }
    
    # 将生成的exposure_to_AF添加到列表中
    exposure_to_AF_list <- c(exposure_to_AF_list, list(exposure_to_AF))
  }, error = function(e) {
    cat("Error occurred for gene:", gene_name, "\n", "Error message:", conditionMessage(e), "\n")
    next
  })
}

# 合并所有的exposure_to_AF数据框为最终的表格
final_table <- rbindlist(exposure_to_AF_list, fill = TRUE)

final_table_subset_0.05<-final_table[final_table$pval<0.05,]

saveRDS(final_table, file = "~/房颤与炎症性疾病/Data/GCST006414验证队列循环蛋白最终MR结果.RDS")
saveRDS(final_table_subset_0.05, file = "~/房颤与炎症性疾病/Data/GCST006414验证队列循环蛋白最终MR结果_0.05.RDS")


Discovery_queue_MR<-readRDS("~/房颤与炎症性疾病/Data/芬兰发现队列循环蛋白最终MR结果0.05.RDS")

Validation_queue_MR<-readRDS("~/房颤与炎症性疾病/Data/GCST006414验证队列循环蛋白最终MR结果_0.05.RDS")



















