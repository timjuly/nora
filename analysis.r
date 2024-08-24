# 设置工作路径，所有的输出结果会放到这个文件夹里，没有请新建
setwd(system("pwd", intern = T)) # or setwd('~/Documents/Nora')

# 安装依赖
install.packages("meta", repos = "https://cran.rstudio.com")
install.packages("metafor", repos = "https://cran.rstudio.com")

library(meta)
library(metafor)

# 准备数据
data <- data.frame(
  Study = c("Schneuer et al. (2014)", "Flood-Nichols et al. (2015)", 
            "Yue et al. (2021)", "Achkar et al. (2015)", "Zhao et al. (2017)", 
            "Arisoy et al. (2016)", "Benachi et al. (2020)", 
            "Benachi et al. (2020)", "Zeng et al. (2020)"),
  Stage = c("Early", "Early", "Mid", "Mid", "Mid", "Late", "Early", "Late", "Mid"),
  SampleSize = c(5109, 235, 7976, 2144, 11151, 257, 402, 302, 370),
  OR = c(0.98, 1.36, 1.48, 2.23, 3.16, 3.27, 1.75, 2.32, 1.28),
  LowerCI = c(0.73, 0.48, 1.16, 1.29, 1.77, 1.32, 0.99, 1.25, 1.15),
  UpperCI = c(1.33, 3.88, 1.90, 3.83, 5.65, 8.10, 3.33, 4.35, 1.79)
)

# 进行 Meta 分析
meta_analysis <- metagen(
  TE = log(data$OR), # 对数效应量
  seTE = (log(data$UpperCI) - log(data$LowerCI)) / (2 * 1.96), # 效应量的标准误
  studlab = data$Study,
  data = data,
  sm = "OR", 
  method.tau = "REML" # 使用REML方法估计异质性
)

print(summary(meta_analysis))
# 保存结果
sink("元分析结果.txt")
print(summary(meta_analysis))
sink() # 恢复控制台输出

# 森林图
pdf("森林图.pdf", width = 11.69, height = 8.27) # A4 大小
forest(meta_analysis,
  print.subgroup.labels = TRUE,
  print.tau2 = TRUE,
  print.I2 = TRUE,
  print.Q = TRUE
)
dev.off()

# 查看异质性统计量
meta_analysis$I2
meta_analysis$Q
meta_analysis$pval.Q
# 保存结果
sink("异质性统计量.txt")
print("$I2")
meta_analysis$I2
print("$Q")
meta_analysis$Q
print("$pval.Q")
meta_analysis$pval.Q
sink() # 恢复控制台输出

# 进行逐项排除敏感性分析
sensitivity_result <- metainf(metagen(TE = meta_analysis$TE, seTE = meta_analysis$seTE, studlab = data$Study, data = data, comb.random = TRUE))
print(sensitivity_result)
# 保存结果
sink("敏感性分析.txt")
print(sensitivity_result)
sink() # 恢复控制台输出

# 进行一比一替换分析
replace_one_sensitivity <- metainf(meta_analysis, method.incr = "replace")
print(replace_one_sensitivity)
# 保存结果
sink("逐项替换灵敏度分析.txt")
print(replace_one_sensitivity)
sink() # 恢复控制台输出

# 绘制漏斗图
pdf("漏斗图.pdf", width = 11.69, height = 8.27) # A4 大小
funnel(meta::metagen(TE = meta_analysis$TE, seTE = meta_analysis$seTE), xlab="OR")
dev.off()

# 拟合随机效应模型
rma_model <- metafor::rma(yi = meta_analysis$TE, sei = meta_analysis$seTE, method = "REML")
# 查看模型结果
print(rma_model)
# 保存结果
sink("拟合随机效应模型.txt")
print(rma_model)
sink() # 恢复控制台输出

# 应用Trim and Fill方法
trimfill_result <- trimfill(rma_model)

# 绘制调整后的漏斗图
pdf("漏斗图_调整后.pdf", width = 11.69, height = 8.27) # A4 大小
funnel(trimfill_result)
dev.off()

# 拟合随机效应模型 DerSimonian-Laird (DL)
rma_model_DL <- metafor::rma(yi = meta_analysis$TE, sei = meta_analysis$seTE, method = "DL")
# 查看模型结果
print(rma_model_DL)
# 保存结果
sink("拟合随机效应模型_DL.txt")
print(rma_model_DL)
sink() # 恢复控制台输出

# 应用Trim and Fill方法
trimfill_result_DL <- trimfill(rma_model_DL)

# 绘制调整后的漏斗图
pdf("漏斗图_调整后_DL.pdf", width = 11.69, height = 8.27) # A4 大小
funnel(trimfill_result_DL)
dev.off()

# 进行影响分析
influence_analysis <- influence(rma_model)
print(influence_analysis)
# 保存结果
sink("影响分析.txt")
print(influence_analysis)
sink() # 恢复控制台输出

# 影响图
pdf("影响图.pdf", width = 11.69, height = 8.27) # A4 大小
plot(influence_analysis)
dev.off()

# 进行Egger’s测试
egger_test <- metafor::regtest(rma_model, model = "lm")
# 输出Egger’s测试结果
print(egger_test)
# 保存结果
sink("Egger测试结果.txt")
print(egger_test)
sink() # 恢复控制台输出

# 进行Begg’s 测试
begg_test <- metafor::ranktest(rma_model)
# 输出Begg’s测试结果
print(begg_test)
sink("Begg测试结果.txt")
print(begg_test)
sink() # 恢复控制台输出

# 生成Galbraith图
pdf("Galbraith图.pdf", width = 11.69, height = 8.27) # A4 大小
galbraith(meta_analysis)
dev.off()

# 进行 subgroup meta 分析
subgroup_meta_analysis <- metagen(
  TE = log(data$OR),       # Log-transformed Odds Ratio
  seTE = (log(data$UpperCI) - log(data$LowerCI)) / (2 * 1.96),  # Standard Error
  studlab = data$Study,
  data = data,
  subgroup = data$Stage,  # Subgrouping by Stage
  sm = "OR",
  method.tau = "REML" # 使用REML方法估计异质性
)

print(subgroup_meta_analysis)
# 保存结果
sink("元分析结果_亚组.txt")
print(summary(subgroup_meta_analysis))
sink() # 恢复控制台输出

# 亚组森林图
pdf("亚组森林图.pdf", width = 11.69, height = 8.27) # A4 大小
forest(subgroup_meta_analysis,
  print.subgroup.labels = TRUE,
  print.tau2 = TRUE,
  print.I2 = TRUE,
  print.Q = TRUE
)
dev.off()