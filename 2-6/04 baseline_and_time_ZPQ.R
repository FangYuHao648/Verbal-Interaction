# --- 1. 加载包 ---
library(lme4)
library(lmerTest)  # 用于为lmer模型提供p值
library(dplyr)
library(readxl)  # 添加 readxl 包用于读取 Excel 文件

# 读取你的 Excel 数据
df <- read_excel("D:\\桌面\\Matlab Working path\\residual effect\\Re_ISC_end_filled_end.xlsx")

# --- 2. 确保数据格式正确 ---
df$subject <- factor(df$subject)
df$Feedback <- factor(df$Feedback)
df$Information <- factor(df$InformationType)
df$EyeContact <- factor(df$EyeContact)

df$Pre_Feedback <- factor(df$Pre_Feedback)
df$Pre_Information <- factor(df$Pre_InformationType)
df$Pre_EyeContact <- factor(df$Pre_EyeContact)

df$trial_number <- df$trial_number # 假设 trial_number 就是试验顺序
df$trial_number_z <- scale(df$trial_number) # Z-标准化

cat("数据准备完毕。\n")

# --- 3. 模型构建函数 ---
build_model <- function(formula, data) {
  return(lmer(formula, data = data, REML = FALSE))
}

# --- 4. 层次化模型构建 ---
model_baseline <- build_model(
  ISC ~ Feedback * Information * EyeContact + (1 | subject),
  df
)

# 模型1：加入时间效应 (trial_number_z)
model_time <- build_model(
  ISC ~ Feedback * Information * EyeContact + trial_number_z + (1 | subject),
  df
)

# 模型2：加入一个滞留效应 (Pre_Feedback)
model_pre_F <- build_model(
  ISC ~ Feedback * Information * EyeContact + trial_number_z + Pre_Feedback + (1 | subject),
  df
)

# 模型3：加入一个滞留效应 (Pre_InformationType)
model_pre_I<- build_model(
  ISC ~ Feedback * Information * EyeContact + trial_number_z + Pre_Information+ (1 | subject),
  df
)

# 模型4：加入一个滞留效应 (Pre_Eyecontact)
model_pre_E <- build_model(
  ISC ~ Feedback * Information * EyeContact + trial_number_z + Pre_EyeContact+ (1 | subject),
  df
)


# 模型5：(Pre_Feedback+Pre_Information)
model_pre_FI <- build_model(
  ISC ~ Feedback * Information * EyeContact + trial_number_z + Pre_Feedback + Pre_Information + (1 | subject),
  df
)

# 模型6：加入两个滞留效应的交互 (Pre_Feedback * Pre_Information)
model_pre_FI_int <- build_model(
  ISC ~ Feedback * Information * EyeContact + trial_number_z + Pre_Feedback * Pre_Information + (1 | subject),
  df
)

# 模型7：加入两个滞留效应的交互 (Pre_Feedback + Pre_Eyecontact)
model_pre_FE<- build_model(
  ISC ~ Feedback * Information * EyeContact + trial_number_z + Pre_Feedback + Pre_EyeContact + (1 | subject),
  df
)

# 模型8：加入两个滞留效应的交互 (Pre_Feedback *Pre_Eyecontact)
model_pre_FE_int <- build_model(
  ISC ~ Feedback * Information * EyeContact + trial_number_z + Pre_Feedback *Pre_EyeContact + (1 | subject),
  df
)

# 模型9：加入两个滞留效应的交互 (Pre_Information +Pre_Eyecontact)
model_pre_IE <- build_model(
  ISC ~ Feedback * Information * EyeContact + trial_number_z + Pre_Information+Pre_EyeContact + (1 | subject),
  df
)

# 模型10：加入两个滞留效应的交互 (Pre_Information*Pre_Eyecontact)
model_pre_IE_int <- build_model(
  ISC ~ Feedback * Information * EyeContact + trial_number_z + Pre_Information*Pre_EyeContact + (1 | subject),
  df
)
# 模型11：加入所有三个滞留效应的交互 (Pre_Feedback + Pre_Information + Pre_EyeContact)
model_pre_FIE <- build_model(
  ISC ~ Feedback * Information * EyeContact + trial_number_z + Pre_Feedback + Pre_Information + Pre_EyeContact + (1 | subject),
  df
)

# 模型12：加入所有三个滞留效应的交互 (Pre_Feedback * Pre_Information * Pre_EyeContact)
model_pre_FIE_int <- build_model(
  ISC ~ Feedback * Information * EyeContact + trial_number_z + Pre_Feedback * Pre_Information * Pre_EyeContact + (1 | subject),
  df
)

# 高级滞留效应模型
model_true_carryover <- build_model(
  ISC ~ Feedback * Information * EyeContact + 
    (Feedback:Pre_Feedback) + 
    (Information:Pre_Information) + 
    (EyeContact:Pre_EyeContact) + 
    (1 | subject),
  df
)


# 比较基线模型与其他模型
results <- list()

# 比较基线模型与其他模型
results[[1]] <- anova(model_baseline, model_time)
results[[2]] <- anova(model_baseline, model_pre_F)
results[[3]] <- anova(model_baseline, model_pre_I)
results[[4]] <- anova(model_baseline, model_pre_E)
results[[5]] <- anova(model_baseline, model_pre_FI)
results[[6]] <- anova(model_baseline, model_pre_FI_int)
results[[7]] <- anova(model_baseline, model_pre_FE)
results[[8]] <- anova(model_baseline, model_pre_FE_int)
results[[9]] <- anova(model_baseline, model_pre_IE)
results[[10]] <- anova(model_baseline, model_pre_IE_int)
results[[11]] <- anova(model_baseline, model_pre_FIE)
results[[12]] <- anova(model_baseline, model_pre_FIE_int)
results[[13]] <- anova(model_baseline, model_true_carryover)



# 合并所有显著比较结果
final_results <- bind_rows(results)

# 输出最终结果
cat("\n\n--- 显著增量模型比较结果 ---\n")
print(final_results)


