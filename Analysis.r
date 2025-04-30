#### Cox分析 ####

end_date <- as.Date("2022-12-31")
data_osteoporosis$ost_date <- as.Date(data_osteoporosis$ost_date)
data_osteoporosis$recruitment_date <- as.Date(data_osteoporosis$recruitment_date)
data_osteoporosis$death_date <- as.Date(data_osteoporosis$death_date)
data_osteoporosis <- data_osteoporosis %>%
  mutate(survival_time = case_when(
    osteoporosis == 1 ~ as.numeric(ost_date - recruitment_date), # 患病
    osteoporosis == 0 & !is.na(death_date) ~ as.numeric(death_date - recruitment_date), # 未患病但死亡
    osteoporosis == 0 & is.na(death_date) ~ as.numeric(end_date - recruitment_date) # 未患病且未死亡
  ))
data_osteoporosis$sex <- as.factor(data_osteoporosis$sex)
data_osteoporosis$race <- as.factor(data_osteoporosis$race)
data_osteoporosis$income <- as.factor(data_osteoporosis$income)
data_osteoporosis$smoking <- as.factor(data_osteoporosis$smoking)
data_osteoporosis$drinking <- as.factor(data_osteoporosis$drinking)
cox_model <- coxph(Surv(survival_time, osteoporosis) ~ 
                     pm2.5_10 +  
                     age + sex + race + income + bmi + tdi + smoking + drinking + humidity + temperature,
                     data = data_osteoporosis)
cox_model
summary(cox_model)

#### 4.RCS ####
ps <- c("foreign","ggplot2","rms","survival","Hmisc","splines")
for(i in ps){library(i, character.only = T)}; rm(i)

dd <- datadist(data_osteoporosis)
options(datadist='dd')

aic_values <- sapply(3:5, function(k) {
  model <- cph(Surv(survival_time, osteoporosis) ~ rcs(aps, k) +  
                 age + sex + race + income + bmi + tdi + smoking + drinking + humidity + temperature, 
               data = data_osteoporosis)
  AIC(model)
})
aic_values
names(aic_values) <- c("3", "4", "5")
best_k <- as.numeric(names(which.min(aic_values)))

fit<-cph(Surv(survival_time,osteoporosis) ~ rcs(no2,best_k) +  
           age + sex + race + income + bmi + tdi + smoking + drinking + humidity + temperature, 
         data=data_osteoporosis)
fit
an <- anova(fit)
an

HR <- Predict(fit, no2, fun=exp, ref.zero = TRUE)

p1 <- ggplot(anova=an, pval=T) +
  # 画曲线
  geom_line(data=HR, aes(no2, yhat), linetype=1, linewidth=1, alpha=0.9, colour="red") +
  # 画置信区间
  geom_ribbon(data=HR, aes(no2, ymin=lower, ymax=upper), alpha=0.3, fill="red") +
  # x轴任意刻度：增加一条竖线
  geom_hline(yintercept=1, linetype=2, linewidth=1) +  # 替换 size 为 linewidth
  # 去除背景
  theme_classic() +
  # 增加标签
  labs(title="RCS", x="no2", y="HR (95%CI)")
p1

##### mediation analysis  #####
library(foreign)
library(plyr)
library(survival)
library(survminer)
boot.med <- function(data_olink_3, m_value_name){
  data_olink_3 <- data_olink_3[sample(1:nrow(data_olink_3), replace = T), ]
  lm_formula <- as.formula(paste(m_value_name, "~ pm2.5 + age + sex + race + income + bmi + tdi + smoking + drinking"))
  alpha.temp <- coefficients(lm(lm_formula, data_olink_3))[2]  
  cox_formula1 <- as.formula(paste("Surv(survival_time, osteoporosis == 1) ~", m_value_name, "+ age + sex + race + income + bmi + tdi + smoking + drinking"))
  beta.temp <- coefficients(coxph(cox_formula1, data_olink_3, method = "breslow"))[1]
  cox_formula2 <- as.formula(paste("Surv(survival_time, osteoporosis == 1) ~ pm2.5 +", m_value_name, "+ age + sex + race + income + bmi + tdi + smoking + drinking"))
  c_prime.temp <- coefficients(coxph(cox_formula2, data_olink_3, method = "breslow"))[1]
  cox_formula3 <- as.formula(paste("Surv(survival_time, osteoporosis == 1) ~ pm2.5 + age + sex + race + income + bmi + tdi + smoking + drinking"))
  c.temp <- coefficients(coxph(cox_formula3, data_olink_3, method = "breslow"))[1]
  IE1.l <- alpha.temp * beta.temp
  IE2.l <- c.temp - c_prime.temp
  DE.l <- c_prime.temp
  TOT.l <- c.temp
  results <- c(IE1.l, IE2.l, DE.l, TOT.l)
  return(results)
}
m.values <- variable_data$Protein[1:30]
result.m <- matrix(nrow = length(m.values), ncol = 13)
result.m <- data.frame(result.m)
rownames(result.m) <- m.values
for (m.value in m.values) {
  results <- boot.med(data_olink_3, m.value)
  G <- 1000
  IE1.cox  <- matrix(rep(0,G),G,1)
  IE2.cox  <- matrix(rep(0,G),G,1)  
  DE.cox  <- matrix(rep(0,G),G,1) 
  med.boot.cox  <- lapply(1:G, FUN = function(i) boot.med(data_olink_3, m.value))
  dim(med.boot.cox)
  IE1.cox<- unlist(lapply(med.boot.cox, '[[', 1)) 
  IE2.cox<- unlist(lapply(med.boot.cox, '[[', 2))
  DE.cox <- unlist(lapply(med.boot.cox, '[[', 3)) 
  TOT.cox <- unlist(lapply(med.boot.cox, '[[', 4)) 
  dir.cox <-mean(DE.cox) # IE c-c' method
  ci.dir.cox <- as.table(quantile (DE.cox, c(0.025, 0.975))) #95% CI's c-c'
  tot.cox <-mean(TOT.cox) # IE c-c' method
  ci.tot.cox <- as.table(quantile (TOT.cox, c(0.025, 0.975)))
  ind.ab.cox <- mean(IE1.cox) # IE ab method
  ind.ccp.cox <-mean(IE2.cox) # IE c-c' method
  ci.ab.cox <- as.table(quantile (IE1.cox, c(0.025, 0.975))) #95% CI's ab
  ci.ccp.cox <- as.table(quantile (IE2.cox, c(0.025, 0.975))) #95% CI's c-c'
  PM_1 <- ind.ab.cox / tot.cox  
  result.m[which(rownames(result.m) == m.value), ] <- 
    c(dir.cox, ci.dir.cox[1],ci.dir.cox[2], 
      tot.cox,ci.tot.cox[1],ci.tot.cox[2],
      ind.ab.cox, ind.ccp.cox,
      ci.ab.cox[1],ci.ab.cox[2], ci.ccp.cox[1],ci.ccp.cox[2],PM_1)
}
colnames(result.m) <- c('dir.cox', 'ci.dir.cox_2.5', 'ci.dir.cox_97.5',
                        'tot.cox', 'ci.tot.cox_2.5', 'ci.tot.cox_97.5', 'ind.ab.cox', 'ind.ccp.cox',
                        'ci.ab.cox_2.5', 'ci.ab.cox_97.5', 'ci.ccp.cox_2.5', 'ci.ccp.cox_97.5','PM_1')
