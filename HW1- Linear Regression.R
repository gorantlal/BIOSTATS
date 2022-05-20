library(ggplot2) 

#Exploring explanatory variables as predictive models

#loading dataset
medcost = load("Desktop/Medical_Cost.RData")
df = na.omit(df) #discard data points with missing values
ggplot(df, aes(x=bmi, y=charges, color=smoker)) + geom_point()

#charges = B0 + B1*bmi
fit1 = lm(charges ~ bmi, data=df) #bmi as the sole predictor predictor
fit1_sum = summary(fit1)
fit1_confint = confint(fit1)
fit1_mse = mean(fit1_sum$residuals^2)
predict(fit1, data.frame(bmi=c(32)))

with(df, plot(bmi,charges))
abline(fit1)

p1 = ggplot(df, aes(x=bmi, y=charges, color=smoker)) + geom_point() +
  geom_abline(slope = coef(fit1)[[2]], intercept = coef(fit1)[[1]])
p1 + ggtitle("charge = b[o] + b[1]*bmi")

#The average medical cost is about $1192.94. Additionally, an increase of about 
#$393.87 medical charges is associated with an increase by about 1 unit of BMI.

#The 95% confidence interval for β0 is [-2072.97, 4458.8487] and the 95% 
#confidence interval for β1 is [289.41, 498.34]. Therefore, we can conclude that 
#in the absence of any BMI consideration, medical cost on average, fall somewhere 
#between $-2072.97 and $4458.85. Furthermore, for every 1 unit increase in BMI, 
#there will be an average increase in medical cost of between $289.41 and $498.34. 


#to make smoker into a dummy variable
df$smoker = factor(df$smoker, levels =c("no","yes"), ordered = TRUE)

#charges = B0 + B1*bmi + B2*smoker
fit2 = lm(formula = charges ~ bmi + smoker, data = df) #bmi,smoker as the predictor
fit2_sum = summary(fit2)
fit2_confint = confint(fit2)
fit2_mse = mean(fit2_sum$residuals^2)
predict(fit2, data.frame(bmi=c(32), smoker="yes"))

p2 = ggplot(df, aes(x=bmi, y=charges, color=smoker)) + geom_point() +
  geom_abline(slope = coef(fit2)[[2]], intercept = coef(fit2)[[1]] + coef(fit2)[[3]], colour = "orange") + 
  geom_abline(slope = coef(fit2)[[2]], intercept = coef(fit2)[[1]], colour = "purple")
p2 + ggtitle("charge = b[o] + b[1]*bmi + b[2]*smoker")

with(df, plot(bmi,charges))
abline(fit2)

#The dummy variable for smoking is 1 for smokers and 0 for non-smokers. The 
#average medical charges for a non-smoker is estimated to be $-3459.10 while the 
#average medical charge for a smoker is estimated to be an additional $23593.98 
#making the total $20,134.88. Additionally, an increase of about $388.02 is 
#associated with an increase by about 1 unit of BMI. 

#We can conclude that in the absence of any BMI consideration, medical cost on 
#average for non-smokers, fall somewhere between $-5417.4628 and $-1500.73. 
#The medical cost on average for smokers, fall somewhere between $17,234.53 and 
#$23,035.24. Furthermore, for every 1 unit increase in BMI, there will be an 
#average increase in medical cost of between $325.66 and $450.37.

#charges = B0 + B1*bmi + B2*smoker + B3*(bmi*smoker)
fit3 = lm(formula = charges ~ bmi*smoker, data=df) #bmi,smoker interaction
fit3_sum = summary(fit3)
fit3_confint = confint(fit3)
fit3_mse = mean(fit3_sum$residuals^2)
predict(fit3, data.frame(bmi=c(32), smoker="yes"))
predict(fit3, data.frame(bmi=c(28), smoker="yes"))

p3 = ggplot(df, aes(x=bmi, y=charges, color=smoker)) + geom_point() +
  geom_abline(slope = coef(fit3)[[2]] + coef(fit3)[[4]] , intercept = coef(fit3)[[1]] + coef(fit3)[[3]], colour = "orange") + 
  geom_abline(slope = coef(fit3)[[2]], intercept = coef(fit3)[[1]], colour = "purple")
p3 + ggtitle("charge = b[o] + b[1]*bmi + b[2]*smoker + b[3](bmi*smoker)")

with(df, plot(bmi,charges))
abline(fit3)

#The dummy variable for smoking is 1 for smokers and 0 for non-smokers. The 
#average medical charges for a non-smoker is estimated to be $5879.42 while the 
#average medical charge for a smoker is estimated to be $19,066 less making the 
#total $-13186.58. Additionally, an increase of about $83.35 is associated with 
#an increase by about 1 unit of BMI for non-smokers. There is increase of about 
#$1,473.11 medical cost for every unit of BMI for smokers; this slope value 
#includes the $83.35 and the $1,389.76 from interaction term. 

#We can conclude that in the absence of any BMI consideration, medical cost on 
#average for non-smokers, fall somewhere between $3963.06 and $7795.79. The 
#medical cost on average for smokers, fall somewhere between -$19206.96 and 
#-$7166.19. Furthermore, for a non-smoker for every 1 unit increase in BMI, 
#there will be an average increase in medical cost of between $22.01 and 
#$144.69. For a smoker for every 1 unit increase in BMI, there will be an 
#increase in medical cost between $1280.75 and $1665.46. 

df$smokerpack = ifelse(df$smoker =="yes" & df$bmi>30, "Yes >30",
                       ifelse(df$smoker =='yes',"Yes <31","No"))
fit4 = lm(formula = charges ~ bmi*smokerpack, data=df)
fit4_sum = summary(fit4)

p4 = ggplot(df,aes(x = bmi, y = charges, color = smokerpack))+geom_point()+
  stat_smooth(method="lm",se=FALSE)
p4 + ggtitle("Linear Model in Figure 01")

#Nonsmokers: charges = 5879.2 + 83.35*BMI
#Smokers, BMI <30: charges = 9070.97 + 485.1*BMI
#Smokers, BMI >30: charges = 23617 + 508.53*BMI

