library(ggplot2)
library(glmnet)
library(ggcorrplot)
library(pracma)

##### Lasso Regression #####
set.seed(167)
#vec = as.matrix(rnorm(30, mean = 10, sd = 3))
X = as.matrix(rnorm(30))
#Y = 3 - 2*vec + 3*vec^2 + rnorm(30)
#Y = as.matrix(Y)
Y = 3 - 2*X + 3*X^2 + rnorm(30)

predmatrix = matrix(0, 30, 7)
for (i in 1:7)
{
  predmatrix[,i] = X^i
}

cv.out = glmnet(x=predmatrix, y=Y, family = "gaussian", alpha = 1, nlambda = 1000)
plot(cv.out)

plot(log(cv.out$lambda), cv.out$beta[1,], type="l", col = "red", main = "Coefficient Value Plot", 
     xlab = expression(paste("log(", lambda,")")), ylab = "Coefficient Values", 
     ylim = c(-2, 3))
lines(log(cv.out$lambda), cv.out$beta[2,], col = "green")
lines(log(cv.out$lambda), cv.out$beta[3,], col = "blue")
lines(log(cv.out$lambda), cv.out$beta[4,], col = "purple")
lines(log(cv.out$lambda), cv.out$beta[5,], col = "orange")
lines(log(cv.out$lambda), cv.out$beta[6,], col = "brown")
lines(log(cv.out$lambda), cv.out$beta[7,], col = "black")
legend("topright", legend = c("beta1", "beta2", "beta3", "beta4", "beta5", "beta6", "beta7"), 
       col = c("red", "green", "blue", "purple", "orange", "brown", "black"), lty=1, cex = 0.8)

#cross validation of nfolds = 10
cv.valid.out = cv.glmnet(x=predmatrix, y=Y, alpha=1, family = "gaussian", nfolds = 10)
plot(cv.valid.out)

bestlam = cv.valid.out$lambda.min
coef(cv.valid.out, s = bestlam)

fitted_model = glmnet(x = predmatrix, y = Y, lambda = bestlam)

X_1000 = as.matrix(rnorm(1000))
Y_1000 = 3 - 2*X_1000 + 3*(X_1000)^2 + rnorm(1000)

predmatrix_1000 = matrix(0, 1000, 7)
for (i in 1:7)
{
  predmatrix_1000[,i] = (X_1000)^i
}
#The best tuning parameter was 0.128 because it was the lambda that gave the 
#smallest value of mean-squared error. 

y_1000 = predict(fitted_model, s = bestlam, newx = predmatrix_1000)
MSE = mean((Y_1000 - y_1000)^2)

#RSS + 0.128|-1.678| + 0.128|2.723|
#Y = 2.948 – 1.678*X+ 2.744*X2

#The equation found through the lasso model has 3 coefficients by zeroing 4 of 
#the original 7 predictors. Through this process, we get an equation similar to 
#Y = 3 – 2*X + 3*X2. When the equation above was used to predict the value for 
#1000 new observations, the MSE generated was 1.351, This is a fairly small value. 


##### Simple Logistic Regression #####

wdbc = read.table("wdbc.data", header = F, sep=",") 
#There are 569 samples and 30 predictors (columns V3 to V32). There are 357 
#observations of benign tumors and 212 observations of malignant tumors. 

set.seed(2)
split = sample(nrow(wdbc),400)
train = wdbc[split, ]
test = wdbc[-split, ]

xtrain = train[,3:ncol(train)]
ytrain = train[,2]

xtest = test[,3:ncol(test)]
ytest = test[,2]

normtrain = matrix(0, nrow = 400, ncol = 30)
nxtrain = matrix(0, nrow = 400, ncol = 30)
normtest = matrix(0, nrow = 169, ncol = 30)
nxtest = matrix(0, nrow = 169, ncol = 30)
for (i in 1:30)
{
  #normalize training predictors
  normtrain[,i] = xtrain[,i] - mean(xtrain[,i])
  nxtrain[,i] = normtrain[,i]/sd(normtrain[,i])
  
  #normalize testing predictors
  normtest[,i] = xtest[,i] - mean(xtest[,i])
  nxtest[,i] = normtest[,i]/sd(normtest[,i])
}
nxtrain = data.frame(nxtrain)
nxtest = data.frame(nxtest)

ytrain = factor(ytrain, levels =c("B","M"))
ytest = factor(ytest, levels =c("B","M"))
numtrain = sapply(ytrain, unclass)
cordata = cor(nxtrain)

ggcorrplot(cordata) + labs(title = "Correlation For Breast Cancer Diagnosis Predictors") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
#The correlation matrix is demonstrating that a lot of the predictors have high 
#correlation with other predictors. Predictors 1-4 have a high correlation with 
#predictors 21-24. This can lead to redundancy in the model when fitting a line 
#and not allow true interpretability. 

#According to the correlation matrix, predictors X1 and X3 have a correlation of 
#0.180. This is closer to 0, which represents a lack of correlation. The value 
#for β_1is -862.997 and the value for β_3 is 989.130. These values are some of 
#the highest coefficient estimates present in the data, alongside β_21 = -1030.876. 
#After these 3 coefficient values, the next biggest coefficient value β_11 = 559.966. 
#This means that β_1,β_3, and β_21 have the biggest impact on the outcome compared 
#to other predictors. However due to the lack of predictor eliminating, a lot of 
#the other predictor values are also present but smaller making us wonder if they 
#have an impact or not. Are they redundant? How much are they adding to model 
#complexity? It’s crazy to see how big the p values are for all the coefficients! 
#They are approximately 0.999!

glm.train = glm(formula = ytrain~., family=binomial(link="logit"), data=nxtrain)
summary(glm.train)

glm.prob.train = predict(glm.train, type="response")
glm.label.train = rep("B", nrow(xtrain))
glm.label.train[glm.prob.train>.5] = "M"
table(glm.label.train, ytrain)
mean(glm.label.train == ytrain)

glm.prob.test = predict(glm.train, type="response", newdata = nxtest)
glm.label.test = rep("B", nrow(xtest))
glm.label.test[glm.prob.test>.5] = "M"
table(glm.label.test, ytest)
mean(glm.label.test == ytest)

#The prediction accuracy for the training set is 100% while the prediction 
#accuracy for the testing set is 94.08%. There is overfitting in this model, 
#which is why the prediction accuracy for training set is 100%. 

##### Ridge Logistic Regression #####

#run initializing codes from Question 2...
#use nxtest, nxtrain

grid = 10^seq(5,-18,length=100)
ridge.mod = glmnet(nxtrain, ytrain, family = "binomial", alpha=0, lambda=grid,
                   thresh=1e-8)
plot(ridge.mod)

plot(log(ridge.mod$lambda), ridge.mod$beta[1,], col = "red", main = "Ridge Regression Coefficient Plot", 
     xlab = expression(paste("log(", lambda,")")), ylab = "Coefficient Values")
points(log(ridge.mod$lambda),ridge.mod$beta[3,], col = "green")
legend("bottomright", legend = c("beta1", "beta3"), 
       col = c("red", "green"), pch = 1, cex = 0.8)
#Both β_1and β_3 head towards an asymptote near 0 as log(λ) increases. The values 
#converge around log(λ)=-10. In addition, β_3 has coefficient values that start 
#off at a smaller value.

#default is 10-fold validation model!
cv.out = cv.glmnet(nxtrain, ytrain, family = "binomial", alpha=0, lambda=grid, 
                   type.measure = "class",thresh=1e-8)
plot(cv.out) #The best lambda that minimizes the misclassification error is 0.00215. 
bestlam = cv.out$lambda.min
coef(cv.out, s = bestlam)
#For this model, all the coefficient estimates are non-zero values. This makes 
#sense because ridge regression models do not perform feature selection to 
#eliminate predictors. This can also be visually represented by the x-axis on 
#top of the graph. 

out = glmnet(nxtrain,ytrain, family="binomial",alpha=0, lambda = bestlam)
glm.prob.train = predict(out, type="response", s=bestlam, newx = nxtrain)

glm.label.train = rep("B", nrow(xtrain))
glm.label.train[glm.prob.train>.5] = "M"
table(glm.label.train, ytrain)
mean(glm.label.train == ytrain)

glm.prob.test = predict(out, type="response", s=bestlam, newx = nxtest)
glm.label.test = rep("B", nrow(xtest))
glm.label.test[glm.prob.test>.5] = "M"
table(glm.label.test, ytest)
mean(glm.label.test == ytest)

#The prediction accuracy from training set decreases by 1.03%. We can conclude, 
#there isn’t as much overfitting with the training set. 

n_segm = 21
TPR = replicate(n_segm, 0)
FPR = replicate(n_segm, 0)
p_th = seq(0, 1, length.out = n_segm)

for (i in 1:n_segm) 
{
  glm.label.test = rep("B", nrow(test))
  glm.label.test[glm.prob.test > p_th[i]] = "M"
  
  TPR[i] = mean(glm.label.test[ytest == "M"] == ytest[ytest == "M"])
  FPR[i] = mean(glm.label.test[ytest == "B"] != ytest[ytest == "B"])
}

auc = -1 * trapz(FPR, TPR)
print(auc)
ggplot() + geom_path(aes(x= FPR, y = TPR)) + geom_point(aes(x = FPR, y = TPR)) + 
  ggtitle("ROC for Ridge Logistic Regression Model") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
#area under ROC curve is 0.991

##### Lasso Logistic Regression #####

#run initializing codes from Question 2
#use nxtest, nxtrain


grid = 10^seq(5,-18,length=100)
lasso.mod = glmnet(nxtrain, ytrain, family = "binomial", alpha=1, lambda=grid,
                   thresh=1e-8)
plot(lasso.mod)

plot(log(lasso.mod$lambda), lasso.mod$beta[1,], col = "red", main = "Lasso Regression Coefficient Plot", 
     xlab = expression(paste("log(", lambda,")")), ylab = "Coefficient Values")
points(log(lasso.mod$lambda),lasso.mod$beta[3,], col = "green")
legend("bottomright", legend = c("beta1", "beta3"), 
       col = c("red", "green"), pch = 1, cex = 0.8)
#β_1 appears to have a constant coefficient value around 0 and β_3 head towards 
#0 as log(λ) increases. The values converge around log(λ)=-10. In addition, β_3
#has coefficient values that start off at a smaller value compared to the 
#initial β_3 for ridge regression. 

#default is 10-fold validation model!
cv.out = cv.glmnet(nxtrain, ytrain, family = "binomial", alpha=1, lambda=grid, 
                   type.measure = "class",thresh=1e-8)
plot(cv.out) #optimal lambda that minimizes misclassification error is 0.00215. 
bestlam = cv.out$lambda.min
coef(cv.out, s = bestlam)
#We started off with 30 predictors in our model but through the feature selection 
#in lasso regression, we end up with only 15 non-zero predictors. This can also 
#be seen through the x-axis on top of the misclassification error graph. 

out = glmnet(nxtrain,ytrain, family="binomial",alpha=1, lambda = bestlam)
glm.prob.train = predict(out, type="response", s=bestlam, newx = nxtrain)

glm.label.train = rep("B", nrow(xtrain))
glm.label.train[glm.prob.train>.5] = "M"
table(glm.label.train, ytrain)
mean(glm.label.train == ytrain)

glm.prob.test = predict(out, type="response", s=bestlam, newx = nxtest)
glm.label.test = rep("B", nrow(xtest))
glm.label.test[glm.prob.test>.5] = "M"
table(glm.label.test, ytest)
mean(glm.label.test == ytest)
#The prediction accuracy from training set decreases by 0.78%. We can conclude, 
#there isn’t as much overfitting with the training set. 

n_segm = 21
TPR = replicate(n_segm, 0)
FPR = replicate(n_segm, 0)
p_th = seq(0, 1, length.out = n_segm)

for (i in 1:n_segm) 
{
  glm.label.test = rep("B", nrow(test))
  glm.label.test[glm.prob.test > p_th[i]] = "M"
  
  TPR[i] = mean(glm.label.test[ytest == "M"] == ytest[ytest == "M"])
  FPR[i] = mean(glm.label.test[ytest == "B"] != ytest[ytest == "B"])
}

auc = -1 * trapz(FPR, TPR)
print(auc)
ggplot() + geom_path(aes(x= FPR, y = TPR)) + geom_point(aes(x = FPR, y = TPR)) + 
  ggtitle("ROC for Lasso Regression Model") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
#area under the ROC curve is 0.998