library(ggplot2)
library(MASS)
library(grid)
library(gridExtra)
library(pracma)
library(numpy)

#Loading the dataset
wdbc = read.table("Desktop/wdbc.data", header = F, sep=",") 
sum(is.na(wdbc)) 
#The data has a sample size n of 569 people, number of predictors p is 32, and 
#the number of observations in each class is 212 malignant cases and 357 benign 
#cases. 
table(wdbc$V2)

set.seed(0) #to ensure reproducibility of data
split = sample(nrow(wdbc),400)
train = wdbc[split, ]
test = wdbc[-split, ]

#### Logistic Regression Model ####
#V2- Diagnosis, V3- Average Radius of Cell Nuclei, V4- Average Texture of Cell Nuclei
plot1 = ggplot(train, aes(x=V3, y=V4, color=V2)) + geom_point() + 
  labs(color = "Diagnosis", x = "Average Radius of Cell Nuclei", 
       y = "Average Texture of Cell Nuclei", title = "Diagnosis Outcomes") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

print(plot1)
#Based on this plot, I don’t think it is possible to predict the outcome with a 
#high accuracy because there is an overlap in diagnosis between the two outcome 
#clusters given average texture and average radius of cell nuclei. However, due 
#to the number of points that are showing two clustered outcomes we can assume 
#a decent amount of accuracy.
  
#glm time
train$V2 = factor(train$V2, levels =c("B","M"))
glm.train = glm(formula = V2~V3+V4, family=binomial(link="logit"), data=train)
summary(glm.train)
contrasts(train$V2)

#The coefficient estimates for average radius of cell nuclei and average texture 
#are both positive in value. As the average radius and texture increases for the 
#cell nuclei, the likelihood of the outcome being a malignant tumor increases. 
#In addition, the p-values for coefficient representing the association between 
#the average radius and diagnosis, and the average texture and diagnosis are 
#less than 0.05 (<2e-16 and 1.1e-07, respectively). This means for both the 
#coefficient estimates, we can reject the null hypothesis that there is no 
#association. 

predict(glm.train,type = "response",data.frame(V3=c(10), V4=c(12)))
#get a log-odds ratio of -6.65 and a predicted probability of 0.00129
glm.prob.train = predict(glm.train, type="response")
glm.label.train = rep("B", nrow(train))
glm.label.train[glm.prob.train>.5] = "M"
table(glm.label.train, train$V2) #confusion matrix
mean(glm.label.train == train$V2)
#accuracy of prediction = 89.8% 

glm.prob.test = predict(glm.train, type="response", newdata = test)
glm.label.test = rep("B", nrow(test))
glm.label.test[glm.prob.test>.5] = "M"
table(glm.label.test, test$V2)
mean(glm.label.test == test$V2)
#accuracy of prediction = 87.6%

#Comparing training and testing: the accuracy of the prediction went down by 
#2.2% for the test data. This means there wasn’t an overfitting issue with the 
#training model. 

x1 = seq(5, 30, length = 100)
x2 = seq(5, 40, length = 100)
grid = expand.grid(x1, x2)
colnames(grid) = c("V3", "V4")

cutoff = c(0.25,0.50,0.75)
gridpred = predict(glm.train, type = "response", newdata=grid)

for (i in cutoff)
{
  glm.label.test = rep("B", nrow(grid))
  glm.label.test[gridpred > i] = "M"
  
  plot2 = ggplot(grid, aes(x=V3, y=V4, color=glm.label.test)) + geom_point() +
    geom_point(train, mapping = aes(x=V3, y=V4, fill = V2), color = "black", pch=21, size = 2) + 
    labs(color= "Test Data", fill = "Train Data", x = "Average Radius of Cell Nuclei", 
         y = "Average Texture of Cell Nuclei", 
         title = paste("Logistic Regression Classifier with Cutoff", i, collapse = ",")) + 
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  print(plot2)
  
  
}

#The number of red points representing benign test data increases as the cutoff 
#increases from 0.25 to 0.75. This is the expected result because only points 
#between 0.75 and 1 are classified as malignant with the 0.75 cutoff. While we 
#are representing more of the benign points in the training set with our increased 
#boundary, we are also capturing more of the malignant testing data points. 
#This is a tradeoff. The boundary slope is linear.

n_segm = 21
TPR = replicate(n_segm, 0)
FPR = replicate(n_segm, 0)
p_th = seq(0, 1, length.out = n_segm)

for (i in 1:n_segm) 
{
  glm.label.test = rep("B", nrow(test))
  glm.label.test[glm.prob.test > p_th[i]] = "M"
  
  TPR[i] = mean(glm.label.test[test$V2 == "M"] == test$V2[test$V2 == "M"])
  FPR[i] = mean(glm.label.test[test$V2 == "B"] != test$V2[test$V2 == "B"])
}

auc = -1 * trapz(FPR, TPR)
print(auc)
ggplot() + geom_path(aes(x= FPR, y = TPR)) + geom_point(aes(x = FPR, y = TPR)) + 
  ggtitle("ROC for Logistic Regression Test") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
#area under the curve is 0.949. 

#### Linear Discriminant Analysis model ####
lda.model = lda(formula = V2~V3+V4, data=train, centre = TRUE, scale = TRUE)
lda.model

lda.pred.train = predict(lda.model, train)
lda.pred.test = predict(lda.model, test)

names(lda.pred.train)
#need to specify $class so that the internal Bayes classifier (set to .5)
#can convert the posterior probabilities to a certain outcome
tt.lda.train = table(lda.pred.train$class, train$V2)
mean(lda.pred.train$class == train$V2)
tt.lda.test = table(lda.pred.test$class, test$V2)
mean(lda.pred.test$class == test$V2)
#accuracy of training data prediction is 89.3%
#accuracy of testing data prediction is 86.4%

lda.pred.test_cutoff = rep("B", nrow(test))
lda.pred.test_cutoff[lda.pred.test$posterior[,2]>.5] = "M"
mean(lda.pred.test_cutoff == test$V2)

x1 = seq(5, 30, length = 100)
x2 = seq(5, 40, length = 100)
grid = expand.grid(x1, x2)
colnames(grid) = c("V3", "V4")

cutoff = c(0.25,0.50,0.75)
gridpred = predict(lda.model, newdata=grid)

for (i in cutoff)
{
  
  lda.label.test = rep("B", nrow(grid))
  lda.label.test[gridpred$posterior[,2] > i] = "M"
  
  plot3 = ggplot(grid, aes(x=V3, y=V4, color=lda.label.test)) + geom_point() +
    geom_point(train, mapping = aes(x=V3, y=V4, fill = V2), color = "black", pch=21, size = 2) + 
    labs(color= "Test Data", fill = "Train Data", x = "Average Radius of Cell Nuclei", 
         y = "Average Texture of Cell Nuclei", 
         title = paste("LDA Classifier with Cutoff", i, collapse = ",")) + 
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  print(plot3)
}

print(plot3)

#The number of red points representing benign test data increases as the cutoff 
#increases from 0.25 to 0.75. This is the expected result because only points 
#between 0.75 and 1 are classified as malignant with the 0.75 cutoff. While we 
#are representing more of the benign points in the training set with our 
#increased boundary, we are also capturing more of the malignant testing data 
#points. This is a tradeoff. And this looks like the image boundary for logistic 
#regression.


#ggplot(grid, aes(x=V3, y=V4, color=lda.label.test)) + geom_point() + 
  #geom_point(train, mapping = aes(x=V3, y=V4, shape = V2), color = "black")

n_segm = 21
TPR = replicate(n_segm, 0)
FPR = replicate(n_segm, 0)
p_th = seq(0, 1, length.out = n_segm)

for (i in 1:n_segm) 
{
  lda.pred.test_cutoff = rep("B", nrow(test))
  lda.pred.test_cutoff[lda.pred.test$posterior[,2]> p_th[i]] = "M"
  TPR[i] = mean(lda.pred.test_cutoff[test$V2 =="M"] == test$V2[test$V2 == "M"])
  FPR[i] = mean(lda.pred.test_cutoff[test$V2 =="B"] != test$V2[test$V2 == "B"])
}

auc = -1 * trapz(FPR, TPR)
print(auc)
ggplot() + geom_path(aes(x= FPR, y = TPR)) + geom_point(aes(x = FPR, y = TPR)) + 
  ggtitle("ROC for LDA") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

#area under the curve is 0.949. 

#### Quadratic Discriminant Aanalysis ####
qda.model = qda(formula = V2~V3+V4, data=train, centre = TRUE)
qda.model

qda.pred.train = predict(qda.model, train)
qda.pred.test = predict(qda.model, test)

names(qda.pred.train)
#need to specify $class so that the internal Bayes classifier (set to .5)
#can convert the posterior probabilities to a certain outcome
tt.qda.train = table(qda.pred.train$class, train$V2)
mean(qda.pred.train$class == train$V2)
tt.qda.test = table(qda.pred.test$class, test$V2)
mean(qda.pred.test$class == test$V2)
#accuracy of prediction if 89.00%
#accuracy of prediction is 86.98%

x1 = seq(5, 30, length = 100)
x2 = seq(5, 40, length = 100)
grid = expand.grid(x1, x2)
colnames(grid) = c("V3", "V4")

cutoff = c(0.25,0.50,0.75)
gridpred = predict(qda.model, newdata=grid)

for (i in cutoff)
{
  
  qda.label.test = rep("B", nrow(grid))
  qda.label.test[gridpred$posterior[,2] > i] = "M"
  
  plot4 = ggplot(grid, aes(x=V3, y=V4, color=qda.label.test)) + geom_point() +
    geom_point(train, mapping = aes(x=V3, y=V4, fill = V2), color = "black", pch=21, size = 2) + 
    labs(color= "Test Data", fill = "Train Data", x = "Average Radius of Cell Nuclei", 
         y = "Average Texture of Cell Nuclei", 
         title = paste("QDA Classifier with Cutoff", i, collapse = ",")) + 
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  print(plot4)
}

print(plot4)

#The number of red points representing benign test data increases as the cutoff 
#increases from 0.25 to 0.75. This is the expected result because only points 
#between 0.75 and 1 are classified as malignant with the 0.75 cutoff. While we 
#are representing more of the benign points in the training set with our 
#increased boundary, we are also capturing more of the malignant testing data 
#points. This is a trade-off. The shape is very interesting! The image boundary 
#is not linear. And in the 0.25 and 0.5 cutoff, we can see that there are 3 
#regions. It’s interesting that the QDA classifier predicts malignant tumors on 
#the leftmost section when there are no data points to suggest anything!

n_segm = 21
TPR = replicate(n_segm, 0)
FPR = replicate(n_segm, 0)
p_th = seq(0, 1, length.out = n_segm)

for (i in 1:n_segm) 
{
  qda.pred.test_cutoff = rep("B", nrow(test))
  qda.pred.test_cutoff[qda.pred.test$posterior[,2]> p_th[i]] = "M"
  TPR[i] = mean(qda.pred.test_cutoff[test$V2 =="M"] == test$V2[test$V2 == "M"])
  FPR[i] = mean(qda.pred.test_cutoff[test$V2 =="B"] != test$V2[test$V2 == "B"])
  
}

auc = -1 * trapz(FPR, TPR)
print(auc)
ggplot() + geom_path(aes(x= FPR, y = TPR)) + geom_point(aes(x = FPR, y = TPR)) + 
  ggtitle("ROC for QDA") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
#area under the curve is 0.953

#### KNN Classifier ####
#check to see what's a categorical value
sapply(train, is.factor) #TRUE = categorical
train$V2 = factor(train$V2, levels =c("B","M"))
train.num = train[,!sapply(train, is.factor)]
train.num = scale(train.num)
train.num = train.num[,2:3]

test$V2 = factor(test$V2, levels =c("B","M"))
sapply(test, is.factor)
test.num = test[,!sapply(test, is.factor)]
test.num = scale(test.num)
test.num = test.num[,2:3]

neighbors = c(1, 2, 3, 4, 20)
set.seed(0)
for (i in neighbors)
{
  knn.train.pred <- knn(train = train.num, test = train.num, cl = train$V2, k= i)
  
  tab = table(knn.train.pred, train$V2)
  me = mean(knn.train.pred == train$V2)
  print(tab)
  print(me)
}

set.seed(123)
for (i in neighbors)
{
  knn.test.pred <- knn(train = train.num, test = test.num, cl= train$V2, k= i)
  tab2 = table(knn.test.pred, test$V2)
  me2 = mean(knn.test.pred == test$V2)
  print(i)
  print(tab2)
  print(me2)
}

#In the training set, there is a general decrease in prediction accuracy as the 
#k-number neighbors increases. However, from k=2 to k=20, it looks more random 
#and like a straight line. In the testing set, there is an increase in prediction 
#accuracy for k=1 to k=4, but then around k=20 there is a slight decrease in the 
#prediction accuracy. This pattern is better represented by the figure in 4c. 
#In addition, the prediction accuracy for training set it higher than that for 
#the testing set. It makes sense that the accuracy is 100% for k=1 in the 
#training set because that is the direct outcome of the given dataset. 

x1 = seq(5, 30, length = 100)
x2 = seq(5, 40, length = 100)
grid = expand.grid(x1, x2)
colnames(grid) = c("V3", "V4")
grid.num = scale(grid)

for (i in neighbors)
{
  knn.grid.pred <- knn(train = train.num, test = grid.num, cl= train$V2, k= i)
  plot5 = ggplot(grid, aes(x=V3, y=V4, color=knn.grid.pred)) + geom_point() +
    geom_point(train, mapping = aes(x=V3, y=V4, fill = V2), color = "black", pch=21, size = 1) + 
    labs(color= "Test Data", fill = "Train Data", x = "Average Radius of Cell Nuclei", 
         y = "Average Texture of Cell Nuclei", 
         title = paste("kNN Classifier set to k =", i, collapse = ",")) + 
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  print(plot5)
}

#The classifier is very interesting! From k = 1 to k = 20, we see the formation 
#of two major regions representing benign and malignant tumors. This makes sense 
#because the testing data is using 20 nearest points around it to estimate the 
#outcome of the point. If it relied on only 1 point, then the image would appear 
#more pixelated depending on the scattering of benign and malignant data points 
#(as seen in the images for k=1, 2, 3, 4); it’s too flexible!

set.seed(123)
max_num = 20
k_vector = c(1:max_num)
train.mean = replicate(max_num,0)
test.mean = replicate(max_num,0)
for (i in k_vector)
{
  knn.train.pred <- knn(train = train.num, test = train.num, cl= train$V2, k= i)
  knn.test.pred <- knn(train = train.num, test = test.num, cl= train$V2, k= i)
  train.mean[i] = mean(knn.train.pred == train$V2)
  test.mean[i] = mean(knn.test.pred == test$V2)
}

df = data.frame(k_vector, train.mean, test.mean)

ggplot(data=df, aes(x=k_vector)) + geom_point(aes(y = train.mean, group =1, color="lightslateblue")) + 
  geom_point(aes(y = test.mean, group = 2, color = "mediumvioletred")) + 
  labs(color = "Type of Dataset", title = "Prediction Accuracy for kNN Classifier", x = "k-Nearest Neighbor",
       y = "Prediction Accuracy") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  scale_color_manual(labels=c("Training","Testing"), values = c("lightslateblue","mediumvioletred")) + 
  geom_line(aes(y = train.mean, color = "lightslateblue")) +
  geom_line(aes(y = test.mean, color = "mediumvioletred"))

#The prediction accuracy decreases until it somewhat stabilizes around a certain 
#region for the training set. For the testing set, it increases and then decreases. 
#Our main goal is to increase the prediction accuracy of the testing set. The 
#highest point of accuracy is k=6 for the testing set. I would pick a range 
#around there. Another interesting thing to note is that when I was testing 
#different seeds, I kept getting slightly different results for the best k. 