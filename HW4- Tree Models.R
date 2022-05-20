library(tree)
library(randomForest)
library(gbm)
library(boot)
library(ggplot2)

##### Tree Model for Heart Disease Prediction #####
load("heart.Rdata")

full$Disease = factor(full$Disease, levels =c("No Disease","Heart Disease"))
#The sample size 371 people, there are 13 predictors, and 171 with no heart 
#disease and 200 with heart disease. 

set.seed(2)
train_ids = sample(nrow(full),200)
test_ids = seq(nrow(full))[-train_ids]

train = full[train_ids, ]
test = full[test_ids, ]

tree.med = tree(Disease~., train)
summary(tree.med)

#visual plot of the overgrown tree
plot(tree.med)
text(tree.med)

#overgrown tree
guess_train = predict(tree.med, newdata = train, type = "class")
table(guess_train, train$Disease)
mean(guess_train == train$Disease)

guess_test = predict(tree.med, newdata = test, type = "class")
table(guess_test, test$Disease)
mean(guess_test == test$Disease)
#The misclassification error for the training set is 12% and for the testing set 
#is 22.81%. The error increases by about 10% so we can see that there is some 
#overfitting in the training dataset. 

set.seed(2)
#10 is the default number of k-folds
cv.med = cv.tree(tree.med, FUN = prune.misclass)
plot(cv.med$size, cv.med$dev, type = "b", main = "Subtree Size vs. Misclassification Error", 
     xlab = "Tree Size", ylab = "Misclassification Error") #shows that overgrown tree is a bad fit

min_idx = which.max(cv.med$dev)
#By using the minimum function, I found that the best tree size that minimizes 
#the misclassification error is 7. 
prune.med = prune.tree(tree.med, best = 7)
#visual plot of the pruned tree
plot(prune.med)
text(prune.med, pretty = 0)

seat_true_prune = prune.misclass(tree.med, best = 7)
summary(seat_true_prune)

set.seed(2)
pred = predict(seat_true_prune, train, type = "class")
table(predicted = pred, actual = train$Disease)
mean(pred == train$Disease)

pred_test = predict(seat_true_prune, test, type = "class")
table(predicted = pred_test, actual = test$Disease)
mean(pred_test == test$Disease)
#The misclassification error for the training set is 12.5% and for the testing 
#set is 21.64%. The error increases by about 10% so we can see that there is 
#some overfitting in the training dataset. However, this is still a slight 
#improvement compared to the overgrown tree. The testing set and training set 
#misclassification error for the overgrown tree was 1.27% higher and 0.5% lower, 
#respectively. 

### Bagged trees
set.seed(2)
bag.med = randomForest(Disease~., data = train, mtry = 12, importance = TRUE)

pred_bag_train = predict(bag.med, newdata = train)
table(pred_bag_train, train$Disease)
mean(pred_bag_train == train$Disease)

pred_bag_test = predict(bag.med, newdata = test)
table(pred_bag_test, test$Disease)
mean(pred_bag_test == test$Disease)

#The misclassification error for the training set is 0% and for the testing set 
#is 22.81%. There is overfitting. It is interesting to see the prediction 
#accuracy for training set going to 0%. 

### Random Forests

set.seed(2)
# m = p/3 ... 12/3 ... 4
rf.med = randomForest(Disease~., data=train, mtry = 4, importance = TRUE)
pred_rf_train = predict(rf.med, newdata = train, type = "class")
table(pred_rf_train, train$Disease)
mean(pred_rf_train == train$Disease)

pred_rf_test = predict(rf.med, newdata = test, type = "class")
table(pred_rf_test, test$Disease)
mean(pred_rf_test == test$Disease)

#The misclassification error for the training set is 0% and for the testing set 
#is 21.05%. There is overfitting but it is an improvement from bagged trees by 
#a little bit. 

### Boosted tree
set.seed(2)
boost.med = gbm((unclass(Disease)-1)~., data=train, distribution = "bernoulli", n.trees = 500,
                interaction.depth = 2, shrinkage = 0.1)

pred_boost_train = predict(boost.med, newdata = train, type = "response")
label.train = rep("No Disease", nrow(train))
label.train[pred_boost_train>.5] = "Heart Disease"
table(label.train, train$Disease)
mean(label.train == train$Disease)


pred_boost_test = predict(boost.med, newdata = test, type = "response")
label.test = rep("No Disease", nrow(test))
label.test[pred_boost_test>.5] = "Heart Disease"
table(label.test, test$Disease)
mean(label.test == test$Disease)

#The misclassification error for the training set is 0% and for the testing set 
#is 19.30%. There is overfitting but it is an improvement from all the other 
#models! 