#Akaash Sidhu, 0850326
#Nikos A. Manesis,
#Phyllis Lam, 1068632
#Assignment 3: Supervised Machine Learning for Breast Cancer Data

#Install and Load the Required Packages----

install.packages("pROC")
install.packages("class")
install.packages("glue")
install.packages("gower")
install.packages("e1071")


library(e1071)
library(glmnet)
library(class)
library(caret)
library(ggplot2)
library(GGally)
library(pROC)
library(MASS)
library(rpart)

##************************************
#Load the Dataset and Set Seed----
##************************************

cancerdata <- read.csv("http://archive.ics.uci.edu/ml/machine-learning-databases/breast-cancer-wisconsin/wdbc.data",head=FALSE)
str(cancerdata)
names(cancerdata) <- c("Subject ID", "IC", "V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13", "V14", "V15", "V16", "V17", "V18", "V19", "V20", "V21", "V22", "V23", "V24", "V25", "V26", "V27", "V28", "V29", "V30")
set.seed(1245) 

#Train Lasso-logistic, KNN and Decision Tree models to predict cancer prognosis using a training dataset. Here are some guiding steps for your analysis.

#Split the entire dataset into a training set and a test set by randomly selecting 119  observations/instances in your test set (use set.seed(1245) before partitioning the dataset). You should keep the train and test set unchanged for the rest of your work. 

##************************************
#Logistic Regression via Lasso----
##************************************

cancerdata$IC <- as.numeric(cancerdata$IC)
#Convert the IC into a numeric where B = 1, M = 2

f0 <- lm(IC~.^2, data=cancerdata[,-1])
#Create a model using all of the variables with the exception of the subject ID

summary(f0)
str(model.matrix(f0))
head(model.matrix(f0))[,1:31]
#There is an intercept that needs to be removed

X <- model.matrix(f0)[,-1]
X <- cbind(X,X[,1:30]^2)
#Add the quadratic effects

colnames(X)
colnames(X) <- c(colnames(model.matrix(f0)[,-1]), paste(colnames(X)[1:30],-2,sep=""))
colnames(X)
#There are duplicates for V1-V30 near the end, thus to differential V1 will become V1-2.

table(cancerdata$IC)
table(cancerdata$IC+1)
y <- (cancerdata$IC+1)/2
table(y)

test.index <- sample.int(dim(X)[1],round(dim(X)[1]*0.20913884), replace = FALSE)
#The test index contains 119 observations. Using this index, the training and test set can be created simply by using -test.index for when training the data. 

train.set <- cancerdata[-test.index,][,-1]
test.set <- cancerdata[test.index,][,-1]

#b)Train a Lasso-logistic model via 10-fold cross-validation.
#i.Extract the model corresponding to lambda.1se. 

cv.lasso <- cv.glmnet(X[-test.index,], y[-test.index], nfolds = 10, family="binomial", alpha=1, type.measure = "auc")

plot(cv.lasso)

lambda.1se <- cv.lasso$lambda.1se
#Use lambda.1se for the model

#ii.Use this model to derive a classification cut-off point using probability estimates obtained for the training set observations. Specifically, use max(min(sensitivity, specificity)) criterion covered in the class to select the optimal probability cut-off .

prds.train <- predict(cv.lasso,newx = X[-test.index,], type = "response", s=lambda.1se)[,1]
prds.test <- predict(cv.lasso,newx = X[test.index,], type = "response", s=lambda.1se)[,1]
#Predictions for the training set and test set.

sn.sp <- function(mat){
  sn <- mat[2,2]/sum(mat[2,])
  sp <- mat[1,1]/sum(mat[1,])
  return(unlist(list(sensitivity=sn, specificity=sp)))
}
#This is a function to compute sensitivity and specificity. 

auc.train <- roc(y[-test.index],prds.train)
auc.train
#The AUC of the is closer to the maximum (1) which means this is very good for the classifier.

par(mfrow=c(1,2))
plot(auc.train, main = "AUC of Train")

auc.test <- roc(y[test.index],prds.test)
auc.test
#The AUC of the test is also close to 1 which means this is very good for the classifier. 

plot(auc.test, main = "AUC of Test")
par(mfrow=c(1,1))

snsp.train <- cbind(auc.train$sensitivities,auc.train$specificities)
#Determine the sensitivity and specificity metrics for the training set for each set of instances. 

indx <- which.max(apply(snsp.train,1,min))
indx
#This is the instance that results in the max specificity and sensitivity

snsp.train[indx,]
#This is the max sensitivity and specificity

cutoff <- auc.train$thresholds[indx]
cutoff
#This is the cutoff for the max specificity and sensitivity. 

sn.sp(table(y=y[-test.index],yhat=as.numeric(prds.train>cutoff)))
#Sensitivity and Specificity of the Training Set
sn.sp(table(y=y[ test.index],yhat=as.numeric(prds.test >cutoff)))
#Sensitivity and Specificity of the Test Set

#iii.Evaluate model performance on train and test sets by visualizing the ROC curve(s) and by computing metrics such as sensitivity, specificity and AUC.

par(mfrow=c(1,2))
plot(auc.train)
abline(h=snsp.train[indx,1],v=snsp.train[indx,2], col='blue', lty=2)
plot(auc.test)
abline(h=snsp.train[indx,1],v=snsp.train[indx,2], col='blue', lty=2)
par(mfrow=c(1,1))

##************************************
#Logistic Regression via KNNN----
##************************************

#i.Use the relevant examples and R code covered in class to train your KNN model for breast cancer data using 10-fold cross validation. 

ctrl <- trainControl(method="cv", number=10)

train.set$IC <- as.factor(train.set$IC)

cv.knn <- train(IC~., data = train.set, method = "knn", trControl = ctrl, preProcess = c("center", "scale"), tuneLength = 30)

cv.knn.results <- cv.knn$results

cv.knn.results[which.max(cv.knn.results[,2]),]

plot(cv.knn)

#ii.Compare model performance on the train and test sets by computing various performance metrics.

knn.prds.test <- predict(cv.knn, newdata = test.set, probability = TRUE)
knn.prds.train <- predict(cv.knn, probability = TRUE)
#Predictions for the train and test set


table(cancerdata$IC[-test.index], knn.prds.train)
table(cancerdata$IC[test.index], knn.prds.test)

IC.test <- as.factor(test.set$IC)

confusionMatrix(knn.prds.test, IC.test)


##************************************
#Logistic Regression via Decision Tree----
##************************************

#i.Use the relevant examples and R code covered in class to train a pruned classification tree. Note that rpart performs 10-fold cross-validation by default to compute cv-error for various cp values. So, you do not need to supply number of folds here.

par(mfrow=c(1,1))

tree.model <- rpart(IC~.,data=train.set)
printcp(tree.model)

summary(tree.model)
plotcp(tree.model)
#Summaries and plot

plot(tree.model, branch=0, uniform=TRUE)
text(tree.model,use.n=T)	
#Comment comment

prune.cp <- function(cptable){
  
  CP <- cptable[,1]
  cp <- sqrt(CP * c(Inf, CP[-length(CP)])) 	### inspect the code inside plotcp()!
  xerr <- cptable[,4]
  minrow <- which.min(xerr)
  
  xerr.1se <- sum(cptable[minrow,4:5])
  index <- min(which(xerr < xerr.1se))
  
  cp.1se <- cp[index]
  return(as.numeric(cp.1se) )}
#Function to prune

prune.cp(tree.model$cptable)

treept <- prune(tree.model, cp = prune.cp (tree.model$cptable))
summary(treept)

par(mfrow=c(1,2))

plot(treept,margin=0.1)
text(treept,use.n=T)

plot(train.set$V28,train.set$V21,pch=15*train.set$IC+1,type='n')

polygon(c(-0.7,0.13,0.13,0.7,0.7,-0.7),c(1.74,1.74,0,0,6,6),col='gray')
points(train.set$V28,train.set$V21,pch=15*train.set$IC)
abline(h=1.74,lty=3)
abline(v=0.13,lty=3)
par(mfrow=c(1,1))
#Plot i think this is wrong

prds <- predict(,type='class')
table(, tree.prds)

tree.prds <- predict(treept, type = 'class')

#ii.Compare model performance on the train and test sets by computing various performance metrics.
