################## Install packages ##########################

install.packages("glmnet") # Package for CV analysis
library(glmnet)
install.packages("lars") # Package for data set 
library(lars)
install.packages("caret")
library(caret)
install.packages("plotly")
library(plotly)

---------------------------------- #PART 1 ----------------------------------

################# Load and split data into training and test #########################

data(diabetes) #load dataset 
smp_size <- floor(0.905 * nrow(diabetes)) #fraction of the sample size
set.seed(1212) # set the seed to make your partition reproducible
train_ind <- sample(seq_len(nrow(diabetes)), size = smp_size) #take random sample
train.set <- diabetes[train_ind, ] #index
test.set <- diabetes[-train_ind, ] #index

######################################### RIDGE ##################

#Cross variant ridge (alpha is 0) test on train set using mean square errors and 10 fold validation
cvridge <- cv.glmnet(train.set$x2, train.set$y, type.measure="mse", nfolds=10, alpha=0) 

################# Plotting MSE against lambda
par(mfrow=c(1,1))
lambda.r <- cvridge$lambda #Name lambda for train ridge
bestlamtrainridge = cvridge$lambda.min #Name minimumlambda for train ridge


######################################### LASSO ########################################
#Cross variant lasso, training dataset
cvlasso <- cv.glmnet(train.set$x2, train.set$y, type.measure="mse", nfolds=10, alpha=1)
attributes(cvlasso)
par(mfrow=c(2,2))
plot(cvlasso,xlab='Log of Lambda values',main='Lasso',cex.main = 0.8)
plot(cvridge,xlab='Log of Lambda values',main='Ridge',cex.main = 0.8) 
mtext('Rregression models of Lambda values against MSE with features', outer=TRUE,  cex=1, line=-1.5)

bestlamtrainlasso = cvlasso$lambda.min #Name minimum lambda for lasso
lambda.l <- cvlasso$lambda #Name lambda for train lasso
################## SPS RIDGE #######################

################# TEST 
wh.r <- which(lambda.r >= bestlamtrainridge & lambda.r <= cvridge$lambda.1se) # Identify Interval indices between lambda min and 1se
sps.testset.ridge <- lambda.r[wh.r] 
prd.sps.testset.ridge <- predict(cvridge,newx=test.set$x2,s=sps.testset.ridge ) #predictions using ridge 
mse <- function (i,y) {
  mean ( (y-i)^2 ) }	
mse.sps.testset.ridge <- apply(prd.sps.testset.ridge , 2, mse, test.set$y) #mean square error perdictions

#################### TRAIN 

wh.r <- which(lambda.r >= bestlamtrainridge & lambda.r <= cvridge$lambda.1se) # Identify Interval indices between lambda min and 1se 
sps.trainset.ridge <- lambda.r[wh.r] 
prd.sps.trainset.ridge <- predict(cvridge,newx=train.set$x2,s=sps.trainset.ridge ) #predictions using ridge 

#create a function to take the mean square error 
mse <- function (i,y) {
  mean ( (y-i)^2 ) }	

mse.sps.trainset.ridge <- apply(prd.sps.trainset.ridge , 2, mse, train.set$y) #mean square error perdictions

################# PLOTTING
par(mfrow=c(2,2))
plot(sps.testset.ridge,mse.sps.testset.ridge,xlab='Ridge Prediction',ylab='MSE',main='Test set')
plot(log(sps.testset.ridge),mse.sps.testset.ridge,xlab='Ridge Log Prediction',ylab='MSE',main='Log Plot of the Test set')
plot(sps.trainset.ridge,mse.sps.trainset.ridge,xlab='Ridge Prediction',ylab='MSE',main='Train set')
plot(log(sps.trainset.ridge),mse.sps.trainset.ridge,xlab='Ridge Prediction',ylab='MSE',main='Log Plot of the Train set')
mtext("Ridge regresion models against MSE ", outer=TRUE,  cex=1, line=-1.5)

################## SPS LASSO  ###############

################ TEST

wh.l <- which(lambda.l >= bestlamtrainlasso & lambda.l <= cvlasso$lambda.1se)  # Identify Interval indices between lambda min and 1se 
sps.testset.lasso <- lambda.l[wh.l] 
prd.sps.testset.lasso <- predict(cvlasso,newx=test.set$x2,s=sps.testset.lasso ) #predictions using ridge 
prd.sps.testset.lasso.fl <- predict(cvlasso,newx=test.set$x2,s=c(bestlamtrainlasso,cvlasso$lambda.1se)) 

mse <- function (i,y) {
  mean ( (y-i)^2 ) }	

## against cvx$lambda.1se
apply(prd.sps.testset.lasso.fl, 2, mse, test.set$y)

mse.sps.testset.lasso <- apply(prd.sps.testset.lasso ,2,mse,test.set$y)  #mean square error perdictions

#################### TRAIN 
sps.trainset.lasso <- lambda.l[wh.l] 
prd.sps.trainset.lasso <- predict(cvlasso,newx=train.set$x2,s=sps.trainset.lasso )
prd.sps.trainset.lasso.fl <- predict(cvlasso,newx=train.set$x2,s=c(bestlamtrainlasso,cvlasso$lambda.1se))
mse <- function (i,y) {
  mean ( (y-i)^2 ) }	
apply(prd.sps.trainset.lasso.fl, 2, mse, train.set$y)
mse.sps.trainset.lasso <- apply(prd.sps.trainset.lasso ,2,mse,train.set$y)

################# PLOTTING

par(mfrow=c(2,2))
plot(sps.testset.lasso,mse.sps.testset.lasso,xlab='Lasso Prediction',ylab='MSE',main='Test set')
plot(log(sps.testset.lasso),mse.sps.testset.lasso,xlab='Lasso Log Prediction',ylab='MSE',main='Log Plot of the Test set')	
plot(sps.trainset.lasso,mse.sps.trainset.lasso,xlab='Lasso Prediction',ylab='MSE',main='Train set')
plot(log(sps.trainset.lasso),mse.sps.trainset.lasso,xlab='Lasso Log Prediction',ylab='MSE',main='Log Plot of Train set')	
mtext("Lasso regresion models against MSE ", outer=TRUE,  cex=1, line=-1.5)

################# Comparing models #########################

#ANOVA TEST
aov.out = aov(mse.sps.testset.lasso ~ sps.testset.lasso, data=diabetes)
aov.out

################# PREDICTIONS 
#Predicting Y using Ridge and lasso with the minimum lambda of test set
ridge_y_predictions_test <- predict(cvridge,newx=test.set$x2,s=bestlamtrainridge )
lasso_y_predictions_test <- predict(cvlasso,newx=test.set$x2,s=bestlamtrainlasso )
#Predicting Y using Ridge and lasso with the minimum lambda of train set
lasso_y_predictions_train <- predict(cvlasso,newx=train.set$x2,s=bestlamtrainlasso ) #Y predictions compare to test.set 
lasso_y_predictions_train <- predict(cvlasso,newx=train.set$x2,s=bestlamtrainlasso )

################ MEANS
mean((test.set$y - ridge_y_predictions_test[,1])^2) # MSE is 3160.594 
mean((test.set$y - cvridge$lambda.min)^2) # MSE is 21039.18 bigger than interval
mean((test.set$y - lasso_y_predictions_test[,1])^2) # MSE is 3050.473 
mean((test.set$y - cvlasso$lambda.min)^2) # MSE is 24512.23


################# Normal distrubution approximation ###########################

################## Test ridge lasso

################## Histograms
par(mfrow=c(1,2))
hist(prd.sps.testset.lasso,xlab='Predicted SPS values',main='Lasso',cex.main = 0.8)
hist(prd.sps.testset.ridge,xlab='Predicted SPS values',main='Ridge',cex.main = 0.8 )
mtext('Histograms of Lasso and Ridge models for test set', outer=TRUE,  cex=1, line=-1.5)

################## QQ plots
par(mfrow=c(2,2))
qqnorm(prd.sps.testset.lasso,main = ' Testing set');qqline(prd.sps.testset.lasso, col = 2)
qqnorm(prd.sps.trainset.lasso,main = 'Training set');qqline(prd.sps.trainset.lasso, col = 2)
mtext('Lasso models Normal QQ plots ', outer=TRUE,  cex=1, line=-1.5)


#train ridge lasso

################## Histograms
par(mfrow=c(1,2))
hist(prd.sps.trainset.lasso,xlab='Predicted SPS values',main='Lasso',cex.main = 0.8)
hist(prd.sps.trainset.ridge,xlab='Predicted SPS values',main='Ridge',cex.main = 0.8 )
mtext('Histograms of Lasso and Ridge models for train set', outer=TRUE,  cex=1, line=-1.5)


################## QQ plots
par(mfrow=c(2,2))
qqnorm(prd.sps.testset.ridge,main = ' Testing set');qqline(prd.sps.testset.ridge, col = 2)
qqnorm(prd.sps.trainset.ridge,main = ' Training set' );qqline(prd.sps.trainset.ridge, col = 2)
mtext('Ridge model Normal QQ plots ', outer=TRUE,  cex=1, line=-1.5)



---------------------------------- #PART 2 ----------------------------------


######################################## ELASTIC NET ##############################

enet <- trainControl(method = "cv", number = 10)
set.seed(1212)
db_net1 <- train(train.set$x2,train.set$y,method = "glmnet", trControl = enet)
db_net1$results
w1 <- which.min(db_net1$results[,3]) 
db_net1$results[w1,] 
mod_enet1 <- glmnet(train.set$x2,train.set$y,alpha=db_net1$results[w1,1],
                    lambda = db_net1$results[w1,2])
mod_enet1 <- glmnet(train.set$x2,train.set$y,alpha=db_net1$results[w1,1])
mean((train.set$y-as.vector(predict(mod_enet1,train.set$x2, s= db_net1$results[w1,2])))^2)


################# Grid search ######################
tuneGrid <- expand.grid(.alpha = seq(0.1, 0.9, length = 10),.lambda = seq(.5,7.5,1)) #ommit ridge and lasso queries by elimiating alpha =0 and alpha = 1
dim(tuneGrid)
head(tuneGrid)

set.seed(1212)
db_net2 <- train(train.set$x,train.set$y,method = "glmnet", trControl = enet, tuneGrid=tuneGrid)

################ LOOPS AND RESULTS #######################
########################PLOTTING
par(mfrow=c(1,2))
plot(db_net2$results[,c(1,3)])
plot(db_net2$results[,c(2,3)])
mtext("Alpha and Lambda values against RMSE ", outer=TRUE,  cex=1, line=-1.5)
par(mfrow=c(1,1))
install.packages("plotly")
library(plotly)
alphas <- db_net2$results[,"alpha"]
lambdas <- db_net2$results[,"lambda"]
RMSE <- db_net2$results[,"RMSE"]
plot_ly(x=alphas, y=lambdas, z=RMSE, type="scatter3d", mode="markers", color=RMSE) #scatterplot generated using plotly


############ For loops 

tune.alpha <- tuneGrid$.alpha
list.of.fits <- list(tuneGrid) # convert tuneGrid variables to a list

# Code for loop adpated from https://statquest.org/2018/10/23/statquest-ridge-lasso-and-elastic-net-regression-in-r/#code

for (i in 0.1:0.9) { #ommit lasso and ridge regressions
  fit.name <- paste0("alpha", i/10) 
  
  list.of.fits[[fit.name]] <-
    cv.glmnet(train.set$x2, train.set$y, type.measure="mse",alpha=i/10)
}

results <- data.frame() #Dataframe storing mse results

for (i in 0.1:0.9) {
  fit.name <- paste0("alpha",i/10) 
  
  predicted <- predict(list.of.fits[[fit.name]], s=list.of.fits[[fit.name]]$lambda.min,newx=test.set$x2) #predict y response variables using list created using tunegrid variables
  
  elastic.mse <- mean((test.set$y - predicted)^2)
  
  temp <- data.frame(alpha=i/10, mse=elastic.mse, fit.name=fit.name)
  
  el.test.results <- rbind(results, temp) #results presenting alpha value and mse value for 
}

################## Find best model 

w2 <- which.min(db_net2$results[,3]) #which is the minimum value of alpha
db_net2$results[w2,] 
summary(db_net2$results[,"RMSE"])
summary(db_net2$results[,"RMSESD"])
db_net2$bestTune # can also be done using bestTune

################### MSE of test set - repeating SPS steps

el.sps <- glmnet(train.set$x2, test.set$y, alpha=db_net2$results[w2 , 1]) #predict using glmnet with the alpha value set as the minimum of the elastic net model
mean((test.set$y - as.vector(predict(el.sps, test.set$x2, s=db_net2$results[w2,2])))^2)
cvelastic <- cv.glmnet(train.set$x2, train.set$y, type.measure="mse", nfolds=10, alpha=db_net2$results[w2 , 1]) #using previously identified minimum alpha
bestlamelastic <- cvelastic$lambda.min #define minimum lambda
cvelastic.lambda <- cvelastic$lambda  # define all lambda values
elastic.lambda <- sort(cvelastic$lambda) #sort lambda values 
elastic.wh <- which(cvelastic.lambda >= bestlamelastic & cvelastic.lambda <= cvelastic$lambda.1se) #find the interval of the lambda indices
elastic.sps.set <- cvelastic.lambda[elastic.wh] 
elastic.prd.sps.set <- predict(cvelastic,newx=test.set$x2,s=elastic.sps.set) 
elastic.prd.sps.set.train <- predict(cvelastic,newx=train.set$x2,s=elastic.sps.set)
elastic.mse.sps.set <- apply(elastic.prd.sps.set,2,mse,test.set$y) #mean square error computation for test 
elastic.mse.sps.set.train <- apply(elastic.prd.sps.set.train ,2, mse, train.set$y) #mean square error computation for training

############# PLOTTING 
par(mfrow = c(2,2))
plot(elastic.sps.set,elastic.mse.sps.set, xlab = "Lambda", ylab = "MSE", main = "Test Set")
plot(log(elastic.sps.set),elastic.mse.sps.set, xlab = "Log Lambda", ylab = "MSE", main = "Test Set")
plot(elastic.sps.set,elastic.mse.sps.set.train, xlab = "Lambda", ylab = "MSE", main = "Training Set")
plot(log(elastic.sps.set),elastic.mse.sps.set.train, xlab = "Log Lambda", ylab = "MSE", main = "Training Set")
mtext('Elastic net model against Lambda on Test and Train sets', outer=TRUE,  cex=1, line=-1.5)


################## Histograms
par(mfrow=c(1,2))
hist(elastic.prd.sps.set,xlab='Predicted SPS values',main='Test',cex.main = 0.8)
hist(elastic.prd.sps.set.train,xlab='Predicted SPS values',main='Train',cex.main = 0.8 )
mtext('Histograms ofElastic models for train and test set', outer=TRUE,  cex=1, line=-1.5)


---------------------------------- #PART 3 ----------------------------------
# Which parameters are least influential 

ridge.matrix <- as.matrix(coef(cvridge,s = sps.testset.ridge))
lasso.matrix <- as.matrix(coef(cvlasso,s = sps.testset.lasso))
elastic.matrix <- as.matrix(coef(cvelastic,s = elastic.sps.set))







