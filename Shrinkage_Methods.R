# Shrinkage-methods-in-R
#Analysis of the prostate dataset using shrinkage methods

### Adding needed packages and importing dataset
library(MASS)
library(lars)
library(pls)
A<-load("ps02p3.Rdata")

### Separating data into test and training
test<-as.matrix(prostate.test)
train<-as.matrix(prostate.train)

### computing the effective degrees of freedom for ridge regression
ridge.df<-function(train,lambda){
  I<-diag(ncol(train))
  df = diag(train%*%ginv(t(train)%*%train + lambda * I)%*%t(train))
  return(sum(df))
}

### Computing ridge regression models for different values of lambda
lambda<-c(seq(0,1000,1))
model<-lm.ridge(lpsa~.,prostate.train,lambda = lambda)
coeframe<-data.frame(model$coef)
dof<-function(lambda){
  degree = c()
  for(i in 1:NROW(lambda)){
    degree = c(degree,ridge.df(train,lambda[i]))
  }
  return(degree)
}
dgf<-dof(lambda)
dgframe<-data.frame(dgf)
matplot(dgframe,t(coeframe),xlab = 'df(lambda)', ylab = 'coefficients', main = "fitted coefficients vs dof")
abline(h=0,col=9,lty=3)
abline(v=5,col=2,lty=3)

### ridge regression and degree of freedom for train and test data 
ridgetrain<-model$coef
ridgetest<-lm.ridge(lpsa~.,prostate.test,lambda = lambda)$coef

doftrain<-c()
for(i in 1:NROW(lambda))
{
  doftrain[i]<-ridge.df(train[,1:8],lambda = i)
}

doftest<-c()
for(i in 1:NROW(lambda))
{
  doftest[i]<-ridge.df(test[,1:8],lambda = i)
}

#### computing the train and test error for ridge regression 

rsstrain<-c()
for(i in 0:1001){
  tra<-lm.ridge(lpsa~.,prostate.train,lambda = i)$coef
  rsstrain[i]<-t(train[,9]-train[,1:8]%*%tra)%*%(train[,9]-train[,1:8]%*%tra) + i * (t(tra)%*%tra)
}

dtrain<-as.matrix(doftrain)
rtrain<-as.matrix(rsstrain)
matplot(dtrain,rtrain,main = "train vs degrees of freedom")

rsstest<-c()
for(j in 0:1001){
  tes<-lm.ridge(lpsa~.,prostate.test,lambda = i)$coef
  rsstest[j]<-t(test[,9]-test[,1:8]%*%tes)%*%(test[,9]-test[,1:8]%*%tes) + j * (t(tes)%*%tes)
}

dtest<-as.matrix(doftest)
rtest<-as.matrix(rsstest)
matplot(dtest,rtest, main = "test vs degrees of freedom")

####computing lasso regression model
train<-prostate.train[,1:8]
resp<-prostate.train[,9]
trainn<-as.matrix(train)
B<-lars(trainn,resp,type='lasso',trace = FALSE, normalize = TRUE, intercept = TRUE,)
plot.lars(B,main = "lassotrain")

test<-prostate.test[,1:8]
output<-prostate.test[,9]
testt<-as.matrix(test)
C<-lars(testt,output,type='lasso',trace = FALSE, normalize = TRUE, intercept = TRUE)
plot.lars(C, main = "lassotest")

testlasso<-predict(C,test)
summary(testlasso)
####Ploting of test data
testlasso$fraction
C$RSS
plot(testlasso$fraction,C$RSS, main = "test vs shrinkage factor")
####Ploting of train data
trainnlasso<-predict(B,train)
summary(trainnlasso)
plot(trainnlasso$fraction,B$RSS, main = "train vs shrinkage factor")