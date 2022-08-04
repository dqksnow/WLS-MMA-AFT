library(MASS)
library(mvtnorm)
library(parallel)
library(ggplot2)
library(reshape2)
library(dplyr)
setwd("D:/OneDrive/研二论文/model averaging/我的实证/PBC")

# par
wls <- function(q){
  library(MASS)
  require(lars)
  library(survival)
  library(mvtnorm)
  library(quadprog)
  library(ncvreg)
  # prepare data
  set.seed(100)
  pbc.row <- pbc
  pbc.row.nonNA <- na.omit(pbc.row)
  X.all <- pbc.row.nonNA[,4:20]
  logU.all<-log(pbc.row.nonNA[[2]])
  logU.all <- scale(logU.all)
  delta.all <- pbc.row.nonNA[[3]]
  #delta.all[which(delta.all!=0)]=1
  delta.all[which(delta.all==1)]=0
  delta.all[which(delta.all==2)]=1
  # X.all['alk.phos'] <- scale(X.all['alk.phos'])
  # X.all['bili'] <- scale(X.all['bili'])
  # X.all['chol'] <- scale(X.all['chol'])
  # X.all['copper'] <- scale(X.all['copper'])
  # X.all['platelet'] <- scale(X.all['platelet'])
  # X.all['protime'] <- scale(X.all['protime'])
  # X.all['ast'] <- scale(X.all['ast'])
  # X.all['trig'] <- scale(X.all['trig'])
  #X.all['age'] <- scale(X.all['age'])
  #X.all['albumin'] <- scale(X.all['albumin'])
  #
  X.all['alk.phos'] <- scale(log(X.all['alk.phos']))
  X.all['bili'] <- scale(log(X.all['bili']))
  X.all['chol'] <- scale(log(X.all['chol']))
  X.all['copper'] <- scale(log(X.all['copper']))
  X.all['platelet'] <- scale(log(X.all['platelet']))
  X.all['protime'] <- scale(log(X.all['protime']))
  X.all['ast'] <- scale(log(X.all['ast']))
  X.all['trig'] <- scale(log(X.all['trig']))
  
  X.all['age'] <- scale(X.all['age'])
  X.all['albumin'] <- scale(X.all['albumin'])
  #
  # X.all['bili'] <- log(X.all['bili'])
  # X.all['protime'] <- log(X.all['protime'])
  # X.all['age'] <- log(X.all['age'])
  # X.all['albumin'] <- log(X.all['albumin'])
  # X.all <- X.all[,c('bili', 'protime', 'edema', 'age', 'albumin')]
  #X.all <- scale(X.all) 
  
  index.f <- which(X.all$sex=='f')
  X.all$gender[index.f]=1
  X.all$gender[-index.f]=0
  X.all <- subset(X.all, select = -c(sex))
  #X.all <- scale(X.all) 
  
  # use LARS to order the variables
  N=length(logU.all)
  n=sum(delta.all)
  M=ncol(X.all)
  v <- c(1:M)
  ord.all <- order(logU.all)
  order.X.all <- X.all[ord.all,]
  order.logU.all <- logU.all[ord.all]
  order.delta.all <- delta.all[ord.all]
  w <- rep(0,N)
  w[1]<-order.delta.all[1]/N
  for(i in 2:N){
    temp1=1
    for(j in 1:(i-1)){temp1 <- temp1*((N-j)/(N-j+1))^order.delta.all[j]}
    w[i] <- (order.delta.all[i]/(N-i+1)) * temp1
  }
  order.X.bar <- c(apply(w*order.X.all,2,sum))
  order.logU.bar <- w%*%order.logU.all
  order.X.w <- sqrt(w)*(order.X.all-order.X.bar)
  order.logU.w <- sqrt(w)*(order.logU.all-rep(order.logU.bar,N))
  fit.lars <- lars(as.matrix(order.X.w), as.vector(order.logU.w), trace=TRUE,  type='lar', intercept=FALSE)
  fit.lars$actions
  candidate<-c()
  for (i in 1:M){candidate<-c(candidate, unlist(fit.lars$actions[i]))}
  candidate <- as.vector(candidate)
  
  # divide data into train data and test data
  set.seed(q+100)
  index=sort(sample(length(logU.all),length(logU.all)*0.7))
  logU.train<-logU.all[index]
  delta.train<-delta.all[index]
  X.train<-as.matrix(X.all[index,])
  delta.train[which.max(logU.train)]=1
  logU.test<-logU.all[-index]
  delta.test<-delta.all[-index]
  X.test<-as.matrix(X.all[-index,])
  
  # model averaging 
  N=length(logU.train)
  n=sum(delta.train)
  M=ncol(X.train)
  H.wls.uncensored <- matrix(0,n,M)
  H.ksv <-  matrix(0,N,M)
  beta.wls.uncensored.all <- beta.ksv.all <- matrix(0,M,M)
  SAIC <- SBIC <- sigma.sq.m.hat.He <- rep(0,M)
  for(m in 1:M){
    candidate.temp <- sort(candidate[1:m])
    X<-X.train[,candidate.temp]
    ord <- order(logU.train) #order返回向量从小到大排序的index
    order.X.train <- matrix(0,N,1:m)
    if(is.vector(X)) order.X.train[,1] <- X[ord] else order.X.train <- X[ord,]
    order.X.train <- as.matrix(order.X.train)
    order.logU.train <- logU.train[ord]
    order.delta.train <- delta.train[ord]
    # calculate KM estimator
    temp <- rep(1,N)
    w <- G <- n.temp <- d <-  rep(0,N) #G是logC的分布，1-G是logC的生存函数
    for(i in 1:N){
      d[i]=ifelse(order.delta.train[i]==0,1,0)
      n.temp[i]=ifelse(order.delta.train[i]==0,N-i+1,1) #delta=1时就无所谓取几了。
      temp[i] <- 1-d[i]/n.temp[i]
      G[i] <- 1-prod(temp[1:i])
    }
    #只用未删失的数据
    if(is.vector(order.X.train)) order.X.train.uncensored <- order.X.train[order.delta.train!=0] else order.X.train.uncensored <- order.X.train[order.delta.train!=0,]
    order.logT.train.uncensored <- order.logU.train[order.delta.train!=0] #未删失的logU，就是logT
    # ksv
    order.logT.train.ksv <- order.logU.train*order.delta.train/(1-G)
    beta.ksv <- solve(t(order.X.train)%*%order.X.train) %*% t(order.X.train)%*%order.logT.train.ksv
    beta.ksv.all[candidate.temp, m] <- beta.ksv
    miu.m.ksv <- order.X.train%*%beta.ksv
    H.ksv[,m] <- miu.m.ksv-order.logT.train.ksv
    # wls
    W.uncensored <- diag( order.delta.train[order.delta.train!=0]/(1-G[order.delta.train!=0]) ) #n*n
    W <- order.delta.train/(1-G) #N*N
    beta.wls.uncensored <- solve(t(order.X.train.uncensored)%*%W.uncensored%*%order.X.train.uncensored) %*% t(order.X.train.uncensored)%*%W.uncensored%*%order.logT.train.uncensored
    beta.wls.uncensored.all[candidate.temp, m] <- beta.wls.uncensored
    miu.m.wls.uncensored <- order.X.train.uncensored%*%beta.wls.uncensored
    temp <- order.logT.train.uncensored - miu.m.wls.uncensored #W.uncensored %*%
    sigma.sq.m.hat.He[m] <- sum(W*(order.logU.train - order.X.train%*%beta.wls.uncensored)^2)/N #(He, 2003)提的方差的估计
    SAIC[m] <- log(sigma.sq.m.hat.He[m]) +2*m/N
    SBIC[m] <- log(sigma.sq.m.hat.He[m]) +log(N)*m/N
    H.wls.uncensored[,m] <-  temp #sqrt(W.uncensored) %*%temp
  }
  Dmat.wls<- 2*t(H.wls.uncensored) %*% H.wls.uncensored + diag(rep(1e-4,length(candidate)))
  
  #fai.n <- log(n)/n
  #dvec <- c(-fai.n * t(v.s))
  A.1 <- rep(1,length(candidate))
  A.2 <- rep(-1,length(candidate))
  A.3 <- diag(x=1, nrow=length(candidate))
  A.4 <- diag(x=-1, nrow=length(candidate))
  Amat <-t(rbind(A.1, A.2, A.3, A.4))
  b.1 <- 1
  b.2 <- -1
  b.3 <- rep(0,length(candidate))
  b.4 <- rep(-1,length(candidate))
  bvec <- c(b.1, b.2, b.3, b.4)
  
  #sigma.sq.hat.wls <- sigma.sq.m.hat.He
  #sigma.sq.hat.Hansen <- as.vector(H.wls.uncensored[,M]%*%H.wls.uncensored[,M]/(N-M))
  # He's estimator of sigma^2
  sigma.sq.hat.He <- sigma.sq.m.hat.He[length(candidate)]
  dvec = -2 * sigma.sq.hat.He *  v 
  result.wls <- solve.QP(Dmat.wls,dvec,Amat,bvec=bvec)
  #sc <- norm(Dmat.wls,"2")
  #result.wls <- solve.QP(Dmat = Dmat.wls/sc, dvec=dvec/sc, Amat=Amat, bvec=bvec, meq=0, factorized=FALSE )
  weight.wls <- result.wls$solution
  weight.wls[abs(weight.wls) < 1e-5] <- 0;
  beta.wls.MA <-  c(apply(weight.wls*t(beta.wls.uncensored.all),2,sum))
  #Hansen's estimator of sigma^2
  sigma.sq.hat.Hansen <- as.vector(sum((order.logT.train.uncensored - miu.m.wls.uncensored)^2)/(sum(delta.train)-M))
  dvec = -2 * sigma.sq.hat.Hansen *  v
  result.wls.Hansen <- solve.QP(Dmat.wls,dvec,Amat,bvec=bvec)
  #result.wls.Hansen <- solve.QP(Dmat = Dmat.wls/sc, dvec=dvec/sc, Amat=Amat, bvec=bvec, meq=0, factorized=FALSE )
  weight.wls.Hansen <- result.wls.Hansen$solution
  weight.wls.Hansen[abs(weight.wls.Hansen) < 1e-5] <- 0
  beta.wls.MA.Hansen <-  c(apply(weight.wls.Hansen*t(beta.wls.uncensored.all),2,sum))
  #ksv
  Dmat.ksv<- 2*t(H.ksv) %*% H.ksv
  sigma.sq.hat.ksv <- as.vector(H.ksv[,M]%*%H.ksv[,M]/(N-M))
  dvec = -2 *sigma.sq.hat.ksv*  v 
  result.ksv <- solve.QP(Dmat.wls,dvec,Amat,bvec=bvec)
  weight.MMA.ksv <- result.wls$solution
  weight.MMA.ksv[abs(weight.MMA.ksv) < 1e-5] <- 0; weight.MMA.ksv
  beta.ksv.MA <-  c(apply(weight.MMA.ksv*t(beta.ksv.all),2,sum))
  # AIC and BIC
  beta.wls.AIC <- beta.wls.uncensored.all[,which.min(SAIC)]
  beta.wls.BIC <- beta.wls.uncensored.all[,which.min(SBIC)]
  # SAIC and SBIC
  weight.SAIC <- exp(-SAIC/2)/sum( exp(-SAIC/2))
  beta.wls.SAIC <-  c(apply(weight.SAIC*t(beta.wls.uncensored.all),2,sum))
  weight.SBIC <-  exp(-SBIC/2)/sum( exp(-SBIC/2))
  beta.wls.SBIC <-  c(apply(weight.SBIC*t(beta.wls.uncensored.all),2,sum))
  # penalized methods (weighted loss function plus penalty function)
  w <- rep(0,N)
  w[1]<-order.delta.train[1]/N
  for(i in 2:N){
    temp1=1
    for(j in 1:(i-1)){temp1 <- temp1*((N-j)/(N-j+1))^order.delta.train[j]}
    w[i] <- (order.delta.train[i]/(N-i+1)) * temp1
  }
  order.X.train.bar <- c(apply(w*order.X.train,2,sum))
  order.logU.train.bar <- w%*%order.logU.train
  order.X.train.w <- sqrt(w)*(order.X.train-order.X.train.bar)
  order.logU.train.w <- sqrt(w)*(order.logU.train-rep(order.logU.train.bar,N))
  
  fit.lasso <- cv.ncvreg(order.X.train.w, order.logU.train.w, penalty="lasso")
  beta.lasso.ncvreg <- as.vector(coef(fit.lasso))#[-1]
  fit.scad <- cv.ncvreg(order.X.train.w, order.logU.train.w, penalty="SCAD")
  beta.scad <- as.vector(coef(fit.scad))#[-1]
  #
  # model <- cv.glmnet(order.X.train.w, order.logU.train.w, family="gaussian",type.measure="default", intercept=F)
  # model.coef <- coef(model$glmnet.fit,s=model$lambda.1se)
  # beta.lasso.glmnet <- as.vector(model.coef)[-1]
  # measure
  # miu.MSE.wls    <- matrix( delta*(X%*%beta.wls.uncensored-logU.test),1,N) %*% matrix(X%*%beta.wls.uncensored-logU.test,N,1) / sum(order.delta)
  # miu.MSE.wls.MA <- matrix( delta*(X%*%beta.wls.MA-logU.test),1,N) %*% matrix(X%*%beta.wls.MA-logU.test,N,1) / sum(order.delta)
  # miu.MSE.ols    <- matrix( delta*(X%*%beta.ols.uncensored-logU.test),1,N) %*% matrix(X%*%beta.ols.uncensored-logU.test,N,1) / sum(order.delta)
  # miu.MSE.lasso  <- matrix( delta*(X%*%beta.lasso-logU.test),1,N)    %*% matrix(X%*%beta.lasso -logU.test,N,1) / sum(order.delta)
  logU.test.uncensored <- logU.test[delta.test!=0]
  X.test.uncensored <- as.matrix(X.test[delta.test!=0,])
  n.test.uncensored <- length(logU.test.uncensored)
  # MSE.ksv       <- matrix( (X.test.uncensored%*%beta.ksv.MA-logU.test.uncensored),1,n.test.uncensored) %*% matrix(X.test.uncensored%*%beta.ksv.MA-logU.test.uncensored,n.test.uncensored,1) / n.test.uncensored
  # MSE.wls       <- matrix( (X.test.uncensored%*%beta.wls.uncensored-logU.test.uncensored),1,n.test.uncensored) %*% matrix(X.test.uncensored%*%beta.wls.uncensored-logU.test.uncensored,n.test.uncensored,1) / n.test.uncensored
  # MSE.wls.MA    <- matrix( (X.test.uncensored%*%beta.wls.MA-logU.test.uncensored),1,n.test.uncensored)   %*% matrix(X.test.uncensored%*%beta.wls.MA-logU.test.uncensored,n.test.uncensored,1) / n.test.uncensored
  # MSE.wls.MA.Hansen    <- matrix( (X.test.uncensored%*%beta.wls.MA.Hansen -logU.test.uncensored),1,n.test.uncensored)   %*% matrix(X.test.uncensored%*%beta.wls.MA.Hansen -logU.test.uncensored,n.test.uncensored,1) / n.test.uncensored
  # MSE.wls.AIC  <- matrix( (X.test.uncensored%*%beta.wls.AIC-logU.test.uncensored),1,n.test.uncensored) %*% matrix(X.test.uncensored%*%beta.wls.AIC-logU.test.uncensored,n.test.uncensored,1) / n.test.uncensored
  # MSE.wls.BIC  <- matrix( (X.test.uncensored%*%beta.wls.BIC-logU.test.uncensored),1,n.test.uncensored) %*% matrix(X.test.uncensored%*%beta.wls.BIC-logU.test.uncensored,n.test.uncensored,1) / n.test.uncensored
  # MSE.wls.SAIC <- matrix( (X.test.uncensored%*%beta.wls.SAIC-logU.test.uncensored),1,n.test.uncensored) %*% matrix(X.test.uncensored%*%beta.wls.SAIC-logU.test.uncensored,n.test.uncensored,1) / n.test.uncensored
  # MSE.wls.SBIC <- matrix( (X.test.uncensored%*%beta.wls.SBIC-logU.test.uncensored),1,n.test.uncensored) %*% matrix(X.test.uncensored%*%beta.wls.SBIC-logU.test.uncensored,n.test.uncensored,1) / n.test.uncensored
  
  MSE.wls       <- matrix( ((X.test.uncensored%*%beta.wls.uncensored)-(logU.test.uncensored)), 1, n.test.uncensored) %*% matrix((X.test.uncensored%*%beta.wls.uncensored)-(logU.test.uncensored), n.test.uncensored,1) / n.test.uncensored
  MSE.ksv.MA       <- matrix( ((X.test.uncensored%*%beta.ksv.MA)-(logU.test.uncensored)),1,n.test.uncensored) %*% matrix((X.test.uncensored%*%beta.ksv.MA)-(logU.test.uncensored),n.test.uncensored,1) / n.test.uncensored
  MSE.wls.MA.Hansen    <- matrix( ((X.test.uncensored%*%beta.wls.MA.Hansen) -(logU.test.uncensored)),1,n.test.uncensored)   %*% matrix((X.test.uncensored%*%beta.wls.MA.Hansen) -(logU.test.uncensored),n.test.uncensored,1) / n.test.uncensored
  MSE.wls.MA    <- matrix( ((X.test.uncensored%*%beta.wls.MA)-(logU.test.uncensored)),1,n.test.uncensored)   %*% matrix((X.test.uncensored%*%beta.wls.MA)-(logU.test.uncensored),n.test.uncensored,1) / n.test.uncensored
  MSE.wls.AIC <- matrix( ((X.test.uncensored%*%beta.wls.AIC)-(logU.test.uncensored)),1,n.test.uncensored) %*% matrix((X.test.uncensored%*%beta.wls.AIC)-(logU.test.uncensored), n.test.uncensored,1) / n.test.uncensored
  MSE.wls.BIC <- matrix( ((X.test.uncensored%*%beta.wls.BIC)-(logU.test.uncensored)),1,n.test.uncensored) %*% matrix((X.test.uncensored%*%beta.wls.BIC)-(logU.test.uncensored), n.test.uncensored,1) / n.test.uncensored
  MSE.wls.SAIC <- matrix( ((X.test.uncensored%*%beta.wls.SAIC)-(logU.test.uncensored)),1,n.test.uncensored) %*% matrix((X.test.uncensored%*%beta.wls.SAIC)-(logU.test.uncensored), n.test.uncensored,1) / n.test.uncensored
  MSE.wls.SBIC <- matrix( ((X.test.uncensored%*%beta.wls.SBIC)-(logU.test.uncensored)),1,n.test.uncensored) %*% matrix((X.test.uncensored%*%beta.wls.SBIC)-(logU.test.uncensored), n.test.uncensored,1) / n.test.uncensored
  
  intercept<-rep(1,n.test.uncensored)
  X.test.uncensored <- scale(X.test.uncensored)
  X.test.uncensored <- cbind(intercept, X.test.uncensored)
  # MSE.lasso.ncvreg <- matrix( (X.test.uncensored%*%beta.lasso.ncvreg-logU.test.uncensored),1,n.test.uncensored) %*% matrix(X.test.uncensored%*%beta.lasso.ncvreg-logU.test.uncensored,n.test.uncensored,1) / n.test.uncensored
  # MSE.scad   <- matrix( (X.test.uncensored%*%beta.scad-logU.test.uncensored),1,n.test.uncensored)   %*% matrix(X.test.uncensored%*%beta.scad-logU.test.uncensored,n.test.uncensored,1) / n.test.uncensored
  MSE.lasso.ncvreg <- matrix( ((X.test.uncensored%*%beta.lasso.ncvreg)-(logU.test.uncensored)),1,n.test.uncensored) %*% matrix((X.test.uncensored%*%beta.lasso.ncvreg)-(logU.test.uncensored), n.test.uncensored,1) / n.test.uncensored
  MSE.scad   <- sum( ((X.test.uncensored%*%beta.scad)-(logU.test.uncensored)) * ((X.test.uncensored%*%beta.scad)-(logU.test.uncensored)) )/ n.test.uncensored
  
  return(list(
    MSE.wls      ,
    MSE.ksv.MA      ,
    MSE.wls.MA.Hansen    ,
    MSE.wls.MA    ,
    MSE.wls.AIC  ,
    MSE.wls.BIC  ,
    MSE.wls.SAIC ,
    MSE.wls.SBIC ,
    MSE.lasso.ncvreg,
    MSE.scad
  )
  )
}
cl <- makeCluster(4)
replication=100
t1=proc.time()
results <- parLapply(cl, 1:replication , wls)
t2=proc.time()
time.cost=t2-t1
print(paste0('cost time：',time.cost[3][[1]]/60,'mins'))

MSE.wls<-MSE.ksv.MA <- MSE.wls.MA  <- MSE.wls.MA.Hansen<- MSE.wls.AIC <- MSE.wls.BIC <-MSE.wls.SAIC <-MSE.wls.SBIC <-
  MSE.lasso.ncvreg <-MSE.scad  <- c()

for(i in 1:replication){
  MSE.wls[i] = results[i][[1]][[1]]
  MSE.ksv.MA[i] = results[i][[1]][[2]]
  MSE.wls.MA.Hansen[i] = results[i][[1]][[3]]
  MSE.wls.MA[i] = results[i][[1]][[4]]
  MSE.wls.AIC[i] <- results[i][[1]][[5]]
  MSE.wls.BIC[i] <-results[i][[1]][[6]]
  MSE.wls.SAIC[i] <-results[i][[1]][[7]]
  MSE.wls.SBIC[i] <-results[i][[1]][[8]]
  MSE.lasso.ncvreg[i] <-results[i][[1]][[9]]
  MSE.scad[i] <- results[i][[1]][[10]]
}
df.MSE <- data.frame(
  MSE.wls=mean(MSE.wls),
  MSE.ksv.MA=mean(MSE.ksv.MA ),
  MSE.wls.MA.Hansen=mean(MSE.wls.MA.Hansen),
  MSE.wls.MA=mean(MSE.wls.MA),
  MSE.wls.AIC=mean(MSE.wls.AIC),
  MSE.wls.BIC=mean(MSE.wls.BIC),
  MSE.wls.SAIC=mean(MSE.wls.SAIC),
  MSE.wls.SBIC=mean(MSE.wls.SBIC),
  MSE.lasso.ncvreg=mean(MSE.lasso.ncvreg),
  MSE.scad=mean(MSE.scad)
)

# box plots
df.MSE.boxplot <- data.frame(
  WLS=MSE.wls,
  KSVMA=MSE.ksv.MA ,
  WLSMA.1=MSE.wls.MA.Hansen,
  WLSMA.2=MSE.wls.MA,
  AIC=MSE.wls.AIC,
  BIC=MSE.wls.BIC,
  SAIC=MSE.wls.SAIC,
  SBIC=MSE.wls.SBIC,
  LASSO=MSE.lasso.ncvreg,
  SCAD=MSE.scad
)
df.MSE.boxplot<-melt( #dataframe是宽数据，需要转换成长数据后再放入boxplot
  df.MSE.boxplot,                       #
  variable.name="Methods",         #转换后的分类字段名称
  value.name="MSE"             #转换后的度量值名称
)

#save(df.MSE.boxplot, file="E:/onedrive/研二论文/model averaging/我的latex/figures0512/DataForPlot/pbc_nest.Rdata")

MSE.mean <-  c(mean(MSE.wls),mean(MSE.ksv.MA),mean(MSE.wls.MA.Hansen),mean(MSE.wls.MA),mean(MSE.wls.AIC),
               mean(MSE.wls.BIC),mean(MSE.wls.SAIC),mean(MSE.wls.SBIC),mean(MSE.lasso.ncvreg),mean(MSE.scad))

plot1 <- ggplot(df.MSE.boxplot,aes(x=Methods, y=MSE))+geom_boxplot()+ylim(0.25, 2 ) + theme_bw()+
  theme(axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18),
        axis.text = element_text(size = 15), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 15))

df.MSE.boxplot %>% 
  group_by(Methods) %>% 
  summarise(mean_value=mean(MSE)) %>%
  rename("Methods_1"="Methods") %>% 
  cbind(ggplot_build(plot1)$data[[1]]) -> df.MSE.boxplot.1

plot1+
  geom_segment(data=df.MSE.boxplot.1,
               aes(x=xmin,xend=xmax,
                   y=mean_value,
                   yend=mean_value),
               color="red",linetype ='dotted')

ggsave("pbc_nested.png",width = 12, height = 6, dpi = 300)

# optimal rate
MSE.all <- cbind(MSE.wls, MSE.ksv.MA, MSE.wls.MA.Hansen, MSE.wls.MA, MSE.wls.AIC, MSE.wls.BIC, MSE.wls.SAIC, MSE.wls.SBIC, MSE.lasso.ncvreg, MSE.scad)
minMSE <- function(X){which(X== min(X), arr.ind = TRUE)}
minMSE.index <- apply(MSE.all, 1, minMSE )


# show the results
df.MSE
#table(unlist(minMSE.index))
#table(unlist(minMSE.un.index))
