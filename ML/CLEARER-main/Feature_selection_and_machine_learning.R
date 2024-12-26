# combine all features from the individual tools before running this code
# or just use the assembled features from /Features directory

# load required packages
library(data.table)
library(dplyr)
library(caret)
library(glmnet)
library(doParallel)
library(PRROC)

setwd("/to/CLEARER_directory/")
# load essential and non-gene labels 
EsInfo <- read.csv("Class_labels/Sc.csv", sep = ",", header = T, stringsAsFactors = F)
head(EsInfo)
table(EsInfo$Essential.CEG)
# generate class labels suitable for R
EsInfo$Essential.CEG <- make.names(EsInfo$Essential.CEG)
table(EsInfo$Essential.CEG)
# load combined features
Data <- fread("Features/Sc_features.csv.gz", stringsAsFactors = F, sep = ",", data.table = F)

rownames(Data) <- Data$genes

# assign class labels
Data$label <- EsInfo$Essential.CEG[match(Data$genes, EsInfo$Gene)]

# randomize Data
set.seed(69)
Data <- Data[sample(1:nrow(Data), size = nrow(Data), replace = F),]

# split Data 
seq <- seq(0,nrow(Data),nrow(Data)/5)
seq <- round(seq,digits = 0)
val1 <- 1:seq[2]
val2 <- (seq[2]+1):(seq[3])
val3 <- (seq[3]+1):(seq[4])
val4 <- (seq[4]+1):(seq[5])
val5 <- (seq[5]+1):(seq[6])

vali1 <- Data[val1,]
vali2 <- Data[val2,]
vali3 <- Data[val3,]
vali4 <- Data[val4,]
vali5 <- Data[val5,]

train1 <- Data[-val1,]
train2 <- Data[-val2,]
train3 <- Data[-val3,]
train4 <- Data[-val4,]
train5 <- Data[-val5,]

train <- list(train1,train2,train3,train4,train5)
vali <- list(vali1,vali2,vali3,vali4,vali5)

# Feature selection
i <- 1
head(train[[i]][1:10,1:10])
for (i in 1:5){
  Dat <- as.matrix(train[[i]][,2:(ncol(train[[i]])-1)])
  labels_ <- as.character(train[[i]][,ncol(train[[i]])])
  cl <- makePSOCKcluster(10)
  registerDoParallel(cl)
  for (j in 10) {  assign(paste("fit", j, sep=""), cv.glmnet(Dat, labels_, type.measure="auc", alpha=j/10,family="binomial",parallel = TRUE))}
  stopCluster(cl)
  # select lasso features 
  Lasso1=as.matrix(coef(fit10,fit10$lambda.min))
  inds<-which(Lasso1!=0)
  VARSELECT2_temp=row.names(Lasso1)[inds]
  Lasso=VARSELECT2_temp[!(VARSELECT2_temp %in% '(Intercept)')]
  Pos <- sapply(as.character(Lasso),function(x) which(x == colnames(train[[i]])))
  train[[i]] <- train[[i]][,c(unlist(Pos),ncol(train[[i]]))]
  # remove correlating features
  correlationMatrix <- cor(train[[i]][,-(ncol(train[[i]]))])
  # find attributes that are highly corrected (ideally >0.75)
  highlyCorrelated <- findCorrelation(correlationMatrix, cutoff=0.7)
  train[[i]] <- train[[i]][,-highlyCorrelated]
  # get the same columns in the validation data 
  Pos <- sapply(colnames(train[[i]]),function(x) which(x == colnames(vali[[i]])))
  vali[[i]] <- vali[[i]][,unlist(Pos)]
}

# machine learning
Control <- trainControl(method = "repeatedcv",
                        number = 5,
                        repeats = 3,
                        verboseIter = FALSE,
                        sampling = "smote",
                        summaryFunction=multiClassSummary,
                        classProbs=T,
                        savePredictions = F)## ML algorithms ##

rf.list = list("rf1"="rf1","rf2"="rf2","rf3"="rf3","rf4"="rf4","rf5"="rf5")
for(i in 1:5) {
  cl <- makePSOCKcluster(30)
  registerDoParallel(cl)
  set.seed(7)
  rf.list[[i]] <- train(x = train[[i]][,1:(ncol(train[[i]])-1)], 
                        y =  factor(train[[i]][,ncol(train[[i]])]),
                        method = "rf",
                        tuneLength = 3,
                        trControl = Control,
                        metric = "Kappa",
                        nthread=30, ntree= 500)
  stopCluster(cl)
}

# performance evaluation test set
CM <- list()
for (i in 1:5){
  pred1 <- predict(rf.list[[i]],vali[[i]][,1:(ncol(vali[[i]])-1)])
  CM[[i]] <- confusionMatrix(pred1,factor(vali[[i]][,ncol(vali[[i]])]),mode = "everything")
}
perf.test <- setNames(vector('list', 5), 1:5)

for (i in 1:5) {
  perf.test[[i]]  <- c(CM[[i]][[3]],CM[[i]][[4]])
}
perf.test <- do.call(rbind, perf.test)

RocAUC <-c()
PrAUC  <- c()

for (i in 1:5){
  pred_prob <-  predict(rf.list[[i]],vali[[i]][,1:(ncol(vali[[i]])-1)],type = "prob")
  pred_prob$label <- vali[[i]][,ncol(vali[[i]])]
  roc<-roc.curve(scores.class0 = pred_prob$Essential[which(pred_prob$label == "Essential")],
                 scores.class1 = pred_prob$Essential[which(pred_prob$label == "NonEssential")],
                 curve = T)
  pr<-pr.curve(scores.class0 = pred_prob$Essential[which(pred_prob$label == "Essential")],
               scores.class1 = pred_prob$Essential[which(pred_prob$label == "NonEssential")],
               curve = T)
  RocAUC[i] <- roc$auc
  PrAUC [i] <- pr$auc.integral
}

metrics <- cbind(perf.test,RocAUC,PrAUC)
mean <- apply(metrics,2,mean, na.rm = T)
SD <- apply(metrics,2,sd, na.rm = T)
metrics <- cbind(t(metrics),mean,SD)
metrics
write.csv(metrics,"test_rf.csv",row.names = T)

# performance evaluation training set

one <- rf.list[[1]][[4]]
two <- rf.list[[2]][[4]]
three <- rf.list[[3]][[4]]
four <- rf.list[[4]][[4]]
five <- rf.list[[5]][[4]]

df <- rbind(one[which(one$Kappa == max(one$Kappa,na.rm = T)),],
            two[which(two$Kappa == max(two$Kappa,na.rm = T)),],
            three[which(three$Kappa == max(three$Kappa,na.rm = T)),],
            four[which(four$Kappa == max(four$Kappa,na.rm = T)),],
            five[which(five$Kappa == max(five$Kappa,na.rm = T)),]
            
)
df
write.csv(df,"train_rf.csv")




