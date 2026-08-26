library(tidyverse)
library(rsample)
library(FNN)
library(igraph)
library(datasets)
library(maps)
library(usmap)
library(GGally)
library(paletteer)
library(tree)
library(pROC)
library(yardstick)
library(rpart)
library(ggplot2)
library(tidyr)
library(RColorBrewer)
library(MASS)
library(glmnet)
library(plotmo)
datas <- read.csv("DataNoNa.csv")
# modellazione -------------------------------------------------
set.seed(123)
#questo è il dataset da usare per fare previsioni e modellazioni
data_split<-datas[,-c(1:3, 93:100, 102:107)]
split <- initial_split(data_split)
dim(split)

train <- training(split)
test <- testing(split)

# modello lineare classico ------------------------------------------------

lm1 <- lm(ViolentCrimesPerPop~., data = train)
summary(lm1)

pred.lm <- predict(lm1, newdata = test)

# Regressione stepwise backward -------------------------------------------

step_back <- lm1 %>% stats:::step(direction = "backward",)
summary(step_back)
pred.step_back <- predict(step_back, newdata = test)

# Regressione stepwise forward --------------------------------------------

lm0 <- lm(ViolentCrimesPerPop~1, data = train)
full_model <- formula(lm(ViolentCrimesPerPop ~ ., train))
step_forw <- lm0 %>% stats::step(direction = "forward", scope = full_model)
summary(step_forw)
pred.step_forw <- predict(step_forw, newdata = test)

# Regressione stepwise ibrida (both) --------------------------------------

step_both1 <- lm1 %>% stats::step(direction = "both")
step_both0 <- lm0 %>% stats::step(direction = "both", scope = full_model)

summary(step_both1)
summary(step_both0)

pred.step_both1 <- predict(step_both1, newdata = test)
pred.step_both0 <- predict(step_both0, newdata = test)

#confronto AIC 
extractAIC(lm1)[2]; extractAIC(step_back)[2]; extractAIC(step_forw)[2]; extractAIC(step_both1)[2]; extractAIC(step_both0)[2] 

# Ridge e Lasso -----------------------------------------------------------
#data l'esplorativa, si sospetta multicollinearità tra gran parte delle variabili
#perciò si applica Ridge e Lasso
library(glmnet)
ridge1 <- glmnet(y = train %>% pull(ViolentCrimesPerPop),
                 x = train %>% dplyr::select(-ViolentCrimesPerPop) %>% as.matrix(),
                 alpha = 0,
                 lambda.min.ratio = 1e-10,
                 nlambda = 200)
previsione_ridge <- predict(ridge1, newx = test %>% 
                              dplyr::select(-ViolentCrimesPerPop) %>% 
                              as.matrix())
errori_ridge <- (test %>% pull(-ViolentCrimesPerPop) - previsione_ridge)^2 %>% 
  apply(2, mean)
plot(y=errori_ridge, x=log(ridge1$lambda))
#lambda ottimale
lambda_opt_ridge <- errori_ridge %>% which.min()
lambda_opt_ridge 
ridge1$lambda[lambda_opt_ridge]
coef(ridge1, s = ridge1$lambda[lambda_opt_ridge])
pred.ridge <- previsione_ridge[, lambda_opt_ridge]
plot_glmnet(ridge1,xvar = "rlambda")

# Lasso -------------------------------------------------------------------

lasso1 <- glmnet(y = train %>% pull(ViolentCrimesPerPop),
                 x = train %>% dplyr::select(-ViolentCrimesPerPop) %>% as.matrix(),
                 alpha = 1,
                 lambda.min.ratio = 1e-10,
                 nlambda = 200)

previsione_lasso <- predict(lasso1, newx = test %>% 
                              dplyr::select(-ViolentCrimesPerPop) %>% 
                              as.matrix())
errori_lasso <- (test %>% pull(-ViolentCrimesPerPop) - previsione_lasso)^2 %>% 
  apply(2, mean)
plot(y=errori_lasso, x=log(lasso1$lambda))  
#lambda ottimale
lambda_opt_lasso <- errori_lasso %>% which.min()
lambda_opt_lasso 

lasso1$lambda[lambda_opt_lasso]
round(coef(lasso1, s = lasso1$lambda[lambda_opt_lasso]),4)
pred.lasso <- previsione_lasso[, lambda_opt_lasso]
plotmo::plot_glmnet(ridge1,xvar = "rlambda")

# Elastic Net -------------------------------------------------------------

elas1 <- glmnet(y = train %>% pull(ViolentCrimesPerPop),
                x = train %>% dplyr::select(-ViolentCrimesPerPop) %>% as.matrix(),
                alpha = 0.5,
                lambda.min.ratio = 1e-10,
                nlambda = 200)
plotmo::plot_glmnet(elas1,xvar = "rlambda")

# Cross Validation --------------------------------------------------------
y <- train %>% pull(ViolentCrimesPerPop) 
x <- train %>% dplyr::select(-ViolentCrimesPerPop) %>% as.matrix() 
xtest <- test %>% dplyr::select(-ViolentCrimesPerPop) %>% as.matrix() 
set.seed(123)
elas_cv <- cv.glmnet(y = y,
                     x = x,
                     alpha = 0.95,
                     nfolds = 50,
                     lambda.min.ratio = 1e-10,
                     nlambda = 100)
plot(elas_cv, sign.lambda = 1)

elas_cv$lambda.1se
elas_cv$lambda.min

p.min <- predict(elas_cv, newx = xtest, s = "lambda.min")
p.1se <- predict(elas_cv, newx = xtest, s = "lambda.1se")

# Scores della PCA come esplicative per lm ------------------------------------------------------------
pca <- prcomp(data_split,scale.=T)
summary(pca)
#teniamo 14 componenti per previsione
#impossibile fare interpretazione
scores <- data.frame(pca$x[,1:14],ViolentCrimesPerPop=datas$ViolentCrimesPerPop)
set.seed(123)
scores_split <- initial_split(scores)
scores_train <- training(scores_split)
scores_test <- testing(scores_split)
lm_pca <- lm(ViolentCrimesPerPop~.,data=scores_train)
summary(lm_pca)
pred.lm_pca <- predict(lm_pca, newdata = scores_test)

#confronto tra modelli

err <- (cbind(pred.lm,
              pred.step_back,
              pred.step_forw,
              pred.step_both1,
              pred.step_both0,
              pred.ridge,
              pred.lasso,
              p.min,
              p.1se,
              pred.lm_pca
              
) - test$ViolentCrimesPerPop)^2 %>% colMeans()
sqrt(err)

tibble(
  model = c(
    "Minimi Quadrati",
    "Stepwise Backward",
    "Stepwise Forward",
    "Stepwise Both (lm1)",
    "Stepwise Both (lm0)",
    "Ridge",
    "Lasso",
    "Elastic Net (lambda minimo)",
    "Elastic Net (lambda 1se)",
    "Modello lineare con PCA"
  ),
  MSE = err) %>% arrange(MSE) %>% mutate(MSE = round(MSE, 3)) %>%
  knitr::kable()


# Interpretazione Elastic Net lambda 1se -------------------------------------------

#come detto, non verranno interpretati i coefficienti del modello lm con PCA
#si sceglie dunque di interpretare il secondo migliore

round(coef(elas_cv, s = elas_cv$lambda.1se),4)

# CLASSIFICAZIONE: LDA E QDA-----------------------------------------------------
# usiamo modelli generativi
# Dividiamo la var. ViolentCrimesPerPop in var. binaria:
train2 <- train
test2 <- test
train = train2 %>% mutate(ViolentCrimesPerPop = ifelse(ViolentCrimesPerPop<250,0,ifelse(ViolentCrimesPerPop>1000,2,1)))
test = test2 %>% mutate(ViolentCrimesPerPop = ifelse(ViolentCrimesPerPop<250,0,ifelse(ViolentCrimesPerPop>1000,2,1)))
# LDA---------
# modello con train:
lmod = lda(ViolentCrimesPerPop~., data = train)
# prediction con test:
p.lda <- predict(lmod, newdata = test, type='response')
# posterior: matrice di probabilità a posteriori restituita da
# predict() su un modello LDA, contiene le probabilità stimate che ogni
# osservazione appartenga a ciascuna classe (0, 1)
# In questo caso estraggo solo la seconda colonna, cioè la probabilità di
# appartenere alla classe 1.
# QDA--------------------------
qmod = qda(ViolentCrimesPerPop~., data = train)
p.qda = predict(qmod, newdata = test, type='response')
# CART--------------------------
tree <- rpart(ViolentCrimesPerPop~., data = train,method="class",control=rpart.control(cp=0.005))
p.tree <- predict(tree,newdata=test)
p.class <- predict(tree,newdata=test,type='class')
# CONFRONTO ---------------
# Valutiamo tutti i modelli adattati
table(test$ViolentCrimesPerPop, p.lda$class, dnn=c("Reale", "Predetto"))
table(test$ViolentCrimesPerPop, p.qda$class, dnn=c("Reale", "Predetto"))
table(test$ViolentCrimesPerPop, p.class, dnn=c("Reale", "Predetto"))
# CURVA ROC 1: ViolentCrimesPerPop è stato diviso in 2 classi: Basso e Alto.
roc_lda = multiclass.roc(test$ViolentCrimesPerPop, p.lda$posterior)
roc_qda = multiclass.roc(test$ViolentCrimesPerPop, p.qda$posterior)
roc_tree = multiclass.roc(test$ViolentCrimesPerPop, p.tree)
aucs <- sapply(roc_lda$rocs, function(r) auc(r[[1]]))
aucs <- rbind(aucs,sapply(roc_qda$rocs, function(r) auc(r[[1]])))
aucs <- rbind(aucs,sapply(roc_tree$rocs, function(r) auc(r[[1]])))
aucs
par(mfrow = c(1, 3),pty="s")
coppie <- sapply(roc_lda$rocs, function(r) paste(r[[1]]$levels, collapse = " vs "))
for(i in 1:length(roc_lda$rocs)) {
  plot.roc(roc_lda$rocs[[i]][[1]], col = "#9cdda5", main = coppie[i],legacy.axes=T)
  plot.roc(roc_qda$rocs[[i]][[1]], col = "#b3caff", add = TRUE,legacy.axes=T)
  plot.roc(roc_tree$rocs[[i]][[1]], col = "#f69593", add = TRUE,legacy.axes=T)
  legend("bottomright", legend = c("LDA", "QDA","CART"), col = c("#9cdda5", "#b3caff","#f69593"), lwd = 2, cex = 0.7)
}
table(train$ViolentCrimesPerPop)
table(test$ViolentCrimesPerPop)
#fortemente sbilanciata
dev.off()
par(mfrow=c(1,1))
# CLASSIFICAZIONE: LDA, QDA e CART binaria -----------------------------------------------------
#Decidiamo di fare una classificazione binaria, visti i precedenti risultati e
#una ricerca.
train = train2 %>% mutate(ViolentCrimesPerPop = ifelse(ViolentCrimesPerPop<500,0,1))
test = test2 %>% mutate(ViolentCrimesPerPop = ifelse(ViolentCrimesPerPop<500,0,1))
# LDA-------------------------
# modello con train:
lmod = lda(ViolentCrimesPerPop~., data = train)
# prediction con test:
p.lda <- predict(lmod, newdata = test, type='response')
# QDA--------------------------
qmod = qda(ViolentCrimesPerPop~., data = train)
p.qda = predict(qmod, newdata = test, type='response')
# CART--------------------------
tree <- rpart(ViolentCrimesPerPop~., data = train,method="class",control=rpart.control(cp=0.005))
p.tree <- predict(tree,newdata=test)
p.class <- predict(tree,newdata=test,type='class')
# CONFRONTO---------------
# Valutiamo tutti i modelli adattati
table(test$ViolentCrimesPerPop, p.lda$class, dnn=c("Reale", "Predetto"))
table(test$ViolentCrimesPerPop, p.qda$class, dnn=c("Reale", "Predetto"))
table(test$ViolentCrimesPerPop, p.class, dnn=c("Reale", "Predetto"))
Q1 = data.frame(y = factor(test$ViolentCrimesPerPop, levels = c(0,1)),
                pred = p.lda$posterior[,1], method = 'LDA')
Q2 = data.frame(y = factor(test$ViolentCrimesPerPop, levels = c(0,1)),
                pred = p.qda$posterior[,1], method = 'QDA')
Q3 = data.frame(y = factor(test$ViolentCrimesPerPop, levels = c(0,1)),
                pred = p.tree[,1], method = 'CART')
QQ = bind_rows(Q1, Q2,Q3) %>% group_by(method)
autoplot(roc_curve(QQ, truth = y,pred,event_level = "first"))
roc_auc(QQ, truth = y, pred,event_level = 'first')