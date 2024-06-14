####rough of machine learning:####
library(tidymodels)
library(themis)
library(skimr)
library(tidyverse)
#library(parsnip) don't thin I actually need this cause getsloaded in with tidymodels
library(randomForest)
library(vip)

nmfMatrices <- readRDS("processed/nmf_wMatrices.rds")
geneExpress <- readRDS("processed/geneExpression.rds")
#inputForML <- basis_matrices
#View(Wmatrices[[1]])

inputForML <- list()
for (each in nmfMatrices) {
  #  dataframe <- data.frame(unlist(each))
  each$significant <- (geneExpress$significant)
  each$significant <- as.factor(each$significant)
  inputForML <- append(inputForML, list(each))
}


####real thing####
#function for preprocessing
split_processingData_fctn <- function(data, proportion){
  #set.seed(234)
  dataSplit <- initial_split(data, prop = proportion)
  trainData <- training(dataSplit)
  testData <- testing(dataSplit)
  #create recipe
  rec <-recipe(significant~., data = trainData) %>%
    step_downsample(significant, under_ratio = 1)
  #return(list(rec, trainData, testData))
  return(list("train" = as.data.frame(trainData), "test"=as.data.frame(testData), "recipe"=rec))
}

#apply split data function to all matrices
processedData <- lapply(X=inputForML, FUN = split_processingData_fctn, proportion =0.8)

#get training data
# trainData <- lapply(processedData, function(w){
#   w$train
# })
# 
#get testing data
testData <- lapply(processedData, function(w){
  w$test
})
# 
# #get the recipe
# theRecipe <- lapply(processedData, function(w){
#   w$recipe
# })

# list_fold <- list()

randForest_fctn <- function(mylist){
  model_RF <- rand_forest(trees = 500, mtry = tune(), min_n = tune(), 
                          mode = "classification") %>% set_engine("randomForest", importance = TRUE)
  ffolds <- vfold_cv(data = mylist$train, v=5)
  limit <- (ncol(mylist$train))-1
  tuningGrid <- grid_regular(
    mtry(range = c(1,limit)),
    min_n(range = c(1,limit)),
    levels = limit
  )
  #build workflow
  wkflow <- workflow() %>%
    add_recipe(mylist$recipe) %>%
    add_model(model_RF)
  #tune model
  res <- tune_grid(
    wkflow,
    resamples = ffolds,
    grid = tuningGrid,
    control = control_grid(save_pred = TRUE),
    metrics = metric_set(roc_auc)
  )
  paramsPlot <- autoplot(res)
  # paramsPlot2 <- res %>%
  #   collect_metrics() %>%
  #   ggplot(aes(x=mtry, y=mean))+
  #   geom_point()+
  #   geom_line()
  
  final_model <- res %>% select_best(metric = "roc_auc")
  
  final_fit <- finalize_workflow(wkflow, final_model) %>%
    fit(data = mylist$train) 
  
  importancePlot <- final_fit %>% extract_fit_parsnip()%>%
    vip(geom='col',aes = list(colour="black", fill='lightblue', alpha=0.7))+
    theme_classic()
  
  aug <- augment(final_fit, mylist$test)
  roc_auc <- roc_auc(aug, significant, .pred_TRUE)
  two_classCurve <- roc_curve(aug, truth = significant,
                              .pred_TRUE)
  rocCurve <- autoplot(two_classCurve)
  
  return(list("tuningPlots" =paramsPlot, "importance"=importancePlot,
              "finalFit" = final_fit, "AUC"= roc_auc, "roc_curve" = rocCurve))
}

# letsTEsthardfctn3 <- lapply(processedData, FUN = hardFunctin2)
mlResList <- lapply(randForest_fctn, FUN = hardFunctin2)

# res_tune <- list()
# for (i in 1:length(processedData)) {
#   current <- processedData[[i]]
#   folds <- vfold_cv(data = current$train, v=5)
#   #list_fold <- append(list_fold, folds)
#   model_RF <- rand_forest(trees = 500, mtry = tune(), min_n = tune(), 
#                           mode = "classification") %>% set_engine("randomForest", importance = TRUE)
#   
#   upper <- ncol(current$train)-1
#   tune_grid_rd <- grid_regular(
#     mtry(range = c(1,upper)),
#     min_n(range = c(1,upper)),
#     levels = upper
#   )
#   theWorkflw <- workflow()%>%
#     add_recipe(current$recipe)%>%
#     add_model(model_RF)
#   
#   res_tun <- tune_grid(theWorkflw,
#                        resamples = folds,
#                        grid= tune_grid_rd,
#                        control=control_grid(save_pred = T),
#                        metrics= metric_set(roc_auc))
#   res_tune[[i]] <- append(res_tune, res_tun)
#   pPlot <- res_tun %>%
#     collect_metrics() %>%
#     ggplot(aes(x=mtry, y=mean))+
#     geom_point()+
#     geom_line()
#   pPlot
# }

# ctrl+shift+C
# preprocessData_fctn <- function(Trainingdata){
#   rec <- recipe(significant ~., data = Trainingdata) %>%
#     step_downsample(significant, under_ratio = 1)
#   return("recipe"=rec)
# }
# 
# processingRec <- lapply(X=inputForML, FUN=preprocessData_fctn)
# #processingRec$


#### some testing on just one dataset####



str(inputForML[[2]])
temp <- inputForML[[5]]
#temp <- temp %>% select(1,4,5,6,7,10)

#convert outcome to a factor
#temp$significant <- as.factor(temp$significant)
str(temp)
skim(temp)
#set seed when creating split
dataSplit <- initial_split(temp, prop = 0.8)
tempTrain <- training(dataSplit)
tempTest <- testing(dataSplit)

#look at before sampling
ggplot(tempTrain, aes(factor(significant)))+
  geom_bar(aes(y=after_stat(count)/sum(after_stat(count))), colour = "black", fill = "lightgrey")+
  scale_y_continuous(labels = percent)+
  xlab("") + ylab("% of genes") +
  theme_minimal(base_size = 20)


temp_rec <- recipe(significant ~ ., data = tempTrain) %>%
  step_downsample(significant, under_ratio =1)
  #step_select_vip()
temp_rec %>% prep() %>% bake(NULL)



#how does undersampling look after
temp_rec %>% 
  prep() %>% 
  juice() %>%
  ggplot(aes(factor(significant))) +
  geom_bar(aes(y = (after_stat(count))/sum(after_stat(count))), colour="black",fill="lightgrey") +
  scale_y_continuous(labels = percent) +
  xlab("") + ylab("% of genes") +
  theme_minimal(base_size = 20)


# set seed for CV:
?yardstick::metric_set
cv_folds <- vfold_cv(data = tempTrain, v = 5)

#model specification: RF

show_engines("rand_forest")
rf_mod <- rand_forest(trees = 500, 
                      mtry = tune(),min_n = tune(),
                      mode = "classification") %>%
  set_engine("randomForest", importance = TRUE)

rf_tune_grid <- grid_regular(
  mtry(range = c(1,6)),
  min_n(range = c(1,2)),
  levels = 6
)
##build workflow:
temp_wf <- workflow() %>%
  add_recipe(temp_rec) %>%
  add_model(rf_mod)

#doParallel::registerDoParallel()

#tune model:
tune_rf <- tune_grid(
  temp_wf,
  resamples = cv_folds,
  grid = rf_tune_grid,
  control = control_grid(save_pred = TRUE),
  metrics = metric_set(roc_auc)
) #had to install package randomForest for this to work
##now fit model on 
#can look at the roc  for each combination of parameters
autoplot(tune_rf)
#the exact same plot above:
plot <- tune_rf%>%
  collect_metrics()%>%
  ggplot(aes(x=mtry, y=mean))+
  geom_point()+
  geom_line()
plot


#best <- tune_rf %>% collect_metrics() %>%
 # arrange(.metric)
#collect_metrics(tune_rf)
#select_best(tune_rf, metric = "roc_auc")
#############################fit_mod <- best_fit(temp_wf, split = dataSplit, metric_set(roc_auc))
best_mod <- tune_rf %>% select_best(metric = "roc_auc")

#fit the final model with best params on training
final_fit <- finalize_workflow(temp_wf, best_mod) %>%
  fit(data = tempTrain)



#final_fit%>% extract_fit_parsnip()%>%tidy()
#library(broom)
#ff <- finalize_workflow(temp_wf, best_mod)
ffff <- best_mod %>% last_fit(split = dataSplit)



final_fit %>% extract_fit_parsnip() %>%
  vip(geom = 'col', aesthetics = list(colour="black", fill ="lightblue", alpha=0.7)) +
  theme_classic()


#get the test performance
temp_aug <- augment(final_fit, tempTest)
head(temp_aug)
roc_auc(temp_aug, significant, .pred_TRUE)
conf_mat(temp_aug, significant, .pred_class)

#plot ROC curve:
two_classCurve <- roc_curve(temp_aug, truth = significant, .pred_TRUE)
autoplot(two_classCurve)
roc_auc(temp_aug, significant, .pred_FALSE)

