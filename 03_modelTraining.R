####rough of machine learning:####
library(tidymodels)
library(themis)
library(skimr)
library(tidyverse)
#library(parsnip) don't thin I actually need this
library(randomForest)

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
processedData <- lapply(X=inputForML, FUN = splitData_fctn, proportion =0.8)

#get training data
trainData <- lapply(processedData, function(w){
  w$train
})

#get testing data
testData <- lapply(processedData, function(w){
  w$test
})

#get the recipe
theRecipe <- lapply(processedData, function(w){
  w$recipe
})

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
                      mtry = tune(),#min_n = tune(),
                      mode = "classification") %>%
  set_engine("randomForest", importance = TRUE)

rf_tune_grid <- grid_regular(
  mtry(range = c(1,6)),
  #min_n(range = c(1,2)),
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


library(vip)

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

