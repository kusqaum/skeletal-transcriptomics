####rough of machine learning:####
library(tidymodels)
library(themis)
library(skimr)
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
splitData_fctn <- function(data, proportion){
  #set.seed(234)
  dataSplit <- initial_split(data, prop = proportion)
  trainData <- training(dataSplit)
  testData <- testing(dataSplit)
  #create recipe
  #rec <-recipe(significant~., data = trainData) %>%
  # step_downsample(significant, under_ratio = 1)
  #return(list(rec, trainData, testData))
  return(list("train" = as.data.frame(trainData), "test"=as.data.frame(testData)))
}

trainTest <- lapply(X=inputForML, FUN = splitData_fctn, proportion =0.8)

trainData <- lapply(trainTest, function(w){
  w$train
})

testData <- lapply(trainTest, function(w){
  w$test
})

preprocessData_fctn <- function(Trainingdata){
  rec <- recipe(significant ~., data = Trainingdata) %>%
    step_downsample(significant, under_ratio = 1)
  return("recipe"=rec)
}

processingRec <- lapply(X=inputForML, FUN=preprocessData_fctn)
#processingRec$


#### some testing on just one dataset####



str(inputForML[[2]])
temp <- inputForML[[2]]


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
  step_downsample(significant, under_ratio =1) %>%
  step_select_vip()
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

cv_folds <- vfold_cv(
  data = tempTrain,
  v = 5
)

#model specification: RF

show_engines("rand_forest")
rf_mod <- rand_forest(trees = 100, 
                      mtry = tune(), min_n = tune(),
                      mode = "classification") %>%
  set_engine("randomForest", importance = TRUE)

rf_tune_grid <- grid_regular(
  mtry(range = c(1,4)),
  levels = 4
)
##workflow:
temp_wf <- workflow() %>%
  add_recipe(temp_rec) %>%
  add_model(rf_mod)

tune_rf <- tune_grid(
  temp_wf,
  resamples = cv_folds,
  grid = 4
) #had to install package randomForest for this to work
##now fit model on 
