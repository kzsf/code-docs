# Capstone Analysis
## 1. Load required libraries (R)
```
library(table1)
library(survival)  
library(cvTools)   
library(tidyverse)
library(risksetROC)
library(pROC)
library(timeROC)

```

## 2. Multiple Imputation (R)
```

#print num of NA in each var
x <- df %>% summarize_all(funs(sum(is.na(.))))
x
length(df$FTV3)
library(dplyr)
library(tidyr)

result <- df %>% summarize_all(funs(sum(is.na(.)))) %>%
  gather(key = "variable", value = "missing count") %>%
  arrange(variable)

result

#imputation easy way: mean and median
#df$FTV0[is.na[df$FTV0]] <- mean(df$FTV0, na.rm = T)

#visualize NA
colnames(df[,25:28])
md.pattern(df[,25:28])
md.pattern(df[,1:35])
png("/path/NA_plot.png")
md.pattern(df)
#dev.off()

#flag variables with NAs
#columns 15, 16, 25:28
#time to local, time to distant, ftv 0-3. columns 16 and 17 are mostly NAs.
colnames(df)

#subset variables
df_copy <- df[, c(1,7,8,13,17:36)]
colnames(df_copy)

imputed_data_copy <- mice(df_copy, m = 10, maxit = 30, method = "pmm")
imputed_data_plot <- plot(imputed_data_copy)

#save plot
pdf("imputed_data_plot.pdf")
print(imputed_data_plot)
dev.off()

#check values of imputed data
imputed_data_copy$imp$FTV0
str(imputed_data_copy)

#pool imputed dataset
df <- complete(imputed_data_copy)
df %>% summarize_all(funs(sum(is.na(.))))
str(df)
colnames(df)

```

## Make table 1 (R)
```

df <- read.csv("/path/cleaned_imputed_df.csv")
colnames(df)

#change efs time to years
df$efs.time.years <- (df$efs.time)/365.25 

#log transformed ftv3
df$log_FTV3 <- log(df$FTV3+1)

#cube root ftv3
df$cube_FTV3 <- df$FTV3^(1/3)

#combined status
df$combined_status <- ifelse(df$HR == 1 & df$HER2 == 1, "HR+/HER2+",
                          ifelse(df$HR == 1, "HR+/HER2-",
                          ifelse(df$HER2 == 1, "HR-/HER2+", "HR-/HER2-")))

#make factor vars
df$HR <- as.factor(df$HR)
df$HER2 <- as.factor(df$HER2)
df$combined_status <- as.factor(df$combined_status)
df$Race <- as.factor(df$Race)
df$efs.ind <- as.factor(df$efs.ind)
df$menopausal_status <- as.factor(df$menopausal_status)
df$Arm <- as.factor(df$Arm)
df$pCR <- as.factor(df$pCR)
df$ethnicity <- as.factor(df$ethnicity)

colnames(df)

tab1 <- table1(~Age_at_Screening + Race + ethnicity + Arm + pCR + 
        menopausal_status + combined_status + efs.time + efs.ind, data = df)

tab1
```

## Read in data and set up variables (R)
```
df <- read.csv("/path/cleaned_imputed_df_2.csv")
colnames(df)

#change efs time to years
df$efs.time.years <- (df$efs.time)/365.25 

#log transformed ftv3
df$log_FTV3 <- log(df$FTV3+1)

#cube root ftv3
df$cube_FTV3 <- df$FTV3^(1/3)

#combined status
df$combined_status <- ifelse(df$HR == 1 & df$HER2 == 1, "HR+/HER2+",
                          ifelse(df$HR == 1, "HR+/HER2-",
                          ifelse(df$HER2 == 1, "HR-/HER2+", "HR-/HER2-")))

#make factor vars
df$HR <- as.factor(df$HR)
df$HER2 <- as.factor(df$HER2)
df$combined_status <- as.factor(df$combined_status)

```

## Run K-fold (5) cross validated cox proportional hazards (R)
We use 3 models with differing forms of FTV3, with concordance index as the evaluation metric. Efs.time is the time variable, and Efs.ind is the indicator variable. We choose the best model based on the highest cross validated concordance index.
```
set.seed(123456)
#make lists of predictor sets
predictor_lists <- list(
  c("FTV3", "combined_status"),
  c("log_FTV3", "combined_status"),
  c("cube_FTV3", "combined_status")
)

#set number of folds for cross-validation
num_folds <- 5

#create folds
folds <- cvFolds(nrow(df), K = num_folds, type = "random")

#initialize a list to store concordance indices for each model
c_indices_list <- vector("list", length(predictor_lists))

#cross-validation func
for(j in seq_along(predictor_lists)){
  c_indices <- numeric(num_folds)
  for(i in 1:num_folds){
    #split the data into train/test sets for this fold
    train_indices <- folds$subsets[folds$which != i]
    test_indices <- folds$subsets[folds$which == i]
    df_train_fold <- df[train_indices, ]
    df_test_fold <- df[test_indices, ]
    
    #fit cox model on training
    cox <- coxph(as.formula(paste("Surv(efs.time, efs.ind) ~", 
                 paste(predictor_lists[[j]], collapse = " + "))), data = df_train_fold)
    
    #calculate concordance index
    c_indices[i] <- concordance(cox, newdata = df_test_fold)
  }
  #store concordance indices for this model
  c_indices_list[[j]] <- list(c_indices)
}

#calculate mean concordance index across all folds for each model
mean_concordance_list <- sapply(c_indices_list, function(x) mean(unlist(x)))
mean_concordance_list

#find the best model based on cv
best_model_index <- which.max(mean_concordance_list)
best_model_predictors <- predictor_lists[[best_model_index]]

best_model_predictors
```

## Sensitivity analysis excluding missing values (R)

```
#read in original dataset
df_original <- read.csv("/path/ISPY2_n985_TCIA_FTVdata_with_EFS_clinical_forN899.csv")
colnames(df_original)

#subset relevant vars (efs.time, efs.ind, FTV3)
df <- df_original[,c(18, 19, 27, 29, 30)]

#summarize NA
x <- df %>% summarize_all(funs(sum(is.na(.))))
#x
result <- df %>% summarize_all(funs(sum(is.na(.)))) %>%
  gather(key = "variable", value = "missing count") %>%
  arrange(variable)
#result

#exclude missing from FTV3
df <- na.exclude(df)

#log transformed ftv3
df$log_FTV3 <- log(df$FTV3+1)

#cube root ftv3
df$cube_FTV3 <- df$FTV3^(1/3)

#combined status
df$combined_status <- ifelse(df$HR == 1 & df$HER2 == 1, "HR+/HER2+",
                          ifelse(df$HR == 1, "HR+/HER2-",
                          ifelse(df$HER2 == 1, "HR-/HER2+", "HR-/HER2-")))

#make factor vars
df$HR <- as.factor(df$HR)
df$HER2 <- as.factor(df$HER2)
df$combined_status <- as.factor(df$combined_status)

table(df$combined_status, df$efs.ind)
table(df$combined_status)
table(df$HR, df$HER2)
set.seed(123456)
#make lists of predictor sets
predictor_lists <- list(
  c("FTV3", "combined_status"),
  c("log_FTV3", "combined_status"),
  c("cube_FTV3", "combined_status")
)

#set number of folds for cross-validation
num_folds <- 5

#create folds
folds <- cvFolds(nrow(df), K = num_folds, type = "random")

#initialize a list to store concordance indices for each model
c_indices_list <- vector("list", length(predictor_lists))

#cross-validation func
for(j in seq_along(predictor_lists)){
  c_indices <- numeric(num_folds)
  for(i in 1:num_folds){
    #split the data into train/test sets for this fold
    train_indices <- folds$subsets[folds$which != i]
    test_indices <- folds$subsets[folds$which == i]
    df_train_fold <- df[train_indices, ]
    df_test_fold <- df[test_indices, ]
    
    #fit cox model on training
    cox <- coxph(as.formula(paste("Surv(efs.time, efs.ind) ~", paste(predictor_lists[[j]], collapse = " + "))), data = df_train_fold)
    
    #calculate concordance index
    c_indices[i] <- concordance(cox, newdata = df_test_fold)
  }
  #store concordance indices for this model
  c_indices_list[[j]] <- list(c_indices)
}

#calculate mean concordance index across all folds for each model
mean_concordance_list <- sapply(c_indices_list, function(x) mean(unlist(x)))
mean_concordance_list

#find the best model based on cv
best_model_index <- which.max(mean_concordance_list)
best_model_predictors <- predictor_lists[[best_model_index]]

best_model_predictors

```

## Time-Dependent ROC (SAS PROC PHREG procedure)

```
options nodate nocenter orientation=landscape missing=' ';
 
%let odate=06072024;
*/libname d1 "C:\Users\kevin\OneDrive\Documents\Data Science\Projects\CAPSTONE\EXCEL" ;
libname db  "C:\Users\kevin\OneDrive\Documents\Data Science\Projects\CAPSTONE\SASfiles" ;
%let results= C:\Users\kevin\OneDrive\Documents\Data Science\Projects\CAPSTONE\SASresults;
 
%let pname=%sysget(SAS_EXECFILEPATH );  /*** capture name of program for footnote on output ***/
%put &pname;
footnote1 "&pname";
  
 
/*** Standard formats ***/
proc format;
	 value meanblnk . = ' ' other = [10.3];
	 value missint  . = ' ' other = [10.];
	 value nblnk . = ' '    other = [10.];
	 value missblnk . = ' ' other = [pvalue6.4];
	 value yesno .= ' ' 0='NO' 1='YES';
	 value HER2_HR 
		.=' '
		0='HR-/HER2-' 
		1='HR+/HER2-' 
		2='HR-/HER2+'
		3='HR+/HER2+';
run; 
/*  import dataset   */
proc import datafile="C:\Users\kevin\OneDrive\Documents\Data Science\Projects\CAPSTONE\cleaned_imputed_df_2.csv" 
	out=imputed dbms=csv replace;
   getnames=yes;
   GUESSINGROWS= 2147483647;
run;
/*
data imputed;
  set d1.Cleaned_imputed_df_cj;
run;
*/

ods listing;
proc freq data=imputed;
table HER2*HR/list missing;
run;

/* feature engineering   */
data imputed2;
 set imputed (drop = combined_status);
	combined_status=.;*set missing to numeric variable;
	if HER2=0 and HR=0 then combined_status=0;
	if HER2=0 and HR=1 then combined_status=1;
	if HER2=1 and HR=0 then combined_status=2;
	if HER2=1 and HR=1 then combined_status=3;
 	 log_FTV3 =log( FTV3+1);
 	cube_FTV3 = FTV3**(1/3);
	efs_time_years=efs_time /365.25;
	format combined_status HER2_HR.;
run; 

ods listing;
proc freq data=imputed2;
table HER2*HR*combined_status/list missing;
run;
/*                                                          Cumulative    Cumulative

The FREQ Procedure

                                                          Cumulative    Cumulative
HER2    HR    combined_status    Frequency     Percent     Frequency      Percent

   0     0    HR-/HER2-               323       35.93           323        35.93
   0     1    HR+/HER2-               357       39.71           680        75.64
   1     0    HR-/HER2+                79        8.79           759        84.43
   1     1    HR+/HER2+               140       15.57           899       100.00


*/


ods listing;
proc contents data=imputed2 short;
run;
/*The SAS System                                                                                                                     6

The CONTENTS Procedure

                                           Alphabetic List of Variables for WORK.IMPUTED2

Age_at_Screening Arm FTV0 FTV1 FTV2 FTV3 HER2 HR MP RESEARCH_ID 
Race RowID combined_status cube_FTV3 drfs_ind drfs_time efs_ind
efs_time efs_time_years ethnicity log_FTV3 menopausal_status 
os_ind os_time pCR patient_diagnosed_distant_progre
patient_diagnosed_local_progress survival_status__c 
time_to_last_followup
  
*/


proc phreg data=imputed2  concordance plots=roc rocoptions(at=1, 2,  5 method = km);
	class combined_status;
   model efs_time_years*efs_ind(0)=log_FTV3 combined_status ;
   output out=imputed3  xbeta=Y_log;
run;

proc phreg data=imputed3  concordance plots=roc rocoptions(at=1, 2,  5 method = km);
	class combined_status;
   model efs_time_years*efs_ind(0)=cube_FTV3  combined_status ;
   output out=imputed4 xbeta=Y_cube;
run;

ods graphics on;
ods rtf file="&results\Time_Dependent_ROC3k_&odate..doc" style=listing bodytitle;
title "&results\Time_Dependent_ROC3k_&odate..doc";
title2 "1a. Time-Dependent ROC curves";
title3 "from the Cox model of efs_time_years";
title4 "with FTV3 and combined_status as predictors";
title5 "Using imputed data";
ods trace on;
 proc phreg data=imputed4   plots(overlay=individual)=roc   rocoptions(at=1, 2, 5 method = km);
 	class combined_status;
   model efs_time_years*efs_ind(0)= FTV3  combined_status / roclabel='FTV3+combined_status';
   roc  'FTV3+combined_status' FTV3  combined_status;
   roc 'log_FTV3+combined_status' pred=Y_log;
    roc 'cube_FTV3+combined_status' pred=Y_cube;
	*roccontrast reference ('FTV3+combined_status') / estimate e;
run;

ods rtf close;


Ods listing;

Ods trace on;

   proc phreg data = imputed2 concordance NAMELEN=200;

     class combined_status/param=ref;

     model efs_time_years*efs_ind(0) = combined_status FTV3/risklimits=PL;

    ods output ConvergenceStatus = _cs;

     ods output FitStatistics = _fs;

     ods output ParameterEstimates = _pe;

       ods output ModelANOVA= _t3; 

     ods output CensoredSummary = _cens;

     ods output  Concordance= _Cd;
   run;

   ods rtf file="&results\Cox_Test.doc" style=listing bodytitle;

title "&results\Cox_Test.doc" ;

title2 "Cox results";

proc phreg data = imputed2 concordance NAMELEN=200;

     class combined_status/param=ref;

     model efs_time_years*efs_ind(0) = combined_status FTV3/risklimits=PL;

    ods output ConvergenceStatus = _cs;

     ods output FitStatistics = _fs;

     ods output ParameterEstimates = _pe;

       ods output ModelANOVA= _t3; 

     ods output CensoredSummary = _cens;

     ods output  Concordance= _Cd;

     ods select ParameterEstimates ModelANOVA Concordance ;

  run;

 

ods rtf close;


```

## 6 Pairwise comparisons (R)
```
#"HR+/HER2+", "HR+/HER2-", "HR-/HER2+", "HR-/HER2-"

  #make HR-/HER2- reference group. 
df$combined_status <- relevel(df$combined_status, ref = "HR-/HER2-")
cox_model <- coxph(Surv(efs.time.years, efs.ind) ~ cube_FTV3 + combined_status, data = df)
summary(cox_model) 

  #make HR+/HER2+ reference group
df$combined_status <- relevel(df$combined_status, ref = "HR+/HER2+")
cox_model <- coxph(Surv(efs.time.years, efs.ind) ~ cube_FTV3 + combined_status, data = df)
summary(cox_model) 

  #make HR+/HER2- reference group
df$combined_status <- relevel(df$combined_status, ref = "HR+/HER2-")
cox_model <- coxph(Surv(efs.time.years, efs.ind) ~ cube_FTV3 + combined_status, data = df)
summary(cox_model) 

  #make HR-/HER2+ reference group
df$combined_status <- relevel(df$combined_status, ref = "HR-/HER2+")
cox_model <- coxph(Surv(efs.time.years, efs.ind) ~ cube_FTV3 + combined_status, data = df)
summary(cox_model) 

```

## ANOVA test for interactions (R)
```
#anova combined and combined w/ interaction
fit1 <- coxph(Surv(efs.time.years, efs.ind) ~ cube_FTV3 + combined_status, data = df)
fit2 <- coxph(Surv(efs.time.years, efs.ind) ~ cube_FTV3 + combined_status + cube_FTV3:combined_status, data = df) 
anova(fit1, fit2)

```

## Subgroup analysis (R)
```

#split data into subgroups based on categories in combined status
splits <- split(df, f = df$combined_status)

  # $`HR-/HER2-`
df_1 <- splits[[2]]
  # $`HR+/HER2-`
df_2 <- splits[[1]]
  # $`HR-/HER2+`
df_3 <- splits[[3]]
  # $`HR+/HER2+`
df_4 <- splits[[4]]

#fit for HR-/HER2- subgroup
cox_1 <- coxph(Surv(efs.time.years, efs.ind) ~ cube_FTV3, data = df_1)
summary(cox_1)

#fit for HR+/HER2- subgroup
cox_2 <- coxph(Surv(efs.time.years, efs.ind) ~ cube_FTV3, data = df_2)
summary(cox_2)

#fit for HR-/HER2+ subgroup
cox_3 <- coxph(Surv(efs.time.years, efs.ind) ~ cube_FTV3, data = df_3)
summary(cox_3)

# fit for HR+/HER2+ subgroup
cox_4 <- coxph(Surv(efs.time.years, efs.ind) ~ cube_FTV3, data = df_4)
summary(cox_4)


```