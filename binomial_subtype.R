#######
##      Binomial Classification (resistant/susceptible)
##
#######

## load all libraries to be used

library(caret)
library(nnet)
library(XLConnect)

## provide all filepaths

input_file = "C:/Users/hp-pc/Documents/India2016/Dissertation/Input.xlsx"
output_file = "C:/Users/hp-pc/Documents/India2016/Dissertation/Report.xlsx"
pi_dataset = "C:/Users/hp-pc/Documents/India2016/Dissertation/PI.csv"
nrti_dataset = "C:/Users/hp-pc/Documents/India2016/Dissertation/NRTI.csv"
nnrti_dataset = "C:/Users/hp-pc/Documents/India2016/Dissertation/NNRTI.csv"

## read from input file.

## select subtype based on which model is to be create from input sheet

wb = loadWorkbook(input_file)

subtype = readWorksheet(wb, sheet = "Input", startRow = 1, endRow = 2,
                        startCol = 1, endCol = 2)

subtype = subtype$Col2

## create input vector (mutation sequence) based upon input sheet

mutation = readWorksheet(wb, sheet = "Input", startRow = 5, endRow = 244,
                         startCol = 1, endCol = 2)

## initialize mutations and features

input = NULL
no_mut_input = NULL
max_mix_input = NULL
no_mut_input = 0
max_mix_input = 0

for (i in 1:240)
{
  for (j in 1:nrow(mutation))
  {
    if(i==as.integer(mutation[j,1]) )
    {
      input[i]=mutation[j,2]
      no_mut_input = no_mut_input+1
      break
    }
    else
    {
      input[i]=0
      
    }
  }
  
}

input = as.matrix(input)
input = t(input)

## now we have the mutation list in input.

## convert mutation to integer equivalent.

for(i in 1:ncol(input))
{
  ## coverts utf of alphabet (mutation) to integer value
  input[1,i] = utf8ToInt(input[1,i]) 
}

input = as.numeric(input)

## initialize vectors and values and cross validation parameter

final = NULL
cv = 10

#################
####### PI FPV - Fosamprenavir workspace
#################
## read the data from csv file

data <- read.csv(pi_dataset,header=TRUE,stringsAsFactors = FALSE)    

## subset data to construct model for FPV

data_drug = subset(data, select = -c(SeqID,ATV,IDV,LPV,NFV,SQV,TPV,DRV))
data_drug = na.omit(data_drug)

## select data as per subtype

if(is.null(subtype))
{
  subtype = "X"
}

## select data correspondin to subtype

if(subtype != "X")
{
  data_drug = data_drug[data_drug$SubType == subtype,]
}

## remove the subtype column now

data_drug = subset(data_drug, select = -c(SubType))

## replace non-mutant positions with 0

data_drug[data_drug == "-"] = 0

## replace the mutant positions with thier integer codes

index = NULL
no_mut = NULL
max_mix = NULL

for (i in 1:nrow(data_drug))
{
  no_mut[i]=0
  max_mix[i]=0
  for (j in 2:NCOL(data_drug))
  {
    x=utf8ToInt(data_drug[i,j])
    sum = 0
    if(length(x)>1)     #### change to 0
    {
      no_mut[i]=no_mut[i]+1
    }
    if(length(x)>max_mix[i])
    {
      max_mix[i]=length(x) 
    }
    for(k in 1:1)
    {
      sum = sum + x[k]
    }
    data_drug[i,j] = as.numeric(sum)
  }
}

## set cut-off manually

cutoff = 4 #11
data_drug$bin = ifelse(data_drug$FPV >= cutoff,1,0)
data_drug = subset(data_drug, select = -c(FPV))

## add features - no of mutations in sequence, maximum mixture length

data_drug$no_mut = no_mut
data_drug$max_mix = max_mix


## convert to numeric  

data_drug$P1 = as.numeric(data_drug$P1)
data_drug$P2 = as.numeric(data_drug$P2)
data_drug$P3 = as.numeric(data_drug$P3)
data_drug$P4 = as.numeric(data_drug$P4)
data_drug$P5 = as.numeric(data_drug$P5)
data_drug$P6 = as.numeric(data_drug$P6)
data_drug$P7 = as.numeric(data_drug$P7)
data_drug$P8 = as.numeric(data_drug$P8)
data_drug$P9 = as.numeric(data_drug$P9)
data_drug$P10 = as.numeric(data_drug$P10)
data_drug$P11 = as.numeric(data_drug$P11)
data_drug$P12 = as.numeric(data_drug$P12)
data_drug$P13 = as.numeric(data_drug$P13)
data_drug$P14 = as.numeric(data_drug$P14)
data_drug$P15 = as.numeric(data_drug$P15)
data_drug$P16 = as.numeric(data_drug$P16)
data_drug$P17 = as.numeric(data_drug$P17)
data_drug$P18 = as.numeric(data_drug$P18)
data_drug$P19 = as.numeric(data_drug$P19)
data_drug$P20 = as.numeric(data_drug$P20)
data_drug$P21 = as.numeric(data_drug$P21)
data_drug$P22 = as.numeric(data_drug$P22)
data_drug$P23 = as.numeric(data_drug$P23)
data_drug$P24 = as.numeric(data_drug$P24)
data_drug$P25 = as.numeric(data_drug$P25)
data_drug$P26 = as.numeric(data_drug$P26)
data_drug$P27 = as.numeric(data_drug$P27)
data_drug$P28 = as.numeric(data_drug$P28)
data_drug$P29 = as.numeric(data_drug$P29)
data_drug$P30 = as.numeric(data_drug$P30)
data_drug$P31 = as.numeric(data_drug$P31)
data_drug$P32 = as.numeric(data_drug$P32)
data_drug$P33 = as.numeric(data_drug$P33)
data_drug$P34 = as.numeric(data_drug$P34)
data_drug$P35 = as.numeric(data_drug$P35)
data_drug$P36 = as.numeric(data_drug$P36)
data_drug$P37 = as.numeric(data_drug$P37)
data_drug$P38 = as.numeric(data_drug$P38)
data_drug$P39 = as.numeric(data_drug$P39)
data_drug$P40 = as.numeric(data_drug$P40)
data_drug$P41 = as.numeric(data_drug$P41)
data_drug$P42 = as.numeric(data_drug$P42)
data_drug$P43 = as.numeric(data_drug$P43)
data_drug$P44 = as.numeric(data_drug$P44)
data_drug$P45 = as.numeric(data_drug$P45)
data_drug$P46 = as.numeric(data_drug$P46)
data_drug$P47 = as.numeric(data_drug$P47)
data_drug$P48 = as.numeric(data_drug$P48)
data_drug$P49 = as.numeric(data_drug$P49)
data_drug$P50 = as.numeric(data_drug$P50)
data_drug$P51 = as.numeric(data_drug$P51)
data_drug$P52 = as.numeric(data_drug$P52)
data_drug$P53 = as.numeric(data_drug$P53)
data_drug$P54 = as.numeric(data_drug$P54)
data_drug$P55 = as.numeric(data_drug$P55)
data_drug$P56 = as.numeric(data_drug$P56)
data_drug$P57 = as.numeric(data_drug$P57)
data_drug$P58 = as.numeric(data_drug$P58)
data_drug$P59 = as.numeric(data_drug$P59)
data_drug$P60 = as.numeric(data_drug$P60)
data_drug$P61 = as.numeric(data_drug$P61)
data_drug$P62 = as.numeric(data_drug$P62)
data_drug$P63 = as.numeric(data_drug$P63)
data_drug$P64 = as.numeric(data_drug$P64)
data_drug$P65 = as.numeric(data_drug$P65)
data_drug$P66 = as.numeric(data_drug$P66)
data_drug$P67 = as.numeric(data_drug$P67)
data_drug$P68 = as.numeric(data_drug$P68)
data_drug$P69 = as.numeric(data_drug$P69)
data_drug$P70 = as.numeric(data_drug$P70)
data_drug$P71 = as.numeric(data_drug$P71)
data_drug$P72 = as.numeric(data_drug$P72)
data_drug$P73 = as.numeric(data_drug$P73)
data_drug$P74 = as.numeric(data_drug$P74)
data_drug$P75 = as.numeric(data_drug$P75)
data_drug$P76 = as.numeric(data_drug$P76)
data_drug$P77 = as.numeric(data_drug$P77)
data_drug$P78 = as.numeric(data_drug$P78)
data_drug$P79 = as.numeric(data_drug$P79)
data_drug$P80 = as.numeric(data_drug$P80)
data_drug$P81 = as.numeric(data_drug$P81)
data_drug$P82 = as.numeric(data_drug$P82)
data_drug$P83 = as.numeric(data_drug$P83)
data_drug$P84 = as.numeric(data_drug$P84)
data_drug$P85 = as.numeric(data_drug$P85)
data_drug$P86 = as.numeric(data_drug$P86)
data_drug$P87 = as.numeric(data_drug$P87)
data_drug$P88 = as.numeric(data_drug$P88)
data_drug$P89 = as.numeric(data_drug$P89)
data_drug$P90 = as.numeric(data_drug$P90)
data_drug$P91 = as.numeric(data_drug$P91)
data_drug$P92 = as.numeric(data_drug$P92)
data_drug$P93 = as.numeric(data_drug$P93)
data_drug$P94 = as.numeric(data_drug$P94)
data_drug$P95 = as.numeric(data_drug$P95)
data_drug$P96 = as.numeric(data_drug$P96)
data_drug$P97 = as.numeric(data_drug$P97)
data_drug$P98 = as.numeric(data_drug$P98)
data_drug$P99 = as.numeric(data_drug$P99)

#### sample CV

totalAccuracy <- c()
totalSensitivity = c()
totalSpecificity = c()
cvDivider <- floor(nrow(data_drug) / (cv+1))

for (cv in seq(1:cv)) {
  # assign chunk to data test
  dataTestIndex <- c((cv * cvDivider):(cv * cvDivider + cvDivider))
  dataTest <- data_drug[dataTestIndex,]
  # everything else to train
  dataTrain <- data_drug[-dataTestIndex,]
  
  cylModel <- multinom(bin~., data=dataTrain, maxit=1000, trace=T) 
  
  pred <- predict(cylModel, newdata=dataTest, type="class")
  
  #  classification error
  cv_ac <- postResample(dataTest$bin, pred)[[1]]
  cv_sen = sensitivity(as.factor(dataTest$bin), pred)
  cv_spe = specificity(as.factor(dataTest$bin), pred)
  print(paste('Current Accuracy:',cv_ac,'for CV:',cv))
  totalAccuracy <- c(totalAccuracy, cv_ac)
  totalSensitivity = c(totalSensitivity, cv_sen)
  totalSpecificity = c(totalSpecificity,cv_spe)
}

mean(totalAccuracy)
mean(totalSensitivity,na.rm = TRUE)
mean(totalSpecificity,na.rm = TRUE)

## check if data available for both classes

comment = "Model Valid"

if(is.na(mean(totalAccuracy)))
{
  comment = "Model Invalid as data available for only one class during cross validation."
  decision = NA
}


if(!is.na(mean(totalAccuracy)))
{  

#################
## main prediction of input
#################

## create test seqence with 99 postitons from input, dummy binary, number of mutations and mixtures length

test_seq = as.matrix(c(input[1:99],0,no_mut_input,max_mix_input))
test_seq = t(test_seq)
colnames(test_seq) = colnames(data_drug)
test_seq = data.frame(test_seq)

## create model on whole data

final_model <- multinom(bin~., data=data_drug, maxit=1000, trace=T)

## predict for input sequence using model

main_pred = predict(final_model, newdata=test_seq, type="class")

## formulate string for prediction / remarks

  if(main_pred == 1)
  {
    decision = "Resistant"
  }
  if(main_pred == 0)
  {
    decision = "Susceptible"
  }

}
## write to result dataframe

final = data.frame(INHIBITOR = "PI",DRUG = "FPV - Fosamprenavir", NO_OF_SEQUENCES = nrow(data_drug),
                   CUTOFF = cutoff,MODEL_ACCURACY = round(mean(totalAccuracy),4)*100,
                   MODEL_SENSITIVITY = round(mean(totalSensitivity,na.rm = TRUE),4)*100,
                   MODEL_SPECIFICITY = round(mean(totalSpecificity,na.rm = TRUE),4)*100,
                   MODEL_PREDICTION = decision, REMARKS = comment, stringsAsFactors = FALSE)



##############
#######
####### PI NFV - Nelfinavir workspace
#######
##############

## read the data from csv file

data <- read.csv(pi_dataset,header=TRUE,stringsAsFactors = FALSE)    

## subset data to construct model for FPV

data_drug = subset(data, select = -c(SeqID,ATV,IDV,FPV,LPV,SQV,TPV,DRV))
data_drug = na.omit(data_drug)

## select data as per subtype

if(is.null(subtype))
{
  subtype = "X"
}

## select data correspondin to subtype

if(subtype != "X")
{
  data_drug = data_drug[data_drug$SubType == subtype,]
}

## remove the subtype column now

data_drug = subset(data_drug, select = -c(SubType))

## replace non-mutant positions with 0

data_drug[data_drug == "-"] = 0

## replace the mutant positions with thier integer codes

index = NULL
no_mut = NULL
max_mix = NULL

for (i in 1:nrow(data_drug))
{
  no_mut[i]=0
  max_mix[i]=0
  for (j in 2:NCOL(data_drug))
  {
    x=utf8ToInt(data_drug[i,j])
    sum = 0
    if(length(x)>1)     #### change to 0
    {
      no_mut[i]=no_mut[i]+1
    }
    if(length(x)>max_mix[i])
    {
      max_mix[i]=length(x) 
    }
    for(k in 1:1)
    {
      sum = sum + x[k]
    }
    data_drug[i,j] = as.numeric(sum)
  }
}

## set cut-off manually

cutoff = 3 #3
data_drug$bin = ifelse(data_drug$NFV >= cutoff,1,0)
data_drug = subset(data_drug, select = -c(NFV))


## add features - no of mutations in sequence, maximum mixture length

data_drug$no_mut = no_mut
data_drug$max_mix = max_mix


## convert to numeric  

data_drug$P1 = as.numeric(data_drug$P1)
data_drug$P2 = as.numeric(data_drug$P2)
data_drug$P3 = as.numeric(data_drug$P3)
data_drug$P4 = as.numeric(data_drug$P4)
data_drug$P5 = as.numeric(data_drug$P5)
data_drug$P6 = as.numeric(data_drug$P6)
data_drug$P7 = as.numeric(data_drug$P7)
data_drug$P8 = as.numeric(data_drug$P8)
data_drug$P9 = as.numeric(data_drug$P9)
data_drug$P10 = as.numeric(data_drug$P10)
data_drug$P11 = as.numeric(data_drug$P11)
data_drug$P12 = as.numeric(data_drug$P12)
data_drug$P13 = as.numeric(data_drug$P13)
data_drug$P14 = as.numeric(data_drug$P14)
data_drug$P15 = as.numeric(data_drug$P15)
data_drug$P16 = as.numeric(data_drug$P16)
data_drug$P17 = as.numeric(data_drug$P17)
data_drug$P18 = as.numeric(data_drug$P18)
data_drug$P19 = as.numeric(data_drug$P19)
data_drug$P20 = as.numeric(data_drug$P20)
data_drug$P21 = as.numeric(data_drug$P21)
data_drug$P22 = as.numeric(data_drug$P22)
data_drug$P23 = as.numeric(data_drug$P23)
data_drug$P24 = as.numeric(data_drug$P24)
data_drug$P25 = as.numeric(data_drug$P25)
data_drug$P26 = as.numeric(data_drug$P26)
data_drug$P27 = as.numeric(data_drug$P27)
data_drug$P28 = as.numeric(data_drug$P28)
data_drug$P29 = as.numeric(data_drug$P29)
data_drug$P30 = as.numeric(data_drug$P30)
data_drug$P31 = as.numeric(data_drug$P31)
data_drug$P32 = as.numeric(data_drug$P32)
data_drug$P33 = as.numeric(data_drug$P33)
data_drug$P34 = as.numeric(data_drug$P34)
data_drug$P35 = as.numeric(data_drug$P35)
data_drug$P36 = as.numeric(data_drug$P36)
data_drug$P37 = as.numeric(data_drug$P37)
data_drug$P38 = as.numeric(data_drug$P38)
data_drug$P39 = as.numeric(data_drug$P39)
data_drug$P40 = as.numeric(data_drug$P40)
data_drug$P41 = as.numeric(data_drug$P41)
data_drug$P42 = as.numeric(data_drug$P42)
data_drug$P43 = as.numeric(data_drug$P43)
data_drug$P44 = as.numeric(data_drug$P44)
data_drug$P45 = as.numeric(data_drug$P45)
data_drug$P46 = as.numeric(data_drug$P46)
data_drug$P47 = as.numeric(data_drug$P47)
data_drug$P48 = as.numeric(data_drug$P48)
data_drug$P49 = as.numeric(data_drug$P49)
data_drug$P50 = as.numeric(data_drug$P50)
data_drug$P51 = as.numeric(data_drug$P51)
data_drug$P52 = as.numeric(data_drug$P52)
data_drug$P53 = as.numeric(data_drug$P53)
data_drug$P54 = as.numeric(data_drug$P54)
data_drug$P55 = as.numeric(data_drug$P55)
data_drug$P56 = as.numeric(data_drug$P56)
data_drug$P57 = as.numeric(data_drug$P57)
data_drug$P58 = as.numeric(data_drug$P58)
data_drug$P59 = as.numeric(data_drug$P59)
data_drug$P60 = as.numeric(data_drug$P60)
data_drug$P61 = as.numeric(data_drug$P61)
data_drug$P62 = as.numeric(data_drug$P62)
data_drug$P63 = as.numeric(data_drug$P63)
data_drug$P64 = as.numeric(data_drug$P64)
data_drug$P65 = as.numeric(data_drug$P65)
data_drug$P66 = as.numeric(data_drug$P66)
data_drug$P67 = as.numeric(data_drug$P67)
data_drug$P68 = as.numeric(data_drug$P68)
data_drug$P69 = as.numeric(data_drug$P69)
data_drug$P70 = as.numeric(data_drug$P70)
data_drug$P71 = as.numeric(data_drug$P71)
data_drug$P72 = as.numeric(data_drug$P72)
data_drug$P73 = as.numeric(data_drug$P73)
data_drug$P74 = as.numeric(data_drug$P74)
data_drug$P75 = as.numeric(data_drug$P75)
data_drug$P76 = as.numeric(data_drug$P76)
data_drug$P77 = as.numeric(data_drug$P77)
data_drug$P78 = as.numeric(data_drug$P78)
data_drug$P79 = as.numeric(data_drug$P79)
data_drug$P80 = as.numeric(data_drug$P80)
data_drug$P81 = as.numeric(data_drug$P81)
data_drug$P82 = as.numeric(data_drug$P82)
data_drug$P83 = as.numeric(data_drug$P83)
data_drug$P84 = as.numeric(data_drug$P84)
data_drug$P85 = as.numeric(data_drug$P85)
data_drug$P86 = as.numeric(data_drug$P86)
data_drug$P87 = as.numeric(data_drug$P87)
data_drug$P88 = as.numeric(data_drug$P88)
data_drug$P89 = as.numeric(data_drug$P89)
data_drug$P90 = as.numeric(data_drug$P90)
data_drug$P91 = as.numeric(data_drug$P91)
data_drug$P92 = as.numeric(data_drug$P92)
data_drug$P93 = as.numeric(data_drug$P93)
data_drug$P94 = as.numeric(data_drug$P94)
data_drug$P95 = as.numeric(data_drug$P95)
data_drug$P96 = as.numeric(data_drug$P96)
data_drug$P97 = as.numeric(data_drug$P97)
data_drug$P98 = as.numeric(data_drug$P98)
data_drug$P99 = as.numeric(data_drug$P99)

#### sample CV

totalAccuracy <- c()
totalSensitivity = c()
totalSpecificity = c()
cvDivider <- floor(nrow(data_drug) / (cv+1))

for (cv in seq(1:cv)) {
  # assign chunk to data test
  dataTestIndex <- c((cv * cvDivider):(cv * cvDivider + cvDivider))
  dataTest <- data_drug[dataTestIndex,]
  # everything else to train
  dataTrain <- data_drug[-dataTestIndex,]
  
  cylModel <- multinom(bin~., data=dataTrain, maxit=1000, trace=T) 
  
  pred <- predict(cylModel, newdata=dataTest, type="class")
  
  #  classification error
  cv_ac <- postResample(dataTest$bin, pred)[[1]]
  cv_sen = sensitivity(as.factor(dataTest$bin), pred)
  cv_spe = specificity(as.factor(dataTest$bin), pred)
  print(paste('Current Accuracy:',cv_ac,'for CV:',cv))
  totalAccuracy <- c(totalAccuracy, cv_ac)
  totalSensitivity = c(totalSensitivity, cv_sen)
  totalSpecificity = c(totalSpecificity,cv_spe)
}

mean(totalAccuracy)
mean(totalSensitivity,na.rm = TRUE)
mean(totalSpecificity,na.rm = TRUE)

## check if data available for both classes

comment = "Model Valid"

if(is.na(mean(totalAccuracy)))
{
  comment = "Model Invalid as data available for only one class during cross validation."
  decision = NA
}


if(!is.na(mean(totalAccuracy)))
{  
  
  #################
  ## main prediction of input
  #################
  
  ## create test seqence with 99 postitons from input, dummy binary, number of mutations and mixtures length
  
  test_seq = as.matrix(c(input[1:99],0,no_mut_input,max_mix_input))
  test_seq = t(test_seq)
  colnames(test_seq) = colnames(data_drug)
  test_seq = data.frame(test_seq)
  
  ## create model on whole data
  
  final_model <- multinom(bin~., data=data_drug, maxit=1000, trace=T)
  
  ## predict for input sequence using model
  
  main_pred = predict(final_model, newdata=test_seq, type="class")
  
  ## formulate string for prediction / remarks
  
  if(main_pred == 1)
  {
    decision = "Resistant"
  }
  if(main_pred == 0)
  {
    decision = "Susceptible"
  }
  
}
## write to result dataframe

final = rbind(final,c("PI","NFV - Nelfinavir",nrow(data_drug),cutoff,round(mean(totalAccuracy),4)*100,
                      round(mean(totalSensitivity,na.rm = TRUE),4)*100,
                      round(mean(totalSpecificity,na.rm = TRUE),4)*100,decision,comment))


##############
#######
####### PI ATV - Atazanavir workspace
#######
##############

## read the data from csv file

data <- read.csv(pi_dataset,header=TRUE,stringsAsFactors = FALSE)    

## subset data to construct model for FPV

data_drug = subset(data, select = -c(SeqID,NFV,IDV,FPV,LPV,SQV,TPV,DRV))
data_drug = na.omit(data_drug)

## select data as per subtype

if(is.null(subtype))
{
  subtype = "X"
}

## select data correspondin to subtype

if(subtype != "X")
{
  data_drug = data_drug[data_drug$SubType == subtype,]
}

## remove the subtype column now

data_drug = subset(data_drug, select = -c(SubType))

## replace non-mutant positions with 0

data_drug[data_drug == "-"] = 0

## replace the mutant positions with thier integer codes

index = NULL
no_mut = NULL
max_mix = NULL

for (i in 1:nrow(data_drug))
{
  no_mut[i]=0
  max_mix[i]=0
  for (j in 2:NCOL(data_drug))
  {
    x=utf8ToInt(data_drug[i,j])
    sum = 0
    if(length(x)>1)     #### change to 0
    {
      no_mut[i]=no_mut[i]+1
    }
    if(length(x)>max_mix[i])
    {
      max_mix[i]=length(x) 
    }
    for(k in 1:1)
    {
      sum = sum + x[k]
    }
    data_drug[i,j] = as.numeric(sum)
  }
}

## set cut-off manually

cutoff = 3 #15
data_drug$bin = ifelse(data_drug$ATV >= cutoff,1,0)
data_drug = subset(data_drug, select = -c(ATV))

## add features - no of mutations in sequence, maximum mixture length

data_drug$no_mut = no_mut
data_drug$max_mix = max_mix


## convert to numeric  

data_drug$P1 = as.numeric(data_drug$P1)
data_drug$P2 = as.numeric(data_drug$P2)
data_drug$P3 = as.numeric(data_drug$P3)
data_drug$P4 = as.numeric(data_drug$P4)
data_drug$P5 = as.numeric(data_drug$P5)
data_drug$P6 = as.numeric(data_drug$P6)
data_drug$P7 = as.numeric(data_drug$P7)
data_drug$P8 = as.numeric(data_drug$P8)
data_drug$P9 = as.numeric(data_drug$P9)
data_drug$P10 = as.numeric(data_drug$P10)
data_drug$P11 = as.numeric(data_drug$P11)
data_drug$P12 = as.numeric(data_drug$P12)
data_drug$P13 = as.numeric(data_drug$P13)
data_drug$P14 = as.numeric(data_drug$P14)
data_drug$P15 = as.numeric(data_drug$P15)
data_drug$P16 = as.numeric(data_drug$P16)
data_drug$P17 = as.numeric(data_drug$P17)
data_drug$P18 = as.numeric(data_drug$P18)
data_drug$P19 = as.numeric(data_drug$P19)
data_drug$P20 = as.numeric(data_drug$P20)
data_drug$P21 = as.numeric(data_drug$P21)
data_drug$P22 = as.numeric(data_drug$P22)
data_drug$P23 = as.numeric(data_drug$P23)
data_drug$P24 = as.numeric(data_drug$P24)
data_drug$P25 = as.numeric(data_drug$P25)
data_drug$P26 = as.numeric(data_drug$P26)
data_drug$P27 = as.numeric(data_drug$P27)
data_drug$P28 = as.numeric(data_drug$P28)
data_drug$P29 = as.numeric(data_drug$P29)
data_drug$P30 = as.numeric(data_drug$P30)
data_drug$P31 = as.numeric(data_drug$P31)
data_drug$P32 = as.numeric(data_drug$P32)
data_drug$P33 = as.numeric(data_drug$P33)
data_drug$P34 = as.numeric(data_drug$P34)
data_drug$P35 = as.numeric(data_drug$P35)
data_drug$P36 = as.numeric(data_drug$P36)
data_drug$P37 = as.numeric(data_drug$P37)
data_drug$P38 = as.numeric(data_drug$P38)
data_drug$P39 = as.numeric(data_drug$P39)
data_drug$P40 = as.numeric(data_drug$P40)
data_drug$P41 = as.numeric(data_drug$P41)
data_drug$P42 = as.numeric(data_drug$P42)
data_drug$P43 = as.numeric(data_drug$P43)
data_drug$P44 = as.numeric(data_drug$P44)
data_drug$P45 = as.numeric(data_drug$P45)
data_drug$P46 = as.numeric(data_drug$P46)
data_drug$P47 = as.numeric(data_drug$P47)
data_drug$P48 = as.numeric(data_drug$P48)
data_drug$P49 = as.numeric(data_drug$P49)
data_drug$P50 = as.numeric(data_drug$P50)
data_drug$P51 = as.numeric(data_drug$P51)
data_drug$P52 = as.numeric(data_drug$P52)
data_drug$P53 = as.numeric(data_drug$P53)
data_drug$P54 = as.numeric(data_drug$P54)
data_drug$P55 = as.numeric(data_drug$P55)
data_drug$P56 = as.numeric(data_drug$P56)
data_drug$P57 = as.numeric(data_drug$P57)
data_drug$P58 = as.numeric(data_drug$P58)
data_drug$P59 = as.numeric(data_drug$P59)
data_drug$P60 = as.numeric(data_drug$P60)
data_drug$P61 = as.numeric(data_drug$P61)
data_drug$P62 = as.numeric(data_drug$P62)
data_drug$P63 = as.numeric(data_drug$P63)
data_drug$P64 = as.numeric(data_drug$P64)
data_drug$P65 = as.numeric(data_drug$P65)
data_drug$P66 = as.numeric(data_drug$P66)
data_drug$P67 = as.numeric(data_drug$P67)
data_drug$P68 = as.numeric(data_drug$P68)
data_drug$P69 = as.numeric(data_drug$P69)
data_drug$P70 = as.numeric(data_drug$P70)
data_drug$P71 = as.numeric(data_drug$P71)
data_drug$P72 = as.numeric(data_drug$P72)
data_drug$P73 = as.numeric(data_drug$P73)
data_drug$P74 = as.numeric(data_drug$P74)
data_drug$P75 = as.numeric(data_drug$P75)
data_drug$P76 = as.numeric(data_drug$P76)
data_drug$P77 = as.numeric(data_drug$P77)
data_drug$P78 = as.numeric(data_drug$P78)
data_drug$P79 = as.numeric(data_drug$P79)
data_drug$P80 = as.numeric(data_drug$P80)
data_drug$P81 = as.numeric(data_drug$P81)
data_drug$P82 = as.numeric(data_drug$P82)
data_drug$P83 = as.numeric(data_drug$P83)
data_drug$P84 = as.numeric(data_drug$P84)
data_drug$P85 = as.numeric(data_drug$P85)
data_drug$P86 = as.numeric(data_drug$P86)
data_drug$P87 = as.numeric(data_drug$P87)
data_drug$P88 = as.numeric(data_drug$P88)
data_drug$P89 = as.numeric(data_drug$P89)
data_drug$P90 = as.numeric(data_drug$P90)
data_drug$P91 = as.numeric(data_drug$P91)
data_drug$P92 = as.numeric(data_drug$P92)
data_drug$P93 = as.numeric(data_drug$P93)
data_drug$P94 = as.numeric(data_drug$P94)
data_drug$P95 = as.numeric(data_drug$P95)
data_drug$P96 = as.numeric(data_drug$P96)
data_drug$P97 = as.numeric(data_drug$P97)
data_drug$P98 = as.numeric(data_drug$P98)
data_drug$P99 = as.numeric(data_drug$P99)

#### sample CV

totalAccuracy <- c()
totalSensitivity = c()
totalSpecificity = c()
cvDivider <- floor(nrow(data_drug) / (cv+1))

for (cv in seq(1:cv)) {
  # assign chunk to data test
  dataTestIndex <- c((cv * cvDivider):(cv * cvDivider + cvDivider))
  dataTest <- data_drug[dataTestIndex,]
  # everything else to train
  dataTrain <- data_drug[-dataTestIndex,]
  
  cylModel <- multinom(bin~., data=dataTrain, maxit=1000, trace=T) 
  
  pred <- predict(cylModel, newdata=dataTest, type="class")
  
  #  classification error
  cv_ac <- postResample(dataTest$bin, pred)[[1]]
  cv_sen = sensitivity(as.factor(dataTest$bin), pred)
  cv_spe = specificity(as.factor(dataTest$bin), pred)
  print(paste('Current Accuracy:',cv_ac,'for CV:',cv))
  totalAccuracy <- c(totalAccuracy, cv_ac)
  totalSensitivity = c(totalSensitivity, cv_sen)
  totalSpecificity = c(totalSpecificity,cv_spe)
}

mean(totalAccuracy)
mean(totalSensitivity)
mean(totalSpecificity)

## check if data available for both classes

comment = "Model Valid"

if(is.na(mean(totalAccuracy)))
{
  comment = "Model Invalid as data available for only one class during cross validation."
  decision = NA
}


if(!is.na(mean(totalAccuracy)))
{  
  
  #################
  ## main prediction of input
  #################
  
  ## create test seqence with 99 postitons from input, dummy binary, number of mutations and mixtures length
  
  test_seq = as.matrix(c(input[1:99],0,no_mut_input,max_mix_input))
  test_seq = t(test_seq)
  colnames(test_seq) = colnames(data_drug)
  test_seq = data.frame(test_seq)
  
  ## create model on whole data
  
  final_model <- multinom(bin~., data=data_drug, maxit=1000, trace=T)
  
  ## predict for input sequence using model
  
  main_pred = predict(final_model, newdata=test_seq, type="class")
  
  ## formulate string for prediction / remarks
  
  if(main_pred == 1)
  {
    decision = "Resistant"
  }
  if(main_pred == 0)
  {
    decision = "Susceptible"
  }
  
}
## write to result dataframe

final = rbind(final,c("PI","ATV - Atazanavir",nrow(data_drug),cutoff,round(mean(totalAccuracy),4)*100,
                      round(mean(totalSensitivity,na.rm = TRUE),4)*100,
                      round(mean(totalSpecificity,na.rm = TRUE),4)*100,decision,comment))



##############
#######
####### PI DRV - Darunavir workspace
#######
##############

## read the data from csv file

data <- read.csv(pi_dataset,header=TRUE,stringsAsFactors = FALSE)    

## subset data to construct model for FPV

data_drug = subset(data, select = -c(SeqID,NFV,IDV,FPV,LPV,SQV,TPV,ATV))
data_drug = na.omit(data_drug)

## select data as per subtype

if(is.null(subtype))
{
  subtype = "X"
}

## select data correspondin to subtype

if(subtype != "X")
{
  data_drug = data_drug[data_drug$SubType == subtype,]
}

## remove the subtype column now

data_drug = subset(data_drug, select = -c(SubType))

## replace non-mutant positions with 0

data_drug[data_drug == "-"] = 0

## replace the mutant positions with thier integer codes

index = NULL
no_mut = NULL
max_mix = NULL

for (i in 1:nrow(data_drug))
{
  no_mut[i]=0
  max_mix[i]=0
  for (j in 2:NCOL(data_drug))
  {
    x=utf8ToInt(data_drug[i,j])
    sum = 0
    if(length(x)>1)     #### change to 0
    {
      no_mut[i]=no_mut[i]+1
    }
    if(length(x)>max_mix[i])
    {
      max_mix[i]=length(x) 
    }
    for(k in 1:1)
    {
      sum = sum + x[k]
    }
    data_drug[i,j] = as.numeric(sum)
  }
}

## set cut-off manually

cutoff = 10 ##90
data_drug$bin = ifelse(data_drug$DRV >= cutoff,1,0)
data_drug = subset(data_drug, select = -c(DRV))

## add features - no of mutations in sequence, maximum mixture length

data_drug$no_mut = no_mut
data_drug$max_mix = max_mix


## convert to numeric  

data_drug$P1 = as.numeric(data_drug$P1)
data_drug$P2 = as.numeric(data_drug$P2)
data_drug$P3 = as.numeric(data_drug$P3)
data_drug$P4 = as.numeric(data_drug$P4)
data_drug$P5 = as.numeric(data_drug$P5)
data_drug$P6 = as.numeric(data_drug$P6)
data_drug$P7 = as.numeric(data_drug$P7)
data_drug$P8 = as.numeric(data_drug$P8)
data_drug$P9 = as.numeric(data_drug$P9)
data_drug$P10 = as.numeric(data_drug$P10)
data_drug$P11 = as.numeric(data_drug$P11)
data_drug$P12 = as.numeric(data_drug$P12)
data_drug$P13 = as.numeric(data_drug$P13)
data_drug$P14 = as.numeric(data_drug$P14)
data_drug$P15 = as.numeric(data_drug$P15)
data_drug$P16 = as.numeric(data_drug$P16)
data_drug$P17 = as.numeric(data_drug$P17)
data_drug$P18 = as.numeric(data_drug$P18)
data_drug$P19 = as.numeric(data_drug$P19)
data_drug$P20 = as.numeric(data_drug$P20)
data_drug$P21 = as.numeric(data_drug$P21)
data_drug$P22 = as.numeric(data_drug$P22)
data_drug$P23 = as.numeric(data_drug$P23)
data_drug$P24 = as.numeric(data_drug$P24)
data_drug$P25 = as.numeric(data_drug$P25)
data_drug$P26 = as.numeric(data_drug$P26)
data_drug$P27 = as.numeric(data_drug$P27)
data_drug$P28 = as.numeric(data_drug$P28)
data_drug$P29 = as.numeric(data_drug$P29)
data_drug$P30 = as.numeric(data_drug$P30)
data_drug$P31 = as.numeric(data_drug$P31)
data_drug$P32 = as.numeric(data_drug$P32)
data_drug$P33 = as.numeric(data_drug$P33)
data_drug$P34 = as.numeric(data_drug$P34)
data_drug$P35 = as.numeric(data_drug$P35)
data_drug$P36 = as.numeric(data_drug$P36)
data_drug$P37 = as.numeric(data_drug$P37)
data_drug$P38 = as.numeric(data_drug$P38)
data_drug$P39 = as.numeric(data_drug$P39)
data_drug$P40 = as.numeric(data_drug$P40)
data_drug$P41 = as.numeric(data_drug$P41)
data_drug$P42 = as.numeric(data_drug$P42)
data_drug$P43 = as.numeric(data_drug$P43)
data_drug$P44 = as.numeric(data_drug$P44)
data_drug$P45 = as.numeric(data_drug$P45)
data_drug$P46 = as.numeric(data_drug$P46)
data_drug$P47 = as.numeric(data_drug$P47)
data_drug$P48 = as.numeric(data_drug$P48)
data_drug$P49 = as.numeric(data_drug$P49)
data_drug$P50 = as.numeric(data_drug$P50)
data_drug$P51 = as.numeric(data_drug$P51)
data_drug$P52 = as.numeric(data_drug$P52)
data_drug$P53 = as.numeric(data_drug$P53)
data_drug$P54 = as.numeric(data_drug$P54)
data_drug$P55 = as.numeric(data_drug$P55)
data_drug$P56 = as.numeric(data_drug$P56)
data_drug$P57 = as.numeric(data_drug$P57)
data_drug$P58 = as.numeric(data_drug$P58)
data_drug$P59 = as.numeric(data_drug$P59)
data_drug$P60 = as.numeric(data_drug$P60)
data_drug$P61 = as.numeric(data_drug$P61)
data_drug$P62 = as.numeric(data_drug$P62)
data_drug$P63 = as.numeric(data_drug$P63)
data_drug$P64 = as.numeric(data_drug$P64)
data_drug$P65 = as.numeric(data_drug$P65)
data_drug$P66 = as.numeric(data_drug$P66)
data_drug$P67 = as.numeric(data_drug$P67)
data_drug$P68 = as.numeric(data_drug$P68)
data_drug$P69 = as.numeric(data_drug$P69)
data_drug$P70 = as.numeric(data_drug$P70)
data_drug$P71 = as.numeric(data_drug$P71)
data_drug$P72 = as.numeric(data_drug$P72)
data_drug$P73 = as.numeric(data_drug$P73)
data_drug$P74 = as.numeric(data_drug$P74)
data_drug$P75 = as.numeric(data_drug$P75)
data_drug$P76 = as.numeric(data_drug$P76)
data_drug$P77 = as.numeric(data_drug$P77)
data_drug$P78 = as.numeric(data_drug$P78)
data_drug$P79 = as.numeric(data_drug$P79)
data_drug$P80 = as.numeric(data_drug$P80)
data_drug$P81 = as.numeric(data_drug$P81)
data_drug$P82 = as.numeric(data_drug$P82)
data_drug$P83 = as.numeric(data_drug$P83)
data_drug$P84 = as.numeric(data_drug$P84)
data_drug$P85 = as.numeric(data_drug$P85)
data_drug$P86 = as.numeric(data_drug$P86)
data_drug$P87 = as.numeric(data_drug$P87)
data_drug$P88 = as.numeric(data_drug$P88)
data_drug$P89 = as.numeric(data_drug$P89)
data_drug$P90 = as.numeric(data_drug$P90)
data_drug$P91 = as.numeric(data_drug$P91)
data_drug$P92 = as.numeric(data_drug$P92)
data_drug$P93 = as.numeric(data_drug$P93)
data_drug$P94 = as.numeric(data_drug$P94)
data_drug$P95 = as.numeric(data_drug$P95)
data_drug$P96 = as.numeric(data_drug$P96)
data_drug$P97 = as.numeric(data_drug$P97)
data_drug$P98 = as.numeric(data_drug$P98)
data_drug$P99 = as.numeric(data_drug$P99)

#### sample CV

totalAccuracy <- c()
totalSensitivity = c()
totalSpecificity = c()
cvDivider <- floor(nrow(data_drug) / (cv+1))

for (cv in seq(1:cv)) {
  # assign chunk to data test
  dataTestIndex <- c((cv * cvDivider):(cv * cvDivider + cvDivider))
  dataTest <- data_drug[dataTestIndex,]
  # everything else to train
  dataTrain <- data_drug[-dataTestIndex,]
  
  cylModel <- multinom(bin~., data=dataTrain, maxit=1000, trace=T) 
  
  pred <- predict(cylModel, newdata=dataTest, type="class")
  
  #  classification error
  cv_ac <- postResample(dataTest$bin, pred)[[1]]
  cv_sen = sensitivity(as.factor(dataTest$bin), pred)
  cv_spe = specificity(as.factor(dataTest$bin), pred)
  print(paste('Current Accuracy:',cv_ac,'for CV:',cv))
  totalAccuracy <- c(totalAccuracy, cv_ac)
  totalSensitivity = c(totalSensitivity, cv_sen)
  totalSpecificity = c(totalSpecificity,cv_spe)
}

mean(totalAccuracy)
mean(totalSensitivity,na.rm = TRUE)
mean(totalSpecificity,na.rm = TRUE)

## check if data available for both classes

comment = "Model Valid"

if(is.na(mean(totalAccuracy)))
{
  comment = "Model Invalid as data available for only one class during cross validation."
  decision = NA
}


if(!is.na(mean(totalAccuracy)))
{  
  
  #################
  ## main prediction of input
  #################
  
  ## create test seqence with 99 postitons from input, dummy binary, number of mutations and mixtures length
  
  test_seq = as.matrix(c(input[1:99],0,no_mut_input,max_mix_input))
  test_seq = t(test_seq)
  colnames(test_seq) = colnames(data_drug)
  test_seq = data.frame(test_seq)
  
  ## create model on whole data
  
  final_model <- multinom(bin~., data=data_drug, maxit=1000, trace=T)
  
  ## predict for input sequence using model
  
  main_pred = predict(final_model, newdata=test_seq, type="class")
  
  ## formulate string for prediction / remarks
  
  if(main_pred == 1)
  {
    decision = "Resistant"
  }
  if(main_pred == 0)
  {
    decision = "Susceptible"
  }
  
}
## write to result dataframe

final = rbind(final,c("PI","DRV - Darunavir",nrow(data_drug),cutoff,round(mean(totalAccuracy),4)*100,
                      round(mean(totalSensitivity,na.rm = TRUE),4)*100,
                      round(mean(totalSpecificity,na.rm = TRUE),4)*100,decision,comment))



##############
#######
####### PI SQV - Saquinavir workspace
#######
##############

## read the data from csv file

data <- read.csv(pi_dataset,header=TRUE,stringsAsFactors = FALSE)    

## subset data to construct model for FPV

data_drug = subset(data, select = -c(SeqID,NFV,IDV,FPV,LPV,DRV,TPV,ATV))
data_drug = na.omit(data_drug)

## select data as per subtype

if(is.null(subtype))
{
  subtype = "X"
}

## select data correspondin to subtype

if(subtype != "X")
{
  data_drug = data_drug[data_drug$SubType == subtype,]
}

## remove the subtype column now

data_drug = subset(data_drug, select = -c(SubType))

## replace non-mutant positions with 0

data_drug[data_drug == "-"] = 0

## replace the mutant positions with thier integer codes

index = NULL
no_mut = NULL
max_mix = NULL

for (i in 1:nrow(data_drug))
{
  no_mut[i]=0
  max_mix[i]=0
  for (j in 2:NCOL(data_drug))
  {
    x=utf8ToInt(data_drug[i,j])
    sum = 0
    if(length(x)>1)     #### change to 0
    {
      no_mut[i]=no_mut[i]+1
    }
    if(length(x)>max_mix[i])
    {
      max_mix[i]=length(x) 
    }
    for(k in 1:1)
    {
      sum = sum + x[k]
    }
    data_drug[i,j] = as.numeric(sum)
  }
}

## set cut-off manually

cutoff = 3 ##15
data_drug$bin = ifelse(data_drug$SQV >= cutoff,1,0)
data_drug = subset(data_drug, select = -c(SQV))

## add features - no of mutations in sequence, maximum mixture length

data_drug$no_mut = no_mut
data_drug$max_mix = max_mix


## convert to numeric  

data_drug$P1 = as.numeric(data_drug$P1)
data_drug$P2 = as.numeric(data_drug$P2)
data_drug$P3 = as.numeric(data_drug$P3)
data_drug$P4 = as.numeric(data_drug$P4)
data_drug$P5 = as.numeric(data_drug$P5)
data_drug$P6 = as.numeric(data_drug$P6)
data_drug$P7 = as.numeric(data_drug$P7)
data_drug$P8 = as.numeric(data_drug$P8)
data_drug$P9 = as.numeric(data_drug$P9)
data_drug$P10 = as.numeric(data_drug$P10)
data_drug$P11 = as.numeric(data_drug$P11)
data_drug$P12 = as.numeric(data_drug$P12)
data_drug$P13 = as.numeric(data_drug$P13)
data_drug$P14 = as.numeric(data_drug$P14)
data_drug$P15 = as.numeric(data_drug$P15)
data_drug$P16 = as.numeric(data_drug$P16)
data_drug$P17 = as.numeric(data_drug$P17)
data_drug$P18 = as.numeric(data_drug$P18)
data_drug$P19 = as.numeric(data_drug$P19)
data_drug$P20 = as.numeric(data_drug$P20)
data_drug$P21 = as.numeric(data_drug$P21)
data_drug$P22 = as.numeric(data_drug$P22)
data_drug$P23 = as.numeric(data_drug$P23)
data_drug$P24 = as.numeric(data_drug$P24)
data_drug$P25 = as.numeric(data_drug$P25)
data_drug$P26 = as.numeric(data_drug$P26)
data_drug$P27 = as.numeric(data_drug$P27)
data_drug$P28 = as.numeric(data_drug$P28)
data_drug$P29 = as.numeric(data_drug$P29)
data_drug$P30 = as.numeric(data_drug$P30)
data_drug$P31 = as.numeric(data_drug$P31)
data_drug$P32 = as.numeric(data_drug$P32)
data_drug$P33 = as.numeric(data_drug$P33)
data_drug$P34 = as.numeric(data_drug$P34)
data_drug$P35 = as.numeric(data_drug$P35)
data_drug$P36 = as.numeric(data_drug$P36)
data_drug$P37 = as.numeric(data_drug$P37)
data_drug$P38 = as.numeric(data_drug$P38)
data_drug$P39 = as.numeric(data_drug$P39)
data_drug$P40 = as.numeric(data_drug$P40)
data_drug$P41 = as.numeric(data_drug$P41)
data_drug$P42 = as.numeric(data_drug$P42)
data_drug$P43 = as.numeric(data_drug$P43)
data_drug$P44 = as.numeric(data_drug$P44)
data_drug$P45 = as.numeric(data_drug$P45)
data_drug$P46 = as.numeric(data_drug$P46)
data_drug$P47 = as.numeric(data_drug$P47)
data_drug$P48 = as.numeric(data_drug$P48)
data_drug$P49 = as.numeric(data_drug$P49)
data_drug$P50 = as.numeric(data_drug$P50)
data_drug$P51 = as.numeric(data_drug$P51)
data_drug$P52 = as.numeric(data_drug$P52)
data_drug$P53 = as.numeric(data_drug$P53)
data_drug$P54 = as.numeric(data_drug$P54)
data_drug$P55 = as.numeric(data_drug$P55)
data_drug$P56 = as.numeric(data_drug$P56)
data_drug$P57 = as.numeric(data_drug$P57)
data_drug$P58 = as.numeric(data_drug$P58)
data_drug$P59 = as.numeric(data_drug$P59)
data_drug$P60 = as.numeric(data_drug$P60)
data_drug$P61 = as.numeric(data_drug$P61)
data_drug$P62 = as.numeric(data_drug$P62)
data_drug$P63 = as.numeric(data_drug$P63)
data_drug$P64 = as.numeric(data_drug$P64)
data_drug$P65 = as.numeric(data_drug$P65)
data_drug$P66 = as.numeric(data_drug$P66)
data_drug$P67 = as.numeric(data_drug$P67)
data_drug$P68 = as.numeric(data_drug$P68)
data_drug$P69 = as.numeric(data_drug$P69)
data_drug$P70 = as.numeric(data_drug$P70)
data_drug$P71 = as.numeric(data_drug$P71)
data_drug$P72 = as.numeric(data_drug$P72)
data_drug$P73 = as.numeric(data_drug$P73)
data_drug$P74 = as.numeric(data_drug$P74)
data_drug$P75 = as.numeric(data_drug$P75)
data_drug$P76 = as.numeric(data_drug$P76)
data_drug$P77 = as.numeric(data_drug$P77)
data_drug$P78 = as.numeric(data_drug$P78)
data_drug$P79 = as.numeric(data_drug$P79)
data_drug$P80 = as.numeric(data_drug$P80)
data_drug$P81 = as.numeric(data_drug$P81)
data_drug$P82 = as.numeric(data_drug$P82)
data_drug$P83 = as.numeric(data_drug$P83)
data_drug$P84 = as.numeric(data_drug$P84)
data_drug$P85 = as.numeric(data_drug$P85)
data_drug$P86 = as.numeric(data_drug$P86)
data_drug$P87 = as.numeric(data_drug$P87)
data_drug$P88 = as.numeric(data_drug$P88)
data_drug$P89 = as.numeric(data_drug$P89)
data_drug$P90 = as.numeric(data_drug$P90)
data_drug$P91 = as.numeric(data_drug$P91)
data_drug$P92 = as.numeric(data_drug$P92)
data_drug$P93 = as.numeric(data_drug$P93)
data_drug$P94 = as.numeric(data_drug$P94)
data_drug$P95 = as.numeric(data_drug$P95)
data_drug$P96 = as.numeric(data_drug$P96)
data_drug$P97 = as.numeric(data_drug$P97)
data_drug$P98 = as.numeric(data_drug$P98)
data_drug$P99 = as.numeric(data_drug$P99)

#### sample CV

totalAccuracy <- c()
totalSensitivity = c()
totalSpecificity = c()
cvDivider <- floor(nrow(data_drug) / (cv+1))

for (cv in seq(1:cv)) {
  # assign chunk to data test
  dataTestIndex <- c((cv * cvDivider):(cv * cvDivider + cvDivider))
  dataTest <- data_drug[dataTestIndex,]
  # everything else to train
  dataTrain <- data_drug[-dataTestIndex,]
  
  cylModel <- multinom(bin~., data=dataTrain, maxit=1000, trace=T) 
  
  pred <- predict(cylModel, newdata=dataTest, type="class")
  
  #  classification error
  cv_ac <- postResample(dataTest$bin, pred)[[1]]
  cv_sen = sensitivity(as.factor(dataTest$bin), pred)
  cv_spe = specificity(as.factor(dataTest$bin), pred)
  print(paste('Current Accuracy:',cv_ac,'for CV:',cv))
  totalAccuracy <- c(totalAccuracy, cv_ac)
  totalSensitivity = c(totalSensitivity, cv_sen)
  totalSpecificity = c(totalSpecificity,cv_spe)
}

mean(totalAccuracy)
mean(totalSensitivity,na.rm = TRUE)
mean(totalSpecificity,na.rm = TRUE)

## check if data available for both classes

comment = "Model Valid"

if(is.na(mean(totalAccuracy)))
{
  comment = "Model Invalid as data available for only one class during cross validation."
  decision = NA
}


if(!is.na(mean(totalAccuracy)))
{  
  
  #################
  ## main prediction of input
  #################
  
  ## create test seqence with 99 postitons from input, dummy binary, number of mutations and mixtures length
  
  test_seq = as.matrix(c(input[1:99],0,no_mut_input,max_mix_input))
  test_seq = t(test_seq)
  colnames(test_seq) = colnames(data_drug)
  test_seq = data.frame(test_seq)
  
  ## create model on whole data
  
  final_model <- multinom(bin~., data=data_drug, maxit=1000, trace=T)
  
  ## predict for input sequence using model
  
  main_pred = predict(final_model, newdata=test_seq, type="class")
  
  ## formulate string for prediction / remarks
  
  if(main_pred == 1)
  {
    decision = "Resistant"
  }
  if(main_pred == 0)
  {
    decision = "Susceptible"
  }
  
}
## write to result dataframe

final = rbind(final,c("PI","SQV - Saquinavir",nrow(data_drug),cutoff,round(mean(totalAccuracy),4)*100,
                      round(mean(totalSensitivity,na.rm = TRUE),4)*100,
                      round(mean(totalSpecificity,na.rm = TRUE),4)*100,decision,comment))


##############
#######
####### PI IDV - Indinavir workspace
#######
##############


## read the data from csv file

data <- read.csv(pi_dataset,header=TRUE,stringsAsFactors = FALSE)    

## subset data to construct model for FPV

data_drug = subset(data, select = -c(SeqID,NFV,SQV,FPV,LPV,DRV,TPV,ATV))
data_drug = na.omit(data_drug)

## select data as per subtype

if(is.null(subtype))
{
  subtype = "X"
}

## select data correspondin to subtype

if(subtype != "X")
{
  data_drug = data_drug[data_drug$SubType == subtype,]
}

## remove the subtype column now

data_drug = subset(data_drug, select = -c(SubType))

## replace non-mutant positions with 0

data_drug[data_drug == "-"] = 0

## replace the mutant positions with thier integer codes

index = NULL
no_mut = NULL
max_mix = NULL

for (i in 1:nrow(data_drug))
{
  no_mut[i]=0
  max_mix[i]=0
  for (j in 2:NCOL(data_drug))
  {
    x=utf8ToInt(data_drug[i,j])
    sum = 0
    if(length(x)>1)     #### change to 0
    {
      no_mut[i]=no_mut[i]+1
    }
    if(length(x)>max_mix[i])
    {
      max_mix[i]=length(x) 
    }
    for(k in 1:1)
    {
      sum = sum + x[k]
    }
    data_drug[i,j] = as.numeric(sum)
  }
}

## set cut-off manually

cutoff = 3 ##15
data_drug$bin = ifelse(data_drug$IDV >= cutoff,1,0)
data_drug = subset(data_drug, select = -c(IDV))

## add features - no of mutations in sequence, maximum mixture length

data_drug$no_mut = no_mut
data_drug$max_mix = max_mix


## convert to numeric  

data_drug$P1 = as.numeric(data_drug$P1)
data_drug$P2 = as.numeric(data_drug$P2)
data_drug$P3 = as.numeric(data_drug$P3)
data_drug$P4 = as.numeric(data_drug$P4)
data_drug$P5 = as.numeric(data_drug$P5)
data_drug$P6 = as.numeric(data_drug$P6)
data_drug$P7 = as.numeric(data_drug$P7)
data_drug$P8 = as.numeric(data_drug$P8)
data_drug$P9 = as.numeric(data_drug$P9)
data_drug$P10 = as.numeric(data_drug$P10)
data_drug$P11 = as.numeric(data_drug$P11)
data_drug$P12 = as.numeric(data_drug$P12)
data_drug$P13 = as.numeric(data_drug$P13)
data_drug$P14 = as.numeric(data_drug$P14)
data_drug$P15 = as.numeric(data_drug$P15)
data_drug$P16 = as.numeric(data_drug$P16)
data_drug$P17 = as.numeric(data_drug$P17)
data_drug$P18 = as.numeric(data_drug$P18)
data_drug$P19 = as.numeric(data_drug$P19)
data_drug$P20 = as.numeric(data_drug$P20)
data_drug$P21 = as.numeric(data_drug$P21)
data_drug$P22 = as.numeric(data_drug$P22)
data_drug$P23 = as.numeric(data_drug$P23)
data_drug$P24 = as.numeric(data_drug$P24)
data_drug$P25 = as.numeric(data_drug$P25)
data_drug$P26 = as.numeric(data_drug$P26)
data_drug$P27 = as.numeric(data_drug$P27)
data_drug$P28 = as.numeric(data_drug$P28)
data_drug$P29 = as.numeric(data_drug$P29)
data_drug$P30 = as.numeric(data_drug$P30)
data_drug$P31 = as.numeric(data_drug$P31)
data_drug$P32 = as.numeric(data_drug$P32)
data_drug$P33 = as.numeric(data_drug$P33)
data_drug$P34 = as.numeric(data_drug$P34)
data_drug$P35 = as.numeric(data_drug$P35)
data_drug$P36 = as.numeric(data_drug$P36)
data_drug$P37 = as.numeric(data_drug$P37)
data_drug$P38 = as.numeric(data_drug$P38)
data_drug$P39 = as.numeric(data_drug$P39)
data_drug$P40 = as.numeric(data_drug$P40)
data_drug$P41 = as.numeric(data_drug$P41)
data_drug$P42 = as.numeric(data_drug$P42)
data_drug$P43 = as.numeric(data_drug$P43)
data_drug$P44 = as.numeric(data_drug$P44)
data_drug$P45 = as.numeric(data_drug$P45)
data_drug$P46 = as.numeric(data_drug$P46)
data_drug$P47 = as.numeric(data_drug$P47)
data_drug$P48 = as.numeric(data_drug$P48)
data_drug$P49 = as.numeric(data_drug$P49)
data_drug$P50 = as.numeric(data_drug$P50)
data_drug$P51 = as.numeric(data_drug$P51)
data_drug$P52 = as.numeric(data_drug$P52)
data_drug$P53 = as.numeric(data_drug$P53)
data_drug$P54 = as.numeric(data_drug$P54)
data_drug$P55 = as.numeric(data_drug$P55)
data_drug$P56 = as.numeric(data_drug$P56)
data_drug$P57 = as.numeric(data_drug$P57)
data_drug$P58 = as.numeric(data_drug$P58)
data_drug$P59 = as.numeric(data_drug$P59)
data_drug$P60 = as.numeric(data_drug$P60)
data_drug$P61 = as.numeric(data_drug$P61)
data_drug$P62 = as.numeric(data_drug$P62)
data_drug$P63 = as.numeric(data_drug$P63)
data_drug$P64 = as.numeric(data_drug$P64)
data_drug$P65 = as.numeric(data_drug$P65)
data_drug$P66 = as.numeric(data_drug$P66)
data_drug$P67 = as.numeric(data_drug$P67)
data_drug$P68 = as.numeric(data_drug$P68)
data_drug$P69 = as.numeric(data_drug$P69)
data_drug$P70 = as.numeric(data_drug$P70)
data_drug$P71 = as.numeric(data_drug$P71)
data_drug$P72 = as.numeric(data_drug$P72)
data_drug$P73 = as.numeric(data_drug$P73)
data_drug$P74 = as.numeric(data_drug$P74)
data_drug$P75 = as.numeric(data_drug$P75)
data_drug$P76 = as.numeric(data_drug$P76)
data_drug$P77 = as.numeric(data_drug$P77)
data_drug$P78 = as.numeric(data_drug$P78)
data_drug$P79 = as.numeric(data_drug$P79)
data_drug$P80 = as.numeric(data_drug$P80)
data_drug$P81 = as.numeric(data_drug$P81)
data_drug$P82 = as.numeric(data_drug$P82)
data_drug$P83 = as.numeric(data_drug$P83)
data_drug$P84 = as.numeric(data_drug$P84)
data_drug$P85 = as.numeric(data_drug$P85)
data_drug$P86 = as.numeric(data_drug$P86)
data_drug$P87 = as.numeric(data_drug$P87)
data_drug$P88 = as.numeric(data_drug$P88)
data_drug$P89 = as.numeric(data_drug$P89)
data_drug$P90 = as.numeric(data_drug$P90)
data_drug$P91 = as.numeric(data_drug$P91)
data_drug$P92 = as.numeric(data_drug$P92)
data_drug$P93 = as.numeric(data_drug$P93)
data_drug$P94 = as.numeric(data_drug$P94)
data_drug$P95 = as.numeric(data_drug$P95)
data_drug$P96 = as.numeric(data_drug$P96)
data_drug$P97 = as.numeric(data_drug$P97)
data_drug$P98 = as.numeric(data_drug$P98)
data_drug$P99 = as.numeric(data_drug$P99)

#### sample CV

totalAccuracy <- c()
totalSensitivity = c()
totalSpecificity = c()
cvDivider <- floor(nrow(data_drug) / (cv+1))

for (cv in seq(1:cv)) {
  # assign chunk to data test
  dataTestIndex <- c((cv * cvDivider):(cv * cvDivider + cvDivider))
  dataTest <- data_drug[dataTestIndex,]
  # everything else to train
  dataTrain <- data_drug[-dataTestIndex,]
  
  cylModel <- multinom(bin~., data=dataTrain, maxit=1000, trace=T) 
  
  pred <- predict(cylModel, newdata=dataTest, type="class")
  
  #  classification error
  cv_ac <- postResample(dataTest$bin, pred)[[1]]
  cv_sen = sensitivity(as.factor(dataTest$bin), pred)
  cv_spe = specificity(as.factor(dataTest$bin), pred)
  print(paste('Current Accuracy:',cv_ac,'for CV:',cv))
  totalAccuracy <- c(totalAccuracy, cv_ac)
  totalSensitivity = c(totalSensitivity, cv_sen)
  totalSpecificity = c(totalSpecificity,cv_spe)
}

mean(totalAccuracy)
mean(totalSensitivity,na.rm = TRUE)
mean(totalSpecificity,na.rm = TRUE)

## check if data available for both classes

comment = "Model Valid"

if(is.na(mean(totalAccuracy)))
{
  comment = "Model Invalid as data available for only one class during cross validation."
  decision = NA
}


if(!is.na(mean(totalAccuracy)))
{  
  
  #################
  ## main prediction of input
  #################
  
  ## create test seqence with 99 postitons from input, dummy binary, number of mutations and mixtures length
  
  test_seq = as.matrix(c(input[1:99],0,no_mut_input,max_mix_input))
  test_seq = t(test_seq)
  colnames(test_seq) = colnames(data_drug)
  test_seq = data.frame(test_seq)
  
  ## create model on whole data
  
  final_model <- multinom(bin~., data=data_drug, maxit=1000, trace=T)
  
  ## predict for input sequence using model
  
  main_pred = predict(final_model, newdata=test_seq, type="class")
  
  ## formulate string for prediction / remarks
  
  if(main_pred == 1)
  {
    decision = "Resistant"
  }
  if(main_pred == 0)
  {
    decision = "Susceptible"
  }
  
}
## write to result dataframe

final = rbind(final,c("PI","IDV - Indinavir",nrow(data_drug),cutoff,round(mean(totalAccuracy),4)*100,
                      round(mean(totalSensitivity,na.rm = TRUE),4)*100,
                      round(mean(totalSpecificity,na.rm = TRUE),4)*100,decision,comment))


##############
#######
####### PI LPV - Lopinavir workspace
#######
##############


## read the data from csv file

data <- read.csv(pi_dataset,header=TRUE,stringsAsFactors = FALSE)    

## subset data to construct model for FPV

data_drug = subset(data, select = -c(SeqID,NFV,SQV,FPV,IDV,DRV,TPV,ATV))
data_drug = na.omit(data_drug)

## select data as per subtype

if(is.null(subtype))
{
  subtype = "X"
}

## select data correspondin to subtype

if(subtype != "X")
{
  data_drug = data_drug[data_drug$SubType == subtype,]
}

## remove the subtype column now

data_drug = subset(data_drug, select = -c(SubType))

## replace non-mutant positions with 0

data_drug[data_drug == "-"] = 0

## replace the mutant positions with thier integer codes

index = NULL
no_mut = NULL
max_mix = NULL

for (i in 1:nrow(data_drug))
{
  no_mut[i]=0
  max_mix[i]=0
  for (j in 2:NCOL(data_drug))
  {
    x=utf8ToInt(data_drug[i,j])
    sum = 0
    if(length(x)>1)     #### change to 0
    {
      no_mut[i]=no_mut[i]+1
    }
    if(length(x)>max_mix[i])
    {
      max_mix[i]=length(x) 
    }
    for(k in 1:1)
    {
      sum = sum + x[k]
    }
    data_drug[i,j] = as.numeric(sum)
  }
}

## set cut-off manually

cutoff = 9 ##55
data_drug$bin = ifelse(data_drug$LPV >= cutoff,1,0)
data_drug = subset(data_drug, select = -c(LPV))

## add features - no of mutations in sequence, maximum mixture length

data_drug$no_mut = no_mut
data_drug$max_mix = max_mix


## convert to numeric  

data_drug$P1 = as.numeric(data_drug$P1)
data_drug$P2 = as.numeric(data_drug$P2)
data_drug$P3 = as.numeric(data_drug$P3)
data_drug$P4 = as.numeric(data_drug$P4)
data_drug$P5 = as.numeric(data_drug$P5)
data_drug$P6 = as.numeric(data_drug$P6)
data_drug$P7 = as.numeric(data_drug$P7)
data_drug$P8 = as.numeric(data_drug$P8)
data_drug$P9 = as.numeric(data_drug$P9)
data_drug$P10 = as.numeric(data_drug$P10)
data_drug$P11 = as.numeric(data_drug$P11)
data_drug$P12 = as.numeric(data_drug$P12)
data_drug$P13 = as.numeric(data_drug$P13)
data_drug$P14 = as.numeric(data_drug$P14)
data_drug$P15 = as.numeric(data_drug$P15)
data_drug$P16 = as.numeric(data_drug$P16)
data_drug$P17 = as.numeric(data_drug$P17)
data_drug$P18 = as.numeric(data_drug$P18)
data_drug$P19 = as.numeric(data_drug$P19)
data_drug$P20 = as.numeric(data_drug$P20)
data_drug$P21 = as.numeric(data_drug$P21)
data_drug$P22 = as.numeric(data_drug$P22)
data_drug$P23 = as.numeric(data_drug$P23)
data_drug$P24 = as.numeric(data_drug$P24)
data_drug$P25 = as.numeric(data_drug$P25)
data_drug$P26 = as.numeric(data_drug$P26)
data_drug$P27 = as.numeric(data_drug$P27)
data_drug$P28 = as.numeric(data_drug$P28)
data_drug$P29 = as.numeric(data_drug$P29)
data_drug$P30 = as.numeric(data_drug$P30)
data_drug$P31 = as.numeric(data_drug$P31)
data_drug$P32 = as.numeric(data_drug$P32)
data_drug$P33 = as.numeric(data_drug$P33)
data_drug$P34 = as.numeric(data_drug$P34)
data_drug$P35 = as.numeric(data_drug$P35)
data_drug$P36 = as.numeric(data_drug$P36)
data_drug$P37 = as.numeric(data_drug$P37)
data_drug$P38 = as.numeric(data_drug$P38)
data_drug$P39 = as.numeric(data_drug$P39)
data_drug$P40 = as.numeric(data_drug$P40)
data_drug$P41 = as.numeric(data_drug$P41)
data_drug$P42 = as.numeric(data_drug$P42)
data_drug$P43 = as.numeric(data_drug$P43)
data_drug$P44 = as.numeric(data_drug$P44)
data_drug$P45 = as.numeric(data_drug$P45)
data_drug$P46 = as.numeric(data_drug$P46)
data_drug$P47 = as.numeric(data_drug$P47)
data_drug$P48 = as.numeric(data_drug$P48)
data_drug$P49 = as.numeric(data_drug$P49)
data_drug$P50 = as.numeric(data_drug$P50)
data_drug$P51 = as.numeric(data_drug$P51)
data_drug$P52 = as.numeric(data_drug$P52)
data_drug$P53 = as.numeric(data_drug$P53)
data_drug$P54 = as.numeric(data_drug$P54)
data_drug$P55 = as.numeric(data_drug$P55)
data_drug$P56 = as.numeric(data_drug$P56)
data_drug$P57 = as.numeric(data_drug$P57)
data_drug$P58 = as.numeric(data_drug$P58)
data_drug$P59 = as.numeric(data_drug$P59)
data_drug$P60 = as.numeric(data_drug$P60)
data_drug$P61 = as.numeric(data_drug$P61)
data_drug$P62 = as.numeric(data_drug$P62)
data_drug$P63 = as.numeric(data_drug$P63)
data_drug$P64 = as.numeric(data_drug$P64)
data_drug$P65 = as.numeric(data_drug$P65)
data_drug$P66 = as.numeric(data_drug$P66)
data_drug$P67 = as.numeric(data_drug$P67)
data_drug$P68 = as.numeric(data_drug$P68)
data_drug$P69 = as.numeric(data_drug$P69)
data_drug$P70 = as.numeric(data_drug$P70)
data_drug$P71 = as.numeric(data_drug$P71)
data_drug$P72 = as.numeric(data_drug$P72)
data_drug$P73 = as.numeric(data_drug$P73)
data_drug$P74 = as.numeric(data_drug$P74)
data_drug$P75 = as.numeric(data_drug$P75)
data_drug$P76 = as.numeric(data_drug$P76)
data_drug$P77 = as.numeric(data_drug$P77)
data_drug$P78 = as.numeric(data_drug$P78)
data_drug$P79 = as.numeric(data_drug$P79)
data_drug$P80 = as.numeric(data_drug$P80)
data_drug$P81 = as.numeric(data_drug$P81)
data_drug$P82 = as.numeric(data_drug$P82)
data_drug$P83 = as.numeric(data_drug$P83)
data_drug$P84 = as.numeric(data_drug$P84)
data_drug$P85 = as.numeric(data_drug$P85)
data_drug$P86 = as.numeric(data_drug$P86)
data_drug$P87 = as.numeric(data_drug$P87)
data_drug$P88 = as.numeric(data_drug$P88)
data_drug$P89 = as.numeric(data_drug$P89)
data_drug$P90 = as.numeric(data_drug$P90)
data_drug$P91 = as.numeric(data_drug$P91)
data_drug$P92 = as.numeric(data_drug$P92)
data_drug$P93 = as.numeric(data_drug$P93)
data_drug$P94 = as.numeric(data_drug$P94)
data_drug$P95 = as.numeric(data_drug$P95)
data_drug$P96 = as.numeric(data_drug$P96)
data_drug$P97 = as.numeric(data_drug$P97)
data_drug$P98 = as.numeric(data_drug$P98)
data_drug$P99 = as.numeric(data_drug$P99)

#### sample CV

totalAccuracy <- c()
totalSensitivity = c()
totalSpecificity = c()
cvDivider <- floor(nrow(data_drug) / (cv+1))

for (cv in seq(1:cv)) {
  # assign chunk to data test
  dataTestIndex <- c((cv * cvDivider):(cv * cvDivider + cvDivider))
  dataTest <- data_drug[dataTestIndex,]
  # everything else to train
  dataTrain <- data_drug[-dataTestIndex,]
  
  cylModel <- multinom(bin~., data=dataTrain, maxit=1000, trace=T) 
  
  pred <- predict(cylModel, newdata=dataTest, type="class")
  
  #  classification error
  cv_ac <- postResample(dataTest$bin, pred)[[1]]
  cv_sen = sensitivity(as.factor(dataTest$bin), pred)
  cv_spe = specificity(as.factor(dataTest$bin), pred)
  print(paste('Current Accuracy:',cv_ac,'for CV:',cv))
  totalAccuracy <- c(totalAccuracy, cv_ac)
  totalSensitivity = c(totalSensitivity, cv_sen)
  totalSpecificity = c(totalSpecificity,cv_spe)
}

mean(totalAccuracy)
mean(totalSensitivity,na.rm = TRUE)
mean(totalSpecificity,na.rm = TRUE)

## check if data available for both classes

comment = "Model Valid"

if(is.na(mean(totalAccuracy)))
{
  comment = "Model Invalid as data available for only one class during cross validation."
  decision = NA
}


if(!is.na(mean(totalAccuracy)))
{  
  
  #################
  ## main prediction of input
  #################
  
  ## create test seqence with 99 postitons from input, dummy binary, number of mutations and mixtures length
  
  test_seq = as.matrix(c(input[1:99],0,no_mut_input,max_mix_input))
  test_seq = t(test_seq)
  colnames(test_seq) = colnames(data_drug)
  test_seq = data.frame(test_seq)
  
  ## create model on whole data
  
  final_model <- multinom(bin~., data=data_drug, maxit=1000, trace=T)
  
  ## predict for input sequence using model
  
  main_pred = predict(final_model, newdata=test_seq, type="class")
  
  ## formulate string for prediction / remarks
  
  if(main_pred == 1)
  {
    decision = "Resistant"
  }
  if(main_pred == 0)
  {
    decision = "Susceptible"
  }
  
}
## write to result dataframe

final = rbind(final,c("PI","LPV - Lopinavir",nrow(data_drug),cutoff,round(mean(totalAccuracy),4)*100,
                      round(mean(totalSensitivity,na.rm = TRUE),4)*100,
                      round(mean(totalSpecificity,na.rm = TRUE),4)*100,decision,comment))



##############
#######
####### PI TPV - Tipranavir workspace
#######
##############

## read the data from csv file

data <- read.csv(pi_dataset,header=TRUE,stringsAsFactors = FALSE)    

## subset data to construct model for FPV

data_drug = subset(data, select = -c(SeqID,NFV,SQV,FPV,IDV,DRV,LPV,ATV))
data_drug = na.omit(data_drug)

## select data as per subtype

if(is.null(subtype))
{
  subtype = "X"
}

## select data correspondin to subtype

if(subtype != "X")
{
  data_drug = data_drug[data_drug$SubType == subtype,]
}

## remove the subtype column now

data_drug = subset(data_drug, select = -c(SubType))

## replace non-mutant positions with 0

data_drug[data_drug == "-"] = 0

## replace the mutant positions with thier integer codes

index = NULL
no_mut = NULL
max_mix = NULL

for (i in 1:nrow(data_drug))
{
  no_mut[i]=0
  max_mix[i]=0
  for (j in 2:NCOL(data_drug))
  {
    x=utf8ToInt(data_drug[i,j])
    sum = 0
    if(length(x)>1)     #### change to 0
    {
      no_mut[i]=no_mut[i]+1
    }
    if(length(x)>max_mix[i])
    {
      max_mix[i]=length(x) 
    }
    for(k in 1:1)
    {
      sum = sum + x[k]
    }
    data_drug[i,j] = as.numeric(sum)
  }
}

## set cut-off manually

cutoff = 2 ##8
data_drug$bin = ifelse(data_drug$TPV >= cutoff,1,0)
data_drug = subset(data_drug, select = -c(TPV))

## add features - no of mutations in sequence, maximum mixture length

data_drug$no_mut = no_mut
data_drug$max_mix = max_mix


## convert to numeric  

data_drug$P1 = as.numeric(data_drug$P1)
data_drug$P2 = as.numeric(data_drug$P2)
data_drug$P3 = as.numeric(data_drug$P3)
data_drug$P4 = as.numeric(data_drug$P4)
data_drug$P5 = as.numeric(data_drug$P5)
data_drug$P6 = as.numeric(data_drug$P6)
data_drug$P7 = as.numeric(data_drug$P7)
data_drug$P8 = as.numeric(data_drug$P8)
data_drug$P9 = as.numeric(data_drug$P9)
data_drug$P10 = as.numeric(data_drug$P10)
data_drug$P11 = as.numeric(data_drug$P11)
data_drug$P12 = as.numeric(data_drug$P12)
data_drug$P13 = as.numeric(data_drug$P13)
data_drug$P14 = as.numeric(data_drug$P14)
data_drug$P15 = as.numeric(data_drug$P15)
data_drug$P16 = as.numeric(data_drug$P16)
data_drug$P17 = as.numeric(data_drug$P17)
data_drug$P18 = as.numeric(data_drug$P18)
data_drug$P19 = as.numeric(data_drug$P19)
data_drug$P20 = as.numeric(data_drug$P20)
data_drug$P21 = as.numeric(data_drug$P21)
data_drug$P22 = as.numeric(data_drug$P22)
data_drug$P23 = as.numeric(data_drug$P23)
data_drug$P24 = as.numeric(data_drug$P24)
data_drug$P25 = as.numeric(data_drug$P25)
data_drug$P26 = as.numeric(data_drug$P26)
data_drug$P27 = as.numeric(data_drug$P27)
data_drug$P28 = as.numeric(data_drug$P28)
data_drug$P29 = as.numeric(data_drug$P29)
data_drug$P30 = as.numeric(data_drug$P30)
data_drug$P31 = as.numeric(data_drug$P31)
data_drug$P32 = as.numeric(data_drug$P32)
data_drug$P33 = as.numeric(data_drug$P33)
data_drug$P34 = as.numeric(data_drug$P34)
data_drug$P35 = as.numeric(data_drug$P35)
data_drug$P36 = as.numeric(data_drug$P36)
data_drug$P37 = as.numeric(data_drug$P37)
data_drug$P38 = as.numeric(data_drug$P38)
data_drug$P39 = as.numeric(data_drug$P39)
data_drug$P40 = as.numeric(data_drug$P40)
data_drug$P41 = as.numeric(data_drug$P41)
data_drug$P42 = as.numeric(data_drug$P42)
data_drug$P43 = as.numeric(data_drug$P43)
data_drug$P44 = as.numeric(data_drug$P44)
data_drug$P45 = as.numeric(data_drug$P45)
data_drug$P46 = as.numeric(data_drug$P46)
data_drug$P47 = as.numeric(data_drug$P47)
data_drug$P48 = as.numeric(data_drug$P48)
data_drug$P49 = as.numeric(data_drug$P49)
data_drug$P50 = as.numeric(data_drug$P50)
data_drug$P51 = as.numeric(data_drug$P51)
data_drug$P52 = as.numeric(data_drug$P52)
data_drug$P53 = as.numeric(data_drug$P53)
data_drug$P54 = as.numeric(data_drug$P54)
data_drug$P55 = as.numeric(data_drug$P55)
data_drug$P56 = as.numeric(data_drug$P56)
data_drug$P57 = as.numeric(data_drug$P57)
data_drug$P58 = as.numeric(data_drug$P58)
data_drug$P59 = as.numeric(data_drug$P59)
data_drug$P60 = as.numeric(data_drug$P60)
data_drug$P61 = as.numeric(data_drug$P61)
data_drug$P62 = as.numeric(data_drug$P62)
data_drug$P63 = as.numeric(data_drug$P63)
data_drug$P64 = as.numeric(data_drug$P64)
data_drug$P65 = as.numeric(data_drug$P65)
data_drug$P66 = as.numeric(data_drug$P66)
data_drug$P67 = as.numeric(data_drug$P67)
data_drug$P68 = as.numeric(data_drug$P68)
data_drug$P69 = as.numeric(data_drug$P69)
data_drug$P70 = as.numeric(data_drug$P70)
data_drug$P71 = as.numeric(data_drug$P71)
data_drug$P72 = as.numeric(data_drug$P72)
data_drug$P73 = as.numeric(data_drug$P73)
data_drug$P74 = as.numeric(data_drug$P74)
data_drug$P75 = as.numeric(data_drug$P75)
data_drug$P76 = as.numeric(data_drug$P76)
data_drug$P77 = as.numeric(data_drug$P77)
data_drug$P78 = as.numeric(data_drug$P78)
data_drug$P79 = as.numeric(data_drug$P79)
data_drug$P80 = as.numeric(data_drug$P80)
data_drug$P81 = as.numeric(data_drug$P81)
data_drug$P82 = as.numeric(data_drug$P82)
data_drug$P83 = as.numeric(data_drug$P83)
data_drug$P84 = as.numeric(data_drug$P84)
data_drug$P85 = as.numeric(data_drug$P85)
data_drug$P86 = as.numeric(data_drug$P86)
data_drug$P87 = as.numeric(data_drug$P87)
data_drug$P88 = as.numeric(data_drug$P88)
data_drug$P89 = as.numeric(data_drug$P89)
data_drug$P90 = as.numeric(data_drug$P90)
data_drug$P91 = as.numeric(data_drug$P91)
data_drug$P92 = as.numeric(data_drug$P92)
data_drug$P93 = as.numeric(data_drug$P93)
data_drug$P94 = as.numeric(data_drug$P94)
data_drug$P95 = as.numeric(data_drug$P95)
data_drug$P96 = as.numeric(data_drug$P96)
data_drug$P97 = as.numeric(data_drug$P97)
data_drug$P98 = as.numeric(data_drug$P98)
data_drug$P99 = as.numeric(data_drug$P99)

#### sample CV

totalAccuracy <- c()
totalSensitivity = c()
totalSpecificity = c()
cvDivider <- floor(nrow(data_drug) / (cv+1))

for (cv in seq(1:cv)) {
  # assign chunk to data test
  dataTestIndex <- c((cv * cvDivider):(cv * cvDivider + cvDivider))
  dataTest <- data_drug[dataTestIndex,]
  # everything else to train
  dataTrain <- data_drug[-dataTestIndex,]
  
  cylModel <- multinom(bin~., data=dataTrain, maxit=1000, trace=T) 
  
  pred <- predict(cylModel, newdata=dataTest, type="class")
  
  #  classification error
  cv_ac <- postResample(dataTest$bin, pred)[[1]]
  cv_sen = sensitivity(as.factor(dataTest$bin), pred)
  cv_spe = specificity(as.factor(dataTest$bin), pred)
  print(paste('Current Accuracy:',cv_ac,'for CV:',cv))
  totalAccuracy <- c(totalAccuracy, cv_ac)
  totalSensitivity = c(totalSensitivity, cv_sen)
  totalSpecificity = c(totalSpecificity,cv_spe)
}

mean(totalAccuracy)
mean(totalSensitivity,na.rm = TRUE)
mean(totalSpecificity,na.rm = TRUE)

## check if data available for both classes

comment = "Model Valid"

if(is.na(mean(totalAccuracy)))
{
  comment = "Model Invalid as data available for only one class during cross validation."
  decision = NA
}


if(!is.na(mean(totalAccuracy)))
{  
  
  #################
  ## main prediction of input
  #################
  
  ## create test seqence with 99 postitons from input, dummy binary, number of mutations and mixtures length
  
  test_seq = as.matrix(c(input[1:99],0,no_mut_input,max_mix_input))
  test_seq = t(test_seq)
  colnames(test_seq) = colnames(data_drug)
  test_seq = data.frame(test_seq)
  
  ## create model on whole data
  
  final_model <- multinom(bin~., data=data_drug, maxit=1000, trace=T)
  
  ## predict for input sequence using model
  
  main_pred = predict(final_model, newdata=test_seq, type="class")
  
  ## formulate string for prediction / remarks
  
  if(main_pred == 1)
  {
    decision = "Resistant"
  }
  if(main_pred == 0)
  {
    decision = "Susceptible"
  }
  
}
## write to result dataframe

final = rbind(final,c("PI","TPV - Tipranavir",nrow(data_drug),cutoff,round(mean(totalAccuracy),4)*100,
                      round(mean(totalSensitivity,na.rm = TRUE),4)*100,
                      round(mean(totalSpecificity,na.rm = TRUE),4)*100,decision,comment))



##########################################################################################
####
#####                 NRTI DRUGS
####
##########################################################################################


############
############
####### NRTI AZT - Zidovudine (ZDV) / azidothymidine (AZT) workspace


## read the data from csv file

data <- read.csv(nrti_dataset,header=TRUE,stringsAsFactors = FALSE)    

## subset data to construct model for FPV

data_drug = subset(data, select = -c(SeqID,X3TC,ABC,D4T,DDI,TDF))
data_drug = na.omit(data_drug)

## select data as per subtype

if(is.null(subtype))
{
  subtype = "X"
}

## select data correspondin to subtype

if(subtype != "X")
{
  data_drug = data_drug[data_drug$SubType == subtype,]
}

## remove the subtype column now

data_drug = subset(data_drug, select = -c(SubType))

## replace non-mutant positions with 0

data_drug[data_drug == "-"] = 0

## replace the mutant positions with thier integer codes

index = NULL
no_mut = NULL
max_mix = NULL

for (i in 1:nrow(data_drug))
{
  no_mut[i]=0
  max_mix[i]=0
  for (j in 2:NCOL(data_drug))
  {
    x=utf8ToInt(data_drug[i,j])
    sum = 0
    if(length(x)>1)     #### change to 0
    {
      no_mut[i]=no_mut[i]+1
    }
    if(length(x)>max_mix[i])
    {
      max_mix[i]=length(x) 
    }
    for(k in 1:1)
    {
      sum = sum + x[k]
    }
    data_drug[i,j] = as.numeric(sum)
  }
}

## set cut-off manually

cutoff = 3 ##3
data_drug$bin = ifelse(data_drug$AZT >= cutoff,1,0)
data_drug = subset(data_drug, select = -c(AZT))

## add features - no of mutations in sequence, maximum mixture length

data_drug$no_mut = no_mut
data_drug$max_mix = max_mix


## convert to numeric  

data_drug$P1 = as.numeric(data_drug$P1)
data_drug$P2 = as.numeric(data_drug$P2)
data_drug$P3 = as.numeric(data_drug$P3)
data_drug$P4 = as.numeric(data_drug$P4)
data_drug$P5 = as.numeric(data_drug$P5)
data_drug$P6 = as.numeric(data_drug$P6)
data_drug$P7 = as.numeric(data_drug$P7)
data_drug$P8 = as.numeric(data_drug$P8)
data_drug$P9 = as.numeric(data_drug$P9)
data_drug$P10 = as.numeric(data_drug$P10)
data_drug$P11 = as.numeric(data_drug$P11)
data_drug$P12 = as.numeric(data_drug$P12)
data_drug$P13 = as.numeric(data_drug$P13)
data_drug$P14 = as.numeric(data_drug$P14)
data_drug$P15 = as.numeric(data_drug$P15)
data_drug$P16 = as.numeric(data_drug$P16)
data_drug$P17 = as.numeric(data_drug$P17)
data_drug$P18 = as.numeric(data_drug$P18)
data_drug$P19 = as.numeric(data_drug$P19)
data_drug$P20 = as.numeric(data_drug$P20)
data_drug$P21 = as.numeric(data_drug$P21)
data_drug$P22 = as.numeric(data_drug$P22)
data_drug$P23 = as.numeric(data_drug$P23)
data_drug$P24 = as.numeric(data_drug$P24)
data_drug$P25 = as.numeric(data_drug$P25)
data_drug$P26 = as.numeric(data_drug$P26)
data_drug$P27 = as.numeric(data_drug$P27)
data_drug$P28 = as.numeric(data_drug$P28)
data_drug$P29 = as.numeric(data_drug$P29)
data_drug$P30 = as.numeric(data_drug$P30)
data_drug$P31 = as.numeric(data_drug$P31)
data_drug$P32 = as.numeric(data_drug$P32)
data_drug$P33 = as.numeric(data_drug$P33)
data_drug$P34 = as.numeric(data_drug$P34)
data_drug$P35 = as.numeric(data_drug$P35)
data_drug$P36 = as.numeric(data_drug$P36)
data_drug$P37 = as.numeric(data_drug$P37)
data_drug$P38 = as.numeric(data_drug$P38)
data_drug$P39 = as.numeric(data_drug$P39)
data_drug$P40 = as.numeric(data_drug$P40)
data_drug$P41 = as.numeric(data_drug$P41)
data_drug$P42 = as.numeric(data_drug$P42)
data_drug$P43 = as.numeric(data_drug$P43)
data_drug$P44 = as.numeric(data_drug$P44)
data_drug$P45 = as.numeric(data_drug$P45)
data_drug$P46 = as.numeric(data_drug$P46)
data_drug$P47 = as.numeric(data_drug$P47)
data_drug$P48 = as.numeric(data_drug$P48)
data_drug$P49 = as.numeric(data_drug$P49)
data_drug$P50 = as.numeric(data_drug$P50)
data_drug$P51 = as.numeric(data_drug$P51)
data_drug$P52 = as.numeric(data_drug$P52)
data_drug$P53 = as.numeric(data_drug$P53)
data_drug$P54 = as.numeric(data_drug$P54)
data_drug$P55 = as.numeric(data_drug$P55)
data_drug$P56 = as.numeric(data_drug$P56)
data_drug$P57 = as.numeric(data_drug$P57)
data_drug$P58 = as.numeric(data_drug$P58)
data_drug$P59 = as.numeric(data_drug$P59)
data_drug$P60 = as.numeric(data_drug$P60)
data_drug$P61 = as.numeric(data_drug$P61)
data_drug$P62 = as.numeric(data_drug$P62)
data_drug$P63 = as.numeric(data_drug$P63)
data_drug$P64 = as.numeric(data_drug$P64)
data_drug$P65 = as.numeric(data_drug$P65)
data_drug$P66 = as.numeric(data_drug$P66)
data_drug$P67 = as.numeric(data_drug$P67)
data_drug$P68 = as.numeric(data_drug$P68)
data_drug$P69 = as.numeric(data_drug$P69)
data_drug$P70 = as.numeric(data_drug$P70)
data_drug$P71 = as.numeric(data_drug$P71)
data_drug$P72 = as.numeric(data_drug$P72)
data_drug$P73 = as.numeric(data_drug$P73)
data_drug$P74 = as.numeric(data_drug$P74)
data_drug$P75 = as.numeric(data_drug$P75)
data_drug$P76 = as.numeric(data_drug$P76)
data_drug$P77 = as.numeric(data_drug$P77)
data_drug$P78 = as.numeric(data_drug$P78)
data_drug$P79 = as.numeric(data_drug$P79)
data_drug$P80 = as.numeric(data_drug$P80)
data_drug$P81 = as.numeric(data_drug$P81)
data_drug$P82 = as.numeric(data_drug$P82)
data_drug$P83 = as.numeric(data_drug$P83)
data_drug$P84 = as.numeric(data_drug$P84)
data_drug$P85 = as.numeric(data_drug$P85)
data_drug$P86 = as.numeric(data_drug$P86)
data_drug$P87 = as.numeric(data_drug$P87)
data_drug$P88 = as.numeric(data_drug$P88)
data_drug$P89 = as.numeric(data_drug$P89)
data_drug$P90 = as.numeric(data_drug$P90)
data_drug$P91 = as.numeric(data_drug$P91)
data_drug$P92 = as.numeric(data_drug$P92)
data_drug$P93 = as.numeric(data_drug$P93)
data_drug$P94 = as.numeric(data_drug$P94)
data_drug$P95 = as.numeric(data_drug$P95)
data_drug$P96 = as.numeric(data_drug$P96)
data_drug$P97 = as.numeric(data_drug$P97)
data_drug$P98 = as.numeric(data_drug$P98)
data_drug$P99 = as.numeric(data_drug$P99)
data_drug$P100 = as.numeric(data_drug$P100)
data_drug$P101 = as.numeric(data_drug$P101)
data_drug$P102 = as.numeric(data_drug$P102)
data_drug$P103 = as.numeric(data_drug$P103)
data_drug$P104 = as.numeric(data_drug$P104)
data_drug$P105 = as.numeric(data_drug$P105)
data_drug$P106 = as.numeric(data_drug$P106)
data_drug$P107 = as.numeric(data_drug$P107)
data_drug$P108 = as.numeric(data_drug$P108)
data_drug$P109 = as.numeric(data_drug$P109)
data_drug$P110 = as.numeric(data_drug$P110)
data_drug$P111 = as.numeric(data_drug$P111)
data_drug$P112 = as.numeric(data_drug$P112)
data_drug$P113 = as.numeric(data_drug$P113)
data_drug$P114 = as.numeric(data_drug$P114)
data_drug$P115 = as.numeric(data_drug$P115)
data_drug$P116 = as.numeric(data_drug$P116)
data_drug$P117 = as.numeric(data_drug$P117)
data_drug$P118 = as.numeric(data_drug$P118)
data_drug$P119 = as.numeric(data_drug$P119)
data_drug$P120 = as.numeric(data_drug$P120)
data_drug$P121 = as.numeric(data_drug$P121)
data_drug$P122 = as.numeric(data_drug$P122)
data_drug$P123 = as.numeric(data_drug$P123)
data_drug$P124 = as.numeric(data_drug$P124)
data_drug$P125 = as.numeric(data_drug$P125)
data_drug$P126 = as.numeric(data_drug$P126)
data_drug$P127 = as.numeric(data_drug$P127)
data_drug$P128 = as.numeric(data_drug$P128)
data_drug$P129 = as.numeric(data_drug$P129)
data_drug$P130 = as.numeric(data_drug$P130)
data_drug$P131 = as.numeric(data_drug$P131)
data_drug$P132 = as.numeric(data_drug$P132)
data_drug$P133 = as.numeric(data_drug$P133)
data_drug$P134 = as.numeric(data_drug$P134)
data_drug$P135 = as.numeric(data_drug$P135)
data_drug$P136 = as.numeric(data_drug$P136)
data_drug$P137 = as.numeric(data_drug$P137)
data_drug$P138 = as.numeric(data_drug$P138)
data_drug$P139 = as.numeric(data_drug$P139)
data_drug$P140 = as.numeric(data_drug$P140)
data_drug$P141 = as.numeric(data_drug$P141)
data_drug$P142 = as.numeric(data_drug$P142)
data_drug$P143 = as.numeric(data_drug$P143)
data_drug$P144 = as.numeric(data_drug$P144)
data_drug$P145 = as.numeric(data_drug$P145)
data_drug$P146 = as.numeric(data_drug$P146)
data_drug$P147 = as.numeric(data_drug$P147)
data_drug$P148 = as.numeric(data_drug$P148)
data_drug$P149 = as.numeric(data_drug$P149)
data_drug$P150 = as.numeric(data_drug$P150)
data_drug$P151 = as.numeric(data_drug$P151)
data_drug$P152 = as.numeric(data_drug$P152)
data_drug$P153 = as.numeric(data_drug$P153)
data_drug$P154 = as.numeric(data_drug$P154)
data_drug$P155 = as.numeric(data_drug$P155)
data_drug$P156 = as.numeric(data_drug$P156)
data_drug$P157 = as.numeric(data_drug$P157)
data_drug$P158 = as.numeric(data_drug$P158)
data_drug$P159 = as.numeric(data_drug$P159)
data_drug$P160 = as.numeric(data_drug$P160)
data_drug$P161 = as.numeric(data_drug$P161)
data_drug$P162 = as.numeric(data_drug$P162)
data_drug$P163 = as.numeric(data_drug$P163)
data_drug$P164 = as.numeric(data_drug$P164)
data_drug$P165 = as.numeric(data_drug$P165)
data_drug$P166 = as.numeric(data_drug$P166)
data_drug$P167 = as.numeric(data_drug$P167)
data_drug$P168 = as.numeric(data_drug$P168)
data_drug$P169 = as.numeric(data_drug$P169)
data_drug$P170 = as.numeric(data_drug$P170)
data_drug$P171 = as.numeric(data_drug$P171)
data_drug$P172 = as.numeric(data_drug$P172)
data_drug$P173 = as.numeric(data_drug$P173)
data_drug$P174 = as.numeric(data_drug$P174)
data_drug$P175 = as.numeric(data_drug$P175)
data_drug$P176 = as.numeric(data_drug$P176)
data_drug$P177 = as.numeric(data_drug$P177)
data_drug$P178 = as.numeric(data_drug$P178)
data_drug$P179 = as.numeric(data_drug$P179)
data_drug$P180 = as.numeric(data_drug$P180)
data_drug$P181 = as.numeric(data_drug$P181)
data_drug$P182 = as.numeric(data_drug$P182)
data_drug$P183 = as.numeric(data_drug$P183)
data_drug$P184 = as.numeric(data_drug$P184)
data_drug$P185 = as.numeric(data_drug$P185)
data_drug$P186 = as.numeric(data_drug$P186)
data_drug$P187 = as.numeric(data_drug$P187)
data_drug$P188 = as.numeric(data_drug$P188)
data_drug$P189 = as.numeric(data_drug$P189)
data_drug$P190 = as.numeric(data_drug$P190)
data_drug$P191 = as.numeric(data_drug$P191)
data_drug$P192 = as.numeric(data_drug$P192)
data_drug$P193 = as.numeric(data_drug$P193)
data_drug$P194 = as.numeric(data_drug$P194)
data_drug$P195 = as.numeric(data_drug$P195)
data_drug$P196 = as.numeric(data_drug$P196)
data_drug$P197 = as.numeric(data_drug$P197)
data_drug$P198 = as.numeric(data_drug$P198)
data_drug$P199 = as.numeric(data_drug$P199)
data_drug$P200 = as.numeric(data_drug$P200)
data_drug$P201 = as.numeric(data_drug$P201)
data_drug$P202 = as.numeric(data_drug$P202)
data_drug$P203 = as.numeric(data_drug$P203)
data_drug$P204 = as.numeric(data_drug$P204)
data_drug$P205 = as.numeric(data_drug$P205)
data_drug$P206 = as.numeric(data_drug$P206)
data_drug$P207 = as.numeric(data_drug$P207)
data_drug$P208 = as.numeric(data_drug$P208)
data_drug$P209 = as.numeric(data_drug$P209)
data_drug$P210 = as.numeric(data_drug$P210)
data_drug$P211 = as.numeric(data_drug$P211)
data_drug$P212 = as.numeric(data_drug$P212)
data_drug$P213 = as.numeric(data_drug$P213)
data_drug$P214 = as.numeric(data_drug$P214)
data_drug$P215 = as.numeric(data_drug$P215)
data_drug$P216 = as.numeric(data_drug$P216)
data_drug$P217 = as.numeric(data_drug$P217)
data_drug$P218 = as.numeric(data_drug$P218)
data_drug$P219 = as.numeric(data_drug$P219)
data_drug$P220 = as.numeric(data_drug$P220)
data_drug$P221 = as.numeric(data_drug$P221)
data_drug$P222 = as.numeric(data_drug$P222)
data_drug$P223 = as.numeric(data_drug$P223)
data_drug$P224 = as.numeric(data_drug$P224)
data_drug$P225 = as.numeric(data_drug$P225)
data_drug$P226 = as.numeric(data_drug$P226)
data_drug$P227 = as.numeric(data_drug$P227)
data_drug$P228 = as.numeric(data_drug$P228)
data_drug$P229 = as.numeric(data_drug$P229)
data_drug$P230 = as.numeric(data_drug$P230)
data_drug$P231 = as.numeric(data_drug$P231)
data_drug$P232 = as.numeric(data_drug$P232)
data_drug$P233 = as.numeric(data_drug$P233)
data_drug$P234 = as.numeric(data_drug$P234)
data_drug$P235 = as.numeric(data_drug$P235)
data_drug$P236 = as.numeric(data_drug$P236)
data_drug$P237 = as.numeric(data_drug$P237)
data_drug$P238 = as.numeric(data_drug$P238)
data_drug$P239 = as.numeric(data_drug$P239)
data_drug$P240 = as.numeric(data_drug$P240)

#### sample CV

totalAccuracy <- c()
totalSensitivity = c()
totalSpecificity = c()
cvDivider <- floor(nrow(data_drug) / (cv+1))

for (cv in seq(1:cv)) {
  # assign chunk to data test
  dataTestIndex <- c((cv * cvDivider):(cv * cvDivider + cvDivider))
  dataTest <- data_drug[dataTestIndex,]
  # everything else to train
  dataTrain <- data_drug[-dataTestIndex,]
  
  cylModel <- multinom(bin~., data=dataTrain, maxit=1000, trace=T) 
  
  pred <- predict(cylModel, newdata=dataTest, type="class")
  
  #  classification error
  cv_ac <- postResample(dataTest$bin, pred)[[1]]
  cv_sen = sensitivity(as.factor(dataTest$bin), pred)
  cv_spe = specificity(as.factor(dataTest$bin), pred)
  print(paste('Current Accuracy:',cv_ac,'for CV:',cv))
  totalAccuracy <- c(totalAccuracy, cv_ac)
  totalSensitivity = c(totalSensitivity, cv_sen)
  totalSpecificity = c(totalSpecificity,cv_spe)
}

mean(totalAccuracy)
mean(totalSensitivity,na.rm = TRUE)
mean(totalSpecificity,na.rm = TRUE)

## check if data available for both classes

comment = "Model Valid"

if(is.na(mean(totalAccuracy)))
{
  comment = "Model Invalid as data available for only one class during cross validation."
  decision = NA
}


if(!is.na(mean(totalAccuracy)))
{  
  
  #################
  ## main prediction of input
  #################
  
  ## create test seqence with 240 postitons from input, dummy binary, number of mutations and mixtures length
  
  test_seq = as.matrix(c(input[1:240],0,no_mut_input,max_mix_input))
  test_seq = t(test_seq)
  colnames(test_seq) = colnames(data_drug)
  test_seq = data.frame(test_seq)
  
  ## create model on whole data
  
  final_model <- multinom(bin~., data=data_drug, maxit=1000, trace=T)
  
  ## predict for input sequence using model
  
  main_pred = predict(final_model, newdata=test_seq, type="class")
  
  ## formulate string for prediction / remarks
  
  if(main_pred == 1)
  {
    decision = "Resistant"
  }
  if(main_pred == 0)
  {
    decision = "Susceptible"
  }
  
}
## write to result dataframe

final = rbind(final,c("NRTI","AZT - Zidovudine",nrow(data_drug),cutoff,round(mean(totalAccuracy),4)*100,
                      round(mean(totalSensitivity,na.rm = TRUE),4)*100,
                      round(mean(totalSpecificity,na.rm = TRUE),4)*100,decision,comment))



############
############
####### NRTI D4T - Stavudine workspace


## read the data from csv file

data <- read.csv(nrti_dataset,header=TRUE,stringsAsFactors = FALSE)    

## subset data to construct model for FPV

data_drug = subset(data, select = -c(SeqID,X3TC,ABC,AZT,DDI,TDF))
data_drug = na.omit(data_drug)

## select data as per subtype

if(is.null(subtype))
{
  subtype = "X"
}

## select data correspondin to subtype

if(subtype != "X")
{
  data_drug = data_drug[data_drug$SubType == subtype,]
}

## remove the subtype column now

data_drug = subset(data_drug, select = -c(SubType))

## replace non-mutant positions with 0

data_drug[data_drug == "-"] = 0

## replace the mutant positions with thier integer codes

index = NULL
no_mut = NULL
max_mix = NULL

for (i in 1:nrow(data_drug))
{
  no_mut[i]=0
  max_mix[i]=0
  for (j in 2:NCOL(data_drug))
  {
    x=utf8ToInt(data_drug[i,j])
    sum = 0
    if(length(x)>1)     #### change to 0
    {
      no_mut[i]=no_mut[i]+1
    }
    if(length(x)>max_mix[i])
    {
      max_mix[i]=length(x) 
    }
    for(k in 1:1)
    {
      sum = sum + x[k]
    }
    data_drug[i,j] = as.numeric(sum)
  }
}

## set cut-off manually

cutoff = 1.5 ##1.5
data_drug$bin = ifelse(data_drug$D4T >= cutoff,1,0)
data_drug = subset(data_drug, select = -c(D4T))


## add features - no of mutations in sequence, maximum mixture length

data_drug$no_mut = no_mut
data_drug$max_mix = max_mix


## convert to numeric  

data_drug$P1 = as.numeric(data_drug$P1)
data_drug$P2 = as.numeric(data_drug$P2)
data_drug$P3 = as.numeric(data_drug$P3)
data_drug$P4 = as.numeric(data_drug$P4)
data_drug$P5 = as.numeric(data_drug$P5)
data_drug$P6 = as.numeric(data_drug$P6)
data_drug$P7 = as.numeric(data_drug$P7)
data_drug$P8 = as.numeric(data_drug$P8)
data_drug$P9 = as.numeric(data_drug$P9)
data_drug$P10 = as.numeric(data_drug$P10)
data_drug$P11 = as.numeric(data_drug$P11)
data_drug$P12 = as.numeric(data_drug$P12)
data_drug$P13 = as.numeric(data_drug$P13)
data_drug$P14 = as.numeric(data_drug$P14)
data_drug$P15 = as.numeric(data_drug$P15)
data_drug$P16 = as.numeric(data_drug$P16)
data_drug$P17 = as.numeric(data_drug$P17)
data_drug$P18 = as.numeric(data_drug$P18)
data_drug$P19 = as.numeric(data_drug$P19)
data_drug$P20 = as.numeric(data_drug$P20)
data_drug$P21 = as.numeric(data_drug$P21)
data_drug$P22 = as.numeric(data_drug$P22)
data_drug$P23 = as.numeric(data_drug$P23)
data_drug$P24 = as.numeric(data_drug$P24)
data_drug$P25 = as.numeric(data_drug$P25)
data_drug$P26 = as.numeric(data_drug$P26)
data_drug$P27 = as.numeric(data_drug$P27)
data_drug$P28 = as.numeric(data_drug$P28)
data_drug$P29 = as.numeric(data_drug$P29)
data_drug$P30 = as.numeric(data_drug$P30)
data_drug$P31 = as.numeric(data_drug$P31)
data_drug$P32 = as.numeric(data_drug$P32)
data_drug$P33 = as.numeric(data_drug$P33)
data_drug$P34 = as.numeric(data_drug$P34)
data_drug$P35 = as.numeric(data_drug$P35)
data_drug$P36 = as.numeric(data_drug$P36)
data_drug$P37 = as.numeric(data_drug$P37)
data_drug$P38 = as.numeric(data_drug$P38)
data_drug$P39 = as.numeric(data_drug$P39)
data_drug$P40 = as.numeric(data_drug$P40)
data_drug$P41 = as.numeric(data_drug$P41)
data_drug$P42 = as.numeric(data_drug$P42)
data_drug$P43 = as.numeric(data_drug$P43)
data_drug$P44 = as.numeric(data_drug$P44)
data_drug$P45 = as.numeric(data_drug$P45)
data_drug$P46 = as.numeric(data_drug$P46)
data_drug$P47 = as.numeric(data_drug$P47)
data_drug$P48 = as.numeric(data_drug$P48)
data_drug$P49 = as.numeric(data_drug$P49)
data_drug$P50 = as.numeric(data_drug$P50)
data_drug$P51 = as.numeric(data_drug$P51)
data_drug$P52 = as.numeric(data_drug$P52)
data_drug$P53 = as.numeric(data_drug$P53)
data_drug$P54 = as.numeric(data_drug$P54)
data_drug$P55 = as.numeric(data_drug$P55)
data_drug$P56 = as.numeric(data_drug$P56)
data_drug$P57 = as.numeric(data_drug$P57)
data_drug$P58 = as.numeric(data_drug$P58)
data_drug$P59 = as.numeric(data_drug$P59)
data_drug$P60 = as.numeric(data_drug$P60)
data_drug$P61 = as.numeric(data_drug$P61)
data_drug$P62 = as.numeric(data_drug$P62)
data_drug$P63 = as.numeric(data_drug$P63)
data_drug$P64 = as.numeric(data_drug$P64)
data_drug$P65 = as.numeric(data_drug$P65)
data_drug$P66 = as.numeric(data_drug$P66)
data_drug$P67 = as.numeric(data_drug$P67)
data_drug$P68 = as.numeric(data_drug$P68)
data_drug$P69 = as.numeric(data_drug$P69)
data_drug$P70 = as.numeric(data_drug$P70)
data_drug$P71 = as.numeric(data_drug$P71)
data_drug$P72 = as.numeric(data_drug$P72)
data_drug$P73 = as.numeric(data_drug$P73)
data_drug$P74 = as.numeric(data_drug$P74)
data_drug$P75 = as.numeric(data_drug$P75)
data_drug$P76 = as.numeric(data_drug$P76)
data_drug$P77 = as.numeric(data_drug$P77)
data_drug$P78 = as.numeric(data_drug$P78)
data_drug$P79 = as.numeric(data_drug$P79)
data_drug$P80 = as.numeric(data_drug$P80)
data_drug$P81 = as.numeric(data_drug$P81)
data_drug$P82 = as.numeric(data_drug$P82)
data_drug$P83 = as.numeric(data_drug$P83)
data_drug$P84 = as.numeric(data_drug$P84)
data_drug$P85 = as.numeric(data_drug$P85)
data_drug$P86 = as.numeric(data_drug$P86)
data_drug$P87 = as.numeric(data_drug$P87)
data_drug$P88 = as.numeric(data_drug$P88)
data_drug$P89 = as.numeric(data_drug$P89)
data_drug$P90 = as.numeric(data_drug$P90)
data_drug$P91 = as.numeric(data_drug$P91)
data_drug$P92 = as.numeric(data_drug$P92)
data_drug$P93 = as.numeric(data_drug$P93)
data_drug$P94 = as.numeric(data_drug$P94)
data_drug$P95 = as.numeric(data_drug$P95)
data_drug$P96 = as.numeric(data_drug$P96)
data_drug$P97 = as.numeric(data_drug$P97)
data_drug$P98 = as.numeric(data_drug$P98)
data_drug$P99 = as.numeric(data_drug$P99)
data_drug$P100 = as.numeric(data_drug$P100)
data_drug$P101 = as.numeric(data_drug$P101)
data_drug$P102 = as.numeric(data_drug$P102)
data_drug$P103 = as.numeric(data_drug$P103)
data_drug$P104 = as.numeric(data_drug$P104)
data_drug$P105 = as.numeric(data_drug$P105)
data_drug$P106 = as.numeric(data_drug$P106)
data_drug$P107 = as.numeric(data_drug$P107)
data_drug$P108 = as.numeric(data_drug$P108)
data_drug$P109 = as.numeric(data_drug$P109)
data_drug$P110 = as.numeric(data_drug$P110)
data_drug$P111 = as.numeric(data_drug$P111)
data_drug$P112 = as.numeric(data_drug$P112)
data_drug$P113 = as.numeric(data_drug$P113)
data_drug$P114 = as.numeric(data_drug$P114)
data_drug$P115 = as.numeric(data_drug$P115)
data_drug$P116 = as.numeric(data_drug$P116)
data_drug$P117 = as.numeric(data_drug$P117)
data_drug$P118 = as.numeric(data_drug$P118)
data_drug$P119 = as.numeric(data_drug$P119)
data_drug$P120 = as.numeric(data_drug$P120)
data_drug$P121 = as.numeric(data_drug$P121)
data_drug$P122 = as.numeric(data_drug$P122)
data_drug$P123 = as.numeric(data_drug$P123)
data_drug$P124 = as.numeric(data_drug$P124)
data_drug$P125 = as.numeric(data_drug$P125)
data_drug$P126 = as.numeric(data_drug$P126)
data_drug$P127 = as.numeric(data_drug$P127)
data_drug$P128 = as.numeric(data_drug$P128)
data_drug$P129 = as.numeric(data_drug$P129)
data_drug$P130 = as.numeric(data_drug$P130)
data_drug$P131 = as.numeric(data_drug$P131)
data_drug$P132 = as.numeric(data_drug$P132)
data_drug$P133 = as.numeric(data_drug$P133)
data_drug$P134 = as.numeric(data_drug$P134)
data_drug$P135 = as.numeric(data_drug$P135)
data_drug$P136 = as.numeric(data_drug$P136)
data_drug$P137 = as.numeric(data_drug$P137)
data_drug$P138 = as.numeric(data_drug$P138)
data_drug$P139 = as.numeric(data_drug$P139)
data_drug$P140 = as.numeric(data_drug$P140)
data_drug$P141 = as.numeric(data_drug$P141)
data_drug$P142 = as.numeric(data_drug$P142)
data_drug$P143 = as.numeric(data_drug$P143)
data_drug$P144 = as.numeric(data_drug$P144)
data_drug$P145 = as.numeric(data_drug$P145)
data_drug$P146 = as.numeric(data_drug$P146)
data_drug$P147 = as.numeric(data_drug$P147)
data_drug$P148 = as.numeric(data_drug$P148)
data_drug$P149 = as.numeric(data_drug$P149)
data_drug$P150 = as.numeric(data_drug$P150)
data_drug$P151 = as.numeric(data_drug$P151)
data_drug$P152 = as.numeric(data_drug$P152)
data_drug$P153 = as.numeric(data_drug$P153)
data_drug$P154 = as.numeric(data_drug$P154)
data_drug$P155 = as.numeric(data_drug$P155)
data_drug$P156 = as.numeric(data_drug$P156)
data_drug$P157 = as.numeric(data_drug$P157)
data_drug$P158 = as.numeric(data_drug$P158)
data_drug$P159 = as.numeric(data_drug$P159)
data_drug$P160 = as.numeric(data_drug$P160)
data_drug$P161 = as.numeric(data_drug$P161)
data_drug$P162 = as.numeric(data_drug$P162)
data_drug$P163 = as.numeric(data_drug$P163)
data_drug$P164 = as.numeric(data_drug$P164)
data_drug$P165 = as.numeric(data_drug$P165)
data_drug$P166 = as.numeric(data_drug$P166)
data_drug$P167 = as.numeric(data_drug$P167)
data_drug$P168 = as.numeric(data_drug$P168)
data_drug$P169 = as.numeric(data_drug$P169)
data_drug$P170 = as.numeric(data_drug$P170)
data_drug$P171 = as.numeric(data_drug$P171)
data_drug$P172 = as.numeric(data_drug$P172)
data_drug$P173 = as.numeric(data_drug$P173)
data_drug$P174 = as.numeric(data_drug$P174)
data_drug$P175 = as.numeric(data_drug$P175)
data_drug$P176 = as.numeric(data_drug$P176)
data_drug$P177 = as.numeric(data_drug$P177)
data_drug$P178 = as.numeric(data_drug$P178)
data_drug$P179 = as.numeric(data_drug$P179)
data_drug$P180 = as.numeric(data_drug$P180)
data_drug$P181 = as.numeric(data_drug$P181)
data_drug$P182 = as.numeric(data_drug$P182)
data_drug$P183 = as.numeric(data_drug$P183)
data_drug$P184 = as.numeric(data_drug$P184)
data_drug$P185 = as.numeric(data_drug$P185)
data_drug$P186 = as.numeric(data_drug$P186)
data_drug$P187 = as.numeric(data_drug$P187)
data_drug$P188 = as.numeric(data_drug$P188)
data_drug$P189 = as.numeric(data_drug$P189)
data_drug$P190 = as.numeric(data_drug$P190)
data_drug$P191 = as.numeric(data_drug$P191)
data_drug$P192 = as.numeric(data_drug$P192)
data_drug$P193 = as.numeric(data_drug$P193)
data_drug$P194 = as.numeric(data_drug$P194)
data_drug$P195 = as.numeric(data_drug$P195)
data_drug$P196 = as.numeric(data_drug$P196)
data_drug$P197 = as.numeric(data_drug$P197)
data_drug$P198 = as.numeric(data_drug$P198)
data_drug$P199 = as.numeric(data_drug$P199)
data_drug$P200 = as.numeric(data_drug$P200)
data_drug$P201 = as.numeric(data_drug$P201)
data_drug$P202 = as.numeric(data_drug$P202)
data_drug$P203 = as.numeric(data_drug$P203)
data_drug$P204 = as.numeric(data_drug$P204)
data_drug$P205 = as.numeric(data_drug$P205)
data_drug$P206 = as.numeric(data_drug$P206)
data_drug$P207 = as.numeric(data_drug$P207)
data_drug$P208 = as.numeric(data_drug$P208)
data_drug$P209 = as.numeric(data_drug$P209)
data_drug$P210 = as.numeric(data_drug$P210)
data_drug$P211 = as.numeric(data_drug$P211)
data_drug$P212 = as.numeric(data_drug$P212)
data_drug$P213 = as.numeric(data_drug$P213)
data_drug$P214 = as.numeric(data_drug$P214)
data_drug$P215 = as.numeric(data_drug$P215)
data_drug$P216 = as.numeric(data_drug$P216)
data_drug$P217 = as.numeric(data_drug$P217)
data_drug$P218 = as.numeric(data_drug$P218)
data_drug$P219 = as.numeric(data_drug$P219)
data_drug$P220 = as.numeric(data_drug$P220)
data_drug$P221 = as.numeric(data_drug$P221)
data_drug$P222 = as.numeric(data_drug$P222)
data_drug$P223 = as.numeric(data_drug$P223)
data_drug$P224 = as.numeric(data_drug$P224)
data_drug$P225 = as.numeric(data_drug$P225)
data_drug$P226 = as.numeric(data_drug$P226)
data_drug$P227 = as.numeric(data_drug$P227)
data_drug$P228 = as.numeric(data_drug$P228)
data_drug$P229 = as.numeric(data_drug$P229)
data_drug$P230 = as.numeric(data_drug$P230)
data_drug$P231 = as.numeric(data_drug$P231)
data_drug$P232 = as.numeric(data_drug$P232)
data_drug$P233 = as.numeric(data_drug$P233)
data_drug$P234 = as.numeric(data_drug$P234)
data_drug$P235 = as.numeric(data_drug$P235)
data_drug$P236 = as.numeric(data_drug$P236)
data_drug$P237 = as.numeric(data_drug$P237)
data_drug$P238 = as.numeric(data_drug$P238)
data_drug$P239 = as.numeric(data_drug$P239)
data_drug$P240 = as.numeric(data_drug$P240)

#### sample CV

totalAccuracy <- c()
totalSensitivity = c()
totalSpecificity = c()
cvDivider <- floor(nrow(data_drug) / (cv+1))

for (cv in seq(1:cv)) {
  # assign chunk to data test
  dataTestIndex <- c((cv * cvDivider):(cv * cvDivider + cvDivider))
  dataTest <- data_drug[dataTestIndex,]
  # everything else to train
  dataTrain <- data_drug[-dataTestIndex,]
  
  cylModel <- multinom(bin~., data=dataTrain, maxit=1000, trace=T) 
  
  pred <- predict(cylModel, newdata=dataTest, type="class")
  
  #  classification error
  cv_ac <- postResample(dataTest$bin, pred)[[1]]
  cv_sen = sensitivity(as.factor(dataTest$bin), pred)
  cv_spe = specificity(as.factor(dataTest$bin), pred)
  print(paste('Current Accuracy:',cv_ac,'for CV:',cv))
  totalAccuracy <- c(totalAccuracy, cv_ac)
  totalSensitivity = c(totalSensitivity, cv_sen)
  totalSpecificity = c(totalSpecificity,cv_spe)
}

mean(totalAccuracy)
mean(totalSensitivity,na.rm = TRUE)
mean(totalSpecificity,na.rm = TRUE)

## check if data available for both classes

comment = "Model Valid"

if(is.na(mean(totalAccuracy)))
{
  comment = "Model Invalid as data available for only one class during cross validation."
  decision = NA
}


if(!is.na(mean(totalAccuracy)))
{  
  
  #################
  ## main prediction of input
  #################
  
  ## create test seqence with 240 postitons from input, dummy binary, number of mutations and mixtures length
  
  test_seq = as.matrix(c(input[1:240],0,no_mut_input,max_mix_input))
  test_seq = t(test_seq)
  colnames(test_seq) = colnames(data_drug)
  test_seq = data.frame(test_seq)
  
  ## create model on whole data
  
  final_model <- multinom(bin~., data=data_drug, maxit=1000, trace=T)
  
  ## predict for input sequence using model
  
  main_pred = predict(final_model, newdata=test_seq, type="class")
  
  ## formulate string for prediction / remarks
  
  if(main_pred == 1)
  {
    decision = "Resistant"
  }
  if(main_pred == 0)
  {
    decision = "Susceptible"
  }
  
}
## write to result dataframe

final = rbind(final,c("NRTI","D4T - Stavudine",nrow(data_drug),cutoff,round(mean(totalAccuracy),4)*100,
                      round(mean(totalSensitivity,na.rm = TRUE),4)*100,
                      round(mean(totalSpecificity,na.rm = TRUE),4)*100,decision,comment))



############
############
####### NRTI TDF - Tenofovir workspace 


## read the data from csv file

data <- read.csv(nrti_dataset,header=TRUE,stringsAsFactors = FALSE)    

## subset data to construct model for FPV

data_drug = subset(data, select = -c(SeqID,X3TC,ABC,AZT,DDI,D4T))
data_drug = na.omit(data_drug)

## select data as per subtype

if(is.null(subtype))
{
  subtype = "X"
}

## select data correspondin to subtype

if(subtype != "X")
{
  data_drug = data_drug[data_drug$SubType == subtype,]
}

## remove the subtype column now

data_drug = subset(data_drug, select = -c(SubType))

## replace non-mutant positions with 0

data_drug[data_drug == "-"] = 0

## replace the mutant positions with thier integer codes

index = NULL
no_mut = NULL
max_mix = NULL

for (i in 1:nrow(data_drug))
{
  no_mut[i]=0
  max_mix[i]=0
  for (j in 2:NCOL(data_drug))
  {
    x=utf8ToInt(data_drug[i,j])
    sum = 0
    if(length(x)>1)     #### change to 0
    {
      no_mut[i]=no_mut[i]+1
    }
    if(length(x)>max_mix[i])
    {
      max_mix[i]=length(x) 
    }
    for(k in 1:1)
    {
      sum = sum + x[k]
    }
    data_drug[i,j] = as.numeric(sum)
  }
}

## set cut-off manually

cutoff = 1.5 ##1.5
data_drug$bin = ifelse(data_drug$TDF >= cutoff,1,0)
data_drug = subset(data_drug, select = -c(TDF))

## add features - no of mutations in sequence, maximum mixture length

data_drug$no_mut = no_mut
data_drug$max_mix = max_mix


## convert to numeric  

data_drug$P1 = as.numeric(data_drug$P1)
data_drug$P2 = as.numeric(data_drug$P2)
data_drug$P3 = as.numeric(data_drug$P3)
data_drug$P4 = as.numeric(data_drug$P4)
data_drug$P5 = as.numeric(data_drug$P5)
data_drug$P6 = as.numeric(data_drug$P6)
data_drug$P7 = as.numeric(data_drug$P7)
data_drug$P8 = as.numeric(data_drug$P8)
data_drug$P9 = as.numeric(data_drug$P9)
data_drug$P10 = as.numeric(data_drug$P10)
data_drug$P11 = as.numeric(data_drug$P11)
data_drug$P12 = as.numeric(data_drug$P12)
data_drug$P13 = as.numeric(data_drug$P13)
data_drug$P14 = as.numeric(data_drug$P14)
data_drug$P15 = as.numeric(data_drug$P15)
data_drug$P16 = as.numeric(data_drug$P16)
data_drug$P17 = as.numeric(data_drug$P17)
data_drug$P18 = as.numeric(data_drug$P18)
data_drug$P19 = as.numeric(data_drug$P19)
data_drug$P20 = as.numeric(data_drug$P20)
data_drug$P21 = as.numeric(data_drug$P21)
data_drug$P22 = as.numeric(data_drug$P22)
data_drug$P23 = as.numeric(data_drug$P23)
data_drug$P24 = as.numeric(data_drug$P24)
data_drug$P25 = as.numeric(data_drug$P25)
data_drug$P26 = as.numeric(data_drug$P26)
data_drug$P27 = as.numeric(data_drug$P27)
data_drug$P28 = as.numeric(data_drug$P28)
data_drug$P29 = as.numeric(data_drug$P29)
data_drug$P30 = as.numeric(data_drug$P30)
data_drug$P31 = as.numeric(data_drug$P31)
data_drug$P32 = as.numeric(data_drug$P32)
data_drug$P33 = as.numeric(data_drug$P33)
data_drug$P34 = as.numeric(data_drug$P34)
data_drug$P35 = as.numeric(data_drug$P35)
data_drug$P36 = as.numeric(data_drug$P36)
data_drug$P37 = as.numeric(data_drug$P37)
data_drug$P38 = as.numeric(data_drug$P38)
data_drug$P39 = as.numeric(data_drug$P39)
data_drug$P40 = as.numeric(data_drug$P40)
data_drug$P41 = as.numeric(data_drug$P41)
data_drug$P42 = as.numeric(data_drug$P42)
data_drug$P43 = as.numeric(data_drug$P43)
data_drug$P44 = as.numeric(data_drug$P44)
data_drug$P45 = as.numeric(data_drug$P45)
data_drug$P46 = as.numeric(data_drug$P46)
data_drug$P47 = as.numeric(data_drug$P47)
data_drug$P48 = as.numeric(data_drug$P48)
data_drug$P49 = as.numeric(data_drug$P49)
data_drug$P50 = as.numeric(data_drug$P50)
data_drug$P51 = as.numeric(data_drug$P51)
data_drug$P52 = as.numeric(data_drug$P52)
data_drug$P53 = as.numeric(data_drug$P53)
data_drug$P54 = as.numeric(data_drug$P54)
data_drug$P55 = as.numeric(data_drug$P55)
data_drug$P56 = as.numeric(data_drug$P56)
data_drug$P57 = as.numeric(data_drug$P57)
data_drug$P58 = as.numeric(data_drug$P58)
data_drug$P59 = as.numeric(data_drug$P59)
data_drug$P60 = as.numeric(data_drug$P60)
data_drug$P61 = as.numeric(data_drug$P61)
data_drug$P62 = as.numeric(data_drug$P62)
data_drug$P63 = as.numeric(data_drug$P63)
data_drug$P64 = as.numeric(data_drug$P64)
data_drug$P65 = as.numeric(data_drug$P65)
data_drug$P66 = as.numeric(data_drug$P66)
data_drug$P67 = as.numeric(data_drug$P67)
data_drug$P68 = as.numeric(data_drug$P68)
data_drug$P69 = as.numeric(data_drug$P69)
data_drug$P70 = as.numeric(data_drug$P70)
data_drug$P71 = as.numeric(data_drug$P71)
data_drug$P72 = as.numeric(data_drug$P72)
data_drug$P73 = as.numeric(data_drug$P73)
data_drug$P74 = as.numeric(data_drug$P74)
data_drug$P75 = as.numeric(data_drug$P75)
data_drug$P76 = as.numeric(data_drug$P76)
data_drug$P77 = as.numeric(data_drug$P77)
data_drug$P78 = as.numeric(data_drug$P78)
data_drug$P79 = as.numeric(data_drug$P79)
data_drug$P80 = as.numeric(data_drug$P80)
data_drug$P81 = as.numeric(data_drug$P81)
data_drug$P82 = as.numeric(data_drug$P82)
data_drug$P83 = as.numeric(data_drug$P83)
data_drug$P84 = as.numeric(data_drug$P84)
data_drug$P85 = as.numeric(data_drug$P85)
data_drug$P86 = as.numeric(data_drug$P86)
data_drug$P87 = as.numeric(data_drug$P87)
data_drug$P88 = as.numeric(data_drug$P88)
data_drug$P89 = as.numeric(data_drug$P89)
data_drug$P90 = as.numeric(data_drug$P90)
data_drug$P91 = as.numeric(data_drug$P91)
data_drug$P92 = as.numeric(data_drug$P92)
data_drug$P93 = as.numeric(data_drug$P93)
data_drug$P94 = as.numeric(data_drug$P94)
data_drug$P95 = as.numeric(data_drug$P95)
data_drug$P96 = as.numeric(data_drug$P96)
data_drug$P97 = as.numeric(data_drug$P97)
data_drug$P98 = as.numeric(data_drug$P98)
data_drug$P99 = as.numeric(data_drug$P99)
data_drug$P100 = as.numeric(data_drug$P100)
data_drug$P101 = as.numeric(data_drug$P101)
data_drug$P102 = as.numeric(data_drug$P102)
data_drug$P103 = as.numeric(data_drug$P103)
data_drug$P104 = as.numeric(data_drug$P104)
data_drug$P105 = as.numeric(data_drug$P105)
data_drug$P106 = as.numeric(data_drug$P106)
data_drug$P107 = as.numeric(data_drug$P107)
data_drug$P108 = as.numeric(data_drug$P108)
data_drug$P109 = as.numeric(data_drug$P109)
data_drug$P110 = as.numeric(data_drug$P110)
data_drug$P111 = as.numeric(data_drug$P111)
data_drug$P112 = as.numeric(data_drug$P112)
data_drug$P113 = as.numeric(data_drug$P113)
data_drug$P114 = as.numeric(data_drug$P114)
data_drug$P115 = as.numeric(data_drug$P115)
data_drug$P116 = as.numeric(data_drug$P116)
data_drug$P117 = as.numeric(data_drug$P117)
data_drug$P118 = as.numeric(data_drug$P118)
data_drug$P119 = as.numeric(data_drug$P119)
data_drug$P120 = as.numeric(data_drug$P120)
data_drug$P121 = as.numeric(data_drug$P121)
data_drug$P122 = as.numeric(data_drug$P122)
data_drug$P123 = as.numeric(data_drug$P123)
data_drug$P124 = as.numeric(data_drug$P124)
data_drug$P125 = as.numeric(data_drug$P125)
data_drug$P126 = as.numeric(data_drug$P126)
data_drug$P127 = as.numeric(data_drug$P127)
data_drug$P128 = as.numeric(data_drug$P128)
data_drug$P129 = as.numeric(data_drug$P129)
data_drug$P130 = as.numeric(data_drug$P130)
data_drug$P131 = as.numeric(data_drug$P131)
data_drug$P132 = as.numeric(data_drug$P132)
data_drug$P133 = as.numeric(data_drug$P133)
data_drug$P134 = as.numeric(data_drug$P134)
data_drug$P135 = as.numeric(data_drug$P135)
data_drug$P136 = as.numeric(data_drug$P136)
data_drug$P137 = as.numeric(data_drug$P137)
data_drug$P138 = as.numeric(data_drug$P138)
data_drug$P139 = as.numeric(data_drug$P139)
data_drug$P140 = as.numeric(data_drug$P140)
data_drug$P141 = as.numeric(data_drug$P141)
data_drug$P142 = as.numeric(data_drug$P142)
data_drug$P143 = as.numeric(data_drug$P143)
data_drug$P144 = as.numeric(data_drug$P144)
data_drug$P145 = as.numeric(data_drug$P145)
data_drug$P146 = as.numeric(data_drug$P146)
data_drug$P147 = as.numeric(data_drug$P147)
data_drug$P148 = as.numeric(data_drug$P148)
data_drug$P149 = as.numeric(data_drug$P149)
data_drug$P150 = as.numeric(data_drug$P150)
data_drug$P151 = as.numeric(data_drug$P151)
data_drug$P152 = as.numeric(data_drug$P152)
data_drug$P153 = as.numeric(data_drug$P153)
data_drug$P154 = as.numeric(data_drug$P154)
data_drug$P155 = as.numeric(data_drug$P155)
data_drug$P156 = as.numeric(data_drug$P156)
data_drug$P157 = as.numeric(data_drug$P157)
data_drug$P158 = as.numeric(data_drug$P158)
data_drug$P159 = as.numeric(data_drug$P159)
data_drug$P160 = as.numeric(data_drug$P160)
data_drug$P161 = as.numeric(data_drug$P161)
data_drug$P162 = as.numeric(data_drug$P162)
data_drug$P163 = as.numeric(data_drug$P163)
data_drug$P164 = as.numeric(data_drug$P164)
data_drug$P165 = as.numeric(data_drug$P165)
data_drug$P166 = as.numeric(data_drug$P166)
data_drug$P167 = as.numeric(data_drug$P167)
data_drug$P168 = as.numeric(data_drug$P168)
data_drug$P169 = as.numeric(data_drug$P169)
data_drug$P170 = as.numeric(data_drug$P170)
data_drug$P171 = as.numeric(data_drug$P171)
data_drug$P172 = as.numeric(data_drug$P172)
data_drug$P173 = as.numeric(data_drug$P173)
data_drug$P174 = as.numeric(data_drug$P174)
data_drug$P175 = as.numeric(data_drug$P175)
data_drug$P176 = as.numeric(data_drug$P176)
data_drug$P177 = as.numeric(data_drug$P177)
data_drug$P178 = as.numeric(data_drug$P178)
data_drug$P179 = as.numeric(data_drug$P179)
data_drug$P180 = as.numeric(data_drug$P180)
data_drug$P181 = as.numeric(data_drug$P181)
data_drug$P182 = as.numeric(data_drug$P182)
data_drug$P183 = as.numeric(data_drug$P183)
data_drug$P184 = as.numeric(data_drug$P184)
data_drug$P185 = as.numeric(data_drug$P185)
data_drug$P186 = as.numeric(data_drug$P186)
data_drug$P187 = as.numeric(data_drug$P187)
data_drug$P188 = as.numeric(data_drug$P188)
data_drug$P189 = as.numeric(data_drug$P189)
data_drug$P190 = as.numeric(data_drug$P190)
data_drug$P191 = as.numeric(data_drug$P191)
data_drug$P192 = as.numeric(data_drug$P192)
data_drug$P193 = as.numeric(data_drug$P193)
data_drug$P194 = as.numeric(data_drug$P194)
data_drug$P195 = as.numeric(data_drug$P195)
data_drug$P196 = as.numeric(data_drug$P196)
data_drug$P197 = as.numeric(data_drug$P197)
data_drug$P198 = as.numeric(data_drug$P198)
data_drug$P199 = as.numeric(data_drug$P199)
data_drug$P200 = as.numeric(data_drug$P200)
data_drug$P201 = as.numeric(data_drug$P201)
data_drug$P202 = as.numeric(data_drug$P202)
data_drug$P203 = as.numeric(data_drug$P203)
data_drug$P204 = as.numeric(data_drug$P204)
data_drug$P205 = as.numeric(data_drug$P205)
data_drug$P206 = as.numeric(data_drug$P206)
data_drug$P207 = as.numeric(data_drug$P207)
data_drug$P208 = as.numeric(data_drug$P208)
data_drug$P209 = as.numeric(data_drug$P209)
data_drug$P210 = as.numeric(data_drug$P210)
data_drug$P211 = as.numeric(data_drug$P211)
data_drug$P212 = as.numeric(data_drug$P212)
data_drug$P213 = as.numeric(data_drug$P213)
data_drug$P214 = as.numeric(data_drug$P214)
data_drug$P215 = as.numeric(data_drug$P215)
data_drug$P216 = as.numeric(data_drug$P216)
data_drug$P217 = as.numeric(data_drug$P217)
data_drug$P218 = as.numeric(data_drug$P218)
data_drug$P219 = as.numeric(data_drug$P219)
data_drug$P220 = as.numeric(data_drug$P220)
data_drug$P221 = as.numeric(data_drug$P221)
data_drug$P222 = as.numeric(data_drug$P222)
data_drug$P223 = as.numeric(data_drug$P223)
data_drug$P224 = as.numeric(data_drug$P224)
data_drug$P225 = as.numeric(data_drug$P225)
data_drug$P226 = as.numeric(data_drug$P226)
data_drug$P227 = as.numeric(data_drug$P227)
data_drug$P228 = as.numeric(data_drug$P228)
data_drug$P229 = as.numeric(data_drug$P229)
data_drug$P230 = as.numeric(data_drug$P230)
data_drug$P231 = as.numeric(data_drug$P231)
data_drug$P232 = as.numeric(data_drug$P232)
data_drug$P233 = as.numeric(data_drug$P233)
data_drug$P234 = as.numeric(data_drug$P234)
data_drug$P235 = as.numeric(data_drug$P235)
data_drug$P236 = as.numeric(data_drug$P236)
data_drug$P237 = as.numeric(data_drug$P237)
data_drug$P238 = as.numeric(data_drug$P238)
data_drug$P239 = as.numeric(data_drug$P239)
data_drug$P240 = as.numeric(data_drug$P240)

#### sample CV

totalAccuracy <- c()
totalSensitivity = c()
totalSpecificity = c()
cvDivider <- floor(nrow(data_drug) / (cv+1))

for (cv in seq(1:cv)) {
  # assign chunk to data test
  dataTestIndex <- c((cv * cvDivider):(cv * cvDivider + cvDivider))
  dataTest <- data_drug[dataTestIndex,]
  # everything else to train
  dataTrain <- data_drug[-dataTestIndex,]
  
  cylModel <- multinom(bin~., data=dataTrain, maxit=1000, trace=T) 
  
  pred <- predict(cylModel, newdata=dataTest, type="class")
  
  #  classification error
  cv_ac <- postResample(dataTest$bin, pred)[[1]]
  cv_sen = sensitivity(as.factor(dataTest$bin), pred)
  cv_spe = specificity(as.factor(dataTest$bin), pred)
  print(paste('Current Accuracy:',cv_ac,'for CV:',cv))
  totalAccuracy <- c(totalAccuracy, cv_ac)
  totalSensitivity = c(totalSensitivity, cv_sen)
  totalSpecificity = c(totalSpecificity,cv_spe)
}

mean(totalAccuracy)
mean(totalSensitivity,na.rm = TRUE)
mean(totalSpecificity,na.rm = TRUE)

## check if data available for both classes

comment = "Model Valid"

if(is.na(mean(totalAccuracy)))
{
  comment = "Model Invalid as data available for only one class during cross validation."
  decision = NA
}


if(!is.na(mean(totalAccuracy)))
{  
  
  #################
  ## main prediction of input
  #################
  
  ## create test seqence with 240 postitons from input, dummy binary, number of mutations and mixtures length
  
  test_seq = as.matrix(c(input[1:240],0,no_mut_input,max_mix_input))
  test_seq = t(test_seq)
  colnames(test_seq) = colnames(data_drug)
  test_seq = data.frame(test_seq)
  
  ## create model on whole data
  
  final_model <- multinom(bin~., data=data_drug, maxit=1000, trace=T)
  
  ## predict for input sequence using model
  
  main_pred = predict(final_model, newdata=test_seq, type="class")
  
  ## formulate string for prediction / remarks
  
  if(main_pred == 1)
  {
    decision = "Resistant"
  }
  if(main_pred == 0)
  {
    decision = "Susceptible"
  }
  
}
## write to result dataframe

final = rbind(final,c("NRTI","TDF - Tenofovir",nrow(data_drug),cutoff,round(mean(totalAccuracy),4)*100,
                      round(mean(totalSensitivity,na.rm = TRUE),4)*100,
                      round(mean(totalSpecificity,na.rm = TRUE),4)*100,decision,comment))



############
############
####### NRTI ABC - Abacavir workspace 


## read the data from csv file

data <- read.csv(nrti_dataset,header=TRUE,stringsAsFactors = FALSE)    

## subset data to construct model for FPV

data_drug = subset(data, select = -c(SeqID,X3TC,TDF,AZT,DDI,D4T))
data_drug = na.omit(data_drug)

## select data as per subtype

if(is.null(subtype))
{
  subtype = "X"
}

## select data correspondin to subtype

if(subtype != "X")
{
  data_drug = data_drug[data_drug$SubType == subtype,]
}

## remove the subtype column now

data_drug = subset(data_drug, select = -c(SubType))

## replace non-mutant positions with 0

data_drug[data_drug == "-"] = 0

## replace the mutant positions with thier integer codes

index = NULL
no_mut = NULL
max_mix = NULL

for (i in 1:nrow(data_drug))
{
  no_mut[i]=0
  max_mix[i]=0
  for (j in 2:NCOL(data_drug))
  {
    x=utf8ToInt(data_drug[i,j])
    sum = 0
    if(length(x)>1)     #### change to 0
    {
      no_mut[i]=no_mut[i]+1
    }
    if(length(x)>max_mix[i])
    {
      max_mix[i]=length(x) 
    }
    for(k in 1:1)
    {
      sum = sum + x[k]
    }
    data_drug[i,j] = as.numeric(sum)
  }
}

## set cut-off manually

cutoff = 3 ##6
data_drug$bin = ifelse(data_drug$ABC >= cutoff,1,0)
data_drug = subset(data_drug, select = -c(ABC))

## add features - no of mutations in sequence, maximum mixture length

data_drug$no_mut = no_mut
data_drug$max_mix = max_mix


## convert to numeric  

data_drug$P1 = as.numeric(data_drug$P1)
data_drug$P2 = as.numeric(data_drug$P2)
data_drug$P3 = as.numeric(data_drug$P3)
data_drug$P4 = as.numeric(data_drug$P4)
data_drug$P5 = as.numeric(data_drug$P5)
data_drug$P6 = as.numeric(data_drug$P6)
data_drug$P7 = as.numeric(data_drug$P7)
data_drug$P8 = as.numeric(data_drug$P8)
data_drug$P9 = as.numeric(data_drug$P9)
data_drug$P10 = as.numeric(data_drug$P10)
data_drug$P11 = as.numeric(data_drug$P11)
data_drug$P12 = as.numeric(data_drug$P12)
data_drug$P13 = as.numeric(data_drug$P13)
data_drug$P14 = as.numeric(data_drug$P14)
data_drug$P15 = as.numeric(data_drug$P15)
data_drug$P16 = as.numeric(data_drug$P16)
data_drug$P17 = as.numeric(data_drug$P17)
data_drug$P18 = as.numeric(data_drug$P18)
data_drug$P19 = as.numeric(data_drug$P19)
data_drug$P20 = as.numeric(data_drug$P20)
data_drug$P21 = as.numeric(data_drug$P21)
data_drug$P22 = as.numeric(data_drug$P22)
data_drug$P23 = as.numeric(data_drug$P23)
data_drug$P24 = as.numeric(data_drug$P24)
data_drug$P25 = as.numeric(data_drug$P25)
data_drug$P26 = as.numeric(data_drug$P26)
data_drug$P27 = as.numeric(data_drug$P27)
data_drug$P28 = as.numeric(data_drug$P28)
data_drug$P29 = as.numeric(data_drug$P29)
data_drug$P30 = as.numeric(data_drug$P30)
data_drug$P31 = as.numeric(data_drug$P31)
data_drug$P32 = as.numeric(data_drug$P32)
data_drug$P33 = as.numeric(data_drug$P33)
data_drug$P34 = as.numeric(data_drug$P34)
data_drug$P35 = as.numeric(data_drug$P35)
data_drug$P36 = as.numeric(data_drug$P36)
data_drug$P37 = as.numeric(data_drug$P37)
data_drug$P38 = as.numeric(data_drug$P38)
data_drug$P39 = as.numeric(data_drug$P39)
data_drug$P40 = as.numeric(data_drug$P40)
data_drug$P41 = as.numeric(data_drug$P41)
data_drug$P42 = as.numeric(data_drug$P42)
data_drug$P43 = as.numeric(data_drug$P43)
data_drug$P44 = as.numeric(data_drug$P44)
data_drug$P45 = as.numeric(data_drug$P45)
data_drug$P46 = as.numeric(data_drug$P46)
data_drug$P47 = as.numeric(data_drug$P47)
data_drug$P48 = as.numeric(data_drug$P48)
data_drug$P49 = as.numeric(data_drug$P49)
data_drug$P50 = as.numeric(data_drug$P50)
data_drug$P51 = as.numeric(data_drug$P51)
data_drug$P52 = as.numeric(data_drug$P52)
data_drug$P53 = as.numeric(data_drug$P53)
data_drug$P54 = as.numeric(data_drug$P54)
data_drug$P55 = as.numeric(data_drug$P55)
data_drug$P56 = as.numeric(data_drug$P56)
data_drug$P57 = as.numeric(data_drug$P57)
data_drug$P58 = as.numeric(data_drug$P58)
data_drug$P59 = as.numeric(data_drug$P59)
data_drug$P60 = as.numeric(data_drug$P60)
data_drug$P61 = as.numeric(data_drug$P61)
data_drug$P62 = as.numeric(data_drug$P62)
data_drug$P63 = as.numeric(data_drug$P63)
data_drug$P64 = as.numeric(data_drug$P64)
data_drug$P65 = as.numeric(data_drug$P65)
data_drug$P66 = as.numeric(data_drug$P66)
data_drug$P67 = as.numeric(data_drug$P67)
data_drug$P68 = as.numeric(data_drug$P68)
data_drug$P69 = as.numeric(data_drug$P69)
data_drug$P70 = as.numeric(data_drug$P70)
data_drug$P71 = as.numeric(data_drug$P71)
data_drug$P72 = as.numeric(data_drug$P72)
data_drug$P73 = as.numeric(data_drug$P73)
data_drug$P74 = as.numeric(data_drug$P74)
data_drug$P75 = as.numeric(data_drug$P75)
data_drug$P76 = as.numeric(data_drug$P76)
data_drug$P77 = as.numeric(data_drug$P77)
data_drug$P78 = as.numeric(data_drug$P78)
data_drug$P79 = as.numeric(data_drug$P79)
data_drug$P80 = as.numeric(data_drug$P80)
data_drug$P81 = as.numeric(data_drug$P81)
data_drug$P82 = as.numeric(data_drug$P82)
data_drug$P83 = as.numeric(data_drug$P83)
data_drug$P84 = as.numeric(data_drug$P84)
data_drug$P85 = as.numeric(data_drug$P85)
data_drug$P86 = as.numeric(data_drug$P86)
data_drug$P87 = as.numeric(data_drug$P87)
data_drug$P88 = as.numeric(data_drug$P88)
data_drug$P89 = as.numeric(data_drug$P89)
data_drug$P90 = as.numeric(data_drug$P90)
data_drug$P91 = as.numeric(data_drug$P91)
data_drug$P92 = as.numeric(data_drug$P92)
data_drug$P93 = as.numeric(data_drug$P93)
data_drug$P94 = as.numeric(data_drug$P94)
data_drug$P95 = as.numeric(data_drug$P95)
data_drug$P96 = as.numeric(data_drug$P96)
data_drug$P97 = as.numeric(data_drug$P97)
data_drug$P98 = as.numeric(data_drug$P98)
data_drug$P99 = as.numeric(data_drug$P99)
data_drug$P100 = as.numeric(data_drug$P100)
data_drug$P101 = as.numeric(data_drug$P101)
data_drug$P102 = as.numeric(data_drug$P102)
data_drug$P103 = as.numeric(data_drug$P103)
data_drug$P104 = as.numeric(data_drug$P104)
data_drug$P105 = as.numeric(data_drug$P105)
data_drug$P106 = as.numeric(data_drug$P106)
data_drug$P107 = as.numeric(data_drug$P107)
data_drug$P108 = as.numeric(data_drug$P108)
data_drug$P109 = as.numeric(data_drug$P109)
data_drug$P110 = as.numeric(data_drug$P110)
data_drug$P111 = as.numeric(data_drug$P111)
data_drug$P112 = as.numeric(data_drug$P112)
data_drug$P113 = as.numeric(data_drug$P113)
data_drug$P114 = as.numeric(data_drug$P114)
data_drug$P115 = as.numeric(data_drug$P115)
data_drug$P116 = as.numeric(data_drug$P116)
data_drug$P117 = as.numeric(data_drug$P117)
data_drug$P118 = as.numeric(data_drug$P118)
data_drug$P119 = as.numeric(data_drug$P119)
data_drug$P120 = as.numeric(data_drug$P120)
data_drug$P121 = as.numeric(data_drug$P121)
data_drug$P122 = as.numeric(data_drug$P122)
data_drug$P123 = as.numeric(data_drug$P123)
data_drug$P124 = as.numeric(data_drug$P124)
data_drug$P125 = as.numeric(data_drug$P125)
data_drug$P126 = as.numeric(data_drug$P126)
data_drug$P127 = as.numeric(data_drug$P127)
data_drug$P128 = as.numeric(data_drug$P128)
data_drug$P129 = as.numeric(data_drug$P129)
data_drug$P130 = as.numeric(data_drug$P130)
data_drug$P131 = as.numeric(data_drug$P131)
data_drug$P132 = as.numeric(data_drug$P132)
data_drug$P133 = as.numeric(data_drug$P133)
data_drug$P134 = as.numeric(data_drug$P134)
data_drug$P135 = as.numeric(data_drug$P135)
data_drug$P136 = as.numeric(data_drug$P136)
data_drug$P137 = as.numeric(data_drug$P137)
data_drug$P138 = as.numeric(data_drug$P138)
data_drug$P139 = as.numeric(data_drug$P139)
data_drug$P140 = as.numeric(data_drug$P140)
data_drug$P141 = as.numeric(data_drug$P141)
data_drug$P142 = as.numeric(data_drug$P142)
data_drug$P143 = as.numeric(data_drug$P143)
data_drug$P144 = as.numeric(data_drug$P144)
data_drug$P145 = as.numeric(data_drug$P145)
data_drug$P146 = as.numeric(data_drug$P146)
data_drug$P147 = as.numeric(data_drug$P147)
data_drug$P148 = as.numeric(data_drug$P148)
data_drug$P149 = as.numeric(data_drug$P149)
data_drug$P150 = as.numeric(data_drug$P150)
data_drug$P151 = as.numeric(data_drug$P151)
data_drug$P152 = as.numeric(data_drug$P152)
data_drug$P153 = as.numeric(data_drug$P153)
data_drug$P154 = as.numeric(data_drug$P154)
data_drug$P155 = as.numeric(data_drug$P155)
data_drug$P156 = as.numeric(data_drug$P156)
data_drug$P157 = as.numeric(data_drug$P157)
data_drug$P158 = as.numeric(data_drug$P158)
data_drug$P159 = as.numeric(data_drug$P159)
data_drug$P160 = as.numeric(data_drug$P160)
data_drug$P161 = as.numeric(data_drug$P161)
data_drug$P162 = as.numeric(data_drug$P162)
data_drug$P163 = as.numeric(data_drug$P163)
data_drug$P164 = as.numeric(data_drug$P164)
data_drug$P165 = as.numeric(data_drug$P165)
data_drug$P166 = as.numeric(data_drug$P166)
data_drug$P167 = as.numeric(data_drug$P167)
data_drug$P168 = as.numeric(data_drug$P168)
data_drug$P169 = as.numeric(data_drug$P169)
data_drug$P170 = as.numeric(data_drug$P170)
data_drug$P171 = as.numeric(data_drug$P171)
data_drug$P172 = as.numeric(data_drug$P172)
data_drug$P173 = as.numeric(data_drug$P173)
data_drug$P174 = as.numeric(data_drug$P174)
data_drug$P175 = as.numeric(data_drug$P175)
data_drug$P176 = as.numeric(data_drug$P176)
data_drug$P177 = as.numeric(data_drug$P177)
data_drug$P178 = as.numeric(data_drug$P178)
data_drug$P179 = as.numeric(data_drug$P179)
data_drug$P180 = as.numeric(data_drug$P180)
data_drug$P181 = as.numeric(data_drug$P181)
data_drug$P182 = as.numeric(data_drug$P182)
data_drug$P183 = as.numeric(data_drug$P183)
data_drug$P184 = as.numeric(data_drug$P184)
data_drug$P185 = as.numeric(data_drug$P185)
data_drug$P186 = as.numeric(data_drug$P186)
data_drug$P187 = as.numeric(data_drug$P187)
data_drug$P188 = as.numeric(data_drug$P188)
data_drug$P189 = as.numeric(data_drug$P189)
data_drug$P190 = as.numeric(data_drug$P190)
data_drug$P191 = as.numeric(data_drug$P191)
data_drug$P192 = as.numeric(data_drug$P192)
data_drug$P193 = as.numeric(data_drug$P193)
data_drug$P194 = as.numeric(data_drug$P194)
data_drug$P195 = as.numeric(data_drug$P195)
data_drug$P196 = as.numeric(data_drug$P196)
data_drug$P197 = as.numeric(data_drug$P197)
data_drug$P198 = as.numeric(data_drug$P198)
data_drug$P199 = as.numeric(data_drug$P199)
data_drug$P200 = as.numeric(data_drug$P200)
data_drug$P201 = as.numeric(data_drug$P201)
data_drug$P202 = as.numeric(data_drug$P202)
data_drug$P203 = as.numeric(data_drug$P203)
data_drug$P204 = as.numeric(data_drug$P204)
data_drug$P205 = as.numeric(data_drug$P205)
data_drug$P206 = as.numeric(data_drug$P206)
data_drug$P207 = as.numeric(data_drug$P207)
data_drug$P208 = as.numeric(data_drug$P208)
data_drug$P209 = as.numeric(data_drug$P209)
data_drug$P210 = as.numeric(data_drug$P210)
data_drug$P211 = as.numeric(data_drug$P211)
data_drug$P212 = as.numeric(data_drug$P212)
data_drug$P213 = as.numeric(data_drug$P213)
data_drug$P214 = as.numeric(data_drug$P214)
data_drug$P215 = as.numeric(data_drug$P215)
data_drug$P216 = as.numeric(data_drug$P216)
data_drug$P217 = as.numeric(data_drug$P217)
data_drug$P218 = as.numeric(data_drug$P218)
data_drug$P219 = as.numeric(data_drug$P219)
data_drug$P220 = as.numeric(data_drug$P220)
data_drug$P221 = as.numeric(data_drug$P221)
data_drug$P222 = as.numeric(data_drug$P222)
data_drug$P223 = as.numeric(data_drug$P223)
data_drug$P224 = as.numeric(data_drug$P224)
data_drug$P225 = as.numeric(data_drug$P225)
data_drug$P226 = as.numeric(data_drug$P226)
data_drug$P227 = as.numeric(data_drug$P227)
data_drug$P228 = as.numeric(data_drug$P228)
data_drug$P229 = as.numeric(data_drug$P229)
data_drug$P230 = as.numeric(data_drug$P230)
data_drug$P231 = as.numeric(data_drug$P231)
data_drug$P232 = as.numeric(data_drug$P232)
data_drug$P233 = as.numeric(data_drug$P233)
data_drug$P234 = as.numeric(data_drug$P234)
data_drug$P235 = as.numeric(data_drug$P235)
data_drug$P236 = as.numeric(data_drug$P236)
data_drug$P237 = as.numeric(data_drug$P237)
data_drug$P238 = as.numeric(data_drug$P238)
data_drug$P239 = as.numeric(data_drug$P239)
data_drug$P240 = as.numeric(data_drug$P240)

#### sample CV

totalAccuracy <- c()
totalSensitivity = c()
totalSpecificity = c()
cvDivider <- floor(nrow(data_drug) / (cv+1))

for (cv in seq(1:cv)) {
  # assign chunk to data test
  dataTestIndex <- c((cv * cvDivider):(cv * cvDivider + cvDivider))
  dataTest <- data_drug[dataTestIndex,]
  # everything else to train
  dataTrain <- data_drug[-dataTestIndex,]
  
  cylModel <- multinom(bin~., data=dataTrain, maxit=1000, trace=T) 
  
  pred <- predict(cylModel, newdata=dataTest, type="class")
  
  #  classification error
  cv_ac <- postResample(dataTest$bin, pred)[[1]]
  cv_sen = sensitivity(as.factor(dataTest$bin), pred)
  cv_spe = specificity(as.factor(dataTest$bin), pred)
  print(paste('Current Accuracy:',cv_ac,'for CV:',cv))
  totalAccuracy <- c(totalAccuracy, cv_ac)
  totalSensitivity = c(totalSensitivity, cv_sen)
  totalSpecificity = c(totalSpecificity,cv_spe)
}

mean(totalAccuracy)
mean(totalSensitivity,na.rm = TRUE)
mean(totalSpecificity,na.rm = TRUE)

## check if data available for both classes

comment = "Model Valid"

if(is.na(mean(totalAccuracy)))
{
  comment = "Model Invalid as data available for only one class during cross validation."
  decision = NA
}


if(!is.na(mean(totalAccuracy)))
{  
  
  #################
  ## main prediction of input
  #################
  
  ## create test seqence with 240 postitons from input, dummy binary, number of mutations and mixtures length
  
  test_seq = as.matrix(c(input[1:240],0,no_mut_input,max_mix_input))
  test_seq = t(test_seq)
  colnames(test_seq) = colnames(data_drug)
  test_seq = data.frame(test_seq)
  
  ## create model on whole data
  
  final_model <- multinom(bin~., data=data_drug, maxit=1000, trace=)
  
  ## predict for input sequence using model
  
  main_pred = predict(final_model, newdata=test_seq, type="class")
  
  ## formulate string for prediction / remarks
  
  if(main_pred == 1)
  {
    decision = "Resistant"
  }
  if(main_pred == 0)
  {
    decision = "Susceptible"
  }
  
}
## write to result dataframe

final = rbind(final,c("NRTI","ABC - Abacavir",nrow(data_drug),cutoff,round(mean(totalAccuracy),4)*100,
                      round(mean(totalSensitivity,na.rm = TRUE),4)*100,
                      round(mean(totalSpecificity,na.rm = TRUE),4)*100,decision,comment))


############
############
####### NRTI 3TC - Lamivudine workspace 


## read the data from csv file

data <- read.csv(nrti_dataset,header=TRUE,stringsAsFactors = FALSE)    

## subset data to construct model for FPV

data_drug = subset(data, select = -c(SeqID,DDI,TDF,AZT,ABC,D4T))
data_drug = na.omit(data_drug)

## select data as per subtype

if(is.null(subtype))
{
  subtype = "X"
}

## select data correspondin to subtype

if(subtype != "X")
{
  data_drug = data_drug[data_drug$SubType == subtype,]
}

## remove the subtype column now

data_drug = subset(data_drug, select = -c(SubType))

## replace non-mutant positions with 0

data_drug[data_drug == "-"] = 0

## replace the mutant positions with thier integer codes

index = NULL
no_mut = NULL
max_mix = NULL

for (i in 1:nrow(data_drug))
{
  no_mut[i]=0
  max_mix[i]=0
  for (j in 2:NCOL(data_drug))
  {
    x=utf8ToInt(data_drug[i,j])
    sum = 0
    if(length(x)>1)     #### change to 0
    {
      no_mut[i]=no_mut[i]+1
    }
    if(length(x)>max_mix[i])
    {
      max_mix[i]=length(x) 
    }
    for(k in 1:1)
    {
      sum = sum + x[k]
    }
    data_drug[i,j] = as.numeric(sum)
  }
}

## set cut-off manually

cutoff = 3 ##20
data_drug$bin = ifelse(data_drug$X3TC >= cutoff,1,0)
data_drug = subset(data_drug, select = -c(X3TC))

## add features - no of mutations in sequence, maximum mixture length

data_drug$no_mut = no_mut
data_drug$max_mix = max_mix


## convert to numeric  

data_drug$P1 = as.numeric(data_drug$P1)
data_drug$P2 = as.numeric(data_drug$P2)
data_drug$P3 = as.numeric(data_drug$P3)
data_drug$P4 = as.numeric(data_drug$P4)
data_drug$P5 = as.numeric(data_drug$P5)
data_drug$P6 = as.numeric(data_drug$P6)
data_drug$P7 = as.numeric(data_drug$P7)
data_drug$P8 = as.numeric(data_drug$P8)
data_drug$P9 = as.numeric(data_drug$P9)
data_drug$P10 = as.numeric(data_drug$P10)
data_drug$P11 = as.numeric(data_drug$P11)
data_drug$P12 = as.numeric(data_drug$P12)
data_drug$P13 = as.numeric(data_drug$P13)
data_drug$P14 = as.numeric(data_drug$P14)
data_drug$P15 = as.numeric(data_drug$P15)
data_drug$P16 = as.numeric(data_drug$P16)
data_drug$P17 = as.numeric(data_drug$P17)
data_drug$P18 = as.numeric(data_drug$P18)
data_drug$P19 = as.numeric(data_drug$P19)
data_drug$P20 = as.numeric(data_drug$P20)
data_drug$P21 = as.numeric(data_drug$P21)
data_drug$P22 = as.numeric(data_drug$P22)
data_drug$P23 = as.numeric(data_drug$P23)
data_drug$P24 = as.numeric(data_drug$P24)
data_drug$P25 = as.numeric(data_drug$P25)
data_drug$P26 = as.numeric(data_drug$P26)
data_drug$P27 = as.numeric(data_drug$P27)
data_drug$P28 = as.numeric(data_drug$P28)
data_drug$P29 = as.numeric(data_drug$P29)
data_drug$P30 = as.numeric(data_drug$P30)
data_drug$P31 = as.numeric(data_drug$P31)
data_drug$P32 = as.numeric(data_drug$P32)
data_drug$P33 = as.numeric(data_drug$P33)
data_drug$P34 = as.numeric(data_drug$P34)
data_drug$P35 = as.numeric(data_drug$P35)
data_drug$P36 = as.numeric(data_drug$P36)
data_drug$P37 = as.numeric(data_drug$P37)
data_drug$P38 = as.numeric(data_drug$P38)
data_drug$P39 = as.numeric(data_drug$P39)
data_drug$P40 = as.numeric(data_drug$P40)
data_drug$P41 = as.numeric(data_drug$P41)
data_drug$P42 = as.numeric(data_drug$P42)
data_drug$P43 = as.numeric(data_drug$P43)
data_drug$P44 = as.numeric(data_drug$P44)
data_drug$P45 = as.numeric(data_drug$P45)
data_drug$P46 = as.numeric(data_drug$P46)
data_drug$P47 = as.numeric(data_drug$P47)
data_drug$P48 = as.numeric(data_drug$P48)
data_drug$P49 = as.numeric(data_drug$P49)
data_drug$P50 = as.numeric(data_drug$P50)
data_drug$P51 = as.numeric(data_drug$P51)
data_drug$P52 = as.numeric(data_drug$P52)
data_drug$P53 = as.numeric(data_drug$P53)
data_drug$P54 = as.numeric(data_drug$P54)
data_drug$P55 = as.numeric(data_drug$P55)
data_drug$P56 = as.numeric(data_drug$P56)
data_drug$P57 = as.numeric(data_drug$P57)
data_drug$P58 = as.numeric(data_drug$P58)
data_drug$P59 = as.numeric(data_drug$P59)
data_drug$P60 = as.numeric(data_drug$P60)
data_drug$P61 = as.numeric(data_drug$P61)
data_drug$P62 = as.numeric(data_drug$P62)
data_drug$P63 = as.numeric(data_drug$P63)
data_drug$P64 = as.numeric(data_drug$P64)
data_drug$P65 = as.numeric(data_drug$P65)
data_drug$P66 = as.numeric(data_drug$P66)
data_drug$P67 = as.numeric(data_drug$P67)
data_drug$P68 = as.numeric(data_drug$P68)
data_drug$P69 = as.numeric(data_drug$P69)
data_drug$P70 = as.numeric(data_drug$P70)
data_drug$P71 = as.numeric(data_drug$P71)
data_drug$P72 = as.numeric(data_drug$P72)
data_drug$P73 = as.numeric(data_drug$P73)
data_drug$P74 = as.numeric(data_drug$P74)
data_drug$P75 = as.numeric(data_drug$P75)
data_drug$P76 = as.numeric(data_drug$P76)
data_drug$P77 = as.numeric(data_drug$P77)
data_drug$P78 = as.numeric(data_drug$P78)
data_drug$P79 = as.numeric(data_drug$P79)
data_drug$P80 = as.numeric(data_drug$P80)
data_drug$P81 = as.numeric(data_drug$P81)
data_drug$P82 = as.numeric(data_drug$P82)
data_drug$P83 = as.numeric(data_drug$P83)
data_drug$P84 = as.numeric(data_drug$P84)
data_drug$P85 = as.numeric(data_drug$P85)
data_drug$P86 = as.numeric(data_drug$P86)
data_drug$P87 = as.numeric(data_drug$P87)
data_drug$P88 = as.numeric(data_drug$P88)
data_drug$P89 = as.numeric(data_drug$P89)
data_drug$P90 = as.numeric(data_drug$P90)
data_drug$P91 = as.numeric(data_drug$P91)
data_drug$P92 = as.numeric(data_drug$P92)
data_drug$P93 = as.numeric(data_drug$P93)
data_drug$P94 = as.numeric(data_drug$P94)
data_drug$P95 = as.numeric(data_drug$P95)
data_drug$P96 = as.numeric(data_drug$P96)
data_drug$P97 = as.numeric(data_drug$P97)
data_drug$P98 = as.numeric(data_drug$P98)
data_drug$P99 = as.numeric(data_drug$P99)
data_drug$P100 = as.numeric(data_drug$P100)
data_drug$P101 = as.numeric(data_drug$P101)
data_drug$P102 = as.numeric(data_drug$P102)
data_drug$P103 = as.numeric(data_drug$P103)
data_drug$P104 = as.numeric(data_drug$P104)
data_drug$P105 = as.numeric(data_drug$P105)
data_drug$P106 = as.numeric(data_drug$P106)
data_drug$P107 = as.numeric(data_drug$P107)
data_drug$P108 = as.numeric(data_drug$P108)
data_drug$P109 = as.numeric(data_drug$P109)
data_drug$P110 = as.numeric(data_drug$P110)
data_drug$P111 = as.numeric(data_drug$P111)
data_drug$P112 = as.numeric(data_drug$P112)
data_drug$P113 = as.numeric(data_drug$P113)
data_drug$P114 = as.numeric(data_drug$P114)
data_drug$P115 = as.numeric(data_drug$P115)
data_drug$P116 = as.numeric(data_drug$P116)
data_drug$P117 = as.numeric(data_drug$P117)
data_drug$P118 = as.numeric(data_drug$P118)
data_drug$P119 = as.numeric(data_drug$P119)
data_drug$P120 = as.numeric(data_drug$P120)
data_drug$P121 = as.numeric(data_drug$P121)
data_drug$P122 = as.numeric(data_drug$P122)
data_drug$P123 = as.numeric(data_drug$P123)
data_drug$P124 = as.numeric(data_drug$P124)
data_drug$P125 = as.numeric(data_drug$P125)
data_drug$P126 = as.numeric(data_drug$P126)
data_drug$P127 = as.numeric(data_drug$P127)
data_drug$P128 = as.numeric(data_drug$P128)
data_drug$P129 = as.numeric(data_drug$P129)
data_drug$P130 = as.numeric(data_drug$P130)
data_drug$P131 = as.numeric(data_drug$P131)
data_drug$P132 = as.numeric(data_drug$P132)
data_drug$P133 = as.numeric(data_drug$P133)
data_drug$P134 = as.numeric(data_drug$P134)
data_drug$P135 = as.numeric(data_drug$P135)
data_drug$P136 = as.numeric(data_drug$P136)
data_drug$P137 = as.numeric(data_drug$P137)
data_drug$P138 = as.numeric(data_drug$P138)
data_drug$P139 = as.numeric(data_drug$P139)
data_drug$P140 = as.numeric(data_drug$P140)
data_drug$P141 = as.numeric(data_drug$P141)
data_drug$P142 = as.numeric(data_drug$P142)
data_drug$P143 = as.numeric(data_drug$P143)
data_drug$P144 = as.numeric(data_drug$P144)
data_drug$P145 = as.numeric(data_drug$P145)
data_drug$P146 = as.numeric(data_drug$P146)
data_drug$P147 = as.numeric(data_drug$P147)
data_drug$P148 = as.numeric(data_drug$P148)
data_drug$P149 = as.numeric(data_drug$P149)
data_drug$P150 = as.numeric(data_drug$P150)
data_drug$P151 = as.numeric(data_drug$P151)
data_drug$P152 = as.numeric(data_drug$P152)
data_drug$P153 = as.numeric(data_drug$P153)
data_drug$P154 = as.numeric(data_drug$P154)
data_drug$P155 = as.numeric(data_drug$P155)
data_drug$P156 = as.numeric(data_drug$P156)
data_drug$P157 = as.numeric(data_drug$P157)
data_drug$P158 = as.numeric(data_drug$P158)
data_drug$P159 = as.numeric(data_drug$P159)
data_drug$P160 = as.numeric(data_drug$P160)
data_drug$P161 = as.numeric(data_drug$P161)
data_drug$P162 = as.numeric(data_drug$P162)
data_drug$P163 = as.numeric(data_drug$P163)
data_drug$P164 = as.numeric(data_drug$P164)
data_drug$P165 = as.numeric(data_drug$P165)
data_drug$P166 = as.numeric(data_drug$P166)
data_drug$P167 = as.numeric(data_drug$P167)
data_drug$P168 = as.numeric(data_drug$P168)
data_drug$P169 = as.numeric(data_drug$P169)
data_drug$P170 = as.numeric(data_drug$P170)
data_drug$P171 = as.numeric(data_drug$P171)
data_drug$P172 = as.numeric(data_drug$P172)
data_drug$P173 = as.numeric(data_drug$P173)
data_drug$P174 = as.numeric(data_drug$P174)
data_drug$P175 = as.numeric(data_drug$P175)
data_drug$P176 = as.numeric(data_drug$P176)
data_drug$P177 = as.numeric(data_drug$P177)
data_drug$P178 = as.numeric(data_drug$P178)
data_drug$P179 = as.numeric(data_drug$P179)
data_drug$P180 = as.numeric(data_drug$P180)
data_drug$P181 = as.numeric(data_drug$P181)
data_drug$P182 = as.numeric(data_drug$P182)
data_drug$P183 = as.numeric(data_drug$P183)
data_drug$P184 = as.numeric(data_drug$P184)
data_drug$P185 = as.numeric(data_drug$P185)
data_drug$P186 = as.numeric(data_drug$P186)
data_drug$P187 = as.numeric(data_drug$P187)
data_drug$P188 = as.numeric(data_drug$P188)
data_drug$P189 = as.numeric(data_drug$P189)
data_drug$P190 = as.numeric(data_drug$P190)
data_drug$P191 = as.numeric(data_drug$P191)
data_drug$P192 = as.numeric(data_drug$P192)
data_drug$P193 = as.numeric(data_drug$P193)
data_drug$P194 = as.numeric(data_drug$P194)
data_drug$P195 = as.numeric(data_drug$P195)
data_drug$P196 = as.numeric(data_drug$P196)
data_drug$P197 = as.numeric(data_drug$P197)
data_drug$P198 = as.numeric(data_drug$P198)
data_drug$P199 = as.numeric(data_drug$P199)
data_drug$P200 = as.numeric(data_drug$P200)
data_drug$P201 = as.numeric(data_drug$P201)
data_drug$P202 = as.numeric(data_drug$P202)
data_drug$P203 = as.numeric(data_drug$P203)
data_drug$P204 = as.numeric(data_drug$P204)
data_drug$P205 = as.numeric(data_drug$P205)
data_drug$P206 = as.numeric(data_drug$P206)
data_drug$P207 = as.numeric(data_drug$P207)
data_drug$P208 = as.numeric(data_drug$P208)
data_drug$P209 = as.numeric(data_drug$P209)
data_drug$P210 = as.numeric(data_drug$P210)
data_drug$P211 = as.numeric(data_drug$P211)
data_drug$P212 = as.numeric(data_drug$P212)
data_drug$P213 = as.numeric(data_drug$P213)
data_drug$P214 = as.numeric(data_drug$P214)
data_drug$P215 = as.numeric(data_drug$P215)
data_drug$P216 = as.numeric(data_drug$P216)
data_drug$P217 = as.numeric(data_drug$P217)
data_drug$P218 = as.numeric(data_drug$P218)
data_drug$P219 = as.numeric(data_drug$P219)
data_drug$P220 = as.numeric(data_drug$P220)
data_drug$P221 = as.numeric(data_drug$P221)
data_drug$P222 = as.numeric(data_drug$P222)
data_drug$P223 = as.numeric(data_drug$P223)
data_drug$P224 = as.numeric(data_drug$P224)
data_drug$P225 = as.numeric(data_drug$P225)
data_drug$P226 = as.numeric(data_drug$P226)
data_drug$P227 = as.numeric(data_drug$P227)
data_drug$P228 = as.numeric(data_drug$P228)
data_drug$P229 = as.numeric(data_drug$P229)
data_drug$P230 = as.numeric(data_drug$P230)
data_drug$P231 = as.numeric(data_drug$P231)
data_drug$P232 = as.numeric(data_drug$P232)
data_drug$P233 = as.numeric(data_drug$P233)
data_drug$P234 = as.numeric(data_drug$P234)
data_drug$P235 = as.numeric(data_drug$P235)
data_drug$P236 = as.numeric(data_drug$P236)
data_drug$P237 = as.numeric(data_drug$P237)
data_drug$P238 = as.numeric(data_drug$P238)
data_drug$P239 = as.numeric(data_drug$P239)
data_drug$P240 = as.numeric(data_drug$P240)

#### sample CV

totalAccuracy <- c()
totalSensitivity = c()
totalSpecificity = c()
cvDivider <- floor(nrow(data_drug) / (cv+1))

for (cv in seq(1:cv)) {
  # assign chunk to data test
  dataTestIndex <- c((cv * cvDivider):(cv * cvDivider + cvDivider))
  dataTest <- data_drug[dataTestIndex,]
  # everything else to train
  dataTrain <- data_drug[-dataTestIndex,]
  
  cylModel <- multinom(bin~., data=dataTrain, maxit=1000, trace=T) 
  
  pred <- predict(cylModel, newdata=dataTest, type="class")
  
  #  classification error
  cv_ac <- postResample(dataTest$bin, pred)[[1]]
  cv_sen = sensitivity(as.factor(dataTest$bin), pred)
  cv_spe = specificity(as.factor(dataTest$bin), pred)
  print(paste('Current Accuracy:',cv_ac,'for CV:',cv))
  totalAccuracy <- c(totalAccuracy, cv_ac)
  totalSensitivity = c(totalSensitivity, cv_sen)
  totalSpecificity = c(totalSpecificity,cv_spe)
}

mean(totalAccuracy)
mean(totalSensitivity,na.rm = TRUE)
mean(totalSpecificity,na.rm = TRUE)

## check if data available for both classes

comment = "Model Valid"

if(is.na(mean(totalAccuracy)))
{
  comment = "Model Invalid as data available for only one class during cross validation."
  decision = NA
}


if(!is.na(mean(totalAccuracy)))
{  
  
  #################
  ## main prediction of input
  #################
  
  ## create test seqence with 240 postitons from input, dummy binary, number of mutations and mixtures length
  
  test_seq = as.matrix(c(input[1:240],0,no_mut_input,max_mix_input))
  test_seq = t(test_seq)
  colnames(test_seq) = colnames(data_drug)
  test_seq = data.frame(test_seq)
  
  ## create model on whole data
  
  final_model <- multinom(bin~., data=data_drug, maxit=1000, trace=T)
  
  ## predict for input sequence using model
  
  main_pred = predict(final_model, newdata=test_seq, type="class")
  
  ## formulate string for prediction / remarks
  
  if(main_pred == 1)
  {
    decision = "Resistant"
  }
  if(main_pred == 0)
  {
    decision = "Susceptible"
  }
  
}
## write to result dataframe

final = rbind(final,c("NRTI","3TC - Lamivudine",nrow(data_drug),cutoff,round(mean(totalAccuracy),4)*100,
                      round(mean(totalSensitivity,na.rm = TRUE),4)*100,
                      round(mean(totalSpecificity,na.rm = TRUE),4)*100,decision,comment))



############
############
####### NRTI DDI - Didanosine workspace 


## read the data from csv file

data <- read.csv(nrti_dataset,header=TRUE,stringsAsFactors = FALSE)    

## subset data to construct model for FPV

data_drug = subset(data, select = -c(SeqID,X3TC,TDF,AZT,ABC,D4T))
data_drug = na.omit(data_drug)

## select data as per subtype

if(is.null(subtype))
{
  subtype = "X"
}

## select data correspondin to subtype

if(subtype != "X")
{
  data_drug = data_drug[data_drug$SubType == subtype,]
}

## remove the subtype column now

data_drug = subset(data_drug, select = -c(SubType))

## replace non-mutant positions with 0

data_drug[data_drug == "-"] = 0

## replace the mutant positions with thier integer codes

index = NULL
no_mut = NULL
max_mix = NULL

for (i in 1:nrow(data_drug))
{
  no_mut[i]=0
  max_mix[i]=0
  for (j in 2:NCOL(data_drug))
  {
    x=utf8ToInt(data_drug[i,j])
    sum = 0
    if(length(x)>1)     #### change to 0
    {
      no_mut[i]=no_mut[i]+1
    }
    if(length(x)>max_mix[i])
    {
      max_mix[i]=length(x) 
    }
    for(k in 1:1)
    {
      sum = sum + x[k]
    }
    data_drug[i,j] = as.numeric(sum)
  }
}

## set cut-off manually

cutoff = 1.5 ##2
data_drug$bin = ifelse(data_drug$DDI >= cutoff,1,0)
data_drug = subset(data_drug, select = -c(DDI))

## add features - no of mutations in sequence, maximum mixture length

data_drug$no_mut = no_mut
data_drug$max_mix = max_mix


## convert to numeric  

data_drug$P1 = as.numeric(data_drug$P1)
data_drug$P2 = as.numeric(data_drug$P2)
data_drug$P3 = as.numeric(data_drug$P3)
data_drug$P4 = as.numeric(data_drug$P4)
data_drug$P5 = as.numeric(data_drug$P5)
data_drug$P6 = as.numeric(data_drug$P6)
data_drug$P7 = as.numeric(data_drug$P7)
data_drug$P8 = as.numeric(data_drug$P8)
data_drug$P9 = as.numeric(data_drug$P9)
data_drug$P10 = as.numeric(data_drug$P10)
data_drug$P11 = as.numeric(data_drug$P11)
data_drug$P12 = as.numeric(data_drug$P12)
data_drug$P13 = as.numeric(data_drug$P13)
data_drug$P14 = as.numeric(data_drug$P14)
data_drug$P15 = as.numeric(data_drug$P15)
data_drug$P16 = as.numeric(data_drug$P16)
data_drug$P17 = as.numeric(data_drug$P17)
data_drug$P18 = as.numeric(data_drug$P18)
data_drug$P19 = as.numeric(data_drug$P19)
data_drug$P20 = as.numeric(data_drug$P20)
data_drug$P21 = as.numeric(data_drug$P21)
data_drug$P22 = as.numeric(data_drug$P22)
data_drug$P23 = as.numeric(data_drug$P23)
data_drug$P24 = as.numeric(data_drug$P24)
data_drug$P25 = as.numeric(data_drug$P25)
data_drug$P26 = as.numeric(data_drug$P26)
data_drug$P27 = as.numeric(data_drug$P27)
data_drug$P28 = as.numeric(data_drug$P28)
data_drug$P29 = as.numeric(data_drug$P29)
data_drug$P30 = as.numeric(data_drug$P30)
data_drug$P31 = as.numeric(data_drug$P31)
data_drug$P32 = as.numeric(data_drug$P32)
data_drug$P33 = as.numeric(data_drug$P33)
data_drug$P34 = as.numeric(data_drug$P34)
data_drug$P35 = as.numeric(data_drug$P35)
data_drug$P36 = as.numeric(data_drug$P36)
data_drug$P37 = as.numeric(data_drug$P37)
data_drug$P38 = as.numeric(data_drug$P38)
data_drug$P39 = as.numeric(data_drug$P39)
data_drug$P40 = as.numeric(data_drug$P40)
data_drug$P41 = as.numeric(data_drug$P41)
data_drug$P42 = as.numeric(data_drug$P42)
data_drug$P43 = as.numeric(data_drug$P43)
data_drug$P44 = as.numeric(data_drug$P44)
data_drug$P45 = as.numeric(data_drug$P45)
data_drug$P46 = as.numeric(data_drug$P46)
data_drug$P47 = as.numeric(data_drug$P47)
data_drug$P48 = as.numeric(data_drug$P48)
data_drug$P49 = as.numeric(data_drug$P49)
data_drug$P50 = as.numeric(data_drug$P50)
data_drug$P51 = as.numeric(data_drug$P51)
data_drug$P52 = as.numeric(data_drug$P52)
data_drug$P53 = as.numeric(data_drug$P53)
data_drug$P54 = as.numeric(data_drug$P54)
data_drug$P55 = as.numeric(data_drug$P55)
data_drug$P56 = as.numeric(data_drug$P56)
data_drug$P57 = as.numeric(data_drug$P57)
data_drug$P58 = as.numeric(data_drug$P58)
data_drug$P59 = as.numeric(data_drug$P59)
data_drug$P60 = as.numeric(data_drug$P60)
data_drug$P61 = as.numeric(data_drug$P61)
data_drug$P62 = as.numeric(data_drug$P62)
data_drug$P63 = as.numeric(data_drug$P63)
data_drug$P64 = as.numeric(data_drug$P64)
data_drug$P65 = as.numeric(data_drug$P65)
data_drug$P66 = as.numeric(data_drug$P66)
data_drug$P67 = as.numeric(data_drug$P67)
data_drug$P68 = as.numeric(data_drug$P68)
data_drug$P69 = as.numeric(data_drug$P69)
data_drug$P70 = as.numeric(data_drug$P70)
data_drug$P71 = as.numeric(data_drug$P71)
data_drug$P72 = as.numeric(data_drug$P72)
data_drug$P73 = as.numeric(data_drug$P73)
data_drug$P74 = as.numeric(data_drug$P74)
data_drug$P75 = as.numeric(data_drug$P75)
data_drug$P76 = as.numeric(data_drug$P76)
data_drug$P77 = as.numeric(data_drug$P77)
data_drug$P78 = as.numeric(data_drug$P78)
data_drug$P79 = as.numeric(data_drug$P79)
data_drug$P80 = as.numeric(data_drug$P80)
data_drug$P81 = as.numeric(data_drug$P81)
data_drug$P82 = as.numeric(data_drug$P82)
data_drug$P83 = as.numeric(data_drug$P83)
data_drug$P84 = as.numeric(data_drug$P84)
data_drug$P85 = as.numeric(data_drug$P85)
data_drug$P86 = as.numeric(data_drug$P86)
data_drug$P87 = as.numeric(data_drug$P87)
data_drug$P88 = as.numeric(data_drug$P88)
data_drug$P89 = as.numeric(data_drug$P89)
data_drug$P90 = as.numeric(data_drug$P90)
data_drug$P91 = as.numeric(data_drug$P91)
data_drug$P92 = as.numeric(data_drug$P92)
data_drug$P93 = as.numeric(data_drug$P93)
data_drug$P94 = as.numeric(data_drug$P94)
data_drug$P95 = as.numeric(data_drug$P95)
data_drug$P96 = as.numeric(data_drug$P96)
data_drug$P97 = as.numeric(data_drug$P97)
data_drug$P98 = as.numeric(data_drug$P98)
data_drug$P99 = as.numeric(data_drug$P99)
data_drug$P100 = as.numeric(data_drug$P100)
data_drug$P101 = as.numeric(data_drug$P101)
data_drug$P102 = as.numeric(data_drug$P102)
data_drug$P103 = as.numeric(data_drug$P103)
data_drug$P104 = as.numeric(data_drug$P104)
data_drug$P105 = as.numeric(data_drug$P105)
data_drug$P106 = as.numeric(data_drug$P106)
data_drug$P107 = as.numeric(data_drug$P107)
data_drug$P108 = as.numeric(data_drug$P108)
data_drug$P109 = as.numeric(data_drug$P109)
data_drug$P110 = as.numeric(data_drug$P110)
data_drug$P111 = as.numeric(data_drug$P111)
data_drug$P112 = as.numeric(data_drug$P112)
data_drug$P113 = as.numeric(data_drug$P113)
data_drug$P114 = as.numeric(data_drug$P114)
data_drug$P115 = as.numeric(data_drug$P115)
data_drug$P116 = as.numeric(data_drug$P116)
data_drug$P117 = as.numeric(data_drug$P117)
data_drug$P118 = as.numeric(data_drug$P118)
data_drug$P119 = as.numeric(data_drug$P119)
data_drug$P120 = as.numeric(data_drug$P120)
data_drug$P121 = as.numeric(data_drug$P121)
data_drug$P122 = as.numeric(data_drug$P122)
data_drug$P123 = as.numeric(data_drug$P123)
data_drug$P124 = as.numeric(data_drug$P124)
data_drug$P125 = as.numeric(data_drug$P125)
data_drug$P126 = as.numeric(data_drug$P126)
data_drug$P127 = as.numeric(data_drug$P127)
data_drug$P128 = as.numeric(data_drug$P128)
data_drug$P129 = as.numeric(data_drug$P129)
data_drug$P130 = as.numeric(data_drug$P130)
data_drug$P131 = as.numeric(data_drug$P131)
data_drug$P132 = as.numeric(data_drug$P132)
data_drug$P133 = as.numeric(data_drug$P133)
data_drug$P134 = as.numeric(data_drug$P134)
data_drug$P135 = as.numeric(data_drug$P135)
data_drug$P136 = as.numeric(data_drug$P136)
data_drug$P137 = as.numeric(data_drug$P137)
data_drug$P138 = as.numeric(data_drug$P138)
data_drug$P139 = as.numeric(data_drug$P139)
data_drug$P140 = as.numeric(data_drug$P140)
data_drug$P141 = as.numeric(data_drug$P141)
data_drug$P142 = as.numeric(data_drug$P142)
data_drug$P143 = as.numeric(data_drug$P143)
data_drug$P144 = as.numeric(data_drug$P144)
data_drug$P145 = as.numeric(data_drug$P145)
data_drug$P146 = as.numeric(data_drug$P146)
data_drug$P147 = as.numeric(data_drug$P147)
data_drug$P148 = as.numeric(data_drug$P148)
data_drug$P149 = as.numeric(data_drug$P149)
data_drug$P150 = as.numeric(data_drug$P150)
data_drug$P151 = as.numeric(data_drug$P151)
data_drug$P152 = as.numeric(data_drug$P152)
data_drug$P153 = as.numeric(data_drug$P153)
data_drug$P154 = as.numeric(data_drug$P154)
data_drug$P155 = as.numeric(data_drug$P155)
data_drug$P156 = as.numeric(data_drug$P156)
data_drug$P157 = as.numeric(data_drug$P157)
data_drug$P158 = as.numeric(data_drug$P158)
data_drug$P159 = as.numeric(data_drug$P159)
data_drug$P160 = as.numeric(data_drug$P160)
data_drug$P161 = as.numeric(data_drug$P161)
data_drug$P162 = as.numeric(data_drug$P162)
data_drug$P163 = as.numeric(data_drug$P163)
data_drug$P164 = as.numeric(data_drug$P164)
data_drug$P165 = as.numeric(data_drug$P165)
data_drug$P166 = as.numeric(data_drug$P166)
data_drug$P167 = as.numeric(data_drug$P167)
data_drug$P168 = as.numeric(data_drug$P168)
data_drug$P169 = as.numeric(data_drug$P169)
data_drug$P170 = as.numeric(data_drug$P170)
data_drug$P171 = as.numeric(data_drug$P171)
data_drug$P172 = as.numeric(data_drug$P172)
data_drug$P173 = as.numeric(data_drug$P173)
data_drug$P174 = as.numeric(data_drug$P174)
data_drug$P175 = as.numeric(data_drug$P175)
data_drug$P176 = as.numeric(data_drug$P176)
data_drug$P177 = as.numeric(data_drug$P177)
data_drug$P178 = as.numeric(data_drug$P178)
data_drug$P179 = as.numeric(data_drug$P179)
data_drug$P180 = as.numeric(data_drug$P180)
data_drug$P181 = as.numeric(data_drug$P181)
data_drug$P182 = as.numeric(data_drug$P182)
data_drug$P183 = as.numeric(data_drug$P183)
data_drug$P184 = as.numeric(data_drug$P184)
data_drug$P185 = as.numeric(data_drug$P185)
data_drug$P186 = as.numeric(data_drug$P186)
data_drug$P187 = as.numeric(data_drug$P187)
data_drug$P188 = as.numeric(data_drug$P188)
data_drug$P189 = as.numeric(data_drug$P189)
data_drug$P190 = as.numeric(data_drug$P190)
data_drug$P191 = as.numeric(data_drug$P191)
data_drug$P192 = as.numeric(data_drug$P192)
data_drug$P193 = as.numeric(data_drug$P193)
data_drug$P194 = as.numeric(data_drug$P194)
data_drug$P195 = as.numeric(data_drug$P195)
data_drug$P196 = as.numeric(data_drug$P196)
data_drug$P197 = as.numeric(data_drug$P197)
data_drug$P198 = as.numeric(data_drug$P198)
data_drug$P199 = as.numeric(data_drug$P199)
data_drug$P200 = as.numeric(data_drug$P200)
data_drug$P201 = as.numeric(data_drug$P201)
data_drug$P202 = as.numeric(data_drug$P202)
data_drug$P203 = as.numeric(data_drug$P203)
data_drug$P204 = as.numeric(data_drug$P204)
data_drug$P205 = as.numeric(data_drug$P205)
data_drug$P206 = as.numeric(data_drug$P206)
data_drug$P207 = as.numeric(data_drug$P207)
data_drug$P208 = as.numeric(data_drug$P208)
data_drug$P209 = as.numeric(data_drug$P209)
data_drug$P210 = as.numeric(data_drug$P210)
data_drug$P211 = as.numeric(data_drug$P211)
data_drug$P212 = as.numeric(data_drug$P212)
data_drug$P213 = as.numeric(data_drug$P213)
data_drug$P214 = as.numeric(data_drug$P214)
data_drug$P215 = as.numeric(data_drug$P215)
data_drug$P216 = as.numeric(data_drug$P216)
data_drug$P217 = as.numeric(data_drug$P217)
data_drug$P218 = as.numeric(data_drug$P218)
data_drug$P219 = as.numeric(data_drug$P219)
data_drug$P220 = as.numeric(data_drug$P220)
data_drug$P221 = as.numeric(data_drug$P221)
data_drug$P222 = as.numeric(data_drug$P222)
data_drug$P223 = as.numeric(data_drug$P223)
data_drug$P224 = as.numeric(data_drug$P224)
data_drug$P225 = as.numeric(data_drug$P225)
data_drug$P226 = as.numeric(data_drug$P226)
data_drug$P227 = as.numeric(data_drug$P227)
data_drug$P228 = as.numeric(data_drug$P228)
data_drug$P229 = as.numeric(data_drug$P229)
data_drug$P230 = as.numeric(data_drug$P230)
data_drug$P231 = as.numeric(data_drug$P231)
data_drug$P232 = as.numeric(data_drug$P232)
data_drug$P233 = as.numeric(data_drug$P233)
data_drug$P234 = as.numeric(data_drug$P234)
data_drug$P235 = as.numeric(data_drug$P235)
data_drug$P236 = as.numeric(data_drug$P236)
data_drug$P237 = as.numeric(data_drug$P237)
data_drug$P238 = as.numeric(data_drug$P238)
data_drug$P239 = as.numeric(data_drug$P239)
data_drug$P240 = as.numeric(data_drug$P240)

#### sample CV

totalAccuracy <- c()
totalSensitivity = c()
totalSpecificity = c()
cvDivider <- floor(nrow(data_drug) / (cv+1))

for (cv in seq(1:cv)) {
  # assign chunk to data test
  dataTestIndex <- c((cv * cvDivider):(cv * cvDivider + cvDivider))
  dataTest <- data_drug[dataTestIndex,]
  # everything else to train
  dataTrain <- data_drug[-dataTestIndex,]
  
  cylModel <- multinom(bin~., data=dataTrain, maxit=1000, trace=T) 
  
  pred <- predict(cylModel, newdata=dataTest, type="class")
  
  #  classification error
  cv_ac <- postResample(dataTest$bin, pred)[[1]]
  cv_sen = sensitivity(as.factor(dataTest$bin), pred)
  cv_spe = specificity(as.factor(dataTest$bin), pred)
  print(paste('Current Accuracy:',cv_ac,'for CV:',cv))
  totalAccuracy <- c(totalAccuracy, cv_ac)
  totalSensitivity = c(totalSensitivity, cv_sen)
  totalSpecificity = c(totalSpecificity,cv_spe)
}

mean(totalAccuracy)
mean(totalSensitivity,na.rm = TRUE)
mean(totalSpecificity,na.rm = TRUE)

## check if data available for both classes

comment = "Model Valid"

if(is.na(mean(totalAccuracy)))
{
  comment = "Model Invalid as data available for only one class during cross validation."
  decision = NA
}


if(!is.na(mean(totalAccuracy)))
{  
  
  #################
  ## main prediction of input
  #################
  
  ## create test seqence with 240 postitons from input, dummy binary, number of mutations and mixtures length
  
  test_seq = as.matrix(c(input[1:240],0,no_mut_input,max_mix_input))
  test_seq = t(test_seq)
  colnames(test_seq) = colnames(data_drug)
  test_seq = data.frame(test_seq)
  
  ## create model on whole data
  
  final_model <- multinom(bin~., data=data_drug, maxit=1000, trace=T)
  
  ## predict for input sequence using model
  
  main_pred = predict(final_model, newdata=test_seq, type="class")
  
  ## formulate string for prediction / remarks
  
  if(main_pred == 1)
  {
    decision = "Resistant"
  }
  if(main_pred == 0)
  {
    decision = "Susceptible"
  }
  
}
## write to result dataframe

final = rbind(final,c("NRTI","DDI - Didanosine",nrow(data_drug),cutoff,round(mean(totalAccuracy),4)*100,
                      round(mean(totalSensitivity,na.rm = TRUE),4)*100,
                      round(mean(totalSpecificity,na.rm = TRUE),4)*100,decision,comment))



##########################################################################################
####
#####                 NNRTI DRUGS
####
##########################################################################################


############
############
####### NNRTI EFV - Efavirenz workspace 


## read the data from csv file

data <- read.csv(nnrti_dataset,header=TRUE,stringsAsFactors = FALSE)    

## subset data to construct model for FPV

data_drug = subset(data, select = -c(SeqID,NVP,ETR,RPV))
data_drug = na.omit(data_drug)

## select data as per subtype

if(is.null(subtype))
{
  subtype = "X"
}

## select data correspondin to subtype

if(subtype != "X")
{
  data_drug = data_drug[data_drug$SubType == subtype,]
}

## remove the subtype column now

data_drug = subset(data_drug, select = -c(SubType))

## replace non-mutant positions with 0

data_drug[data_drug == "-"] = 0

## replace the mutant positions with thier integer codes

index = NULL
no_mut = NULL
max_mix = NULL

for (i in 1:nrow(data_drug))
{
  no_mut[i]=0
  max_mix[i]=0
  for (j in 2:NCOL(data_drug))
  {
    x=utf8ToInt(data_drug[i,j])
    sum = 0
    if(length(x)>1)     #### change to 0
    {
      no_mut[i]=no_mut[i]+1
    }
    if(length(x)>max_mix[i])
    {
      max_mix[i]=length(x) 
    }
    for(k in 1:1)
    {
      sum = sum + x[k]
    }
    data_drug[i,j] = as.numeric(sum)
  }
}

## set cut-off manually

cutoff = 3 ##10
data_drug$bin = ifelse(data_drug$EFV >= cutoff,1,0)
data_drug = subset(data_drug, select = -c(EFV))

## add features - no of mutations in sequence, maximum mixture length

data_drug$no_mut = no_mut
data_drug$max_mix = max_mix


## convert to numeric  

data_drug$P1 = as.numeric(data_drug$P1)
data_drug$P2 = as.numeric(data_drug$P2)
data_drug$P3 = as.numeric(data_drug$P3)
data_drug$P4 = as.numeric(data_drug$P4)
data_drug$P5 = as.numeric(data_drug$P5)
data_drug$P6 = as.numeric(data_drug$P6)
data_drug$P7 = as.numeric(data_drug$P7)
data_drug$P8 = as.numeric(data_drug$P8)
data_drug$P9 = as.numeric(data_drug$P9)
data_drug$P10 = as.numeric(data_drug$P10)
data_drug$P11 = as.numeric(data_drug$P11)
data_drug$P12 = as.numeric(data_drug$P12)
data_drug$P13 = as.numeric(data_drug$P13)
data_drug$P14 = as.numeric(data_drug$P14)
data_drug$P15 = as.numeric(data_drug$P15)
data_drug$P16 = as.numeric(data_drug$P16)
data_drug$P17 = as.numeric(data_drug$P17)
data_drug$P18 = as.numeric(data_drug$P18)
data_drug$P19 = as.numeric(data_drug$P19)
data_drug$P20 = as.numeric(data_drug$P20)
data_drug$P21 = as.numeric(data_drug$P21)
data_drug$P22 = as.numeric(data_drug$P22)
data_drug$P23 = as.numeric(data_drug$P23)
data_drug$P24 = as.numeric(data_drug$P24)
data_drug$P25 = as.numeric(data_drug$P25)
data_drug$P26 = as.numeric(data_drug$P26)
data_drug$P27 = as.numeric(data_drug$P27)
data_drug$P28 = as.numeric(data_drug$P28)
data_drug$P29 = as.numeric(data_drug$P29)
data_drug$P30 = as.numeric(data_drug$P30)
data_drug$P31 = as.numeric(data_drug$P31)
data_drug$P32 = as.numeric(data_drug$P32)
data_drug$P33 = as.numeric(data_drug$P33)
data_drug$P34 = as.numeric(data_drug$P34)
data_drug$P35 = as.numeric(data_drug$P35)
data_drug$P36 = as.numeric(data_drug$P36)
data_drug$P37 = as.numeric(data_drug$P37)
data_drug$P38 = as.numeric(data_drug$P38)
data_drug$P39 = as.numeric(data_drug$P39)
data_drug$P40 = as.numeric(data_drug$P40)
data_drug$P41 = as.numeric(data_drug$P41)
data_drug$P42 = as.numeric(data_drug$P42)
data_drug$P43 = as.numeric(data_drug$P43)
data_drug$P44 = as.numeric(data_drug$P44)
data_drug$P45 = as.numeric(data_drug$P45)
data_drug$P46 = as.numeric(data_drug$P46)
data_drug$P47 = as.numeric(data_drug$P47)
data_drug$P48 = as.numeric(data_drug$P48)
data_drug$P49 = as.numeric(data_drug$P49)
data_drug$P50 = as.numeric(data_drug$P50)
data_drug$P51 = as.numeric(data_drug$P51)
data_drug$P52 = as.numeric(data_drug$P52)
data_drug$P53 = as.numeric(data_drug$P53)
data_drug$P54 = as.numeric(data_drug$P54)
data_drug$P55 = as.numeric(data_drug$P55)
data_drug$P56 = as.numeric(data_drug$P56)
data_drug$P57 = as.numeric(data_drug$P57)
data_drug$P58 = as.numeric(data_drug$P58)
data_drug$P59 = as.numeric(data_drug$P59)
data_drug$P60 = as.numeric(data_drug$P60)
data_drug$P61 = as.numeric(data_drug$P61)
data_drug$P62 = as.numeric(data_drug$P62)
data_drug$P63 = as.numeric(data_drug$P63)
data_drug$P64 = as.numeric(data_drug$P64)
data_drug$P65 = as.numeric(data_drug$P65)
data_drug$P66 = as.numeric(data_drug$P66)
data_drug$P67 = as.numeric(data_drug$P67)
data_drug$P68 = as.numeric(data_drug$P68)
data_drug$P69 = as.numeric(data_drug$P69)
data_drug$P70 = as.numeric(data_drug$P70)
data_drug$P71 = as.numeric(data_drug$P71)
data_drug$P72 = as.numeric(data_drug$P72)
data_drug$P73 = as.numeric(data_drug$P73)
data_drug$P74 = as.numeric(data_drug$P74)
data_drug$P75 = as.numeric(data_drug$P75)
data_drug$P76 = as.numeric(data_drug$P76)
data_drug$P77 = as.numeric(data_drug$P77)
data_drug$P78 = as.numeric(data_drug$P78)
data_drug$P79 = as.numeric(data_drug$P79)
data_drug$P80 = as.numeric(data_drug$P80)
data_drug$P81 = as.numeric(data_drug$P81)
data_drug$P82 = as.numeric(data_drug$P82)
data_drug$P83 = as.numeric(data_drug$P83)
data_drug$P84 = as.numeric(data_drug$P84)
data_drug$P85 = as.numeric(data_drug$P85)
data_drug$P86 = as.numeric(data_drug$P86)
data_drug$P87 = as.numeric(data_drug$P87)
data_drug$P88 = as.numeric(data_drug$P88)
data_drug$P89 = as.numeric(data_drug$P89)
data_drug$P90 = as.numeric(data_drug$P90)
data_drug$P91 = as.numeric(data_drug$P91)
data_drug$P92 = as.numeric(data_drug$P92)
data_drug$P93 = as.numeric(data_drug$P93)
data_drug$P94 = as.numeric(data_drug$P94)
data_drug$P95 = as.numeric(data_drug$P95)
data_drug$P96 = as.numeric(data_drug$P96)
data_drug$P97 = as.numeric(data_drug$P97)
data_drug$P98 = as.numeric(data_drug$P98)
data_drug$P99 = as.numeric(data_drug$P99)
data_drug$P100 = as.numeric(data_drug$P100)
data_drug$P101 = as.numeric(data_drug$P101)
data_drug$P102 = as.numeric(data_drug$P102)
data_drug$P103 = as.numeric(data_drug$P103)
data_drug$P104 = as.numeric(data_drug$P104)
data_drug$P105 = as.numeric(data_drug$P105)
data_drug$P106 = as.numeric(data_drug$P106)
data_drug$P107 = as.numeric(data_drug$P107)
data_drug$P108 = as.numeric(data_drug$P108)
data_drug$P109 = as.numeric(data_drug$P109)
data_drug$P110 = as.numeric(data_drug$P110)
data_drug$P111 = as.numeric(data_drug$P111)
data_drug$P112 = as.numeric(data_drug$P112)
data_drug$P113 = as.numeric(data_drug$P113)
data_drug$P114 = as.numeric(data_drug$P114)
data_drug$P115 = as.numeric(data_drug$P115)
data_drug$P116 = as.numeric(data_drug$P116)
data_drug$P117 = as.numeric(data_drug$P117)
data_drug$P118 = as.numeric(data_drug$P118)
data_drug$P119 = as.numeric(data_drug$P119)
data_drug$P120 = as.numeric(data_drug$P120)
data_drug$P121 = as.numeric(data_drug$P121)
data_drug$P122 = as.numeric(data_drug$P122)
data_drug$P123 = as.numeric(data_drug$P123)
data_drug$P124 = as.numeric(data_drug$P124)
data_drug$P125 = as.numeric(data_drug$P125)
data_drug$P126 = as.numeric(data_drug$P126)
data_drug$P127 = as.numeric(data_drug$P127)
data_drug$P128 = as.numeric(data_drug$P128)
data_drug$P129 = as.numeric(data_drug$P129)
data_drug$P130 = as.numeric(data_drug$P130)
data_drug$P131 = as.numeric(data_drug$P131)
data_drug$P132 = as.numeric(data_drug$P132)
data_drug$P133 = as.numeric(data_drug$P133)
data_drug$P134 = as.numeric(data_drug$P134)
data_drug$P135 = as.numeric(data_drug$P135)
data_drug$P136 = as.numeric(data_drug$P136)
data_drug$P137 = as.numeric(data_drug$P137)
data_drug$P138 = as.numeric(data_drug$P138)
data_drug$P139 = as.numeric(data_drug$P139)
data_drug$P140 = as.numeric(data_drug$P140)
data_drug$P141 = as.numeric(data_drug$P141)
data_drug$P142 = as.numeric(data_drug$P142)
data_drug$P143 = as.numeric(data_drug$P143)
data_drug$P144 = as.numeric(data_drug$P144)
data_drug$P145 = as.numeric(data_drug$P145)
data_drug$P146 = as.numeric(data_drug$P146)
data_drug$P147 = as.numeric(data_drug$P147)
data_drug$P148 = as.numeric(data_drug$P148)
data_drug$P149 = as.numeric(data_drug$P149)
data_drug$P150 = as.numeric(data_drug$P150)
data_drug$P151 = as.numeric(data_drug$P151)
data_drug$P152 = as.numeric(data_drug$P152)
data_drug$P153 = as.numeric(data_drug$P153)
data_drug$P154 = as.numeric(data_drug$P154)
data_drug$P155 = as.numeric(data_drug$P155)
data_drug$P156 = as.numeric(data_drug$P156)
data_drug$P157 = as.numeric(data_drug$P157)
data_drug$P158 = as.numeric(data_drug$P158)
data_drug$P159 = as.numeric(data_drug$P159)
data_drug$P160 = as.numeric(data_drug$P160)
data_drug$P161 = as.numeric(data_drug$P161)
data_drug$P162 = as.numeric(data_drug$P162)
data_drug$P163 = as.numeric(data_drug$P163)
data_drug$P164 = as.numeric(data_drug$P164)
data_drug$P165 = as.numeric(data_drug$P165)
data_drug$P166 = as.numeric(data_drug$P166)
data_drug$P167 = as.numeric(data_drug$P167)
data_drug$P168 = as.numeric(data_drug$P168)
data_drug$P169 = as.numeric(data_drug$P169)
data_drug$P170 = as.numeric(data_drug$P170)
data_drug$P171 = as.numeric(data_drug$P171)
data_drug$P172 = as.numeric(data_drug$P172)
data_drug$P173 = as.numeric(data_drug$P173)
data_drug$P174 = as.numeric(data_drug$P174)
data_drug$P175 = as.numeric(data_drug$P175)
data_drug$P176 = as.numeric(data_drug$P176)
data_drug$P177 = as.numeric(data_drug$P177)
data_drug$P178 = as.numeric(data_drug$P178)
data_drug$P179 = as.numeric(data_drug$P179)
data_drug$P180 = as.numeric(data_drug$P180)
data_drug$P181 = as.numeric(data_drug$P181)
data_drug$P182 = as.numeric(data_drug$P182)
data_drug$P183 = as.numeric(data_drug$P183)
data_drug$P184 = as.numeric(data_drug$P184)
data_drug$P185 = as.numeric(data_drug$P185)
data_drug$P186 = as.numeric(data_drug$P186)
data_drug$P187 = as.numeric(data_drug$P187)
data_drug$P188 = as.numeric(data_drug$P188)
data_drug$P189 = as.numeric(data_drug$P189)
data_drug$P190 = as.numeric(data_drug$P190)
data_drug$P191 = as.numeric(data_drug$P191)
data_drug$P192 = as.numeric(data_drug$P192)
data_drug$P193 = as.numeric(data_drug$P193)
data_drug$P194 = as.numeric(data_drug$P194)
data_drug$P195 = as.numeric(data_drug$P195)
data_drug$P196 = as.numeric(data_drug$P196)
data_drug$P197 = as.numeric(data_drug$P197)
data_drug$P198 = as.numeric(data_drug$P198)
data_drug$P199 = as.numeric(data_drug$P199)
data_drug$P200 = as.numeric(data_drug$P200)
data_drug$P201 = as.numeric(data_drug$P201)
data_drug$P202 = as.numeric(data_drug$P202)
data_drug$P203 = as.numeric(data_drug$P203)
data_drug$P204 = as.numeric(data_drug$P204)
data_drug$P205 = as.numeric(data_drug$P205)
data_drug$P206 = as.numeric(data_drug$P206)
data_drug$P207 = as.numeric(data_drug$P207)
data_drug$P208 = as.numeric(data_drug$P208)
data_drug$P209 = as.numeric(data_drug$P209)
data_drug$P210 = as.numeric(data_drug$P210)
data_drug$P211 = as.numeric(data_drug$P211)
data_drug$P212 = as.numeric(data_drug$P212)
data_drug$P213 = as.numeric(data_drug$P213)
data_drug$P214 = as.numeric(data_drug$P214)
data_drug$P215 = as.numeric(data_drug$P215)
data_drug$P216 = as.numeric(data_drug$P216)
data_drug$P217 = as.numeric(data_drug$P217)
data_drug$P218 = as.numeric(data_drug$P218)
data_drug$P219 = as.numeric(data_drug$P219)
data_drug$P220 = as.numeric(data_drug$P220)
data_drug$P221 = as.numeric(data_drug$P221)
data_drug$P222 = as.numeric(data_drug$P222)
data_drug$P223 = as.numeric(data_drug$P223)
data_drug$P224 = as.numeric(data_drug$P224)
data_drug$P225 = as.numeric(data_drug$P225)
data_drug$P226 = as.numeric(data_drug$P226)
data_drug$P227 = as.numeric(data_drug$P227)
data_drug$P228 = as.numeric(data_drug$P228)
data_drug$P229 = as.numeric(data_drug$P229)
data_drug$P230 = as.numeric(data_drug$P230)
data_drug$P231 = as.numeric(data_drug$P231)
data_drug$P232 = as.numeric(data_drug$P232)
data_drug$P233 = as.numeric(data_drug$P233)
data_drug$P234 = as.numeric(data_drug$P234)
data_drug$P235 = as.numeric(data_drug$P235)
data_drug$P236 = as.numeric(data_drug$P236)
data_drug$P237 = as.numeric(data_drug$P237)
data_drug$P238 = as.numeric(data_drug$P238)
data_drug$P239 = as.numeric(data_drug$P239)
data_drug$P240 = as.numeric(data_drug$P240)

#### sample CV

totalAccuracy <- c()
totalSensitivity = c()
totalSpecificity = c()
cvDivider <- floor(nrow(data_drug) / (cv+1))

for (cv in seq(1:cv)) {
  # assign chunk to data test
  dataTestIndex <- c((cv * cvDivider):(cv * cvDivider + cvDivider))
  dataTest <- data_drug[dataTestIndex,]
  # everything else to train
  dataTrain <- data_drug[-dataTestIndex,]
  
  cylModel <- multinom(bin~., data=dataTrain, maxit=1000, trace=T) 
  
  pred <- predict(cylModel, newdata=dataTest, type="class")
  
  #  classification error
  cv_ac <- postResample(dataTest$bin, pred)[[1]]
  cv_sen = sensitivity(as.factor(dataTest$bin), pred)
  cv_spe = specificity(as.factor(dataTest$bin), pred)
  print(paste('Current Accuracy:',cv_ac,'for CV:',cv))
  totalAccuracy <- c(totalAccuracy, cv_ac)
  totalSensitivity = c(totalSensitivity, cv_sen)
  totalSpecificity = c(totalSpecificity,cv_spe)
}

mean(totalAccuracy)
mean(totalSensitivity,na.rm = TRUE)
mean(totalSpecificity,na.rm = TRUE)

## check if data available for both classes

comment = "Model Valid"

if(is.na(mean(totalAccuracy)))
{
  comment = "Model Invalid as data available for only one class during cross validation."
  decision = NA
}


if(!is.na(mean(totalAccuracy)))
{  
  
  #################
  ## main prediction of input
  #################
  
  ## create test seqence with 240 postitons from input, dummy binary, number of mutations and mixtures length
  
  test_seq = as.matrix(c(input[1:240],0,no_mut_input,max_mix_input))
  test_seq = t(test_seq)
  colnames(test_seq) = colnames(data_drug)
  test_seq = data.frame(test_seq)
  
  ## create model on whole data
  
  final_model <- multinom(bin~., data=data_drug, maxit=1000, trace=T)
  
  ## predict for input sequence using model
  
  main_pred = predict(final_model, newdata=test_seq, type="class")
  
  ## formulate string for prediction / remarks
  
  if(main_pred == 1)
  {
    decision = "Resistant"
  }
  if(main_pred == 0)
  {
    decision = "Susceptible"
  }
  
}
## write to result dataframe

final = rbind(final,c("NNRTI","EFV - Efavirenz",nrow(data_drug),cutoff,round(mean(totalAccuracy),4)*100,
                      round(mean(totalSensitivity,na.rm = TRUE),4)*100,
                      round(mean(totalSpecificity,na.rm = TRUE),4)*100,decision,comment))



############
############
####### NNRTI NVP - Nevirapine workspace


## read the data from csv file

data <- read.csv(nnrti_dataset,header=TRUE,stringsAsFactors = FALSE)    

## subset data to construct model for FPV

data_drug = subset(data, select = -c(SeqID,EFV,ETR,RPV))
data_drug = na.omit(data_drug)

## select data as per subtype

if(is.null(subtype))
{
  subtype = "X"
}

## select data correspondin to subtype

if(subtype != "X")
{
  data_drug = data_drug[data_drug$SubType == subtype,]
}

## remove the subtype column now

data_drug = subset(data_drug, select = -c(SubType))

## replace non-mutant positions with 0

data_drug[data_drug == "-"] = 0

## replace the mutant positions with thier integer codes

index = NULL
no_mut = NULL
max_mix = NULL

for (i in 1:nrow(data_drug))
{
  no_mut[i]=0
  max_mix[i]=0
  for (j in 2:NCOL(data_drug))
  {
    x=utf8ToInt(data_drug[i,j])
    sum = 0
    if(length(x)>1)     #### change to 0
    {
      no_mut[i]=no_mut[i]+1
    }
    if(length(x)>max_mix[i])
    {
      max_mix[i]=length(x) 
    }
    for(k in 1:1)
    {
      sum = sum + x[k]
    }
    data_drug[i,j] = as.numeric(sum)
  }
}

## set cut-off manually

cutoff = 3 ##10
data_drug$bin = ifelse(data_drug$NVP >= cutoff,1,0)
data_drug = subset(data_drug, select = -c(NVP))

## add features - no of mutations in sequence, maximum mixture length

data_drug$no_mut = no_mut
data_drug$max_mix = max_mix


## convert to numeric  

data_drug$P1 = as.numeric(data_drug$P1)
data_drug$P2 = as.numeric(data_drug$P2)
data_drug$P3 = as.numeric(data_drug$P3)
data_drug$P4 = as.numeric(data_drug$P4)
data_drug$P5 = as.numeric(data_drug$P5)
data_drug$P6 = as.numeric(data_drug$P6)
data_drug$P7 = as.numeric(data_drug$P7)
data_drug$P8 = as.numeric(data_drug$P8)
data_drug$P9 = as.numeric(data_drug$P9)
data_drug$P10 = as.numeric(data_drug$P10)
data_drug$P11 = as.numeric(data_drug$P11)
data_drug$P12 = as.numeric(data_drug$P12)
data_drug$P13 = as.numeric(data_drug$P13)
data_drug$P14 = as.numeric(data_drug$P14)
data_drug$P15 = as.numeric(data_drug$P15)
data_drug$P16 = as.numeric(data_drug$P16)
data_drug$P17 = as.numeric(data_drug$P17)
data_drug$P18 = as.numeric(data_drug$P18)
data_drug$P19 = as.numeric(data_drug$P19)
data_drug$P20 = as.numeric(data_drug$P20)
data_drug$P21 = as.numeric(data_drug$P21)
data_drug$P22 = as.numeric(data_drug$P22)
data_drug$P23 = as.numeric(data_drug$P23)
data_drug$P24 = as.numeric(data_drug$P24)
data_drug$P25 = as.numeric(data_drug$P25)
data_drug$P26 = as.numeric(data_drug$P26)
data_drug$P27 = as.numeric(data_drug$P27)
data_drug$P28 = as.numeric(data_drug$P28)
data_drug$P29 = as.numeric(data_drug$P29)
data_drug$P30 = as.numeric(data_drug$P30)
data_drug$P31 = as.numeric(data_drug$P31)
data_drug$P32 = as.numeric(data_drug$P32)
data_drug$P33 = as.numeric(data_drug$P33)
data_drug$P34 = as.numeric(data_drug$P34)
data_drug$P35 = as.numeric(data_drug$P35)
data_drug$P36 = as.numeric(data_drug$P36)
data_drug$P37 = as.numeric(data_drug$P37)
data_drug$P38 = as.numeric(data_drug$P38)
data_drug$P39 = as.numeric(data_drug$P39)
data_drug$P40 = as.numeric(data_drug$P40)
data_drug$P41 = as.numeric(data_drug$P41)
data_drug$P42 = as.numeric(data_drug$P42)
data_drug$P43 = as.numeric(data_drug$P43)
data_drug$P44 = as.numeric(data_drug$P44)
data_drug$P45 = as.numeric(data_drug$P45)
data_drug$P46 = as.numeric(data_drug$P46)
data_drug$P47 = as.numeric(data_drug$P47)
data_drug$P48 = as.numeric(data_drug$P48)
data_drug$P49 = as.numeric(data_drug$P49)
data_drug$P50 = as.numeric(data_drug$P50)
data_drug$P51 = as.numeric(data_drug$P51)
data_drug$P52 = as.numeric(data_drug$P52)
data_drug$P53 = as.numeric(data_drug$P53)
data_drug$P54 = as.numeric(data_drug$P54)
data_drug$P55 = as.numeric(data_drug$P55)
data_drug$P56 = as.numeric(data_drug$P56)
data_drug$P57 = as.numeric(data_drug$P57)
data_drug$P58 = as.numeric(data_drug$P58)
data_drug$P59 = as.numeric(data_drug$P59)
data_drug$P60 = as.numeric(data_drug$P60)
data_drug$P61 = as.numeric(data_drug$P61)
data_drug$P62 = as.numeric(data_drug$P62)
data_drug$P63 = as.numeric(data_drug$P63)
data_drug$P64 = as.numeric(data_drug$P64)
data_drug$P65 = as.numeric(data_drug$P65)
data_drug$P66 = as.numeric(data_drug$P66)
data_drug$P67 = as.numeric(data_drug$P67)
data_drug$P68 = as.numeric(data_drug$P68)
data_drug$P69 = as.numeric(data_drug$P69)
data_drug$P70 = as.numeric(data_drug$P70)
data_drug$P71 = as.numeric(data_drug$P71)
data_drug$P72 = as.numeric(data_drug$P72)
data_drug$P73 = as.numeric(data_drug$P73)
data_drug$P74 = as.numeric(data_drug$P74)
data_drug$P75 = as.numeric(data_drug$P75)
data_drug$P76 = as.numeric(data_drug$P76)
data_drug$P77 = as.numeric(data_drug$P77)
data_drug$P78 = as.numeric(data_drug$P78)
data_drug$P79 = as.numeric(data_drug$P79)
data_drug$P80 = as.numeric(data_drug$P80)
data_drug$P81 = as.numeric(data_drug$P81)
data_drug$P82 = as.numeric(data_drug$P82)
data_drug$P83 = as.numeric(data_drug$P83)
data_drug$P84 = as.numeric(data_drug$P84)
data_drug$P85 = as.numeric(data_drug$P85)
data_drug$P86 = as.numeric(data_drug$P86)
data_drug$P87 = as.numeric(data_drug$P87)
data_drug$P88 = as.numeric(data_drug$P88)
data_drug$P89 = as.numeric(data_drug$P89)
data_drug$P90 = as.numeric(data_drug$P90)
data_drug$P91 = as.numeric(data_drug$P91)
data_drug$P92 = as.numeric(data_drug$P92)
data_drug$P93 = as.numeric(data_drug$P93)
data_drug$P94 = as.numeric(data_drug$P94)
data_drug$P95 = as.numeric(data_drug$P95)
data_drug$P96 = as.numeric(data_drug$P96)
data_drug$P97 = as.numeric(data_drug$P97)
data_drug$P98 = as.numeric(data_drug$P98)
data_drug$P99 = as.numeric(data_drug$P99)
data_drug$P100 = as.numeric(data_drug$P100)
data_drug$P101 = as.numeric(data_drug$P101)
data_drug$P102 = as.numeric(data_drug$P102)
data_drug$P103 = as.numeric(data_drug$P103)
data_drug$P104 = as.numeric(data_drug$P104)
data_drug$P105 = as.numeric(data_drug$P105)
data_drug$P106 = as.numeric(data_drug$P106)
data_drug$P107 = as.numeric(data_drug$P107)
data_drug$P108 = as.numeric(data_drug$P108)
data_drug$P109 = as.numeric(data_drug$P109)
data_drug$P110 = as.numeric(data_drug$P110)
data_drug$P111 = as.numeric(data_drug$P111)
data_drug$P112 = as.numeric(data_drug$P112)
data_drug$P113 = as.numeric(data_drug$P113)
data_drug$P114 = as.numeric(data_drug$P114)
data_drug$P115 = as.numeric(data_drug$P115)
data_drug$P116 = as.numeric(data_drug$P116)
data_drug$P117 = as.numeric(data_drug$P117)
data_drug$P118 = as.numeric(data_drug$P118)
data_drug$P119 = as.numeric(data_drug$P119)
data_drug$P120 = as.numeric(data_drug$P120)
data_drug$P121 = as.numeric(data_drug$P121)
data_drug$P122 = as.numeric(data_drug$P122)
data_drug$P123 = as.numeric(data_drug$P123)
data_drug$P124 = as.numeric(data_drug$P124)
data_drug$P125 = as.numeric(data_drug$P125)
data_drug$P126 = as.numeric(data_drug$P126)
data_drug$P127 = as.numeric(data_drug$P127)
data_drug$P128 = as.numeric(data_drug$P128)
data_drug$P129 = as.numeric(data_drug$P129)
data_drug$P130 = as.numeric(data_drug$P130)
data_drug$P131 = as.numeric(data_drug$P131)
data_drug$P132 = as.numeric(data_drug$P132)
data_drug$P133 = as.numeric(data_drug$P133)
data_drug$P134 = as.numeric(data_drug$P134)
data_drug$P135 = as.numeric(data_drug$P135)
data_drug$P136 = as.numeric(data_drug$P136)
data_drug$P137 = as.numeric(data_drug$P137)
data_drug$P138 = as.numeric(data_drug$P138)
data_drug$P139 = as.numeric(data_drug$P139)
data_drug$P140 = as.numeric(data_drug$P140)
data_drug$P141 = as.numeric(data_drug$P141)
data_drug$P142 = as.numeric(data_drug$P142)
data_drug$P143 = as.numeric(data_drug$P143)
data_drug$P144 = as.numeric(data_drug$P144)
data_drug$P145 = as.numeric(data_drug$P145)
data_drug$P146 = as.numeric(data_drug$P146)
data_drug$P147 = as.numeric(data_drug$P147)
data_drug$P148 = as.numeric(data_drug$P148)
data_drug$P149 = as.numeric(data_drug$P149)
data_drug$P150 = as.numeric(data_drug$P150)
data_drug$P151 = as.numeric(data_drug$P151)
data_drug$P152 = as.numeric(data_drug$P152)
data_drug$P153 = as.numeric(data_drug$P153)
data_drug$P154 = as.numeric(data_drug$P154)
data_drug$P155 = as.numeric(data_drug$P155)
data_drug$P156 = as.numeric(data_drug$P156)
data_drug$P157 = as.numeric(data_drug$P157)
data_drug$P158 = as.numeric(data_drug$P158)
data_drug$P159 = as.numeric(data_drug$P159)
data_drug$P160 = as.numeric(data_drug$P160)
data_drug$P161 = as.numeric(data_drug$P161)
data_drug$P162 = as.numeric(data_drug$P162)
data_drug$P163 = as.numeric(data_drug$P163)
data_drug$P164 = as.numeric(data_drug$P164)
data_drug$P165 = as.numeric(data_drug$P165)
data_drug$P166 = as.numeric(data_drug$P166)
data_drug$P167 = as.numeric(data_drug$P167)
data_drug$P168 = as.numeric(data_drug$P168)
data_drug$P169 = as.numeric(data_drug$P169)
data_drug$P170 = as.numeric(data_drug$P170)
data_drug$P171 = as.numeric(data_drug$P171)
data_drug$P172 = as.numeric(data_drug$P172)
data_drug$P173 = as.numeric(data_drug$P173)
data_drug$P174 = as.numeric(data_drug$P174)
data_drug$P175 = as.numeric(data_drug$P175)
data_drug$P176 = as.numeric(data_drug$P176)
data_drug$P177 = as.numeric(data_drug$P177)
data_drug$P178 = as.numeric(data_drug$P178)
data_drug$P179 = as.numeric(data_drug$P179)
data_drug$P180 = as.numeric(data_drug$P180)
data_drug$P181 = as.numeric(data_drug$P181)
data_drug$P182 = as.numeric(data_drug$P182)
data_drug$P183 = as.numeric(data_drug$P183)
data_drug$P184 = as.numeric(data_drug$P184)
data_drug$P185 = as.numeric(data_drug$P185)
data_drug$P186 = as.numeric(data_drug$P186)
data_drug$P187 = as.numeric(data_drug$P187)
data_drug$P188 = as.numeric(data_drug$P188)
data_drug$P189 = as.numeric(data_drug$P189)
data_drug$P190 = as.numeric(data_drug$P190)
data_drug$P191 = as.numeric(data_drug$P191)
data_drug$P192 = as.numeric(data_drug$P192)
data_drug$P193 = as.numeric(data_drug$P193)
data_drug$P194 = as.numeric(data_drug$P194)
data_drug$P195 = as.numeric(data_drug$P195)
data_drug$P196 = as.numeric(data_drug$P196)
data_drug$P197 = as.numeric(data_drug$P197)
data_drug$P198 = as.numeric(data_drug$P198)
data_drug$P199 = as.numeric(data_drug$P199)
data_drug$P200 = as.numeric(data_drug$P200)
data_drug$P201 = as.numeric(data_drug$P201)
data_drug$P202 = as.numeric(data_drug$P202)
data_drug$P203 = as.numeric(data_drug$P203)
data_drug$P204 = as.numeric(data_drug$P204)
data_drug$P205 = as.numeric(data_drug$P205)
data_drug$P206 = as.numeric(data_drug$P206)
data_drug$P207 = as.numeric(data_drug$P207)
data_drug$P208 = as.numeric(data_drug$P208)
data_drug$P209 = as.numeric(data_drug$P209)
data_drug$P210 = as.numeric(data_drug$P210)
data_drug$P211 = as.numeric(data_drug$P211)
data_drug$P212 = as.numeric(data_drug$P212)
data_drug$P213 = as.numeric(data_drug$P213)
data_drug$P214 = as.numeric(data_drug$P214)
data_drug$P215 = as.numeric(data_drug$P215)
data_drug$P216 = as.numeric(data_drug$P216)
data_drug$P217 = as.numeric(data_drug$P217)
data_drug$P218 = as.numeric(data_drug$P218)
data_drug$P219 = as.numeric(data_drug$P219)
data_drug$P220 = as.numeric(data_drug$P220)
data_drug$P221 = as.numeric(data_drug$P221)
data_drug$P222 = as.numeric(data_drug$P222)
data_drug$P223 = as.numeric(data_drug$P223)
data_drug$P224 = as.numeric(data_drug$P224)
data_drug$P225 = as.numeric(data_drug$P225)
data_drug$P226 = as.numeric(data_drug$P226)
data_drug$P227 = as.numeric(data_drug$P227)
data_drug$P228 = as.numeric(data_drug$P228)
data_drug$P229 = as.numeric(data_drug$P229)
data_drug$P230 = as.numeric(data_drug$P230)
data_drug$P231 = as.numeric(data_drug$P231)
data_drug$P232 = as.numeric(data_drug$P232)
data_drug$P233 = as.numeric(data_drug$P233)
data_drug$P234 = as.numeric(data_drug$P234)
data_drug$P235 = as.numeric(data_drug$P235)
data_drug$P236 = as.numeric(data_drug$P236)
data_drug$P237 = as.numeric(data_drug$P237)
data_drug$P238 = as.numeric(data_drug$P238)
data_drug$P239 = as.numeric(data_drug$P239)
data_drug$P240 = as.numeric(data_drug$P240)

#### sample CV

totalAccuracy <- c()
totalSensitivity = c()
totalSpecificity = c()
cvDivider <- floor(nrow(data_drug) / (cv+1))

for (cv in seq(1:cv)) {
  # assign chunk to data test
  dataTestIndex <- c((cv * cvDivider):(cv * cvDivider + cvDivider))
  dataTest <- data_drug[dataTestIndex,]
  # everything else to train
  dataTrain <- data_drug[-dataTestIndex,]
  
  cylModel <- multinom(bin~., data=dataTrain, maxit=1000, trace=T) 
  
  pred <- predict(cylModel, newdata=dataTest, type="class")
  
  #  classification error
  cv_ac <- postResample(dataTest$bin, pred)[[1]]
  cv_sen = sensitivity(as.factor(dataTest$bin), pred)
  cv_spe = specificity(as.factor(dataTest$bin), pred)
  print(paste('Current Accuracy:',cv_ac,'for CV:',cv))
  totalAccuracy <- c(totalAccuracy, cv_ac)
  totalSensitivity = c(totalSensitivity, cv_sen)
  totalSpecificity = c(totalSpecificity,cv_spe)
}

mean(totalAccuracy)
mean(totalSensitivity,na.rm = TRUE)
mean(totalSpecificity,na.rm = TRUE)

## check if data available for both classes

comment = "Model Valid"

if(is.na(mean(totalAccuracy)))
{
  comment = "Model Invalid as data available for only one class during cross validation."
  decision = NA
}


if(!is.na(mean(totalAccuracy)))
{  
  
  #################
  ## main prediction of input
  #################
  
  ## create test seqence with 240 postitons from input, dummy binary, number of mutations and mixtures length
  
  test_seq = as.matrix(c(input[1:240],0,no_mut_input,max_mix_input))
  test_seq = t(test_seq)
  colnames(test_seq) = colnames(data_drug)
  test_seq = data.frame(test_seq)
  
  ## create model on whole data
  
  final_model <- multinom(bin~., data=data_drug, maxit=1000, trace=T)
  
  ## predict for input sequence using model
  
  main_pred = predict(final_model, newdata=test_seq, type="class")
  
  ## formulate string for prediction / remarks
  
  if(main_pred == 1)
  {
    decision = "Resistant"
  }
  if(main_pred == 0)
  {
    decision = "Susceptible"
  }
  
}
## write to result dataframe

final = rbind(final,c("NNRTI","NVP - Nevirapine",nrow(data_drug),cutoff,round(mean(totalAccuracy),4)*100,
                      round(mean(totalSensitivity,na.rm = TRUE),4)*100,
                      round(mean(totalSpecificity,na.rm = TRUE),4)*100,decision,comment))


############
############
####### NNRTI ETR - Etravirine workspace


## read the data from csv file

data <- read.csv(nnrti_dataset,header=TRUE,stringsAsFactors = FALSE)    

## subset data to construct model for FPV

data_drug = subset(data, select = -c(SeqID,EFV,NVP,RPV))
data_drug = na.omit(data_drug)

## select data as per subtype

if(is.null(subtype))
{
  subtype = "X"
}

## select data correspondin to subtype

if(subtype != "X")
{
  data_drug = data_drug[data_drug$SubType == subtype,]
}

## remove the subtype column now

data_drug = subset(data_drug, select = -c(SubType))

## replace non-mutant positions with 0

data_drug[data_drug == "-"] = 0

## replace the mutant positions with thier integer codes

index = NULL
no_mut = NULL
max_mix = NULL

for (i in 1:nrow(data_drug))
{
  no_mut[i]=0
  max_mix[i]=0
  for (j in 2:NCOL(data_drug))
  {
    x=utf8ToInt(data_drug[i,j])
    sum = 0
    if(length(x)>1)     #### change to 0
    {
      no_mut[i]=no_mut[i]+1
    }
    if(length(x)>max_mix[i])
    {
      max_mix[i]=length(x) 
    }
    for(k in 1:1)
    {
      sum = sum + x[k]
    }
    data_drug[i,j] = as.numeric(sum)
  }
}

## set cut-off manually

cutoff = 3 ##10
data_drug$bin = ifelse(data_drug$ETR >= cutoff,1,0)
data_drug = subset(data_drug, select = -c(ETR))

## add features - no of mutations in sequence, maximum mixture length

data_drug$no_mut = no_mut
data_drug$max_mix = max_mix


## convert to numeric  

data_drug$P1 = as.numeric(data_drug$P1)
data_drug$P2 = as.numeric(data_drug$P2)
data_drug$P3 = as.numeric(data_drug$P3)
data_drug$P4 = as.numeric(data_drug$P4)
data_drug$P5 = as.numeric(data_drug$P5)
data_drug$P6 = as.numeric(data_drug$P6)
data_drug$P7 = as.numeric(data_drug$P7)
data_drug$P8 = as.numeric(data_drug$P8)
data_drug$P9 = as.numeric(data_drug$P9)
data_drug$P10 = as.numeric(data_drug$P10)
data_drug$P11 = as.numeric(data_drug$P11)
data_drug$P12 = as.numeric(data_drug$P12)
data_drug$P13 = as.numeric(data_drug$P13)
data_drug$P14 = as.numeric(data_drug$P14)
data_drug$P15 = as.numeric(data_drug$P15)
data_drug$P16 = as.numeric(data_drug$P16)
data_drug$P17 = as.numeric(data_drug$P17)
data_drug$P18 = as.numeric(data_drug$P18)
data_drug$P19 = as.numeric(data_drug$P19)
data_drug$P20 = as.numeric(data_drug$P20)
data_drug$P21 = as.numeric(data_drug$P21)
data_drug$P22 = as.numeric(data_drug$P22)
data_drug$P23 = as.numeric(data_drug$P23)
data_drug$P24 = as.numeric(data_drug$P24)
data_drug$P25 = as.numeric(data_drug$P25)
data_drug$P26 = as.numeric(data_drug$P26)
data_drug$P27 = as.numeric(data_drug$P27)
data_drug$P28 = as.numeric(data_drug$P28)
data_drug$P29 = as.numeric(data_drug$P29)
data_drug$P30 = as.numeric(data_drug$P30)
data_drug$P31 = as.numeric(data_drug$P31)
data_drug$P32 = as.numeric(data_drug$P32)
data_drug$P33 = as.numeric(data_drug$P33)
data_drug$P34 = as.numeric(data_drug$P34)
data_drug$P35 = as.numeric(data_drug$P35)
data_drug$P36 = as.numeric(data_drug$P36)
data_drug$P37 = as.numeric(data_drug$P37)
data_drug$P38 = as.numeric(data_drug$P38)
data_drug$P39 = as.numeric(data_drug$P39)
data_drug$P40 = as.numeric(data_drug$P40)
data_drug$P41 = as.numeric(data_drug$P41)
data_drug$P42 = as.numeric(data_drug$P42)
data_drug$P43 = as.numeric(data_drug$P43)
data_drug$P44 = as.numeric(data_drug$P44)
data_drug$P45 = as.numeric(data_drug$P45)
data_drug$P46 = as.numeric(data_drug$P46)
data_drug$P47 = as.numeric(data_drug$P47)
data_drug$P48 = as.numeric(data_drug$P48)
data_drug$P49 = as.numeric(data_drug$P49)
data_drug$P50 = as.numeric(data_drug$P50)
data_drug$P51 = as.numeric(data_drug$P51)
data_drug$P52 = as.numeric(data_drug$P52)
data_drug$P53 = as.numeric(data_drug$P53)
data_drug$P54 = as.numeric(data_drug$P54)
data_drug$P55 = as.numeric(data_drug$P55)
data_drug$P56 = as.numeric(data_drug$P56)
data_drug$P57 = as.numeric(data_drug$P57)
data_drug$P58 = as.numeric(data_drug$P58)
data_drug$P59 = as.numeric(data_drug$P59)
data_drug$P60 = as.numeric(data_drug$P60)
data_drug$P61 = as.numeric(data_drug$P61)
data_drug$P62 = as.numeric(data_drug$P62)
data_drug$P63 = as.numeric(data_drug$P63)
data_drug$P64 = as.numeric(data_drug$P64)
data_drug$P65 = as.numeric(data_drug$P65)
data_drug$P66 = as.numeric(data_drug$P66)
data_drug$P67 = as.numeric(data_drug$P67)
data_drug$P68 = as.numeric(data_drug$P68)
data_drug$P69 = as.numeric(data_drug$P69)
data_drug$P70 = as.numeric(data_drug$P70)
data_drug$P71 = as.numeric(data_drug$P71)
data_drug$P72 = as.numeric(data_drug$P72)
data_drug$P73 = as.numeric(data_drug$P73)
data_drug$P74 = as.numeric(data_drug$P74)
data_drug$P75 = as.numeric(data_drug$P75)
data_drug$P76 = as.numeric(data_drug$P76)
data_drug$P77 = as.numeric(data_drug$P77)
data_drug$P78 = as.numeric(data_drug$P78)
data_drug$P79 = as.numeric(data_drug$P79)
data_drug$P80 = as.numeric(data_drug$P80)
data_drug$P81 = as.numeric(data_drug$P81)
data_drug$P82 = as.numeric(data_drug$P82)
data_drug$P83 = as.numeric(data_drug$P83)
data_drug$P84 = as.numeric(data_drug$P84)
data_drug$P85 = as.numeric(data_drug$P85)
data_drug$P86 = as.numeric(data_drug$P86)
data_drug$P87 = as.numeric(data_drug$P87)
data_drug$P88 = as.numeric(data_drug$P88)
data_drug$P89 = as.numeric(data_drug$P89)
data_drug$P90 = as.numeric(data_drug$P90)
data_drug$P91 = as.numeric(data_drug$P91)
data_drug$P92 = as.numeric(data_drug$P92)
data_drug$P93 = as.numeric(data_drug$P93)
data_drug$P94 = as.numeric(data_drug$P94)
data_drug$P95 = as.numeric(data_drug$P95)
data_drug$P96 = as.numeric(data_drug$P96)
data_drug$P97 = as.numeric(data_drug$P97)
data_drug$P98 = as.numeric(data_drug$P98)
data_drug$P99 = as.numeric(data_drug$P99)
data_drug$P100 = as.numeric(data_drug$P100)
data_drug$P101 = as.numeric(data_drug$P101)
data_drug$P102 = as.numeric(data_drug$P102)
data_drug$P103 = as.numeric(data_drug$P103)
data_drug$P104 = as.numeric(data_drug$P104)
data_drug$P105 = as.numeric(data_drug$P105)
data_drug$P106 = as.numeric(data_drug$P106)
data_drug$P107 = as.numeric(data_drug$P107)
data_drug$P108 = as.numeric(data_drug$P108)
data_drug$P109 = as.numeric(data_drug$P109)
data_drug$P110 = as.numeric(data_drug$P110)
data_drug$P111 = as.numeric(data_drug$P111)
data_drug$P112 = as.numeric(data_drug$P112)
data_drug$P113 = as.numeric(data_drug$P113)
data_drug$P114 = as.numeric(data_drug$P114)
data_drug$P115 = as.numeric(data_drug$P115)
data_drug$P116 = as.numeric(data_drug$P116)
data_drug$P117 = as.numeric(data_drug$P117)
data_drug$P118 = as.numeric(data_drug$P118)
data_drug$P119 = as.numeric(data_drug$P119)
data_drug$P120 = as.numeric(data_drug$P120)
data_drug$P121 = as.numeric(data_drug$P121)
data_drug$P122 = as.numeric(data_drug$P122)
data_drug$P123 = as.numeric(data_drug$P123)
data_drug$P124 = as.numeric(data_drug$P124)
data_drug$P125 = as.numeric(data_drug$P125)
data_drug$P126 = as.numeric(data_drug$P126)
data_drug$P127 = as.numeric(data_drug$P127)
data_drug$P128 = as.numeric(data_drug$P128)
data_drug$P129 = as.numeric(data_drug$P129)
data_drug$P130 = as.numeric(data_drug$P130)
data_drug$P131 = as.numeric(data_drug$P131)
data_drug$P132 = as.numeric(data_drug$P132)
data_drug$P133 = as.numeric(data_drug$P133)
data_drug$P134 = as.numeric(data_drug$P134)
data_drug$P135 = as.numeric(data_drug$P135)
data_drug$P136 = as.numeric(data_drug$P136)
data_drug$P137 = as.numeric(data_drug$P137)
data_drug$P138 = as.numeric(data_drug$P138)
data_drug$P139 = as.numeric(data_drug$P139)
data_drug$P140 = as.numeric(data_drug$P140)
data_drug$P141 = as.numeric(data_drug$P141)
data_drug$P142 = as.numeric(data_drug$P142)
data_drug$P143 = as.numeric(data_drug$P143)
data_drug$P144 = as.numeric(data_drug$P144)
data_drug$P145 = as.numeric(data_drug$P145)
data_drug$P146 = as.numeric(data_drug$P146)
data_drug$P147 = as.numeric(data_drug$P147)
data_drug$P148 = as.numeric(data_drug$P148)
data_drug$P149 = as.numeric(data_drug$P149)
data_drug$P150 = as.numeric(data_drug$P150)
data_drug$P151 = as.numeric(data_drug$P151)
data_drug$P152 = as.numeric(data_drug$P152)
data_drug$P153 = as.numeric(data_drug$P153)
data_drug$P154 = as.numeric(data_drug$P154)
data_drug$P155 = as.numeric(data_drug$P155)
data_drug$P156 = as.numeric(data_drug$P156)
data_drug$P157 = as.numeric(data_drug$P157)
data_drug$P158 = as.numeric(data_drug$P158)
data_drug$P159 = as.numeric(data_drug$P159)
data_drug$P160 = as.numeric(data_drug$P160)
data_drug$P161 = as.numeric(data_drug$P161)
data_drug$P162 = as.numeric(data_drug$P162)
data_drug$P163 = as.numeric(data_drug$P163)
data_drug$P164 = as.numeric(data_drug$P164)
data_drug$P165 = as.numeric(data_drug$P165)
data_drug$P166 = as.numeric(data_drug$P166)
data_drug$P167 = as.numeric(data_drug$P167)
data_drug$P168 = as.numeric(data_drug$P168)
data_drug$P169 = as.numeric(data_drug$P169)
data_drug$P170 = as.numeric(data_drug$P170)
data_drug$P171 = as.numeric(data_drug$P171)
data_drug$P172 = as.numeric(data_drug$P172)
data_drug$P173 = as.numeric(data_drug$P173)
data_drug$P174 = as.numeric(data_drug$P174)
data_drug$P175 = as.numeric(data_drug$P175)
data_drug$P176 = as.numeric(data_drug$P176)
data_drug$P177 = as.numeric(data_drug$P177)
data_drug$P178 = as.numeric(data_drug$P178)
data_drug$P179 = as.numeric(data_drug$P179)
data_drug$P180 = as.numeric(data_drug$P180)
data_drug$P181 = as.numeric(data_drug$P181)
data_drug$P182 = as.numeric(data_drug$P182)
data_drug$P183 = as.numeric(data_drug$P183)
data_drug$P184 = as.numeric(data_drug$P184)
data_drug$P185 = as.numeric(data_drug$P185)
data_drug$P186 = as.numeric(data_drug$P186)
data_drug$P187 = as.numeric(data_drug$P187)
data_drug$P188 = as.numeric(data_drug$P188)
data_drug$P189 = as.numeric(data_drug$P189)
data_drug$P190 = as.numeric(data_drug$P190)
data_drug$P191 = as.numeric(data_drug$P191)
data_drug$P192 = as.numeric(data_drug$P192)
data_drug$P193 = as.numeric(data_drug$P193)
data_drug$P194 = as.numeric(data_drug$P194)
data_drug$P195 = as.numeric(data_drug$P195)
data_drug$P196 = as.numeric(data_drug$P196)
data_drug$P197 = as.numeric(data_drug$P197)
data_drug$P198 = as.numeric(data_drug$P198)
data_drug$P199 = as.numeric(data_drug$P199)
data_drug$P200 = as.numeric(data_drug$P200)
data_drug$P201 = as.numeric(data_drug$P201)
data_drug$P202 = as.numeric(data_drug$P202)
data_drug$P203 = as.numeric(data_drug$P203)
data_drug$P204 = as.numeric(data_drug$P204)
data_drug$P205 = as.numeric(data_drug$P205)
data_drug$P206 = as.numeric(data_drug$P206)
data_drug$P207 = as.numeric(data_drug$P207)
data_drug$P208 = as.numeric(data_drug$P208)
data_drug$P209 = as.numeric(data_drug$P209)
data_drug$P210 = as.numeric(data_drug$P210)
data_drug$P211 = as.numeric(data_drug$P211)
data_drug$P212 = as.numeric(data_drug$P212)
data_drug$P213 = as.numeric(data_drug$P213)
data_drug$P214 = as.numeric(data_drug$P214)
data_drug$P215 = as.numeric(data_drug$P215)
data_drug$P216 = as.numeric(data_drug$P216)
data_drug$P217 = as.numeric(data_drug$P217)
data_drug$P218 = as.numeric(data_drug$P218)
data_drug$P219 = as.numeric(data_drug$P219)
data_drug$P220 = as.numeric(data_drug$P220)
data_drug$P221 = as.numeric(data_drug$P221)
data_drug$P222 = as.numeric(data_drug$P222)
data_drug$P223 = as.numeric(data_drug$P223)
data_drug$P224 = as.numeric(data_drug$P224)
data_drug$P225 = as.numeric(data_drug$P225)
data_drug$P226 = as.numeric(data_drug$P226)
data_drug$P227 = as.numeric(data_drug$P227)
data_drug$P228 = as.numeric(data_drug$P228)
data_drug$P229 = as.numeric(data_drug$P229)
data_drug$P230 = as.numeric(data_drug$P230)
data_drug$P231 = as.numeric(data_drug$P231)
data_drug$P232 = as.numeric(data_drug$P232)
data_drug$P233 = as.numeric(data_drug$P233)
data_drug$P234 = as.numeric(data_drug$P234)
data_drug$P235 = as.numeric(data_drug$P235)
data_drug$P236 = as.numeric(data_drug$P236)
data_drug$P237 = as.numeric(data_drug$P237)
data_drug$P238 = as.numeric(data_drug$P238)
data_drug$P239 = as.numeric(data_drug$P239)
data_drug$P240 = as.numeric(data_drug$P240)

#### sample CV

totalAccuracy <- c()
totalSensitivity = c()
totalSpecificity = c()
cvDivider <- floor(nrow(data_drug) / (cv+1))

for (cv in seq(1:cv)) {
  # assign chunk to data test
  dataTestIndex <- c((cv * cvDivider):(cv * cvDivider + cvDivider))
  dataTest <- data_drug[dataTestIndex,]
  # everything else to train
  dataTrain <- data_drug[-dataTestIndex,]
  
  cylModel <- multinom(bin~., data=dataTrain, maxit=1000, trace=T) 
  
  pred <- predict(cylModel, newdata=dataTest, type="class")
  
  #  classification error
  cv_ac <- postResample(dataTest$bin, pred)[[1]]
  cv_sen = sensitivity(as.factor(dataTest$bin), pred)
  cv_spe = specificity(as.factor(dataTest$bin), pred)
  print(paste('Current Accuracy:',cv_ac,'for CV:',cv))
  totalAccuracy <- c(totalAccuracy, cv_ac)
  totalSensitivity = c(totalSensitivity, cv_sen)
  totalSpecificity = c(totalSpecificity,cv_spe)
}

mean(totalAccuracy)
mean(totalSensitivity,na.rm = TRUE)
mean(totalSpecificity,na.rm = TRUE)

## check if data available for both classes

comment = "Model Valid"

if(is.na(mean(totalAccuracy)))
{
  comment = "Model Invalid as data available for only one class during cross validation."
  decision = NA
}


if(!is.na(mean(totalAccuracy)))
{  
  
  #################
  ## main prediction of input
  #################
  
  ## create test seqence with 240 postitons from input, dummy binary, number of mutations and mixtures length
  
  test_seq = as.matrix(c(input[1:240],0,no_mut_input,max_mix_input))
  test_seq = t(test_seq)
  colnames(test_seq) = colnames(data_drug)
  test_seq = data.frame(test_seq)
  
  ## create model on whole data
  
  final_model <- multinom(bin~., data=data_drug, maxit=1500, trace=T)
  
  ## predict for input sequence using model
  
  main_pred = predict(final_model, newdata=test_seq, type="class")
  
  ## formulate string for prediction / remarks
  
  if(main_pred == 1)
  {
    decision = "Resistant"
  }
  if(main_pred == 0)
  {
    decision = "Susceptible"
  }
  
}
## write to result dataframe

final = rbind(final,c("NNRTI","ETR - Etravirine",nrow(data_drug),cutoff,round(mean(totalAccuracy),4)*100,
                      round(mean(totalSensitivity,na.rm = TRUE),4)*100,
                      round(mean(totalSpecificity,na.rm = TRUE),4)*100,decision,comment))


final

## writing to the REPORT XLSX

wb = loadWorkbook(output_file)

## prepare PI dataframe to write and write it to report

pi_output = final[final$INHIBITOR == "PI",]
pi_output = subset(pi_output, select = -(INHIBITOR))

writeWorksheet(wb, pi_output, sheet = "Report", startRow = 5, startCol = 1, header = FALSE)
saveWorkbook(wb)

## prepare NRTI dataframe to write and write it to report

nrti_output = final[final$INHIBITOR == "NRTI",]
nrti_output = subset(nrti_output, select = -(INHIBITOR))

writeWorksheet(wb, nrti_output, sheet = "Report", startRow = 16, startCol = 1, header = FALSE)
saveWorkbook(wb)

## prepare NNRTI dataframe to write and write it to report

nnrti_output = final[final$INHIBITOR == "NNRTI",]
nnrti_output = subset(nnrti_output, select = -(INHIBITOR))

writeWorksheet(wb, nnrti_output, sheet = "Report", startRow = 25, startCol = 1, header = FALSE)
saveWorkbook(wb)

## END OF CODE




