BiocManager::install('pROC', force = TRUE)

install.packages("BiocManager")

library(BiocManager)

library(DESeq2)
library(GenomeInfoDb)
library(edgeR)
library(limma)
library(gplots)
library(RColorBrewer)
library(tibble)
library(dplyr)
library(devtools)
library(ggpubr)
library(EnhancedVolcano)
library(caret)
library(randomForest)
library(pROC)


#Load in the data and format for deseq2 analysis
countdata <- read.table('/Users/rk092/Downloads/GSE263022_main.txt', sep = '\t', header=TRUE, row.names = 1)
count_matrix <- data.matrix(countdata)

## condition is the experimental setup so thats why reordering samples would be important to have all samples of one group together
condition <- factor(c(rep("Control", 36), rep("Treated", 48)))

coldata <- data.frame(row.names=colnames(count_matrix), condition)
dds <- DESeqDataSetFromMatrix(countData=count_matrix, colData=coldata, design=~condition)
dds$condition <- relevel(dds$condition, ref = "Control")

#Run deseq2
dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds)

#get results
res <- res[order(res$padj), ]
ressig <- subset(res, padj < 0.05)
summary(ressig)

# out of 164 with nonzero total read count
adjusted p-value < 0.1
LFC > 0 (up)       : 28, 17%
LFC < 0 (down)     : 136, 83%
outliers [1]       : 0, 0%
low counts [2]     : 0, 0%
(mean count < 1)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results


resdata <- merge(as.data.frame(ressig), as.data.frame(counts(dds, normalized=TRUE)), by= 'row.names', sort=FALSE)
names(resdata)[1] <- "gene"
head(resdata)
summary(resdata)

write.csv(as.data.frame(resdata), file="/Users/rk092/Downloads/GSE263022_main.csv")

# Load your data (replace with your file path)
data <- read.csv("/Users/rk092/Downloads/input_main.csv", row.names = 1)
# Replace with your file
# Ensure data contains DE genes and a 'Response' column

# Define DE genes identified by DESeq2
de_genes <- c('ADRB2', 'CXCL8', 'GOPC', 'KCNN4', 'LDLR', 'PSMD5', 'TNF') 
# Replace with your DE genes
features <- data[, de_genes]
response <- as.factor(data$Response)  
# Response: Binary (e.g., 0 = Non-responder, 1 = Responder)

# Split data into training and testing sets (75/25 split)
set.seed(42)
train_index <- createDataPartition(response, p = 0.75, list = FALSE)
train_data <- features[train_index, ]
train_response <- response[train_index]
test_data <- features[-train_index, ]
test_response <- response[-train_index]

# Normalize the data (center and scale)
preProc <- preProcess(train_data, method = c("center", "scale"))
train_data <- predict(preProc, train_data)
test_data <- predict(preProc, test_data)

# Train a logistic regression model
log_model <- train(x = train_data, y = train_response, method = "glm", family = "binomial")

# Summarize the model
summary(log_model$finalModel)

# Call:
NULL

Coefficients:
  Estimate Std. Error z value Pr(>|z|)  
(Intercept)  -0.3924     0.6574  -0.597   0.5506  
ADRB2        -0.6666     0.6540  -1.019   0.3081  
CXCL8        -2.1913     2.4973  -0.877   0.3802  
GOPC         -0.8335     0.4967  -1.678   0.0933 .
KCNN4         1.1913     0.6205   1.920   0.0549 .
LDLR         -0.7431     0.5713  -1.301   0.1933  
PSMD5        -0.1481     0.4362  -0.340   0.7342  
TNF          -0.9904     1.4242  -0.695   0.4868  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

(Dispersion parameter for binomial family taken to be 1)

Null deviance: 86.046  on 62  degrees of freedom
Residual deviance: 45.236  on 55  degrees of freedom
AIC: 61.236

Number of Fisher Scoring iterations: 7


# Predict on the test set
log_pred <- predict(log_model, newdata = test_data)  # Class predictions
log_prob <- predict(log_model, newdata = test_data, type = "prob")[, 2]  # Probabilities

# Evaluate model performance
# Confusion matrix
conf_matrix <- confusionMatrix(log_pred, test_response)
print(conf_matrix)

# Confusion Matrix and Statistics

Reference
Prediction  0  1
0  5  1
1  4 11

Accuracy : 0.7619          
95% CI : (0.5283, 0.9178)
No Information Rate : 0.5714          
P-Value [Acc > NIR] : 0.05841         

Kappa : 0.4928          

Mcnemar's Test P-Value : 0.37109         
                                          
            Sensitivity : 0.5556          
            Specificity : 0.9167          
         Pos Pred Value : 0.8333          
         Neg Pred Value : 0.7333          
             Prevalence : 0.4286          
         Detection Rate : 0.2381          
   Detection Prevalence : 0.2857          
      Balanced Accuracy : 0.7361          
                                          
       'Positive' Class : 0 '
       
# iteration 1

# ROC Curve and AUC
roc_curve <- roc(test_response, log_prob)
auc_value <- auc(roc_curve)
cat("AUC:", auc_value, "\n")

# Plot ROC Curve
plot(roc_curve, col = "blue", main = "Logistic Regression ROC Curve")
abline(a = 0, b = 1, lty = 2, col = "gray")


