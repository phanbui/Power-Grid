#########################################################################################################
# PART 1
#########################################################################################################


# Use tidyverse package
library(tidyverse)


# Read data
df <- read.table("TermProjectData.txt", header = TRUE, sep = ",")
head(df)
nrow(df)


# Filter out infinite value
filter_infinite_value <- function(df){
  df <- df[is.finite(df$Global_active_power) & is.finite(df$Global_reactive_power) 
           & is.finite(df$Global_reactive_power) & is.finite(df$Global_intensity)
           & is.finite(df$Sub_metering_1) & is.finite(df$Sub_metering_2) & is.finite(df$Sub_metering_3),]
}

df <- filter_infinite_value(df)
head(df)
nrow(df)


# Select numeric fields for PCA
df_numeric_fields <- subset(df, select = -c(Date, Time))


# Picking principle component
pca <- prcomp(df_numeric_fields, scale = TRUE)
summary(pca)


# Finding features with most impact
loading_scores <- pca$rotation[,1]
variable_scores <- abs(loading_scores)
variable_scores_ranked <- sort(variable_scores, decreasing = TRUE)
variable_scores_ranked


#########################################################################################################
# PART 2
#########################################################################################################


# Replace Date and Time in df with POSIX
add_posix_vector <- function(df){
  posix_vector <- vector(mode = "character", nrow(df))
  for (i in 1:nrow(df)){
    posix_vector[i] <- paste(df[i, 1], df[i, 2])
  }
  posix_vector <- as.POSIXlt(posix_vector, format = "%d/%m/%Y %H:%M:%S")
  df$POSIX <- posix_vector
  df <- subset(df, select = -c(Date, Time))
  return(df)
}

df <- add_posix_vector(df)
head(df)


# Pick days for training data
observation <- df[df$POSIX$wday == 2 & df$POSIX$hour >= 9 & df$POSIX$hour < 11,]
nrow(observation)


# Remove incomplete days from training data 
remove_incomplete_day <- function(df){
  iterator <- 1
  while(iterator <= nrow(df)){
    if(iterator + 119 <= nrow(df)){
      if(df[iterator,"POSIX"]$yday != df[iterator+119, "POSIX"]$yday){
        # Filter out incomplete day
        print(paste("Remove day at row", iterator, ", year is", df[iterator,"POSIX"]$year
                    , ", yd is", df[iterator,"POSIX"]$yd))
        df <- df[!(df$POSIX$year == df[iterator, "POSIX"]$year 
                   & df$POSIX$yd == df[iterator, "POSIX"]$yd),]
      }
      else{
        iterator <- iterator + 120
      }
    }
    else{
      df <- df[!(df$POSIX$year == df[iterator, "POSIX"]$year 
                 & df$POSIX$yd == df[iterator, "POSIX"]$yd),]
      break
    }
  }
  return(df)
}

observation <- remove_incomplete_day(observation)
nrow(observation)


# Scale Global_intensity & Global_active_power
observation$Global_intensity_scaled <- scale(observation$Global_intensity, center = TRUE, scale = TRUE)
observation$Global_active_power_scaled <- scale(observation$Global_active_power, center = TRUE, scale = TRUE)
head(observation)


# Choose number of training and testing observation
observation_num <- nrow(observation)/120
observation_num
train_observation_num <- 124
test_observation_num <- observation_num - train_observation_num
train_observation <- observation[1:(train_observation_num*120),]
test_observation <- observation[(train_observation_num*120+1):(train_observation_num*120+test_observation_num*120),]



# Use depmixS4 package
library(depmixS4)
set.seed(1)


# Train ntimes
train_ntimes <- numeric(train_observation_num)
for(i in 1:train_observation_num){
  train_ntimes[i] = 120
}
train_ntimes

# Test ntimes
test_ntimes <- numeric(test_observation_num)
for(i in 1:test_observation_num){
  test_ntimes[i] = 120
}
test_ntimes

# Training model
state_num <- c(4:15)
state_num_len <- length(state_num)

train_log_likelihood <- numeric(state_num_len)
train_bic <- numeric(state_num_len)

test_log_likelihood <- numeric(state_num_len)
test_bic <- numeric(state_num_len)

for (i in 1:length(state_num)){
  writeLines("\n\n\n\n##################################################################################################")
  writeLines(paste("NUMBER OF STATES:", state_num[i]))
  writeLines("##################################################################################################")
  # Training
  writeLines("------------------------- TRAINING -------------------------")
  train_model <- depmix(response = list(Global_intensity_scaled ~ 1, Global_active_power_scaled ~ 1)
                        , data = train_observation, nstates = state_num[i], ntimes = train_ntimes
                        , family = list(gaussian(), gaussian()))
  fit_train_model <- fit(train_model)
  train_log_likelihood[i] = logLik(fit_train_model)
  train_bic[i] = BIC(fit_train_model)
  summary(fit_train_model)
  print(fit_train_model)
  
  # Testing
  writeLines("\n------------------------ TEST RESULT -----------------------")
  test_model <- depmix(response = list(Global_intensity_scaled ~ 1, Global_active_power_scaled ~ 1)
                       , data = test_observation, nstates = state_num[i], ntimes = test_ntimes
                       , family = list(gaussian(), gaussian()))
  test_model <- setpars(test_model, getpars(fit_train_model))
  
  print("Hello")
  print(BIC(test_model))
  
  x <- forwardbackward(test_model)
  test_log_likelihood[i] <- x$logLik
  test_bic[i] <- BIC(test_model)
  writeLines(paste("log Lik:", test_log_likelihood[i]))
  writeLines(paste("BIC:", test_bic[i]))
}


# Plot Log Likelihood
plot(state_num, train_log_likelihood, col = "blue", xlab = "Number of states", ylab = "Log   Likelihood")
lines(state_num, train_log_likelihood, col = "blue")
points(state_num, test_log_likelihood, col = "red")
lines(state_num, test_log_likelihood, col = "red")
abline(h=0)
legend("bottomright", legend = c("train", "test"), lty = c(1,1), col = c("blue", "red"))

# Plot BIC
plot(state_num, train_bic, col = "blue", xlab = "Number of states", ylab = "BIC")
lines(state_num, train_bic, col = "blue")
points(state_num, test_bic, col = "red")
lines(state_num, test_bic, col = "red")
abline(h=0)
legend("topright", legend = c("train", "test"), lty = c(1,1), col = c("blue", "red"))


#########################################################################################################
# PART 3
#########################################################################################################


# Read data
anomaly_data1 <- read.table("DataWithAnomalies1.txt", header = TRUE, sep = ",")
anomaly_data2 <- read.table("DataWithAnomalies2.txt", header = TRUE, sep = ",")
anomaly_data3 <- read.table("DataWithAnomalies3.txt", header = TRUE, sep = ",")


# Filter out infinite value
anomaly_data1 <- filter_infinite_value(anomaly_data1)
anomaly_data2 <- filter_infinite_value(anomaly_data2)
anomaly_data3 <- filter_infinite_value(anomaly_data3)


# Replace Date and Time in df with POSIX
anomaly_data1 <- add_posix_vector(anomaly_data1)
anomaly_data2 <- add_posix_vector(anomaly_data2)
anomaly_data3 <- add_posix_vector(anomaly_data3)


# Pick days for training data
anomaly_observation1 <- anomaly_data1[anomaly_data1$POSIX$wday == 2 & anomaly_data1$POSIX$hour >= 9 
                                      & anomaly_data1$POSIX$hour < 11,]
anomaly_observation2 <- anomaly_data2[anomaly_data2$POSIX$wday == 2 & anomaly_data2$POSIX$hour >= 9 
                                      & anomaly_data2$POSIX$hour < 11,]
anomaly_observation3 <- anomaly_data3[anomaly_data3$POSIX$wday == 2 & anomaly_data3$POSIX$hour >= 9 
                                      & anomaly_data3$POSIX$hour < 11,]


# Remove incomplete days from training data 
anomaly_observation1 <- remove_incomplete_day(anomaly_observation1)
anomaly_observation2 <- remove_incomplete_day(anomaly_observation2)
anomaly_observation3 <- remove_incomplete_day(anomaly_observation3)


# Scale Global_intensity & Global_active_power
anomaly_observation1$Global_intensity_scaled <- scale(anomaly_observation1$Global_intensity, center = TRUE, scale = TRUE)
anomaly_observation1$Global_active_power_scaled <- scale(anomaly_observation1$Global_active_power, center = TRUE, scale = TRUE)

anomaly_observation2$Global_intensity_scaled <- scale(anomaly_observation2$Global_intensity, center = TRUE, scale = TRUE)
anomaly_observation2$Global_active_power_scaled <- scale(anomaly_observation2$Global_active_power, center = TRUE, scale = TRUE)

anomaly_observation3$Global_intensity_scaled <- scale(anomaly_observation3$Global_intensity, center = TRUE, scale = TRUE)
anomaly_observation3$Global_active_power_scaled <- scale(anomaly_observation3$Global_active_power, center = TRUE, scale = TRUE)


# Training
train_model <- depmix(response = list(Global_intensity_scaled ~ 1, Global_active_power_scaled ~ 1)
                      , data = train_observation, nstates = 12, ntimes = train_ntimes
                      , family = list(gaussian(), gaussian()))
fit_train_model <- fit(train_model)


# Anomaly observation test ntimes
anomaly_observation1_test_ntimes <- numeric(nrow(anomaly_observation1)/120)
for(i in 1:nrow(anomaly_observation1)/120){
  anomaly_observation1_test_ntimes[i] = 120
}

anomaly_observation2_test_ntimes <- numeric(nrow(anomaly_observation2)/120)
for(i in 1:nrow(anomaly_observation2)/120){
  anomaly_observation2_test_ntimes[i] = 120
}

anomaly_observation3_test_ntimes <- numeric(nrow(anomaly_observation3)/120)
for(i in 1:nrow(anomaly_observation3)/120){
  anomaly_observation3_test_ntimes[i] = 120
}


# Testing
anomaly_observation1_test_model <- depmix(response = list(Global_intensity_scaled ~ 1, Global_active_power_scaled ~ 1)
                                          , data = anomaly_observation1, nstates = 12
                                          , ntimes = anomaly_observation1_test_ntimes, family = list(gaussian(), gaussian()))

anomaly_observation2_test_model <- depmix(response = list(Global_intensity_scaled ~ 1, Global_active_power_scaled ~ 1)
                                          , data = anomaly_observation2, nstates = 12
                                          , ntimes = anomaly_observation2_test_ntimes, family = list(gaussian(), gaussian()))

anomaly_observation3_test_model <- depmix(response = list(Global_intensity_scaled ~ 1, Global_active_power_scaled ~ 1)
                                          , data = anomaly_observation3, nstates = 12
                                          , ntimes = anomaly_observation3_test_ntimes, family = list(gaussian(), gaussian()))

anomaly_observation1_test_model <- setpars(anomaly_observation1_test_model, getpars(fit_train_model))

anomaly_observation2_test_model <- setpars(anomaly_observation2_test_model, getpars(fit_train_model))

anomaly_observation3_test_model <- setpars(anomaly_observation3_test_model, getpars(fit_train_model))


# Log Likelihood 
forwardbackward(anomaly_observation1_test_model)$logLik
forwardbackward(anomaly_observation2_test_model)$logLik
forwardbackward(anomaly_observation3_test_model)$logLik

head(anomaly_observation1)
anomaly_observation1[nrow(anomaly_observation1),]$POSIX$wday
anomaly_observation1[1,]$POSIX$wday


