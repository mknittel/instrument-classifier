#!/usr/bin/env Rscript
library('tuneR')
library('seewave')
library('data.table')
library('ProjectTemplate')
library('caret')
library('kknn')
library('class')
library('neuralnet')
library('base')
library('e1071')
library('pROC')

knn_predict <- function(train, valid) {
    train_pred <- knn(train[, 1:5], valid[, 1:5], train[, 6], k=11)

    acc <- c()

    for (i in 1:length(train_pred)) {
        acc <- c(acc, train_pred[i] == valid[i, 6])
    }

    return(c(acc, train_pred))
}

find_k <- function(train, valid) {
    k <- c()
    knn_accs <- c()

    for (i in 1:49) {
        k <- c(k, i)
        train_pred <- knn(train[, 1:5], valid[, 1:5], train[, 6], k=i)
        knn_acc <- c()

        for (j in 1:length(train_pred)) {
            knn_acc <- c(knn_acc, train_pred[j] == (valid[, 6])[j])
        }

        knn_accs <- c(knn_accs, 1 - sum(knn_acc) / length(train_pred))
    }

    return(c(k, knn_accs))
}

net_predict <- function(train, valid) {
    net <- neuralnet(instr~r1+r2+r3+r4+fund_freqs, train, hidden=0)
    valid_pred <- compute(net, valid[, 1:5])$net.result

    acc <- c()

    for (i in 1:length(valid_pred)) {
        acc <- c(acc, round(valid_pred[i]) == valid[i, 6])
    }

    return(c(acc, valid_pred))
}

bayes_predict <- function(train, valid) {
    net <- naiveBayes(as.factor(instr)~., data=train)
    valid_pred <- predict(net, valid[, 1:5])

    acc <- c()

    for (i in 1:length(valid_pred)) {
        acc <- c(acc, valid_pred[i] == valid[i, 6])
    }

    return(c(acc, valid_pred))
}

set_feats_to_ratios <- function(data) {
    data$f2 <- scale(data$f2/data$f1, center = FALSE, scale = stddev[1]) 
    data$f3 <- scale(data$f3/data$f1, center = FALSE, scale = stddev[1]) 
    data$f4 <- scale(data$f4/data$f1, center = FALSE, scale = stddev[1]) 
    data$f5 <- scale(data$f5/data$f1, center = FALSE, scale = stddev[1]) 

    return(data)
}

process_data <- function(files) {
    # Data setup
    f1 <- c()
    f2 <- c()
    f3 <- c()
    f4 <- c()
    f5 <- c()
    fund_freqs <- c()
    instr <- c()

    for (i in 1:(sum(complete.cases(files)) - 1)) {
        in_directory <- files[[i, 4]]
        note <- files[[i,3]]
        octave <- files[[i,2]]
        this_instr <- files[[i,1]]

        file <- paste("data/", in_directory, sep="")
        feats <- calc_feats(file, note, octave)

        f1[i] <- feats[1]
        f2[i] <- feats[2]
        f3[i] <- feats[3]
        f4[i] <- feats[4]
        f5[i] <- feats[5]
        fund_freqs[i] <- feats[6]
        instr[i] <- this_instr

        if (i %% 20 == 0) {
            print(i)
        }
        if (i == 2000) {
            break
        }
    }
    r1 <- f2/f1
    r2 <- f3/f1
    r3 <- f4/f1
    r4 <- f5/f1

    data <- data.frame(r1, r2, r3, r4, fund_freqs, instr)

    return(data)    
}

calc_feats <- function(file, note, octave) {
    wave <- readMP3(file)

    5292 / wave@samp.rate
    s1 <- wave@left 

    n <- length(s1)
    p <- fft(s1)

    nUniquePts <- ceiling((n+1)/2)
    p <- p[1:nUniquePts]
    p <- abs(p)
    p <- p / n
    p <- p^2

    if (n %% 2 > 0) {
        p[2:length(p)] <- p[2:length(p)]*2
    } else {
        p[2:(length(p) - 1)] <- p[2:(length(p) - 1)]*2
    }

    freqArray <- (0:(nUniquePts-1)) * (wave@samp.rate / n)
    powers <- 10 * log10(p)

    #length <- length(freqArray) 
    #freqs1 <- head(freqArray, length)
    #p1 <- head(powers, length)

    #plot(freqs1 / 1000, p1, type='l', col='black', xlab='Frequency (kHz)', ylab='Power (dB)')

    fund_freq <- note_from_freq(note, octave)
    peaks <- find_peaks(powers, fund_freq, freqArray, 5)

    #f1 <- powers[peaks[1]]
    #f2 <- powers[peaks[2]]
    #f3 <- powers[peaks[3]]
    #f4 <- powers[peaks[4]]
    #f5 <- powers[peaks[5]]
    #ret <- list(f1, f2, f3, f4, f5, fund_freq)
    #print(append(peaks, fund_freq))

    return(append(peaks, fund_freq))
}

find_peaks <- function(data, fund_freq, freqs, npeaks) {
    if (npeaks == 0) {
        return(c())
    }

    peaks <- c()
    harmonic <- 1

    for (i in 1:npeaks) {
        freq <- fund_freq * harmonic
        index <- get_closest(freq, freqs)
        peaks <- c(index, peaks)

        harmonic <- harmonic + 1
    }

    return(peaks)
}

note_from_freq <- function(note, octave) {
    table <- read.table("note_freqs.txt")
    freq <- 16.35

    note <- as.character(note)
    octave <- as.character(octave)

    for (i in 1:(sum(complete.cases(table)) - 1)) {
        this_note <- table[[i, 1]]
        this_octave <- table[[i, 2]]

        if (note == this_note && octave == this_octave) {
            freq <- table[[i,3]]
            break
        }
    }

    return(freq)
}

get_subset <- function(data, indices) {
    ret <- c()

    for (i in 1:length(indices)) {
        index <- indices[i]
        ret <- c(c(data[index]), ret)
    }

    return(ret)
}

get_closest <- function(element, lst) {
    close_dist <- abs(element - lst[1])
    close_index <- 1

    for (i in 1:(length(lst) - 1)) {
        dist <- abs(element - lst[i])
        
        if (dist < close_dist) {
            close_dist <- dist
            close_index <- i
        }
    }

    return(close_index)
}

files <- read.table('data.txt')
data <- process_data(files)

data[, 1:5] <- scale(data[, 1:5])

samp_size <- floor(.8 * nrow(data))
set.seed(123)
train_ind <- sample(seq_len(nrow(data)), size=samp_size)

test <- data[-train_ind, ]
train_partition <- data[train_ind, ]

samp_size <- floor(.8 * nrow(train_partition))
train_ind <- sample(seq_len(nrow(train_partition)), size=samp_size)

train <- train_partition[train_ind, ]
valid <- train_partition[-train_ind, ]

ret1 <- knn_predict(train, valid)
knn_acc <- ret1[1]
knn_pred <- ret1[2]

ret2 <- net_predict(train, valid)
net_acc <- ret2[1]
net_pred <- ret2[2]

ret3 <- bayes_predict(train, valid)
bayes_acc <- ret3[1]
bayes_pred <- ret3[2]

