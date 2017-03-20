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

main <- function() {
    files <- read.table('data.txt')
    data <- process_data(files)

    stddev = sapply(data, sd)
    data$f1 <- scale(data$f1, center = TRUE, scale = stddev[1])
    data$f2 <- scale(data$f2, center = TRUE, scale = stddev[2])
    data$f3 <- scale(data$f3, center = TRUE, scale = stddev[3])
    data$f4 <- scale(data$f4, center = TRUE, scale = stddev[4])
    data$f5 <- scale(data$f5, center = TRUE, scale = stddev[5])
    data$fund_freqs <- scale(data$fund_freqs, center = TRUE, scale = stddev[6])

    samp_size <- floor(.8 * nrow(data))
    set.seed(123)
    train_ind <- sample(seq_len(nrow(data)), size=samp_size)

    test_data <- data[-train_ind, 1:6]
    test_label <- data[-train_ind, 7]
    test <- data[-train_ind, ]

    train_partition <- data[train_ind, ]

    samp_size <- floor(.8 * nrow(train_partition))
    train_ind <- sample(seq_len(nrow(train_partition)), size=samp_size)

    train_data <- train_partition[train_ind, 1:6]
    valid_data <- train_partition[-train_ind, 1:6]
    train_label <- train_partition[train_ind, 7]
    valid_label <- train_partition[-train_ind, 7]

    train <- train_partition[train_ind, ]
    valid <- train_partition[-train_ind, ]

    knn_acc <- knn_pred(train_data, train_label)
    net_acc <- net_pred(train, valid_data, valid_label)

    print("done")
}

knn_pred <- function(train_data, train_label) {
    train_pred <- knn.cv(train_data, train_label)

    acc <- c()
    
    for (i in 1:length(train_pred)) {
        acc <- c(acc, train_pred[i] == train_label[i])
    }

    return(acc)
}

net_pred <- function(train, valid_data, valid_label) {
    net <- neuralnet(instr~f1+f2+f3+f4+f5+fund_freqs, train, hidden=10)
    valid_pred = compute(net, valid_data)$net.result

    acc <- c()

    for (i in 1:length(valid_pred)) {
        acc <- c(acc, round(valid_pred[i]) == valid_label[i])
    }

    return(acc)
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
        feats <- calc_feats(file)

        f1 <- c(f1, feats[1])
        f2 <- c(f2, feats[2])
        f3 <- c(f3, feats[3])
        f4 <- c(f4, feats[4])
        f5 <- c(f5, feats[5])
        fund_freqs <- c(fund_freqs, feats[6])
        instr <- c(instr, this_instr)

        if (i %% 20 == 0) {
            print(i)
        }
        if (i == 2000) {
            break
        }
    }

    data <- data.frame(f1, f2, f3, f4, f5, fund_freqs, instr)

    return(data)    
}

calc_feats <- function(file) {
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

    f1 <- powers[peaks[1]]
    f2 <- powers[peaks[2]]
    f3 <- powers[peaks[3]]
    f4 <- powers[peaks[4]]
    f5 <- powers[peaks[5]]

    return(list(f1, f2, f3, f4, f5, fund_freq))
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

main()

