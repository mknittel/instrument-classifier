#!/usr/bin/env Rscript
library('tuneR')
library('seewave')
library('data.table')
library('ProjectTemplate')
library('caret')
library('kknn')
library('class')

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

get_subset <- function(data, indices) {
    ret <- c()

    for (i in 1:length(indices)) {
        index <- indices[i]
        ret <- c(c(data[index]), ret)
    }

    return(ret)
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
        if (i == 10) {
            break
        }
    }

    data <- data.frame(f1, f2, f3, f4, f5, fund_freqs, instr)

    return(data)    
}

main <- function() {
    files <- read.table('data.txt')
    data <- process_data(files)

    samp_size <- floor(.8 * nrow(data))
    set.seed(123)
    train_ind <- sample(seq_len(nrow(data)), size=samp_size)

    train_data <- data[train_ind, 1:6]
    test_data <- data[-train_ind, 1:6]
    train_label <- data[train_ind, 7]
    test_label <- data[-train_ind, 7]

    print(typeof(train_label))
    print(typeof(train_data))

    train_pred <- knn.cv(train = train_data, cl = train_label)
    print(train_pred)
    print(train_label)

    acc <- c()

    for (i in 1:length(train_pred)) {
        acc <- c(acc, train_pred[i] == train_label[i])
    }

    print(acc)

    print("done")
}

main()

