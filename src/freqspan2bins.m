function [ nFreq ] = freqspan2bins(fs, nfft, freq_span)
%freqspan2bins finds the number of bins in the given frequency span
nFreq = freq_span .* nfft/fs;