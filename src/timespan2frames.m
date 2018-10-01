function [ nFrames ] = timespan2frames(fs, nwin, ninc, time_span)
%TIMESPAN2FRAMES finds the number of frames in the given time span
nFrames = (((time_span.*fs)-nwin) ./ ninc) + 1;
