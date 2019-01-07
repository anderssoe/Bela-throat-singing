% Real-time throat singing simulation based on the Bela platform
% 
% By Rasmus Tinndahn & Anders Søe during B.Sc in Electrical Engineering at the Technical University of Denmark

%% Setup 
clc,
clear all;

% load test audio
[y, Fs] = audioread('SantaHoHo.wav');
y = y(:,1);

SampleRate = 44100; % sample rate
FFTsize = 1024; % size of FFT
binsize = SampleRate/FFTsize; % size of every bin
window = sqrt(hann(1024));

%% Processing

% process sample-by-sample here:
for i = 1: floor(length(y)/FFTsize)
    y(1 + i*FFTsize - FFTsize: i * FFTsize) 
    
end


%% Playback

% playback test audio
sound(y, Fs);