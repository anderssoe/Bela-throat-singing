% Real-time throat singing simulation based on the Bela platform
% 
% By Rasmus Tinndahn & Anders Søe during B.Sc in Electrical Engineering at the Technical University of Denmark

%% Setup 
clc,
clear all;

% load test audio
[y, Fs] = audioread('SantaHoHo.wav');
y = y(1 : length(y) - mod(length(y), 128));
counter = 1;

%% Processing

% process sample-by-sample here:
for i = 1:128:length(y)
    
    FFT = fft(y(i : i +127));
    FFT(1 : 32) = FFT(33 : 64);
    FFT(97: 128) =FFT(65 : 96);
    FFT(1 : 32)  = zeros(32, 1);
    FFT(97: 128) = zeros(32, 1);
    y(i : i +127) = real(ifft(FFT));
    
end


%% Playback

% playback test audio
sound(y, Fs);