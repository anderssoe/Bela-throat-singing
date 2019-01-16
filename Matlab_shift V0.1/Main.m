% Real-time throat singing simulation based on the Bela platform
% 
% By Rasmus Tinndahn & Anders S�e during B.Sc in Electrical Engineering at the Technical University of Denmark

%% Setup 
clc,
clear all;
close all;

global Fs
% load test audio
[Analysis_signal, Fs] = audioread('Kargyraa_O.wav');
[Voice_signal] = audioread('Voice_O.wav');
[Voice_signal] = Voice_signal(:,1);

%[Voice_signal] = audioread('Kargyraa_O.wav');


% Set constants
Fs = 44100; % sample rate
FFTsize = 44100/10;
Blocksize = 1024*3; % size of FFT
BINsize = Fs/Blocksize; % size of every bin
Hishelf = 1000; %-30dB highshelf in hz;
dryWet = 0.55; %percentage of dry to wet (1 = completely dry)

window = (hann(Blocksize, 'periodic'))';
M = Blocksize/2; % block size

% making the filter for the pitch tracking algorithm
sig_filter = designfilt('bandpassiir', 'FilterOrder', 4, 'PassbandFrequency1', 100, 'PassbandFrequency2', 180, 'PassbandRipple', 1, 'SampleRate', 44100);

% making LP filter for final processing
sig_filter_lp = designfilt('lowpassiir', 'FilterOrder', 2, 'HalfPowerFrequency', 1200, 'SampleRate', 44100);
%% Processing

output = 0;
% Overlap buffer
buffer = zeros(1,M);
% process block-to-block
for i = 1: floor(length(Voice_signal)/M) - 1

    input = Voice_signal(1 + i*M - M: i * M + M)';
    
    % Track pitch per block
    f0 = pitchTrack_zc(input, length(input), sig_filter);
    f0_vec(i) = f0; 
    
    
    %zeropadding input block to Fs Hz for 1Hz resolution
    input = [input zeros(1,FFTsize-length(input))];
    
    %FFT input block
    original_fft = fft(input);
    
    
    input = circshift(original_fft, -round((f0*(FFTsize/Fs)+1)/2));
    
    
    % Hishelf on the pitch shifted version
    fc1 = round(600*(FFTsize/Fs)+1);
    fc2 = round(500*(FFTsize/Fs)+1);
    filter = [ones(1,fc1) ones(1, fc2)*0.7 ones(1,length(input)-fc1-fc2)*0.1];
    original_fft = original_fft.* filter;
    
    % Hishelf from 1k on the input
    input = input.*[ones(1,round(Hishelf*(FFTsize/Fs)+1)) ones(1,round(FFTsize - (Hishelf*(FFTsize/Fs)+1)))*0.05];
    % Mix dry and wet signal
    input = (input)*dryWet + (original_fft)*(1-dryWet);
    
    % Notch filter between 650 and 880Hz
    input(round(650*(FFTsize/Fs)+1):round(880*(FFTsize/Fs)+1)) = input(round(650*(FFTsize/Fs)+1):round(880*(FFTsize/Fs)+1))*0.5;
       
    input(1:round((f0*(FFTsize/Fs)+1)-(f0/4 *(FFTsize/Fs)+1))) =  input(1:round((f0*(FFTsize/Fs)+1)-(f0/4 *(FFTsize/Fs)+1))) * 0.1;

    input(1) = 0;
    input(FFTsize/2 + 1) = 0;
    
    input(FFTsize/2+2 : end) = conj(fliplr(input(2 : FFTsize/2)));
    
    input = ifft(input);
  
    
    input = input(1 : Blocksize);
    input = input .* window;
    
    output = [output (input(1 : 1 * M) + buffer)];

    buffer = input(1+ M : 2 * M);
    
    
end


%% Playback
% output =  0.9 .* output ./ max(abs(output));
% 
% audiowrite('../Analysis/SpectralShift_RT.wav',output,Fs)
% % playback test audio
% sound(output', Fs);


%%
moving = 4;
plot(movmean(mag2db(abs(fft(output,Fs))),moving), 'LineWidth', 2.3);
hold on
plot(movmean(mag2db(abs(fft(Analysis_signal,Fs))),moving),'LineStyle','-.', 'LineWidth', 1.2);
%plot(abs(fft(Voice_signal,Fs)), 'LineWidth', 1.5);
xlim([0 1000]);
legend('output','kargyraa_O','input')
title(sprintf('Output vs targeted, %d Blocksize, MA = %d',Blocksize,moving))

saveas(gcf, sprintf('../Analysis/Plots/SpectralShiftComparison%d',Blocksize), 'epsc')

