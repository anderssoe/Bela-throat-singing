% Real-time throat singing simulation based on the Bela platform
% 
% By Rasmus Tinndahn & Anders S�e during B.Sc in Electrical Engineering at the Technical University of Denmark

%% Setup 
clc,
clear all;
close all;


global FFTsize BINsize Fs;
% load test audio
[Analysis_signal, Fs] = audioread('Kargyraa_O.wav');
[Voice_signal] = audioread('Voice_O.wav');
[Voice_signal] = Voice_signal(:,1);

%[Voice_signal] = audioread('Kargyraa_O.wav');


% Set constants
Fs = 44100; % sample rate
FFTsize = 1024*4; % size of FFT
BINsize = Fs/FFTsize; % size of every bin
phase = ones(1, Fs);

window = (hann(FFTsize, 'periodic'))';
M = FFTsize/2; % block size

% making the filter for the pitch tracking algorithm
%sig_filter = designfilt('lowpassiir', 'FilterOrder', 4, 'HalfPowerFrequency', f0_max, 'SampleRate', SampleRate);
sig_filter = designfilt('bandpassiir', 'FilterOrder', 4, 'PassbandFrequency1', 100, 'PassbandFrequency2', 180, 'PassbandRipple', 1, 'SampleRate', 44100);
sig_filter_lp = designfilt('lowpassiir', 'FilterOrder', 2, 'HalfPowerFrequency', 1200, 'SampleRate', 44100);;
% making w_lp for the spectral envelope
N1 = 75;
w_lp(1) = 1;
w_lp(2 : N1) = 2;
w_lp(N1 + 1 : FFTsize) = 0;

w_hp(1) = 0;
w_hp(2 : N1) = 0;
w_hp(N1 + 1 : FFTsize) = 2;
% %% Making a sinus at the length of 1 sec
% 
% % amplitude
% a = 1;
% 
% % frequency of 1st sinus
% f = 1000;
% 
% % phase
% phi_k = 0;
% 
% % duration in seconds
% T_s = 1;
% 
% % create 1st single sinusoid
% [time_vector Voice_signal] = generate_sinusoid(a, f, phi_k, Fs, T_s);
% Voice_signal = Voice_signal';
%% Processing

output = 0;
buffer = zeros(1,M);
Average_X = 5; %(mod(i, Average_X) + 1,:)
% process sample-by-sample here:
for i = 1: floor(length(Voice_signal)/M) - 1;
    
    % using the same spectral envelope for all iterations
    analysis = Analysis_signal(1: FFTsize)';
    
    % changing spectral envelope for every iteration.
    %analysis = Analysis_signal(1 + i*M - M: i * M + M)';
    
    input = Voice_signal(1 + i*M - M: i * M + M)';
   
    %f0 = 120;
    f0 = pitchTrack_zc(input, length(input), sig_filter);
    f0_vec(i) = f0; 
    input = [input zeros(1,Fs-length(input))];
    original_fft = fft(input);

    input = circshift(original_fft, -round(f0/2));
    
    fc1 = 600;
    fc2 = 500;
    filter = [ones(1,fc1) ones(1, fc2)*0.7 ones(1,length(input)-fc1-fc2)*0.1];
    original_fft = original_fft.* filter;
    input = input.*[ones(1,1000) ones(1,Fs-1000)*0.05];
    input = (input)*0.55 + (original_fft)*0.45;
    
    
    input(650:880) = input(650:880)*0.5;

    %theta = angle(original_fft);
    
    %input = input.*(cos(theta) + 1i*sin(theta));
    
    input(1:round(f0-f0/4)) = input(1:round(f0-f0/4)) * 0.1;
    
    
    [X] = Spectral_envelope(Analysis_signal(1:FFTsize)', w_lp);
    
    %input = input .* X;
    
    input(1) = 0;
    input(Fs/2 + 1) = 0;
    
    input(Fs/2+2 : end) = conj(fliplr(input(2 : Fs/2)));
    
    input = ifft(input);
  
    
    input = input(1 : FFTsize);
    input = input .* window;
    
    %input = input / rms(input);
    output = [output (input(1 : 1 * M) + buffer)];

    buffer = input(1+ M : 2 * M);
    
    
end


%% Playback
output =  0.9 .* output ./ max(output);

audiowrite('Analysis/SpectralShift.wav',output,Fs)
% playback test audio
sound(output', Fs);


%%

plot(mag2db(abs(fft(output,Fs))), 'LineWidth', 2.3);
hold on
plot(mag2db(abs(fft(Analysis_signal,Fs))), 'LineWidth', 1.6);
%plot(abs(fft(Voice_signal,Fs)), 'LineWidth', 1.5);
xlim([0 1400]);
legend('output','kargyraa_O','input')

saveas(gcf, 'Analysis/Plots/SpectralShiftComparison', 'epsc')

