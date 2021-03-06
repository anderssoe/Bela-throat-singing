% Real-time throat singing simulation based on the Bela platform
% 
% By Rasmus Tinndahn & Anders S�e during B.Sc in Electrical Engineering at the Technical University of Denmark

%% Setup 
clc,
clear all;

global FFTsize BINsize Fs;

% load test audio
[signal, Fs] = audioread('SantaHoHo.wav');
signal = signal(:,1);

% Set constants
Fs = 44100; % sample rate
FFTsize = 1024*4; % size of FFT
BINsize = Fs/FFTsize; % size of every bin

window = (hann(FFTsize))';
M = FFTsize/2; % block size

% making the filter for the pitch tracking algorithm
%sig_filter = designfilt('lowpassiir', 'FilterOrder', 4, 'HalfPowerFrequency', f0_max, 'SampleRate', SampleRate);
sig_filter = designfilt('bandpassiir', 'FilterOrder', 4, 'PassbandFrequency1', 80, 'PassbandFrequency2', 180, 'PassbandRipple', 1, 'SampleRate', 44100);

% making w_lp for the spectral envelope
w_lp(1) = 1;
w_lp(2 : FFTsize/2) = 2;
w_lp(FFTsize/2 + 1 : FFTsize) = 0;

%% Making a sinus at the length of 1 sec

% % amplitude
% a = 1;
% 
% % frequency of 1st sinus
% f = 200;
% 
% % phase
% phi_k = 0;
% 
% % duration
% sec = 1
% T_s = sec * SampleRate/Fs;
% 
% % create 1st single sinusoid
% [time_vector y] = generate_sinusoid(a, f, phi_k, Fs, T_s);
% y = y';
%% Processing

output = 0;
buffer = zeros(1,M);
Average_X = 10;
% process sample-by-sample here:
for i = 1: floor(length(signal)/M) - 1;
    
    input = signal(1 + i*M - M: i * M + M)';
    input = input .* window;
    
    f0 = pitchTrack_zc(input, length(input), sig_filter);
    [X(mod(Average_X,i)+1,:)] = Spectral_envelope(input, w_lp);
    [SourceSpectrum] = Source_spectrum(f0);
    
    input = mean(X, 1) .* SourceSpectrum;
    input = real(ifft(input));
    input = input(1 : FFTsize);
    %input = input';
    input = input .* window;
    
    input = input / rms(input);
        
    output = [output (input(1 : 1 * M) + buffer)];

    buffer = input(1+ M : 2 * M);
    
    
end


%% Playback
output =  output ./ max(output);

%audiowrite('Normal+phase.wav',output,Fs)
% playback test audio
sound(output', Fs);