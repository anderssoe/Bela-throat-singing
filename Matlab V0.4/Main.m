% Real-time throat singing simulation based on the Bela platform
% 
% By Rasmus Tinndahn & Anders Søe during B.Sc in Electrical Engineering at the Technical University of Denmark

%% Setup 
clc,
clear all;


global FFTsize BINsize Fs;
% load test audio
[Analysis_signal, Fs] = audioread('Kargyraa_O.wav');
[Voice_signal] = audioread('Voice_O.wav');
[Voice_signal] = Voice_signal(:,1);

%[Voice_signal] = audioread('Kargyraa_O.wav');


% Set constants
Fs = 44100; % sample rate
FFTsize = 1024; % size of FFT
BINsize = Fs/FFTsize; % size of every bin
phase = ones(1, Fs);

window = (hann(FFTsize, 'periodic'))';
M = FFTsize/2; % block size

% making the filter for the pitch tracking algorithm
%sig_filter = designfilt('lowpassiir', 'FilterOrder', 4, 'HalfPowerFrequency', f0_max, 'SampleRate', SampleRate);
sig_filter = designfilt('bandpassiir', 'FilterOrder', 4, 'PassbandFrequency1', 80, 'PassbandFrequency2', 180, 'PassbandRipple', 1, 'SampleRate', 44100);

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
    input = input .* window;
    f0(mod(i, Average_X) + 1) = pitchTrack_zc(input, length(input), sig_filter);
    %[X(mod(i, Average_X) + 1,:)] = Spectral_envelope(input, w_lp);
    
    % spectral envelope
    H = Spectral_envelope(input, w_lp);
    
    % source envelope
    X = Spectral_envelope(input, w_hp);
    %if(i > 2)
    %    f0(i) = (f0(i) + f0(i-1) + f0(i-2))/3;
    %end

    [SourceSpectrum] = Source_spectrum(mean(f0), X);
    
    input = H .* X;
    %input = SourceSpectrum;
    input = real(ifft(input));%.*((cos(phase) + 1i * sin(phase)))));

    %input = ones(1, 1024);
    input = fftshift(input);
    %phase = angle(fft(circshift(input, 512)));

    input = input(1 : FFTsize);
    input = input .* window;
    
    %input = input / rms(input);
    output = [output (input(1 : 1 * M) + buffer)];

    buffer = input(1+ M : 2 * M);
    
    
end


%% Playback
output =  0.9 .* output ./ max(output);

audiowrite('Analysis/VoiceEnv.wav',output,Fs)
% playback test audio
sound(output', Fs);