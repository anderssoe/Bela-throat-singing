% Real-time throat singing simulation based on the Bela platform
% 
% By Rasmus Tinndahn & Anders S�e during B.Sc in Electrical Engineering at the Technical University of Denmark

%% Setup 
clc,
clear all;

% load test audio
[y, Fs] = audioread('SantaHoHo.wav');
y = y(:,1);

SampleRate = 44100; % sample rate
FFTsize = 1024; % size of FFT
binsize = SampleRate/FFTsize; % size of every bin
window = sqrt(hann(1024))';
f0 = 146.83; % fundamental frequency (d)
M = FFTsize/2; % block size

% 
% %% Making a sinus at the length of 256 amples
% 
% % amplitude
% a = 1;
% 
% % frequency of 1st sinus
% f = 500;
% 
% % frequency of 2nd sinus
% f2 = 1100;
% 
% % phase
% phi_k = 0;
% 
% % duration
% T_s = 20000/Fs;
% 
% % create 1st single sinusoid
% [time_vector y] = generate_sinusoid(a, f, phi_k, Fs, T_s);
% y = y'
%% Processing
output = 0;
buffer = zeros(1,M);
% process sample-by-sample here:
for i = 1: floor(length(y)/M) - 1;
    
    input = y(1 + i*M - M: i * M + M)';
    input = input .* window;

    [source_spectrum, source_spectrum_fft, spectral_envelope] = Source_spectrum_synthesis(f0, input);
    input = real(ifft(spectral_envelope .* source_spectrum_fft'));
    
    %input = input';
    input = input .* window;
    
    
    output = [output (input(1 : 1 * M) + buffer)];

    buffer = input(1+ M : 2 * M);
    
    
end


%% Playback
output =  output ./ max(output);
% playback test audio
sound(output', Fs);