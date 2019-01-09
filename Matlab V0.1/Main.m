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
FFTsize = 1024*8; % size of FFT
binsize = SampleRate/FFTsize; % size of every bin
window = (hann(FFTsize))';
f0 = 146.83; % fundamental frequency (d)
M = FFTsize/2; % block size
f0_max = 180;
%sig_filter = designfilt('lowpassiir', 'FilterOrder', 4, 'HalfPowerFrequency', f0_max, 'SampleRate', SampleRate);
sig_filter = designfilt('bandpassiir', 'FilterOrder', 4, 'PassbandFrequency1', 80, 'PassbandFrequency2', 180, 'PassbandRipple', 1, 'SampleRate', 44100);

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
average_envelope = 0;
meanammount = 50;
spectral_envelope = zeros(meanammount, FFTsize);
output = 0;
buffer = zeros(1,M);
% process sample-by-sample here:
for i = 1: floor(length(y)/M) - 1;
    
    input = y(1 + i*M - M: i * M + M)';
    input = input .* window;
    f0 = pitchTrack_zc(input, length(input), Fs, sig_filter);
    %f0 = 200;
    [source_spectrum_fft, spectral_envelope(mod(meanammount, i) + 1, :)] = Source_spectrum_synthesis(f0, input, FFTsize);
    
    source_spectrum_ifft = ifft(source_spectrum_fft);
    source_spectrum_ifft = [source_spectrum_ifft ; zeros(1023, 1)];
    spectral_envelope_ifft = ifft(mean(spectral_envelope(mod(meanammount, i) + 1, :), 1));
    spectral_envelope_ifft = [spectral_envelope_ifft zeros(1, 1023)];
    
    source_spectrum_fft = fft(source_spectrum_ifft);
    spectral_envelope_mean = fft(spectral_envelope_ifft);
    
    input = [spectral_envelope_mean .* source_spectrum_fft']; %zero padding for circular convolution
    input = real(ifft(input));
    input = input(1 : FFTsize);
    %input = input';
    input = input .* window;
    
    
    output = [output (input(1 : 1 * M) + buffer)];

    buffer = input(1+ M : 2 * M);
    
    
    save(i, :) = spectral_envelope_mean;
end


%% Playback
output =  output ./ max(output);

audiowrite('Normal+phase.wav',output,Fs)

% playback test audio
sound(output', Fs);