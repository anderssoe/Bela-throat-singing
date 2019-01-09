clear all
close all
clc

%% Test
[signal, Fs] = audioread('SantaHoHo.wav');
signal = signal(:,1)';


%% Function
%function [out] = FD_SOURCE_SPECTRUM(signal, P, N, f0, Fs)


%TD-PSOLA

%s[k]: input signal of length M
%N: block size
%H: hop size
%P[k] = the glottal pulse train
%T: offset for pulse train¨
sig_filter = designfilt('bandpassiir', 'FilterOrder', 4, 'PassbandFrequency1', 80, 'PassbandFrequency2', 180, 'PassbandRipple', 1, 'SampleRate', 44100);

M = length(signal);
N = 1024;
H = N/2;
win = hann(N)';

w_lp(1) = 1;
w_lp(2 : N/2) = 2;
w_lp(N/2 + 1 : N) = 0;

T = 0;



i = 1;

%% Begin processing
while i + N < M
    %Extract x as windowed segment of input
    x = win .* signal(i: i+N -1);
    
    % Calculate spectral envelope
    Xm =  fft(x);
    theta = angle(Xm);
    Xc = ifft(log(10^(-9) + Xm));
    Xc = Xc .* w_lp;
    spectral_envelope = exp(real(fft(Xc))).*(cos(theta) + 1i * sin(theta));
    
    % Calculate periodicity in samples
    [f0] = pitchTrack_zc(x ,N , Fs, sig_filter);
    %f0 = 150;
    d = f0 * N / Fs;
    
    P = zeros(1,N/2);
    % Calculate number of pulses, np, and offset T for next block
    np = ceil((N/2-T)/d);
    T = np * d - N;
    summen = [];
    for k = 1 : N/2
        for n = 0 : np
             summen = summen + dirac(k - T - n * d);
        end
        P(k) = summen
    end
%     %% Make pulse train
%     P = [];
%     highestHarmonic = Fs/2;
%     harmonicStep = f0;
%     ii = 1;
%     while harmonicStep <= highestHarmonic
%     source_spectrum(ii) = 1 / harmonicStep^2;
%     harmonicStep = harmonicStep + f0;
%     ii = ii+1;
%     end
    
    % Add negative frequencies to P 
    P(N/2 + 1 : N) = fliplr(P(1 : N/2));
    
    temp = length(P);
    P = [P zeros(1, length(spectral_envelope))];
    spectral_envelope = [spectral_envelope zeros(1, temp)];
   
    y(i : i + 2 * N-1) = fftshift(ifft(P.*spectral_envelope));
    
    i = i + H;
end

    
