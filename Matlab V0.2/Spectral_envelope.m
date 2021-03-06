function [X] = Spectral_envelope(y, w_lp)
% y is a input signal of size (1, FFTsize)

% [signal, Fs] = audioread('SantaHoHo.wav');
% signal = signal(:,1);
% y = signal(90001:90001 + 1024)';

global BINsize Fs;
Xm =  fft(y);

theta = angle(Xm);

Xc = ifft(log(10^(-9) + Xm));

Xc = Xc .* w_lp;

X = exp(real(fft(Xc))).*(cos(theta) + 1i * sin(theta));

%%begin interpolation

X  = interp1(1 : BINsize : Fs , X , 1 : 1 : Fs,'linear',0);






