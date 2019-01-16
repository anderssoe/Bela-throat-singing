function [X] = Spectral_envelope(y, w_lp)
% y is a input signal of size (1, FFTsize)

% [signal, Fs] = audioread('SantaHoHo.wav');
% signal = signal(:,1);
% y = signal(90001:90001 + 1024)';

global BINsize Fs FFTsize;
Xm =  fft(y);

Xc = ifft(log(10^(-18) + abs(Xm)));

Xc = Xc .* w_lp;

X = exp(real(fft(Xc)));%.*normrnd(0,pi/2,[1,FFTsize]);%.*(cos(theta) + 1i * sin(theta));

%X = real(fft(Xc));
%%begin interpolation

X  = interp1(1 : BINsize : Fs , X , 1 : 1 : Fs,'linear',0);
%X = smooth(X,20,'moving');
%X = X';
X(Fs/2+1) = 0;
X(Fs/2 + 2 : end) = fliplr(X(2:Fs/2));







