% Outputs a source spectrum. Input fftsize and fundamental frequency.

function [X] = Source_envelope(y, w_hp)
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

%X  = interp1(1 : BINsize : Fs , X , 1 : 1 : Fs,'linear',0);
%X = smooth(X,20,'moving');
%X = X';
X(Fs/2+1) = 0;
X(Fs/2 + 2 : end) = fliplr(X(2:Fs/2));



% source_spectrum = zeros(1, Fs/2);
% 
% for idx = round(f0/2): round(f0/2) : Fs/2
%     source_spectrum(i+1) = 1/(idx^2)
% end
% 
% source_spectrum(Fs/2 + 1 : Fs) = [fliplr(source_spectrum(2 : end))];

% 
% for i = round(f0/2) : round(f0/2) : Fs/2
%     source_spectrum(i+1) = 1/(i^2);
% end
% 
% source_spectrum(Fs/2 + 1 : Fs) = [fliplr(source_spectrum(2 : end)) 0];

% % while loop forming the source spetrum up to the chosen highest harmonic.
% highestHarmonic = 10000;
% ii = 1;
% while harmonicStep <= highestHarmonic
%     source_spectrum(ii) = 1 / harmonicStep^2;
%     harmonicStep = harmonicStep + f0;
%     ii = ii+1;
% end

% %% plot spectrum
% figure('Name', 'Source spectrum')
% stem(S)
% set(gca, 'YScale', 'log')
% xt = get(gca, 'XTick');
% set(gca, 'XTick', xt, 'XTickLabel', xt*f0)
% xtickangle(45)

%% creating source spectrum one octave lower

% harmonicStep = f0/2;
% highestHarmonic = 10000;
% ii = 1;
% while harmonicStep <= highestHarmonic
%     if mod(ii, 2) == 0
%         source_spectrum_octave(ii) = source_spectrum(ii/2);
%     else
%         source_spectrum_octave(ii) = 1 / harmonicStep^2;
%     end
%     harmonicStep = harmonicStep + f0/2;
%     ii = ii+1;
% end
% 


%% adding source spectrum to FFT-bins

% % while loop converting the source spectrum to match the chosen FFTsize.
% harmonicStep = f0/2;
% ii = 1;
% j = 1;
% while harmonicStep <= 10000
%     while (harmonicStep >= binsize * ii - binsize / 2) && (harmonicStep < binsize * ii + binsize / 2)
%         source_spectrum_fft(ii) = source_spectrum_fft(ii) + source_spectrum_octave(j);
%         harmonicStep = harmonicStep + f0;
%         j = j + 1;
%     end
%     ii = ii + 1;
% end
% 
% label('Amplitude')


