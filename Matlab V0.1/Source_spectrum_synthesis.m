% Outputs a source spectrum. Input fftsize and fundamental frequency.

function [source_spectrum_fft, spectral_envelope] = Source_spectrum_synthesis(f0, y, FFTsize)

%% Setup
%close all
%clear all
%clc

SampleRate = 44100; % sample rate
binsize = SampleRate/FFTsize; % size of every bin
source_spectrum_fft = zeros(FFTsize, 1); % source spectrum, fourier transformed
%
harmonicStep = f0; % initialize the step of the harmonics


%% load audio for testing
% loading 1024 samples for this test, the phoneme (o:) with a fundamental
% frequency of 146.83 Hz
%[y] = audioread('SantaHoHo.wav');
%y1 = y(90001: 91024);
%y = y(90001: 91024);
%y_fft = fft(y); %

%% creating source spectrum

% while loop forming the source spetrum up to the chosen highest harmonic.
highestHarmonic = 10000;
ii = 1;
while harmonicStep <= highestHarmonic
    source_spectrum(ii) = 1 / harmonicStep^2;
    harmonicStep = harmonicStep + f0;
    ii = ii+1;
end

% %% plot spectrum
% figure('Name', 'Source spectrum')
% stem(S)
% set(gca, 'YScale', 'log')
% xt = get(gca, 'XTick');
% set(gca, 'XTick', xt, 'XTickLabel', xt*f0)
% xtickangle(45)

%% creating source spectrum one octave lower

harmonicStep = f0/2;
highestHarmonic = 10000;
ii = 1;
while harmonicStep <= highestHarmonic
    if mod(ii, 2) == 0
        source_spectrum_octave(ii) = source_spectrum(ii/2);
    else
        source_spectrum_octave(ii) = 1 / harmonicStep^2;
    end
    harmonicStep = harmonicStep + f0/2;
    ii = ii+1;
end



%% adding source spectrum to FFT-bins

% while loop converting the source spectrum to match the chosen FFTsize.
harmonicStep = f0/2;
ii = 1;
j = 1;
while harmonicStep <= 10000
    while (harmonicStep >= binsize * ii - binsize / 2) && (harmonicStep < binsize * ii + binsize / 2)
        source_spectrum_fft(ii) = source_spectrum_fft(ii) + source_spectrum_octave(j);
        harmonicStep = harmonicStep + f0;
        j = j + 1;
    end
    ii = ii + 1;
end

% flipping spectrum for fft.
%S_fft = [S_fft(1:FFTsize/2) ; flip(S_fft(1:FFTsize/2))];

%
% figure('Name', 'Source spectrum converted to size of FFT')
% stem(S_fft)
% set(gca, 'YScale', 'log')
% xt = get(gca, 'XTick');
% set(gca, 'XTick', xt, 'XTickLabel', round(xt*binsize - binsize/2))
% xtickangle(45)
% xlabel('Frequency [Hz]')
% ylabel('Amplitude')

% while loop converting the source spectrum to match the chosen FFTsize.
% This one distributes the harmonics to bins around the harmonic.
%harmonicStep = f0
% i = 1; j = 1; while harmonicStep <= 10000
%
%    while (harmonicStep >= binsize * i - binsize / 2) && (harmonicStep <
%    binsize * i + binsize / 2)
%
%         distance = harmonicStep - i * binsize factor = (abs(distance) /
%         binsize);
%
%         if i == 1
%             S_fft(i) = S_fft(i) + S(j);
%         elseif distance > 0
%             S_fft(i) = S_fft(i) + S(j) * (1 - factor); S_fft(i) = S_fft(i
%             + 1) + S(j) * factor;
%         else
%             S_fft(i) = S_fft(i) + S(j) * (1 - factor); S_fft(i) = S_fft(i
%             - 1) + S(j) * factor;
%         end harmonicStep = harmonicStep + f0; j = j + 1;
%     end i = i + 1;
% end

%% Calculate the envelope spectrum
% Calculations from special course paper TODO: document this.

Xm =  fft(y);
theta = angle(Xm);
Xc = ifft(log(10^(-9) + Xm));

w_lp(1) = 1;
w_lp(2 : FFTsize/2) = 2;
w_lp(FFTsize/2 + 1 : FFTsize) = 0;

Xc = Xc .* w_lp;

spectral_envelope = exp(real(fft(Xc)));%.*(cos(theta) + 1i * sin(theta));



%Xenv = 10 * log(Xenv ./ max(Xenv));

%figure(4)
%plot(Xenv)

%% Calculate the transfer function (filter / throat)


% calculating the tf output / input = TF

%y = real(ifft(Xenv .* S_fft'));
%y = y / max(y) * 0.6;

%% plotting Source and envelope
% figure(10)
% plot(Xenv / max(Xenv))
% hold on
% plot(S_fft' / (max(S_fft)))
% %plot(TF)
% set(gca, 'YScale', 'log')


