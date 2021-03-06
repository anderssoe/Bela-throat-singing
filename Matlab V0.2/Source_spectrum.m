% Outputs a source spectrum. Input fftsize and fundamental frequency.

function [source_spectrum] = Source_spectrum(f0)

global Fs

%% creating source spectrum

source_spectrum = zeros(1, Fs/2);

for i = round(f0) : round(f0) : Fs/2
    source_spectrum(i) = 1/(i^2);
end

source_spectrum(Fs/2 + 1 : Fs) = fliplr(source_spectrum);

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


