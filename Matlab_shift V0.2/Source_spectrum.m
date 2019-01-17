% Outputs a source spectrum. Input fftsize and fundamental frequency.

function [source_spectrum] = Source_spectrum(f0, Envelope)

global Fs

%% creating source spectrum

y_pos = zeros(Fs/2+1,1);

for kk = round(f0/2) : round(f0/2) : Fs/2
  y_pos(kk+1) = 1;
end

y_pos = y_pos.*Envelope(1 : Fs/2+1)';
y_pos = circshift(y_pos,-round(f0/2));
y_pos(end - round(f0/2):end) = 0;

y_neg = flipud(y_pos(2:end-1));

source_spectrum = [y_pos;y_neg]';


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


