clear all
close all
clc

%% Analysis of kargyraa vs normal voice

[kargyraa, Fs]= audioread('Kargyraa_O.wav');
voice = audioread('Voice_O.wav');

% % normalize in respect to RMS
% kargyraa = 0.1*kargyraa/rms(kargyraa);
% normal = 0.1*normal/rms(normal);

fft_k = fft(kargyraa);
fft_n = fft(voice);

fft_k = fft_k(1: length(fft_k)/2);
fft_n = fft_n(1: length(fft_n)/2);

%% in dB

frequency = 0: 0.5*Fs/length(fft_n) : Fs/2-1/Fs;

figure('Name', 'Fourier comparison_dB', 'position', [300 300 900 500])
plot(frequency, mag2db(abs(fft_k)), 'Linewidth', 2.4)
hold all
plot(frequency, mag2db(abs(fft_n)), 'Linewidth', 1.99, 'LineStyle', ':')
legend('Kargyraa','Voice');
set(gca,'fontsize', 16)
title('Kargyraa vs voice signal in db', 'FontSize', 20)
ylabel('Amplitude [dB]')
xlabel('Frequency [Hz]')
xlim([0 1550])
ax = gca;
xtickangle(45)

saveas(gcf, 'Plots/Fourier_comparison_dB', 'epsc')



%% Amplitude not in dB

norm = max(fft_k);
fft_k = fft_k/norm;
fft_n = fft_n/norm;

frequency = 0: 0.5*Fs/length(fft_n) : Fs/2-1/Fs;

figure('Name', 'Fourier comparison', 'position', [300 300 900 500])
plot(frequency, abs(fft_k), 'Linewidth', 2.4)
hold all
plot(frequency, abs(fft_n), 'Linewidth', 1.99, 'LineStyle', ':')
legend('Kargyraa','Voice');
set(gca,'fontsize', 16)
title('Kargyraa vs voice signal', 'FontSize', 20)
ylabel('Amplitude [dB]')
xlabel('Frequency [Hz]')
ylim([-0.001 1.05])
xlim([0 1550])
ax = gca;
ax.XTick = 0:200:1500
xtickangle(45)

saveas(gcf, 'Plots/Fourier_comparison', 'epsc')