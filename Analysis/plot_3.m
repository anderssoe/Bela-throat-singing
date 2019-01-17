clear all
close all
clc
meanammount = 4;
load('means_for_plot_3.mat')

figure('name', 'Results', 'Position', [800 0 800 400])
plot(movmean(kargyraa_mean_gathered, meanammount), 'LineWidth', 2.6)
hold on
plot(movmean(sourcefilter_mean, meanammount), 'LineWidth', 1.8)
legend('Kargyraa spectrum','Source-filter method')
xlim([0 1500])
set(gca,'fontsize', 16)
xlabel('Frequency [Hz]')
ylabel('Amplitude normalized [dB]')
title('Kargyraa vs source-filter method')
xtickangle(45);
ylim([-60 0])

saveas(gcf, 'Analysis/Plots/Kargyraa_Vs_SourceFilter', 'epsc')

meanammount = 4;
load('means_for_plot_3.mat')

figure('name', 'Results', 'Position', [800 0 800 400])
plot(movmean(kargyraa_mean_gathered, meanammount), 'LineWidth', 2.6)
hold on
plot(movmean(shiftedspectrum_mean, meanammount), 'LineWidth', 1.8)
legend('Kargyraa spectrum','Shifted spectrum ')
xlim([0 1500])
set(gca,'fontsize', 16)
xlabel('Frequency [Hz]')
ylabel('Amplitude normalized [dB]')
title('Kargyraa vs shifted spectrum method')
ylim([-60 0])
xtickangle(45);
saveas(gcf, 'Analysis/Plots/Kargyraa_Vs_ShiftedSpectral', 'epsc')