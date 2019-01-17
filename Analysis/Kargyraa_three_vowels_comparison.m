clear all
close all
clc
% The kargyraa sound file of samplerate = 44100 Hz contains the Phonemes e(everything)),
% æ (last), ? (toast). First 10 different e's followed by 10 æ's and 10 ?'s. Each of them
% has the length of one second which corresponds to 44100 samples/vowel.
xlimmax = 1500; 


[Kargyraa Fs] = audioread('Kargyraa_three_vowels_10_times_each.wav');

vowels = (['e','æ','upsilon']);

for ii = 1 : 3
    for jj = 1 : 10
        temp(jj, :) = Kargyraa(jj*Fs + ii*10*Fs - Fs - 10*Fs +1 : jj*Fs + ii*10*Fs - 10*Fs);
        temp(jj, :) = abs(fft(temp(jj, :)));
    end
    mean_fft(ii,:) = mean(temp, 1);
    mean_fft(ii,:) = mean_fft(ii,:)/max(mean_fft(ii,:));
%     figure(ii)
%     %figure('name', vowels(ii) , 'position', [0 0 1200 600])
%     plot(movmean(mean_fft(ii,:),4), 'LineWidth', 2.1)
%     xlim([0 1500])
%     title(vowels(ii))
end

% modalvoice = audioread('Voice_O.wav');
% modalvoice = modalvoice(1:Fs);
% modalvoice = abs(fft(modalvoice));
% modalvoice = modalvoice/max(modalvoice);
% 
% figure('name', 'Comparison of kargyraa vowels', 'Position', [0 0 1100 600])
% plot(movmean(mean_fft(1,:),4), 'LineWidth', 2.6)
% hold on
% plot(movmean(mean_fft(2,:),4), 'LineWidth', 1.8)
% plot(movmean(mean_fft(3,:),4), 'LineWidth', 1.6)
% legend('e','æ','\upsilon')
% xlim([0 1000])
% set(gca,'fontsize', 16)
% xlabel('Frequency [Hz]')
% ylabel('Amplitude Normalized')
% saveas(gcf, 'Analysis/Plots/Kargyraa_vowel_comparison', 'epsc')

figure('name', 'Comparison of kargyraa vowels in dB', 'Position', [0 0 800 400])
plot(movmean(mag2db(mean_fft(1,:)),4), 'LineWidth', 2.6)
hold on
plot(movmean(mag2db(mean_fft(2,:)),4), 'LineWidth', 1.8)
plot(movmean(mag2db(mean_fft(3,:)),4), 'LineWidth', 1.6)
legend('e','æ','\upsilon')
xlim([0 xlimmax])
set(gca,'fontsize', 16)
xlabel('Frequency [Hz]')
title('Comparison of kargyraa vowels')
xtickangle(45);
ylabel('Amplitude normalized [dB]')
saveas(gcf, 'Analysis/Plots/Kargyraa_Vowel_Comparison', 'epsc')


kargyraa_mean_gathered = movmean(mag2db(mean(mean_fft, 1)),4);

[Modal Fs] = audioread('Modal_three_vowels_10_times_each2.wav');


for ii = 1 : 3
    for jj = 1 : 10
        temp(jj, :) = Modal(jj*Fs + ii*10*Fs - Fs - 10*Fs +1 : jj*Fs + ii*10*Fs - 10*Fs);
        temp(jj, :) = abs(fft(temp(jj, :)));
    end
    mean_fft(ii,:) = mean(temp, 1);
    mean_fft(ii,:) = mean_fft(ii,:)/max(mean_fft(ii,:));
%     figure(ii)
    %figure('name', vowels(ii) , 'position', [0 0 1200 600])
%     hold on
%     plot(movmean(mean_fft(ii,:),4), 'LineWidth', 1.7)
%     legend('Kargyraa', 'Modal voice')
%     xlim([0 1500])
%     title(vowels(ii))
end

% figure('name', 'Comparison of modal', 'Position', [0 0 1100 600])
% plot(movmean(mean_fft(1,:),4), 'LineWidth', 2.6)
% hold on
% plot(movmean(mean_fft(2,:),4), 'LineWidth', 1.8)
% plot(movmean(mean_fft(3,:),4), 'LineWidth', 1.6)
% legend('e','æ','\upsilon')
% xlim([0 1000])
% set(gca,'fontsize', 16)
% xlabel('Frequency [Hz]')
% ylabel('Amplitude Normalized')
% saveas(gcf, 'Analysis/Plots/Modal_vowel_comparison', 'epsc')
% 
% figure('name', 'Comparison of modal vowels in dB', 'Position', [1100 0 1100 600])
% plot(movmean(mag2db(mean_fft(1,:)),4), 'LineWidth', 2.6)
% hold on
% plot(movmean(mag2db(mean_fft(2,:)),4), 'LineWidth', 1.8)
% plot(movmean(mag2db(mean_fft(3,:)),4), 'LineWidth', 1.6)
% legend('e','æ','\upsilon')
% xlim([0 1000])
% set(gca,'fontsize', 16)
% xlabel('Frequency [Hz]')
% saveas(gcf, 'Analysis/Plots/Modal_vowel_comparison_dB', 'epsc')
% 
% ylabel('Amplitude [dB]')

modal_mean_gathered = movmean(mag2db(mean(mean_fft, 1)),4);

%% kargyraa gathered vs modal gathered

figure('name', 'Gathered', 'Position', [800 0 800 400])
plot(kargyraa_mean_gathered, 'LineWidth', 2.6)
hold on
plot(modal_mean_gathered, 'LineWidth', 1.8)
legend('Kargyraa spectrum','Modal spectrum')
xlim([0 xlimmax])
set(gca,'fontsize', 16)
xlabel('Frequency [Hz]')
ylabel('Amplitude normalized [dB]')
title('Kargyraa vs modal voice')
xtickangle(45);
saveas(gcf, 'Analysis/Plots/Kargyraa_Vs_Model', 'epsc')