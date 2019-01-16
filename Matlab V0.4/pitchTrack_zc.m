
%% TEST INIT
% clear
% fs = 44100;
% f0_min = 0;
% f0_max = 180; %-3dB LP frequency
% f = 90; 
% T = 1024/44100;
% time = [0:1/fs:T-1/fs];
% signal = sawtooth(2*pi*f*time); % T length f hz sawtooth
% [y] = audioread('../Hoho_NoAutotune.wav');
% time = [0:1/fs:(length(signal)-1)/fs]; %Create time vector
% N = length(signal);
% %4th order IIR LP filter
% sig_filter = designfilt('lowpassiir', 'FilterOrder', 4, 'HalfPowerFrequency', f0_max, 'SampleRate', 44100);
% 

%% FUNCTION START
function [f0] = pitchTrack_zc(signal,N, sig_filter)

global Fs;

signal_filtered = filter(sig_filter,signal);
start_sample = 1;
end_sample = N-1;
pitch_point = [];
pitch_time = zeros(1,N);

% Mark pitch points
for idx = start_sample:end_sample
    if(signal_filtered(idx)<=0 && signal_filtered(idx+1) > 0)
        pitch_point = [pitch_point idx];
        pitch_time(idx) = 1;
    end
end

%Calculate possible f0 values
f0 = -1; % Default value
%calculate distance between pitch points and convert to frequency
for idx = 2:length(pitch_point)
    f0(idx -1) = 1/((pitch_point(idx) - pitch_point(idx-1))/Fs);
end

%return lowest frequency if multiple possible f0 are detected
if(length(f0) >1)
    f0 = min(f0);
end



end
