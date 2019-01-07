clear
fs = 44100;
f0_min = 60;
f0_max = 200;
f = 347.5;
T = 1;
time = [0:1/fs:T-1/fs];
%signal = sin(2*pi*f*time);
%signal = sawtooth(2*pi*f*time); % T length f hz sawtooth
[y] = audioread('../SantaHoHo.wav');
signal = y;
signal = y(90001: 91024);
time = [0:1/fs:(length(signal)-1)/fs];
N = length(signal);
sig_filter = designfilt('lowpassiir', 'FilterOrder', 4, 'HalfPowerFrequency', f0_max, 'SampleRate', 44100);


%% FUNCTION START
%function [f0] = pitchTrack_zc(signal,N, fs, sig_filter)

signal_filtered = filter(sig_filter,signal);
start_sample = 100;
end_sample = N-1;
pitch_point = [];
% Mark pitch points
for idx = start_sample:end_sample
    if(signal_filtered(idx)<=0 & signal_filtered(idx+1) > 0)
        pitch_point = [pitch_point idx];
    end
end

%Calculate f0 depending on which zero crossing used
for idx = 2:length(pitch_point)
    f0(idx) = 1/((pitch_point(idx) - pitch_point(idx-1))/fs);
end
%f0_rounded = round(f0

%end