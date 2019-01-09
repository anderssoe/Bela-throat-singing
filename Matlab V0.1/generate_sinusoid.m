% generates a sinusoid with given specifications.
function [time_vector signal] = generate_sinusoid(a, f, phi_k, fs, T_s)
% Function call: 
% >> [time_vector signal] = generate_sinusoid(a, f, phi, fs, T_s) 
% 
% INPUT: 
% a : amplitude 
% f : frequency of sinusoid (in Hz) 
% phi_k : phase (in multiples of 2pi) 
% fs : sampling frequency (in Hz) 
% T_s : duration (in seconds) 
% 
% OUTPUT: 
% time_vector : time vector with sampling points 
% signal : the output signal 
% 
% This function generates a sinusoid at amplitde <a>, frequency <f>, phase <phi> 
% sampled at a sampling frequency <fs> and a duration of <T> seconds. 
% 
% Example: 
% 
% >> [t sig] = generate_sinusoid(0.1, 1000, 0, 44100, 1) 
% 
% generates a sinusoid with amplitude 0.1, frequency of 1000 Hz, phase zero 
% sampled at a sampling frequency of 44100 Hz and a duration of 1 second.

time_vector = 0:1/fs:T_s-1/fs; %time vector
signal = a*sin(2*pi*f*time_vector+phi_k); %signal

end

