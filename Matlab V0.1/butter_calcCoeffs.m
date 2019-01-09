cutoff_Hz = 1000;
Fs = 44100;

%function [b,a] = butter_calcCoeffs(cutoff_Hz,Fs)
%Calculates a 2nd order Butterworth IIR filter with a chosen cutoff
%frequency in Hz.
b = NaN(1,3); %3 b values
a = NaN(1,2); %3 a values

cutoff_normalised = (2*pi*cutoff_Hz)/Fs; %normalise cutoff by sample rate and in radians. Nyq = pi
denom = 4+2*sqrt(2)*cutoff_normalised + cutoff_normalised * cutoff_normalised;

b(1) = cutoff_normalised*cutoff_normalised/denom; %b0
b(2) = 2*b(1); %b1
b(3) = b(1); %b2

a(1) = 1;
a(2) = (2*cutoff_normalised*cutoff_normalised - 8) / denom; %a1
a(3) = (cutoff_normalised*cutoff_normalised + 4 - 2*sqrt(2)*cutoff_normalised)/denom;


%[b_matlab, a_matlab] = butter(2,cutoff_Hz/(Fs/2)); %to test if code is similar to matlab