clear
fs = 44100;
f0_min = 80;
f0_max = 500;
c = 0.9;
R = 1; % Downsampling factor

f = 440;
T = 1;
time = [0:1/fs:T-1/fs];
signal = sawtooth(2*pi*f*time); % T length f hz sinusoid
N = length(signal);


tau_min = round(fs/f0_max);
tau_max = round(fs/f0_min);


my_sum = 0;
stdDev = 0;

%Calculate mean
for idx = 1:N
    my_sum = my_sum + signal(idx);
end
my_mean = my_sum/N;

for idx = 1:N
    stdDev = stdDev + (signal(idx) - my_mean)^2;
end
stdDev = sqrt(stdDev/N);


stdDevTau = 0;
for idx = 1:tau_max-tau_min +1
    stdDevTau = stdDevTau + (signal(idx) - my_mean)^2;
end
stdDevTau = sqrt(stdDevTau/(tau_max-tau_min));
 

% for k = 1:(N-tau)
%    ACF = ACF +  (signal(k+tau)-my_mean)*(signal(k)-my_mean);
% end
% NACF = 1/stdDev^2 * ACF;
% 
signal_NCCF = [signal zeros(1,tau_max +1)];
NCCF = zeros(1,tau_max - tau_min + 1);
for tau = 1:tau_max - tau_min + 1
    for k = 1:tau_max-tau_min +1
        NCCF(tau) = NCCF(tau) +  ((signal(k + tau_min + tau-1)-my_mean)*(signal(k) - my_mean));% 
    end
    
    for l = tau:tau+tau_max-1
        
    
    NCCF(tau) = NCCF(tau)/sqrt
end
NCCF =  NCCF./(stdDevTau*stdDev);
% 
% NCCF_max = find(max(NCCF));
% f0 = [];
% for idx = 1:length(NCCF)
%     if NCCF(idx) > c*NCCF_max
%         f0 = [f0 (fs/R*(idx+tau_min))];
%     end
% end

sqrt(
