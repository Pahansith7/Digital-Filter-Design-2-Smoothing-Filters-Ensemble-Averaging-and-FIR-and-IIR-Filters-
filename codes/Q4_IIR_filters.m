% 4. IIR filters

% 4.1. Realising IIR filters

% Clear the workspace and the command window
clear all;
clc;

load("ECG_with_noise.mat")
fs = 500;
t = (0:length(nECG)-1)/fs;  

% i) Obtain filter coefficients of a Butterworth lowpass filter with the same cut off 
% frequency and of the same order that was used to implement the FIR lowpass filter in 
% the previous section. 

n_low = 38;
Wn_low = 0.32;
[b_low, a_low] = butter(n_low,Wn_low, 'low');

% ii) Visualise the magnitude response, phase response and the group delay 
fvtool(b_low, a_low);

% iii) Obtain the filter coefficients of the highpass and the comb filter and visualise them 
% High Pass filter
n_high = 4;
Wn_high = 0.004;
[b_high, a_high] = butter(n_high,Wn_high, 'high');
fvtool(b_high, a_high);

% Comb Filter
fs = 500; 
fo = 50;  
q = 50; 
bw = (fo/(fs/2))/q;
[b_comb,a_comb] = iircomb(fs/fo, bw, 'notch'); % Note type flag 'notch'
fvtool(b_comb,a_comb);


% Combine all three filters
Low_pass = dfilt.df1(b_low, a_low);
High_pass = dfilt.df1(b_high, a_high);
Comb = dfilt.df1(b_comb, a_comb);
Cascaded_filter = dfilt.cascade(Low_pass,High_pass,Comb);

% Plot the frequency response of the cascaded IIR filter
freqz(Cascaded_filter, 1024)





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4.2. Filtering methods using IIR filters

% i) Apply forward filtering using IIR filters

iir_forward_filtered_ecg = filter(Cascaded_filter, nECG);
figure(1)
plot(t(1:2500), nECG(1:2500), t(1:2500), iir_forward_filtered_ecg(1:2500));
legend('Noisy ECG', 'IIR forward filtered ECG');
title('ECG after forward filtering');
xlabel('Time (s)');
ylabel('Amplitude');

%  ii) Apply forward-backward filtering using IIR filters
filtered_ecg_1 = filtfilt(b_low, a_low, nECG);
filtered_ecg_2 = filtfilt(b_high, a_high, filtered_ecg_1);
iir_forward_backward_filtered_ecg = filtfilt(b_comb, a_comb,filtered_ecg_2 );
figure(2)
plot(t(1:2500), nECG(1:2500), t(1:2500), iir_forward_backward_filtered_ecg(1:2500));
legend('Noisy ECG', 'IIR forward-backward filtered ECG');
title('ECG after forward-backward filtering');
xlabel('Time (s)');
ylabel('Amplitude');

% Generate overlapping time domain plots of the FIR filtered ECG, IIR forward filtered 
% ECG and IIR forward-backward filtered ECG. 
load("FIR_filtered_ECG.mat")
figure(3)
plot(t(1:800), nECG(1:800), t(1:800), fir_filtered_ecg(1:800),...
     t(1:800), iir_forward_filtered_ecg(1:800),...
     t(1:800), iir_forward_backward_filtered_ecg(1:800));
legend('Noisy ECG signal', 'FIR filtered ECG signal', 'IIR forward filtered ECG signal',...
    'IIR forward-backward filtered ECG signal');
title('IIR, FIR fIltered ECG signals comparison');
xlabel('Time (s)');
ylabel('Amplitude');

window = rectwin(length(fir_filtered_ecg)); % to avoid spectral leakage
nfft = length(nECG);
[PSD_n, f_f] = periodogram(nECG, window , nfft, fs);
[PSD_f, ~] = periodogram(fir_filtered_ecg, window , nfft, fs);
[PSD_i_f, ~] = periodogram(iir_forward_filtered_ecg, window , nfft, fs);
[PSD_i_f_b, ~] = periodogram(iir_forward_backward_filtered_ecg, window , nfft, fs);

figure(4)
plot(f_f, 10*log10(PSD_n), f_f, 10*log10(PSD_f), f_f, 10*log10(PSD_i_f), f_f, 10*log10(PSD_i_f_b));
legend('PSD of Noisy ECG signal', 'PSD of FIR filtered ECG signal',...
    'PSD of IIR forward filtered ECG signal', 'PSD of IIR forward-backward filtered ECG signal');
title('IIR, FIR FIltered ECG signal PSD comparison');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
