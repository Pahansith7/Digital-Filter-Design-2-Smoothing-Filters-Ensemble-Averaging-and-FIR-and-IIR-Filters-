% Designing FIR filters using windows

% 3.1. Characteristics of window functions (use the fdatool)

% Clear the workspace and the command window
clear all;
clc;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% 3.2. FIR Filter design and application using the Kaiser window

% i) Plot the time domain signal and the power spectral density (PSD) estimate 
%    (or in other words frequency spectrum) of the above signal

load("ECG_with_noise.mat")
fs = 500;

% Plot the time-domain signal
t = (0:length(nECG)-1)/fs;  

figure(1);
plot(t, nECG);
title('Time-domain ECG Signal with Noise');
xlabel('Time (s)');
ylabel('Amplitude');

figure(2);
plot(t(1:1000), nECG(1:1000));
title('Time-domain ECG Signal with Noise');
xlabel('Time (s)');
ylabel('Amplitude');

% Compute and plot the PSD
window = rectwin(length(nECG)); % to avoid spectral leakage
nfft = length(nECG);
[PSD_n, f_n] = periodogram(nECG, window , nfft, fs);

figure(3)
plot(f_n, 10*log10(PSD_n));
title('Power Spectral Density of Noisy ECG');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');

% iii) Calculate the relevant beta and M values for the highpass and lowpass filters.

% Passband and stopband ripple 
delta = 0.01; 

% For highpass filter
f_p_high = 1.5;
f_s_high = 0.5;

% For lowpass filter
f_p_low = 65; 
f_s_low = 95; 

% Calculate the required filter order and beta for Kaiser window
[M_high, WM_high, beta_high, filtype_high] = kaiserord([f_s_high f_p_high], [0 1], [delta delta], fs);
[M_low, Wn_low, beta_low, filtype_low] = kaiserord([f_p_low f_s_low], [1 0], [delta delta], fs);
fprintf('For High pass filter M = %d, Beta = %d ',M_high, beta_high)
fprintf('For Low pass filter M = %d and Beta = %d ',M_low, beta_low)

% Create the filters using fir1 with Kaiser window
b_high = fir1(M_high, WM_high, filtype_high, kaiser(M_high+1, beta_high)); % numerator coeff. (High Pass)
b_low = fir1(M_low, Wn_low, filtype_low, kaiser(M_low+1, beta_low)); % numerator coeff. (Low Pass)
a_high = 1; % denominator coeff. (High Pass)
a_low = 1;% denominator coeff. (Low Pas)

% Visualize the highpass filter
figure(4);
freqz(b_high, a_high, 1024, fs);
title('Highpass Filter Response');


% Visualize the lowpass filter
figure(5);
freqz(b_low, a_low, 1024, fs);
title('Lowpass Filter Response');


% Apply the highpass filter
filtered_ecg_high = filter(b_high, a_high, nECG);

%--------------------------------------------------------------------------
% Not Required

figure(6);
plot(t, filtered_ecg_high);
title('ECG after Highpass Filter Before Compensating the Group Delay');
xlabel('Time (s)');
ylabel('Amplitude');

% Observe the effect of group delay 
figure(7);
plot(t(1:2500), nECG(1:2500), t(1:2500), filtered_ecg_high(1:2500));
legend('Noisy ECG', 'Filtered ECG');
title('ECG after Highpass Filter(Observe the group delay at the begining)');
xlabel('Time (s)');
ylabel('Amplitude');

figure(8);
plot(t(length(nECG)-2500:length(nECG)), nECG(length(nECG)-2500:length(nECG)),...
     t(length(nECG)-2500:length(nECG)), filtered_ecg_high(length(nECG)-2500:length(nECG)));
legend('Noisy ECG', 'Filtered ECG');
title('ECG after Highpass Filter(Observe the group delay at the end)');
xlabel('Time (s)');
ylabel('Amplitude');

%--------------------------------------------------------------------------

% compensate the group delay
filtered_ecg_high = circshift(filtered_ecg_high, -1*ceil(M_high/2));

figure(9);
plot(t, filtered_ecg_high);
title('ECG after Highpass Filter');
xlabel('Time (s)');
ylabel('Amplitude');

figure(10);
plot(t(1:1500), nECG(1:1500), t(1:1500), filtered_ecg_high(1:1500));
legend('Noisy ECG', 'Filtered ECG');
title('ECG after Highpass Filter');
xlabel('Time (s)');
ylabel('Amplitude');


%--------------------------------------------------------------------------
% Not Required

figure(11);
plot(t(1:1000), nECG(1:1000), t(1:1000), filtered_ecg_high(1:1000));
legend('Noisy ECG', 'Filtered ECG');
title('ECG after Highpass Filter(Observe the group delay at the begining-after Compensating)');
xlabel('Time (s)');
ylabel('Amplitude');

figure(12);
plot(t(length(nECG)-1000:length(nECG)), nECG(length(nECG)-1000:length(nECG)),...
     t(length(nECG)-1000:length(nECG)), filtered_ecg_high(length(nECG)-1000:length(nECG)));
legend('Noisy ECG', 'Filtered ECG');
title('ECG after Highpass Filter(Observe the group delay at the end - after Compensating)');
xlabel('Time (s)');
ylabel('Amplitude');

% -------------------------------------------------------------------------

% Apply the lowpass filter
filtered_ecg_low = filter(b_low, a_low, nECG);
filtered_ecg_low = circshift(filtered_ecg_low, -1*ceil((M_low/2)));

% Plot the filtered signal

figure(13)
plot(t, filtered_ecg_low);
title('ECG after Lowpass Filter');
xlabel('Time (s)');
ylabel('Amplitude');

figure(14);
plot(t(1:1500), nECG(1:1500), t(1:1500), filtered_ecg_low(1:1500));
legend('Noisy ECG', 'Filtered ECG');
title('ECG after Lowpass FIlter');
xlabel('Time (s)');
ylabel('Amplitude');


% Design Comb Filter

% Cutoff Frequencies
f1 = 50;
f2 = 100;
f3 = 150;

% Specify zeros in the complex z-plane
zeros_locations = [exp(1i*(2*pi*f1/fs)), exp(-1i*(2*pi*f1/fs)), ...
    exp(1i*(2*pi*f2/fs)), exp(-1i*(2*pi*f2/fs)), ...
    exp(1i*(2*pi*f3/fs)), exp(-1i*(2*pi*f3/fs))];  

% Compute FIR Filter Coefficients
% Convert to (b_0 + b_1*z^-1 +  b_2*z^-2 ...) format
b_comb = poly(zeros_locations);

% Set the gain equals to 1 at oHz.
b_comb = b_comb/sum(b_comb);
a_comb = 1;

figure(15);
freqz(b_comb, a_comb, 2048, fs); % Frequency response of the filter
title('Frequency Response of FIR Filter from Specified Zeros');
fvtool(b_comb, a_comb);

% Apply All thrree filters
% Apply High Pass Filter
filtered_ecg_high = filter(b_high, a_high, nECG);
% Compensate Group Delay
filtered_ecg_high = circshift(filtered_ecg_high, -1*ceil(M_high/2));
% Apply Low Pass Filter
filtered_ecg_high_low = filter(b_low, a_low, filtered_ecg_high); 
% Compensate Group Delay
filtered_ecg_high_low = circshift(filtered_ecg_high_low, -1*ceil(M_low/2));
% Apply the comb Filter
fir_filtered_ecg = filter(b_comb, a_comb, filtered_ecg_high_low);
%Compensate Group Delay
fir_filtered_ecg = circshift(fir_filtered_ecg, -1*(6/2));
save('FIR_filtered_ECG.mat', 'fir_filtered_ecg');

figure(15)
plot(t, fir_filtered_ecg);
title('ECG After Applying All 3 Filters');
xlabel('Time (s)');
ylabel('Amplitude');

figure(16);
plot(t(1:2000), nECG(1:2000), t(1:2000), fir_filtered_ecg(1:2000));
legend('Noisy ECG', 'Filtered ECG');
title('ECG After Applying All 3 Filters');
xlabel('Time (s)');
ylabel('Amplitude');

figure()
[H_high, w] = freqz(b_high, a_high, 1024, fs); % Frequency response of filter 1
[H_low, ~] = freqz(b_low, a_low, 1024, fs); % Frequency response of filter 2
[H_comb, ~] = freqz(b_comb, a_comb, 1024, fs); % Frequency response of filter 3

% Combine the Frequency Responses
H_combined = H_high .* H_low .* H_comb; 
magnitude_combined = 20*log10(abs(H_combined)); 

figure(17);
plot(w, magnitude_combined);
title('Magnitude Response of Combined Filters');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

window = rectwin(length(fir_filtered_ecg)); % to avoid spectral leakage
nfft = length(fir_filtered_ecg);
[PSD_f, f_f] = periodogram(fir_filtered_ecg, window , nfft, fs);

figure(18)
plot(f_n, 10*log10(PSD_n), f_f, 10*log10(PSD_f));
title('Power Spectral Density of Noisy ECG');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');


