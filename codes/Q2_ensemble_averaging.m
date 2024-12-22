% 2.1. Signal with multiple measurements

% Preliminaries

% i) ) Clear the workspace and the command window
clear all;
clc;

% ii) Load ABR_rec.mat
load("ABR_rec.mat");
stimuli_train = ABR_rec(:,1);
ABR_signal = ABR_rec(:,2);

% iii) Plot the train of stimuli and ABRs on a single plot
figure(1);
start_sample = 100000;
no_of_samples_to_plot = 10000;
plot(start_sample:(start_sample + no_of_samples_to_plot - 1), ...
    stimuli_train(start_sample:(start_sample + no_of_samples_to_plot - 1)), ...
    start_sample:(start_sample + no_of_samples_to_plot - 1), ...
    ABR_signal(start_sample:(start_sample + no_of_samples_to_plot - 1)));
legend('Stimuli Train', 'ABR Signal');

% iv) Determine a voltage threshold (e.g. 50) to automatically detect the likely stimuli occurrences.
thresh = find(stimuli_train>100);

% v) Extract actual stimulus points
stim_point = [];
j=1; 
for i=1:length(thresh)-1
 if thresh(i+1)-thresh(i)>1
    stim_point(j,1)=thresh(i+1); 
    j=j+1;
 end 
end

% vi) Window ABR epochs according to the extracted stimulus points 
j = 0;
epochs = [];
for i=1:length(stim_point)
    j = j + 1;
    epochs(:,j) = ABR_signal(stim_point(i)-80:stim_point(i)+399); 
end

% vii) Calculate the ensemble average of all the extracted epochs
ensmbl_avg = mean(epochs(:,(1:length(stim_point))),2);

% viii) Plot the ensemble averaged ABR waveform
figure(2),
plot((-80:399)/40,ensmbl_avg)
xlabel('Time (ms)'), ylabel('Voltage(uV)')
title(['Ensemble averaged ABR from ', num2str(length(epochs)), ' epochs']);


%--------------------------------------------------------------------------

% Improvement of the SNR

% i) Calculate progressive MSEs
MSEs = calculate_progressive_MSE(ensmbl_avg, epochs, length(stim_point));

% ii)Plot a graph of MSEk against k 
figure(3);
plot(1:length(stim_point), MSEs, '-x');
xlabel('k');
ylabel('Mean Squared Error (MSE)');
title('Change of MSE with number of epochs used for averaging')





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 2.2. Signal with repetitive patterns

% Viewing the signal and addition of Gaussian white noise

% i) Load ECG_rec.mat to MATLAB workspace
load("ECG_rec.mat");

% ii) Plot the data
figure(4);
fs_ECG = 128;
t_ECG = (0: length(ECG_rec)-1)/fs_ECG;
plot(t_ECG(30:222), ECG_rec(30:222));
xlabel("Time(s)")
ylabel("Amplitude(mV)")
title('A Small Segment from Recorded ECG Signal')

% iii) Extract a single PQRST waveform
figure(5);
plot(t_ECG(143:215), ECG_rec(143:215));
xlabel("Time(s)");
ylabel("Amplitude(mV)");
title('Template ECG Waveform ')
ECG_template = ECG_rec(143:215); % Define template ECG

% iv) Add Gaussian white noise of 5 dB
snr = 5; % in dB
rng(0); % set random feed to zero
nECG = awgn(ECG_rec, snr, 'measured');



%--------------------------------------------------------------------------

% Segmenting ECG into separate epochs and ensemble averaging

% i) Calculate Normalized Cross-correlation using xcorr() function
nECG_normalized = (nECG-mean(nECG))/std(nECG); % normalise nECG
ECG_template_normalized = (ECG_template-mean(ECG_template))/ std(ECG_template); % normalise ECG_template
[cross_corr, lags] = xcorr(nECG_normalized,ECG_template_normalized, 'none');  

% ii)Plot the normalized cross-correlation values against the adjusted lag axis converted to time axis.
adjusted_lag_axis = lags(length(nECG):length(lags)); % adjust the lag
cross_corr = cross_corr(length(nECG):length(cross_corr)); % adjust the cross_corr
time_axis = adjusted_lag_axis/fs_ECG;
figure(6);
plot(time_axis, cross_corr);
xlabel('Time (seconds)');
ylabel('Normalized Cross-Correlation');
title('Normalized Cross-Correlation of ECG Template and Noisy ECG Signal');
grid on;

% iii) Segment ECG pulses by defining a threshold and store in a separate matrix.
cross_corr_thresh = find(cross_corr>39); % defining a threshold

% Identify starting points of the ECG pulses
pulse_starting_points = [];
j=1; 
for i=1:length(cross_corr_thresh)-1
 if cross_corr_thresh(i+1)-cross_corr_thresh(i)>1
    pulse_starting_points(j,1)=cross_corr_thresh(i-1); 
    j=j+1;
 end 
end

% Segment and store ECG pulses
j = 0;
ECG_pulses = [];
for i=1:length(pulse_starting_points)
    j = j + 1;
    ECG_pulses(:,j) = nECG(pulse_starting_points(i,1):pulse_starting_points(i,1)+72); 
    
end

% iv) Calculate and plot the improvement in SNR as the number of ECG pulses included in 
% the ensemble average is increased.
SNRs = [];
for i = 1:length(pulse_starting_points)  
        y_i = mean(ECG_pulses(:,(1:i)),2);
        SNR_i = calculate_SNR(ECG_template', y_i);
        SNRs(i,1) = SNR_i;
end 
figure(7);
plot(1:length(pulse_starting_points), SNRs, '-x');
xlabel('number of ECG segments used for ensemble averaging');
ylabel('SNR (dB)');
title('SNR change as with the number of ECG pulses included in the ensemble averagung')

% v) Plot and compare (in the one graph), a selected noisy ECG pulse and two arbitrary 
% selected ensemble averaged ECG pulses.
figure(8);
ensmbl_avg_ECG_1 = mean(ECG_pulses(:,(1:20)),2);
ensmbl_avg_ECG_2 = mean(ECG_pulses(:,(1:50)),2);
plot((0:72)/fs_ECG, ECG_pulses(:,1), (0:72)/fs_ECG, ensmbl_avg_ECG_1, (0:72)/fs_ECG, ensmbl_avg_ECG_2);
legend('noisy ECG Pulse', 'Averaged signal using 20 ECG pulses', 'Averaged signal using 50 ECG pulses');
xlabel('Time(S)');
ylabel('Amplitude(mV)');
