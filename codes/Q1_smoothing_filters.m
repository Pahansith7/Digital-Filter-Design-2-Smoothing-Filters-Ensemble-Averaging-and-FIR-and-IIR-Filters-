% 1.1. Moving Average MA(N) Filter

% Prliminaries

% Clear the workspace and the command window
clear all;
clc;

% i) load ECG_template.mat
load("ECG_template.mat");
ECG = ECG_template;

%signal parameters
fs = 500; % Sampling Frequency (Hz)

% ii) Plot the loaded signal with adjusted time scale
t = (0:length(ECG)-1)/fs;
figure(1)
plot(t,ECG);
title('ECG Template');
xlabel('Time (s)');
ylabel('Amplitude (mV)');

% iii) Add white Gaussian noise to the ECG_template,
% so that the SNR of resultig signal is 5 dB
snr = 5; % in dB
rng(0); % set random feed to zero
nECG= awgn(ECG, snr, 'measured'); % 'measured' argument is used to calculate
% the signal power of the input ECG signal automatically and add a noise
% to get the required snr.

%Plot the ECG Template and Noisy ECG in a same figure
figure(2)
plot(t, ECG, t, nECG);
legend('ECG Template', 'Noisy ECG');
title('Noisy ECG SIgnal');
xlabel('Time (s)');
ylabel('Amplitude (mV)');

% iv) Power spectral density (PSD) estimate of the nECG
nfft = length(nECG); 
window = rectwin(length(nECG));
[ps1, f1] = periodogram(nECG, window, nfft, fs);
figure(3)
plot(f1, 10*log10(ps1)); 
title('Power Spectral Density of Noisy ECG signal');
xlabel('frequency(Hz)');
ylabel('Power Spectral Densityh (dB/Hz)');



%--------------------------------------------------------------------------

% MA(3)  filter implementation with a customised script

% i) Write a MATLAB script for a MA(3) 
N = 3; % Length of the moving window(order)
ma3ECG_1 = MA(nECG, N); 
figure(4)

% Observe the effect of group delay
plot(t*fs, nECG, (1:length(ma3ECG_1)), ma3ECG_1);
legend('Noisy ECG signal', 'Filtered ECG signal');
xlabel('Time(S)');
ylabel('Amplitude (mV)');
title('Noisy ECG signal and filtered ECG signal before compansating the group delay')

% ii) Derive group delay
% Group delay(tau) = (N-1)/2 

% iii) Plot the delay compensated ma3ECG_1and nECG on the same plot
tau = (N-1)/2;
ma3ECG_1 = ma3ECG_1(1+tau:length(ma3ECG_1)-tau); % Compansate the group delay
figure(5)
plot(t, nECG, t, ma3ECG_1);
legend('Noisy ECG signal', 'Filtered ECG signal');
xlabel('Time(S)');
ylabel('Amplitude (mV)');
title('Noisy ECG signal and filtered ECG signal after compansating the group delay')

% Compare PSDs of ma3ECG_1and nECG
figure(6)
[ps2, f2] = periodogram(ma3ECG_1, window, nfft, fs);
plot(f1, 10*log10(ps1), f2, 10*log10(ps2)) 
legend('PSD of noisy ECG signal', 'PSD of filtered ECG signal');
xlabel('frequency(Hz)');
ylabel('Power Spectral Densityh (dB/Hz)');
title('Compare power spectral density of noisy ECG signal and the filtered ECG signal')



%--------------------------------------------------------------------------

% MA(3) filter implementation with the MATLAB built-in function

% i) filter nECG using MATLAB built-in function filter(b,a,X)
N = 3;
b_3 = (1/N)*ones(1,N); % Numerator
a_3 = 1; % Denominator 
ma3ECG_2 = filter(b_3, a_3, nECG);
ma3ECG_2= circshift(ma3ECG_2, -1); % group delay compansate

% ii) plot nECG, ECG_template and ma3ECG_2 on the same plot
figure(7)
plot(t, ECG_template, 'r', t, nECG, 'g', t, ma3ECG_2, 'b');
legend('Template ECG', 'Noisy ECG signal', 'Filtered ECG signal');
xlabel('Time(S)');
ylabel('Amplitude (mV)');
title('Template ECG, Noisy ECG signal and filtered ECG signal after compansating the group delay')


% iii) magnitude response, phase response and the polezero plot
% fvtool(b_3, a_3);



%--------------------------------------------------------------------------

% MA(10) filter implementation with the MATLAB built-in function

N = 10;
b_10 = (1/N)*ones(1,N);
a_10 = 1;

% i) Identify the improvement of the MA(10) filter over the MA(3)
% fvtool(b_10, a_10);

% ii) Filter the nECG signal using the MA(10) filter
ma10ECG = filter(b_10, a_10, nECG);
ma10ECG = circshift(ma10ECG, -5); % group delay compansate

% iv) Plot nECG, ECG_template, ma3ECG_2 and ma10ECG on the same plot to compare 
% the improvement.
figure(8)
plot(t, nECG, t, ECG_template, t, ma3ECG_2, t, ma10ECG); 
xlabel('Time (s)');
ylabel('Amplitude');
legend('nECG', 'ECG-template', 'ma3ECG-2', 'ma10ECG', 'Location', 'Best');
title('Comparison of noisy ECG signal and filtered ECG signals')



%--------------------------------------------------------------------------

% Optimum MA(N) filter order

% i) MSE value is calculted in this script itself.

% ii) test for a rage of N values and determine the optimum filter order
% which gives the minimum MSE

% Define the range of N values to test
N_values = 1:50;  
MSE_values = zeros(size(N_values));  % Initialize an array to store MSE values

for i = 1:length(N_values)
    N = N_values(i);  % Current filter order
    filtered_ECG = filter(ones(1, N) / N, 1, nECG);
    % compansate the group delay
    filtered_ECG = circshift(filtered_ECG, -1*ceil((N-1)/2));
    % Calculate the Mean Squared Error (MSE) between the noise-free signal and the filtered signal
    mse = mean((ECG - filtered_ECG).^2);
    MSE_values(i) = mse;
end

% ii) Find the optimal N that gives the minimum MSE
[~, optimal_index] = min(MSE_values);
optimal_N = N_values(optimal_index);

figure(9);
plot(N_values, MSE_values, '-x');
xlabel('Filter Order N');
ylabel('Mean Squared Error (MSE)');
title('Mean Square Error(MSE) with the order of filter')

fprintf('The optimal order for moving average filter is N = %d ', optimal_N);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1.2. Savitzky-Golay SG(N,L) filter

% Application of SG(N,L)

% i) Apply SG(3,11) filter
N = 20;
L = 11;
L_ = 2*L + 1;
sg310ECG = sgolayfilt(nECG, N, L_);
% sg310ECG = circshift(sg310ECG,-1);

% ii) Plot nECG, ECG_template and sg310ECG on the same plot 
figure(10)
plot(t, ECG_template, 'k', t, nECG, 'b', t, sg310ECG, 'r', 'LineWidth', 1);
legend('Original ECG signal', 'Noisy ECG signal', 'filtered ECG signal');
xlabel('Time (s)');
ylabel('Amplitude (mV)');



%--------------------------------------------------------------------------

% Optimum SG(N,L) filter parameters

% i) Calculate the MSE values for a range of parameters of the SG(N,L) filter and determine 
% the optimum filter parameters which gives the minimum MSE.

% Define the range of parameters
L_values = 1:25;  % Window size range (e.g., from 5 to 50 with step 10)
N_values = 0:24;  % Polynomial order range (e.g., from 1 to 5)

% Initialize the MSE matrix
MSE_matrix = ones(length(L_values), length(N_values)).*inf; %initialize to infinity

% Loop through each combination of N and L
for i = 1:length(L_values)
    L = 2*L_values(i)+1;
    for j = 1:length(N_values)
        N = N_values(j);
        if N > L-1 % N>2L
            continue;
        end
        filtered_signal = sgolayfilt(nECG, N, L);
        % Calculate MSE
        MSE_matrix(i, j) = mean((ECG_template - filtered_signal).^2);
    end
end

% Create a grid for the parameters
[N_grid, L_grid] = meshgrid(N_values, L_values);

figure(11);
surf(N_grid, L_grid, MSE_matrix');
shading interp; 
title('MSE vs. Filter Order (N) and Window Size (L)');
xlabel('(L)');
ylabel('(N))');
colorbar; 

% Find the optimal (N, L) pair
[~, minIndex] = min(MSE_matrix(:));
[optimal_L_index, optimal_N_index] = ind2sub(size(MSE_matrix), minIndex);
optimal_N = N_values(optimal_N_index);
optimal_L = L_values(optimal_L_index);
optimal_MSE = MSE_matrix(optimal_L_index, optimal_N_index);

fprintf('The optimal filter parameters are N = %d and L = %d with a minimum MSE of %.6f\n', ...
        optimal_N, optimal_L, optimal_MSE);

sg_optimum_ECG = sgolayfilt(nECG, optimal_N, optimal_L*2+1); 
    
% ii) Plot ECG_template, sg310ECG and the signal obtained from optimum SG(N,L) filter
figure(12);
plot(t, ECG_template, 'k', t, sg310ECG, 'b', t, sg_optimum_ECG, 'r');
legend('Original ECG signal', 'SG(3,11)', 'SG(Optimum N, Optimum L');
xlabel('Time (s)');
ylabel('Amplitude (mV)');

% iii) Plot ECG_template, signal obtained from optimum MA(N) filter the signal obtained from optimum SG(N,L) filter
N = 12; %optimal filter order of MA filter
b_10 = (1/N)*ones(1,N);
a_10 = 1;
ma_optimum_ECG = filter(b_10, a_10, nECG);
ma_optimum_ECG = circshift(ma_optimum_ECG, -5); % group delay compansate

figure(13);
plot(t, ECG_template, 'r', t, ma_optimum_ECG, 'b', t, sg_optimum_ECG, 'g');
legend('Original ECG signal Template', 'MA(Optimum N)', 'SG(Optimum N, Optimum L');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
