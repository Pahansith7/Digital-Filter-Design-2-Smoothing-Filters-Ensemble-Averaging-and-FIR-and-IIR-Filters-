function SNR = calculate_SNR(reference_signal, noisy_signal)
    % calculateSNR computes the Signal-to-Noise Ratio (SNR) in dB
    % between a reference (clean) signal and a noisy signal.
    %
    %
    % Inputs:
    %   reference_signal - Clean signal(vector)
    %   noisy_signal     - The noisy signal(vector)
    %
    % Output:
    %   SNR              - The Signal-to-Noise Ratio in decibels (dB)

    % Calculate the noise (difference between noisy and clean signals)
    noise = noisy_signal - reference_signal;

    % Calculate the power of the signal and noise
    signal_power = sum(reference_signal.^2);  % Sum of squared clean signal values
    noise_power = sum(noise.^2);  % Sum of squared noise values
    % Calculate SNR in dB
    SNR = 10 * log10(signal_power / noise_power);
end