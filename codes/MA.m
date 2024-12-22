function y = MA(x, N)
    % MA applies a moving average filter to the input signal.
    %
    % Inputs:
    %   x - Input signal (vector)
    %   N - Order (integer)
    %
    % Outputs:
    %   y - Filtered signal (vector)

    % Create the moving average filter coefficients
    window = ones(1, N) / N;
    
    % Apply the filter to the input signal using convolution
    y = conv(x, window); % 'same' ensures the output is the same length as the input
end