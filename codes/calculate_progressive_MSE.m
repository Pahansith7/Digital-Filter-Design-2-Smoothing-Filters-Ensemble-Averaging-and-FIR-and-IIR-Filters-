function MSEs = calculate_progressive_MSE(template, epochs, M)
    % calculate_progressive_MSE - Calculate progressive MSE between a template 
    % and the mean of increasing epochs.
    %
    % Inputs:
    %   template - Reference(template) signal(vector)
    %   epochs   - A matrix where each column is an epoch signal (matrix)
    %   M        - The maximum number of epochs to use in calculating the
    %              progressive MSE (integer)
    %
    % Output:
    %   MSEs     - A column vector of MSE values calculated progressively up to M.

    MSEs = []; % Initialize an empty array to store MSE values.

    for i = 1:M    
        % Compute the mean of the first i epochs.
        y_i = mean(epochs(:, (1:i)), 2);
        
        % Calculate the Mean Squared Error (MSE) between the template and the averaged epochs.
        mse_value = sqrt(mean((template - y_i).^2));
        
        % Store the computed MSE in the MSEs array.
        MSEs(i, 1) = mse_value;
    end 
end
