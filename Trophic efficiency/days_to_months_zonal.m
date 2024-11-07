% Example 1x73 vector of daily values
% days_vector = rand(1, 73); % Replace with your actual data
function months_vector=days_to_months_zonal(lambdaHTL_zonal)
days_vector = lambdaHTL_zonal;
% Define days per month assuming each month has 30 days
days_in_months = [30, 30, 13]; % 73 days: January = 30, February = 30, part of March = 13

% Initialize the monthly vector
months_vector = zeros(12,64);

% Accumulate data into months
start_idx = 1;
for i = 1:length(days_in_months)
    end_idx = start_idx + days_in_months(i) - 1;
    months_vector(i,:) = sum(days_vector(start_idx:end_idx,:),1); % Use mean() instead of sum() if averaging
    start_idx = end_idx + 1;
end

% Display the resulting months vector
%disp(months_vector);
