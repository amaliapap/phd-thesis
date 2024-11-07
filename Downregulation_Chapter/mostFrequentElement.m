function most_frequent = mostFrequentElement(A)
    % Initialize the output cell array to hold most frequent elements per row
    most_frequent = cell(size(A, 1), 1);
    
    % Loop through each row
    for i = 1:size(A, 1)
        unique_elements = unique(A(i,:)); % Get unique elements in the row
        element_counts = zeros(1, length(unique_elements)); % Array to store counts

        % Count occurrences of each unique element
        for j = 1:length(unique_elements)
            element = unique_elements(j);
            element_counts(j) = sum(A(i,:) == element); % Count occurrences
        end

        % Find the most frequent element(s)
        max_count = max(element_counts); % Get the highest count
        most_frequent_elements = unique_elements(element_counts == max_count); % Elements with max count

        % Store the result (can be a single number or an array in case of a tie)
        most_frequent{i} = most_frequent_elements;
    end
end
