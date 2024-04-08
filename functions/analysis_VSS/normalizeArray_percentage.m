function normalizedArray = normalizeArray_percentage(yourArray,threshold)
%
%   Detailed explanation goes here
% Assuming your array is named 'yourArray'
maxValue = max(yourArray);
minValue = min(yourArray);

% Check for the special cases you mentioned
if maxValue == 0 && minValue < 0
    % If the max value is 0 and min value is negative, divide by the min value
    normalizedArray = yourArray / minValue;
elseif minValue == 0 && maxValue > 0
    % If the min value is 0 and max value is positive, divide by the max value
    normalizedArray = yourArray / maxValue;
else
    % Otherwise, perform the division based on the sign of each element
    normalizedArray = zeros(size(yourArray)); % Initialize the array

    % Loop through each element of yourArray
    for i = 1:length(yourArray)
        if yourArray(i) > 0
            % If the element is positive, divide by the max value
            normalizedArray(i) = yourArray(i) / maxValue;
        elseif yourArray(i) < 0
            % If the element is negative, divide by the min value
            normalizedArray(i) = yourArray(i) / minValue;
        % If the element is zero, it remains zero after normalization
        end
    end
end

normalizedArray(abs(normalizedArray)<threshold)=nan;
% Make sure the sign of the result matches the sign of the original element
normalizedArray = normalizedArray .* sign(yourArray);

end