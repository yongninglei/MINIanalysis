function normalizedArray = normalizeArray_percentage_posonly(yourArray,threshold)
%
%   Detailed explanation goes here
% Assuming your array is named 'yourArray'

% Otherwise, perform the division based on the sign of each element
normalizedArray = zeros(size(yourArray)); % Initialize the array

% replace all the negative to nan
yourArray(yourArray<0)=nan;

% Loop through each element of yourArray
for i = 1:length(yourArray)

    % If the element is positive, divide by the max value
    normalizedArray(i) = yourArray(i) / maxValue;

end

normalizedArray(normalizedArray<threshold)=nan;

end