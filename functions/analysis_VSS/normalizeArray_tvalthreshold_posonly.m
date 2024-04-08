function normalizedArray = normalizeArray_tvalthreshold_posonly(yourArray, threshold);
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% Otherwise, perform the division based on the sign of each element
normalizedArray = yourArray; % Initialize the array

normalizedArray(abs(normalizedArray)<threshold)=nan;
normalizedArray(normalizedArray<0)=nan;
end