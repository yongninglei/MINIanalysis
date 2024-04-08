function myMultimodalPCA()
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% first test it for 1 vertex
% convert the struct into a matrix

% normalize the matrix
normalize=true;
[LEX_matrix, col_name_LEX]=prepareMatrix(vertex_2489_Struct_LEX ,normalize);

normalize=true;
[PER_matrix, col_name_PER]=prepareMatrix(vertex_2489_Struct_PER ,normalize);

%Do the PCA
[coeff, score, latent, tsquared, explained]= pca(LEX_matrix);

[coeff_per, score_per, latent_per, tsquared_per, explained_per]= pca(PER_matrix);





% so I will normalized it and see the change
% draw the distribution figure to see the different measurement's
% distributiuon

% Create a figure for the scatter plots
figure;

% Number of measurements
numMeasurements = size(LEX_matrix, 2);

% Plot each pair of measurements
for i = 1:numMeasurements
    for j = i+1:numMeasurements
        % Calculate subplot position
        subplotIndex = (i-1) * (numMeasurements-1) + j-1;
        
        % Create subplot
        subplot(numMeasurements-1, numMeasurements-1, subplotIndex);
        
        % Scatter plot of measurement i vs measurement j
        scatter(LEX_matrix(:, i), LEX_matrix(:, j));
        
        % Labeling
        xlabel(sprintf('%s', col_name_LEX(i)));
        ylabel(sprintf('%s', col_name_LEX(j)));
        title(sprintf('Measure %s vs %s', col_name_LEX(i), col_name_LEX(j)));
    end
end

% Adjust layout to prevent label overlap
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 1, 1]);

% draw everything to 1 figure
colors = ['r', 'g', 'b', 'c', 'm', 'y']; % One color for each measurement
%% plot all measurement to 1 
figure;
hold on; % Keep the plot active to overlay multiple scatter plots

% Assuming x-axis is simply the observation index
x = 1:size(LEX_matrix, 1); % Or create a meaningful x-axis if you have one

% Loop through each measurement
for i = 1:size(LEX_matrix, 2)
    scatter(x, LEX_matrix(:, i), colors(i), 'DisplayName', sprintf('Measure %d', i));
end

hold off;
legend('show'); % Show legend to identify each measurement
xlabel('Observation index'); % Adjust as per your x-axis description
ylabel('Measurement value');
title('Scatter Plot of MRI Measurements');


end