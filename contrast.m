%  -*- coding: utf-8 -*-
%
%  Copyright 2024 Technical University of Denmark
%  Authored by: Piotr Budzisz (DTU)
%
%  Licensed under the Apache License, Version 2.0 (the "License");
%  you may not use this file except in compliance with the License.
%  You may obtain a copy of the License at
%
%      http://www.apache.org/licenses/LICENSE-2.0
%
%  Unless required by applicable law or agreed to in writing, software
%  distributed under the License is distributed on an "AS IS" BASIS,
%  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%  See the License for the specific language governing permissions and
%  limitations under the License.
%
%  Note: The following code has been created for the purpose of a master's thesis by Piotr Budzisz titled
%  "Computational reduction of reconstruction artifacts in stitching-based biomedical computed tomography,"
%  written at the Technical University of Denmark. 
%
%  The dataset and files mentioned belong to the CFU (Center for Fast
%  Ultrasound Imaging) at DTU




%% Example line profiles

clear;
clc;
close all;

stitched_image_512 = read_tif('/dtu/cfu/data/userdata18/s220464/storage/PHASE_1/fbp_stitched_3_4_no_clip/stitched_slices_3_4_pad_512.tiff');
stitched_image_8192 = read_tif('/dtu/cfu/data/userdata18/s220464/storage/PHASE_1/fbp_stitched_3_4_no_clip/stitched_slices_3_4_pad_8192.tiff');

line_profile1 = stitched_image_512(:, 512);
line_profile2 = stitched_image_8192(:, 512);


figure; % Create a new figure window
plot(line_profile1, 'b'); %og line profile with clipped values
hold on;
plot(line_profile2, 'r'); %stitched, at the same level as og but without clipping 
hold off;
title('Stitched Image column 1024');
xlabel('Column Index');
ylabel('Intensity Value');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%draw rectangle
figure;
%imshow(imadjust(stitched_image_512));
imagesc(stitched_image_512);
axis on;
hold on;
axis image
rectangle('Position', [150, 100, 50, 30], 'EdgeColor', 'r'); % Draws a rectangle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%combo
figure;
%imshow(imadjust(stitched_image_512));
imagesc(stitched_image_512);
axis on; % Show axis
hold on;
axis image
% Marking a specific pixel
plot(150, 100, 'r+'); % Example pixel
% Highlighting a region
rectangle('Position', [1000, 1500, 500, 300], 'EdgeColor', 'y'); % Example region
% Zoom into a region
xlim([1000, 1500]);
ylim([1500, 1800]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% entropy and std
% Assuming image1 and image2 are your grayscale images
image1 = read_tif('/dtu/cfu/data/userdata18/s220464/storage/PHASE_1/pictures/rec_slice1_512_idx_0000.tiff');
image2 = read_tif('/dtu/cfu/data/userdata18/s220464/storage/PHASE_1/pictures/rec_slice1_8192_idx_0000.tiff');
figure; imshow(imadjust(image1));
figure; imshow(imadjust(image2));
image1 = double(image1);
figure; imshow(imadjust(image1));
image2 = double(image2);
% Calculate entropy
entropy1 = entropy(image1);
entropy2 = entropy(image2);

% Calculate standard deviation
stdDev1 = std2(image1);
stdDev2 = std2(image2);

% Display the results
fprintf('Image 1 - Entropy: %f, Standard Deviation: %f\n', entropy1, stdDev1);
fprintf('Image 2 - Entropy: %f, Standard Deviation: %f\n', entropy2, stdDev2);

% Interpretation
% Higher entropy suggests more information content.
% Higher standard deviation suggests higher contrast.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STD and Entropy measure in ROI
close all;

image1 = read_tif('/dtu/cfu/data/userdata18/s220464/storage/PHASE_1/pictures/rec_slice1_512_idx_0000.tiff');
image2 = read_tif('/dtu/cfu/data/userdata18/s220464/storage/PHASE_1/pictures/rec_slice1_1024_idx_0000.tiff');
image1 = double(image1);
image2 = double(image2);

% Define the ROI coordinates
% x1 = 250;  y1 = 500;  % Top-left corner
% x2 = 1800; y2 = 1500; % Bottom-right corner

x1 = 370; y1 = 1500;
x2 = 1200; y2 = 1800;
% Calculate width and height
width = x2 - x1;
height = y2 - y1;

% Define the position vector for the rectangle [x, y, width, height]
rect = [x1, y1, width, height];

roiImage1 = imcrop(image1, rect);
roiImage2 = imcrop(image2, rect);

% Entropy
entropyRoiImage1 = entropy(roiImage1);
entropyRoiImage2 = entropy(roiImage2);

normRoiImage1 = roiImage1 - mean(roiImage1);
normRoiImage2 = roiImage2 - mean(roiImage2);

entropyNormRoiImage1 = entropy(normRoiImage1);
entropyNormRoiImage2 = entropy(normRoiImage2);

% Standard deviation - ensure the image is in double format for std2
stdRoiImage1 = std2((roiImage1));
stdRoiImage2 = std2((roiImage2));

% Displaying the calculated values
disp(['Entropy of Image 1 ROI: ', num2str(entropyRoiImage1)]);
disp(['Entropy of Image 2 ROI: ', num2str(entropyRoiImage2)]);

disp(['Entropy of Image 1 normalized ROI: ', num2str(entropyNormRoiImage1)]);
disp(['Entropy of Image 2 normalized ROI: ', num2str(entropyNormRoiImage2)]);

disp(['Std Dev of Image 1 ROI: ', num2str(stdRoiImage1)]);
disp(['Std Dev of Image 2 ROI: ', num2str(stdRoiImage2)]);

figure; imshow(imadjust(image1)); hold on;
rectangle('Position', rect, 'EdgeColor', 'r', 'LineWidth', 2); hold off;
title('Image 1 with ROI'); 

figure; imshow(imadjust(image2)); hold on;
rectangle('Position', rect, 'EdgeColor', 'r', 'LineWidth', 2); hold off;
title('Image 2 with ROI');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% std and entropyin ROI from stiched image

close all;

image1 = read_tif('/dtu/cfu/data/userdata18/s220464/storage/PHASE_1/stitched_pictures/stitched_slices_3_4_pad_512_Z2.35.tiff');
image2 = read_tif('/dtu/cfu/data/userdata18/s220464/storage/PHASE_1/stitched_pictures/stitched_slices_3_4_pad_8192_Z2.35.tiff');
%figure; imshow(imadjust(image1));
%figure; imshow(imadjust(image2));
image1 = double(image1);
figure; imshow(imadjust(image1));
image2 = double(image2);
figure; imshow(imadjust(image2));
%limit rectangle
% Define the ROI coordinates
% x1 = 250;  y1 = 500;  % Top-left corner
% x2 = 1800; y2 = 1500; % Bottom-right corner
%datacursormode on;
x1 = 500; y1 = 900;
x2 = 1600; y2 = 1800;
% Calculate width and height
width = x2 - x1;
height = y2 - y1;

% Define the position vector for the rectangle [x, y, width, height]
rect = [x1, y1, width, height];

roiImage1 = imcrop(image1, rect);
roiImage2 = imcrop(image2, rect);

% Entropy
entropyRoiImage1 = entropy(roiImage1);
entropyRoiImage2 = entropy(roiImage2);

normRoiImage1 = roiImage1 - mean2(roiImage1);
normRoiImage2 = roiImage2 - mean2(roiImage2);

entropyNormRoiImage1 = entropy(normRoiImage1);
entropyNormRoiImage2 = entropy(normRoiImage2);

% Standard deviation - ensure the image is in double format for std2
stdRoiImage1 = std2((roiImage1));
stdRoiImage2 = std2((roiImage2));


% Displaying the calculated values
disp(['Entropy of Image 1 ROI: ', num2str(entropyRoiImage1)]);
disp(['Entropy of Image 2 ROI: ', num2str(entropyRoiImage2)]);

disp(['Entropy of Image 1 normalized ROI: ', num2str(entropyNormRoiImage1)]);
disp(['Entropy of Image 2 normalized ROI: ', num2str(entropyNormRoiImage2)]);

disp(['Std Dev of Image 1 ROI: ', num2str(stdRoiImage1)]);
disp(['Std Dev of Image 2 ROI: ', num2str(stdRoiImage2)]);

figure; imshow(imadjust(image1)); hold on;
rectangle('Position', rect, 'EdgeColor', 'r', 'LineWidth', 2); hold off;
title('Image 1 with ROI'); 

figure; imshow(imadjust(image2)); hold on;
rectangle('Position', rect, 'EdgeColor', 'r', 'LineWidth', 2); hold off;
title('Image 2 with ROI');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% vessel line profile in single slice 

clear;
clc;
close all;

image1 = read_tif('/dtu/cfu/data/userdata18/s220464/storage/PHASE_1/pictures/rec_slice2_512_idx_0000.tiff');
image2 = read_tif('/dtu/cfu/data/userdata18/s220464/storage/PHASE_1/pictures/rec_slice2_8192_idx_0000.tiff');
%figure; imshow(imadjust(image1)); 
line_profile1 = image1(1865, 515:583);
line_profile2 = image2(1865, 515:583);

line_profile1 = double(line_profile1);
line_profile2 = double(line_profile2);


figure;
plot(515:583, line_profile1, 'LineWidth', 2);
hold on; % Hold on to plot the second line profile on the same figure
plot(515:583, line_profile2, 'LineWidth', 2);
legend('Image 1', 'Image 2');
xlabel('Column Index');
ylabel('Intensity Value');
title('Line Profile Comparison');
hold off;

croppedImage1 = image1(1835:1895, 515:583);
croppedImage2 = image2(1835:1895, 515:583);

figure;
%imshow(imadjust(image1), []);
subplot(2,2,1)
imagesc(image1, [-0.001 0.002])
%imagesc(image1)
colormap gray;
axis image;
colorbar;
hold on;
% Draw rectangle to mark the ROI
% Note: rectangle position is [X, Y, width, height]
rectangle('Position', [515, 1835, 583-515, 1895-1835], 'EdgeColor', 'r', 'LineWidth', 2);
title('ROI position marked');

% For Image 1
subplot(2,2,2)
imagesc(croppedImage1, [-0.001 0.002]);
%imagesc(croppedImage1);
colormap gray;
colorbar;
hold on;
axis image
line([1, 69], [30, 30], 'Color', 'b', 'LineWidth', 2); % Line from col 515 to 583 at row 30
title('ROI with 512 padding');

%
% For Image 2
subplot(2,2,4)
%imshow(imadjust(croppedImage2), []);
imagesc(croppedImage2, [-0.001 0.002]);
axis image;
colormap gray
colorbar;
hold on;
line([1, 69], [30, 30], 'Color', 'r', 'LineWidth', 2); % Similar to above
title('ROI with 8192 padding');

%scaling
scaled_profile1 = line_profile1 * 1e6;
scaled_profile2 = line_profile2 * 1e6;

% execute function
assess_vessel_contrast(scaled_profile1, scaled_profile2);

% plot
subplot(2,2,3) % Create a new figure window
plot(scaled_profile1, 'b'); % Original line profile with clipped values
hold on;
plot(scaled_profile2, 'r'); % Stitched, at the same level as original but without clipping

% Find and plot peaks for line_profile1 with black triangles
[pks1, locs1] = findpeaks(scaled_profile1);
scatter(locs1, pks1, 'k^'); % Use black triangles to mark peaks on line_profile1

% Find and plot peaks for line_profile2 with black triangles
[pks2, locs2] = findpeaks(scaled_profile2);
scatter(locs2, pks2, 'k^'); % Use black triangles to mark peaks on line_profile2

[troughs1, locs_troughs1] = findpeaks(-scaled_profile1);
troughs1 = -troughs1;
scatter(locs_troughs1, troughs1, 'k^');
[troughs2, locs_troughs2] = findpeaks(-scaled_profile2);
troughs2 = -troughs2;
scatter(locs_troughs2, troughs2, 'k^');

hold off;
title('Vessel line profiles');
xlabel('Column Index');
ylabel('Intensity Value');
legend('512 padding', 'xyz padding', 'Peaks and Troughs');

% Base directory for storing images
savePath = '/dtu/cfu/data/userdata18/s220464/storage/all_figures/';
%% Save the figure
figureFileName = 'vessel_contrast.png';  % Change only the name here
exportgraphics(gcf, fullfile(savePath, figureFileName), 'Resolution', 300);

%% histograms for line profiles of previous section

% Plotting the histogram for both line profiles
figure; % Create a new figure window for histograms

% Histogram for the first line profile
subplot(1, 2, 1); % Divide the figure into a 1x2 grid, and use the first slot
histogram(scaled_profile1, 100,  'FaceColor', 'b');
title('Histogram of Line Profile 1');
xlabel('Intensity Value');
ylabel('Frequency');
xlim([min([scaled_profile1, scaled_profile2]), max([scaled_profile1, scaled_profile2])]); % Align both histograms

% Histogram for the second line profile
subplot(1, 2, 2); % Use the second slot for the second histogram
histogram(scaled_profile2, 100, 'FaceColor', 'r');
title('Histogram of Line Profile 2');
xlabel('Intensity Value');
ylabel('Frequency');
xlim([min([scaled_profile1, scaled_profile2]), max([scaled_profile1, scaled_profile2])]); % Align both histograms

% Adjust the bins for finer granularity or specific analysis needs
% e.g., histogram(scaled_profile1, 40, 'FaceColor', 'b');


%% normalization of profiles

% Normalize both profiles to the range [0, 1]
norm_profile1 = (scaled_profile1 - min(scaled_profile1)) / (max(scaled_profile1) - min(scaled_profile1));
norm_profile2 = (scaled_profile2 - min(scaled_profile2)) / (max(scaled_profile2) - min(scaled_profile2));

% Now proceed with the range calculation and comparison as before
range_norm_profile1 = max(norm_profile1) - min(norm_profile1);
range_norm_profile2 = max(norm_profile2) - min(norm_profile2);

% Display the normalized ranges
disp(['Normalized intensity range for Profile 1: ', num2str(range_norm_profile1)]);
disp(['Normalized intensity range for Profile 2: ', num2str(range_norm_profile2)]);

% Assess which normalized profile has a greater intensity range
if range_norm_profile1 > range_norm_profile2
    disp('Normalized Line Profile 1 has a greater intensity range.');
elseif range_norm_profile2 > range_norm_profile1
    disp('Normalized Line Profile 2 has a greater intensity range.');
else
    disp('Both normalized line profiles have the same intensity range.');
end


%% plot normalized histograms

% Assuming norm_profile1 and norm_profile2 are your normalized data

figure;

% Plot histogram for normalized profile 1
subplot(1, 2, 1); % First subplot in a 1x2 grid
histogram(norm_profile1, 'FaceColor', 'b');
title('Normalized Histogram of Profile 1');
xlabel('Normalized Intensity Value');
ylabel('Frequency');

% Plot histogram for normalized profile 2
subplot(1, 2, 2); % Second subplot in the same grid
histogram(norm_profile2, 'FaceColor', 'r');
title('Normalized Histogram of Profile 2');
xlabel('Normalized Intensity Value');
ylabel('Frequency');

%% plot in one subplot

% Plotting both histograms in the same figure for direct comparison
figure;
histogram(norm_profile1, 0.01, 'FaceColor', 'b', 'FaceAlpha', 0.5, 'EdgeAlpha', 0.5);
hold on; % Keep the figure active to overlay the second histogram
histogram(norm_profile2, 0.01, 'FaceColor', 'r', 'FaceAlpha', 0.5, 'EdgeAlpha', 0.5);
hold off; % Release the figure

title('Overlayed Normalized Histograms of Both Profiles');
xlabel('Normalized Intensity Value');
ylabel('Frequency');
legend('Profile 1', 'Profile 2');

%%  change bin edges

% Generating bin edges from 0 to 0.03 with fine granularity
binEdges = linspace(0, 1, 101); % Creates 100 bins between 0 and 0.03

% Normalized Profile 1 Histogram with custom bin edges
figure;
histogram(norm_profile1, binEdges, 'FaceColor', 'b');
title('Normalized Histogram of Profile 1');
xlabel('Normalized Intensity Value');
ylabel('Frequency');

% Normalized Profile 2 Histogram with custom bin edges
figure; % Use a new figure or subplot as needed
histogram(norm_profile2, binEdges, 'FaceColor', 'r');
title('Normalized Histogram of Profile 2');
xlabel('Normalized Intensity Value');
ylabel('Frequency');

% Plotting both histograms in the same figure for direct comparison
figure;
histogram(norm_profile1, binEdges, 'FaceColor', 'b', 'FaceAlpha', 0.5, 'EdgeAlpha', 0.5);
hold on; % Keep the figure active to overlay the second histogram
histogram(norm_profile2, binEdges, 'FaceColor', 'r', 'FaceAlpha', 0.5, 'EdgeAlpha', 0.5);
hold off; % Release the figure

title('Overlayed Normalized Histograms of Both Profiles');
xlabel('Normalized Intensity Value');
ylabel('Frequency');
legend('Profile 1', 'Profile 2');

%% other section
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lets count the difference
% Step 1: Limit matrices to values from rows 1800 to 2300
limited_matrix1 =scaled_profile1(1800:2300);
limited_matrix2 = scaled_profile2(1800:2300);

% Step 2: Subtract corresponding values
differences = limited_matrix1 - limited_matrix2;

% Step 3: Plot the differences
figure; % Create a new figure window
plot(1800:2300, differences);
title('Differences Between Matrix1 and Matrix2 (Rows 1800 to 2300)');
xlabel('Row Number');
ylabel('Difference (Matrix1 - Matrix2)');

% Step 4: Calculate and display error statistics
mean_difference = mean(differences);
std_difference = std(differences);

% Display the statistics
disp(['Mean of Differences: ', num2str(mean_difference)]);
disp(['Standard Deviation of Differences: ', num2str(std_difference)]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assuming you have matrix1 and matrix2 as before

% Step 1: Calculate the range of values in matrix1 in the specified rows
limited_matrix1_range = [min(scaled_profile1(1800:2300)), max(scaled_profile2(1800:2300))];
disp(['Range of values in matrix1 (Rows 1800 to 2300): [', num2str(limited_matrix1_range(1)), ', ', num2str(limited_matrix1_range(2)), ']']);

% Step 2: Subtract corresponding values (as before)
differences = limited_matrix1 - limited_matrix2;

% Calculate the mean of matrix1 in the specified rows for percentage error calculation
mean_limited_matrix1 = mean(abs(limited_matrix1));

% Calculate percentage error for each difference relative to matrix1
percentage_errors = abs(differences) / mean_limited_matrix1 * 100;

% Step 3: Plot the percentage errors
figure;
plot(1800:2300, percentage_errors);
title('Percentage Errors Between Matrix1 and Matrix2 (Rows 1800 to 2300)');
xlabel('Row Number');
ylabel('Percentage Error (%)');

% Display overall error statistics
mean_percentage_error = mean(percentage_errors);
disp(['Mean Percentage Error: ', num2str(mean_percentage_error), '%']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assuming matrix1 and matrix2 are defined as before
% And you've already limited and subtracted values to get 'differences'

% Calculate the range for plotting
row_range = 1800:2300;

% Plot the two line profiles and their differences
figure; 
hold on;

% Plot line profile 1
plot(row_range, scaled_profile1(row_range), 'b', 'DisplayName', 'Profile 1');

% Plot line profile 2
plot(row_range, scaled_profile2(row_range), 'r', 'DisplayName', 'Profile 2');

% Plot differences as a third line
% Depending on the scale of differences, you might want to adjust the factor for visualization
plot(row_range, differences, 'k', 'DisplayName', 'Differences');

hold off;

title('Comparison of Two Line Profiles and Their Differences');
xlabel('Row Number');
ylabel('Value / Difference');
legend show;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot differences through whole stitched image

% Plot the two line profiles and their differences
figure; 
hold on;

% Plot line profile 1
plot(scaled_profile1, 'b', 'DisplayName', 'Profile 1');

% Plot line profile 2
plot(scaled_profile2, 'r', 'DisplayName', 'Profile 2');

% Plot differences as a third line
% Depending on the scale of differences, you might want to adjust the factor for visualization
minus = scaled_profile1 - scaled_profile2;
plot(minus, 'k', 'DisplayName', 'Differences');

hold off;

title('Comparison of Two Line Profiles and Their Differences');
xlabel('Row Number');
ylabel('Value / Difference');
legend show;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%line_profile1 = cos(linspace(0, 3*pi, 100)) * 20 + 50 + randn(1, 100);
%line_profile2 = cos(linspace(0, 3*pi, 100)) * 20 + 40 + randn(1, 100); % Slightly darker

visualize_contrast(line_profile1, line_profile2);

%% Presentation settings
set(groot,'defaultAxesFontName','times')
set(groot,'defaultAxesFontSize',16)
set(0,'DefaultTextFontSize',12);
set(0,'DefaultLineLineWidth',1.2);


function visualize_contrast(line_profile1, line_profile2)
    % Visualize two line profiles and the contrast between peaks and nearest troughs

    figure;
    
    % Plot first line profile
    subplot(2, 1, 1);
    plot_and_annotate(line_profile1);
    title('Line Profile 1');
    
    % Plot second line profile
    subplot(2, 1, 2);
    plot_and_annotate(line_profile2);
    title('Line Profile 2');
end

function plot_and_annotate(line_profile)
    % Plot line profile and annotate with differences between peaks and nearest troughs
    plot(line_profile, 'LineWidth', 1.5);
    hold on;
    
    % Find peaks
    [pks, locs] = findpeaks(line_profile);
    scatter(locs, pks, 'vr');
    
    % Find troughs by inverting the signal
    [troughs, trough_locs] = findpeaks(-line_profile);
    troughs = -troughs; % Invert back to original values
    scatter(trough_locs, troughs, '^b');
    
%     % Calculate and annotate differences
%     for i = 1:length(locs)
%         % Find the nearest trough for each peak
%         [minDifference, nearestTroughIndex] = min(abs(trough_locs - locs(i)));
%         nearestTroughVal = troughs(nearestTroughIndex);
%         
%         % Calculate difference
%         difference = abs(pks(i) - nearestTroughVal);
%         
%         % Annotate on plot
%         text(locs(i), pks(i) + 0.05 * range(line_profile), sprintf('%.2f', difference), 'HorizontalAlignment', 'center', 'Color', 'black');
%     end
    
    hold off;
    grid on;
    xlabel('Sample');
    ylabel('Intensity');
    legend('Profile', 'Maxima', 'Minima');
end



