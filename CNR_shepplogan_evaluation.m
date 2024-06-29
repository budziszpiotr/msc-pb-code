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



% evaluate CNR code on shepp logan phantom
clear;
clc;
close all;

% Base directory for storing images
savePath = '/dtu/cfu/data/userdata18/s220464/storage/all_figures/';

%% Presentation settings
set(groot,'defaultAxesFontName','times')
set(groot,'defaultAxesFontSize',16)
set(0,'DefaultTextFontSize',12);
set(0,'DefaultLineLineWidth',1.2);

%% import phantom

close all;

phantomImage = phantom(512); 
%figure;
%imshow(phantomImage); 

%roi = [215, 383, 85, 55]; % [X, Y, Width, Height] 3 thingies
roi = [330, 200, 40, 40]; % [X, Y, Width, Height] 3 thingies
mask = createMaskFromROI(phantomImage, roi);
figure;
imshow(mask);

figure; 
imshow(phantomImage); 
hold on; 
axis on ;
title('Shepp Logan phantom with marked ROI');
rectangle('Position', roi, 'EdgeColor', 'r', 'LineWidth', 2); 
hold off; 

figure; 
subplot(1,2,1)
imshow(phantomImage); 
hold on; 
axis on ;
title('Shepp Logan phantom with marked ROI');
rectangle('Position', roi, 'EdgeColor', 'r', 'LineWidth', 2); 
hold off; 
subplot(1,2,2)
imshow(mask)
title('Mask based on selected ROI');
axis on ;
%% Save the figure
figureFileName = 'shepp_logan_roi_mask_subfigure.png';  % Change only the name here
exportgraphics(gcf, fullfile(savePath, figureFileName), 'Resolution', 300);

%% Define noise levels with a finer range and increase to 8 levels
noise_levels = linspace(0, 0.1, 8); % Scale standard deviation from 0 to 0.5

%% Initialize a figure for displaying noisy images and their histograms
close all;
figure;

for i = 1:length(noise_levels)
    std_dev = noise_levels(i);
    noise = std_dev * randn(size(phantomImage)); % Generate noise
    noisy_image = phantomImage + noise; % Add noise to the image

    % Crop to ROI
    noisy_roi = imcrop(noisy_image, roi);

    % Display noisy images
    subplot(2, 4, i); % Arrange subplots in 2 rows and 4 columns
    imshow(noisy_roi);
    colorbar;
    %imshow(noisy_roi, []) % uncomment  to scale the image
    colormap gray;
    title(sprintf('Noise STD %.2f', std_dev));
end

%% Save the figure
figureFileName = 'shepp_logan_noise_roi.png';  % Change only the name here
exportgraphics(gcf, fullfile(savePath, figureFileName), 'Resolution', 300);

%% Display histograms in a separate figure for clarity
figure;
for i = 1:length(noise_levels)
    std_dev = noise_levels(i);
    noise = std_dev * randn(size(phantomImage)); % Generate noise again
    noisy_image = phantomImage + noise; % Add noise to the image
    noisy_roi = imcrop(noisy_image, roi); % Crop to ROI

    % Create histograms
    subplot(2, 4, i); % Maintain the same layout for histograms
    histogram(noisy_roi(:), 50); % Create histogram with 50 bins
    title(sprintf('Histogram STD %.2f', std_dev));
    xlabel('Intensity');
    ylabel('Frequency');
end

%% Save the figure
figureFileName = 'shepp_logan_histograms.png';  % Change only the name here
exportgraphics(gcf, fullfile(savePath, figureFileName), 'Resolution', 300);

%% single loop and two figures

% Initialize figure for displaying noisy images
figureNoisy = figure;
title('Noisy Images');

% Initialize figure for displaying histograms
figureHist = figure;
title('Histograms of Noisy Images');

% Generate noisy images and their histograms
for i = 1:length(noise_levels)
    std_dev = noise_levels(i);
    noise = std_dev * randn(size(phantomImage)); % Generate noise
    noisy_image = phantomImage + noise; % Add noise to the image

    % Crop to ROI
    noisy_roi = imcrop(noisy_image, roi);

    % Display noisy images
    figure(figureNoisy); % Switch to figure for noisy images
    subplot(2, 4, i); % Arrange subplots in 2 rows and 4 columns
    imshow(noisy_roi, []);
    title(sprintf('Noise STD %.2f', std_dev));

    % Create histograms
    figure(figureHist); % Switch to figure for histograms
    subplot(2, 4, i); % Arrange subplots in 2 rows and 4 columns
    histogram(noisy_roi(:), 50); % Create histogram with 50 bins
    title(sprintf('Histogram STD %.2f', std_dev));
    xlabel('Intensity');
    ylabel('Frequency');
end



%% Calculate CNR multiple times and find average CNR

num_repeats = 10;  % Number of times to repeat CNR calculation per noise level
cnr_values = zeros(length(noise_levels), num_repeats);  % Store CNR values

for i = 1:length(noise_levels)
    for j = 1:num_repeats
        noise = noise_levels(i) * randn(size(phantomImage)); % Generate noise
        noisy_image = phantomImage + noise; % Add noise to the image
        [cnr, ~, ~] = computeAndDisplayCNR(noisy_image, roi, sprintf('Phantom Image ROI with Noise Level %.2f', noise_levels(i)), mask);
        
        cnr_values(i, j) = cnr;  % Store each CNR measurement
    end
end

% Calculate average CNR for each noise level
cnr_avg = mean(cnr_values, 2);

% Display average CNR values for each noise level
disp('Average CNR values for each noise level:');
for i = 1:length(noise_levels)
    fprintf('Noise Level %.2f: Average CNR = %.4f\n', noise_levels(i), cnr_avg(i));
end

%% other section
% Assuming phantomImage, roi1, and mask1 are already defined

% Initialize arrays to store the outputs
cnr_values = zeros(1, length(noise_levels));
inverse_denominator_values = zeros(1, length(noise_levels));

% Loop through each noise level and compute CNR
for i = 1:length(noise_levels)
    noisy_image = imnoise(phantomImage, 'gaussian', 0, noise_levels(i)); % Add Gaussian noise
    [cnr, num, denom] = computeAndDisplayCNR(noisy_image, roi, sprintf('Phantom Image ROI 1 with Noise Level %.2f', noise_levels(i)), mask1);
    
    % Store the CNR and 1/denominator values
    cnr_values(i) = cnr;
    if denom ~= 0 % Avoid division by zero
        inverse_denominator_values(i) = 1 / denom;
    else
        inverse_denominator_values(i) = NaN; % Handle potential division by zero
    end
end


%% figure
figure; 
plot(inverse_denominator_values, cnr_values, '-o'); 
xlabel('1/Denominator');
ylabel('CNR');
title('CNR vs. 1/Denominator');
grid on; 





%% ALL IN ONE


clear;
clc;
close all;

% Import phantom
phantomImage = phantom(512); 

% Display original phantom image
figure;
imshow(phantomImage);
title('Original Phantom Image');

% Define the ROI
roi = [330, 200, 40, 40]; % [X, Y, Width, Height]
mask = createMaskFromROI(phantomImage, roi);

% Define noise levels with a finer range
noise_levels = linspace(0, 0.2, 8); % Scale standard deviation from 0 to 0.2

% Initialize figure for displaying noisy images and assign a handle
figureNoisy = figure;
set(figureNoisy, 'Name', 'Noisy Images');

% Initialize figure for displaying histograms and assign a handle
figureHist = figure;
set(figureHist, 'Name', 'Histograms of Noisy Images');
% Variables for CNR calculation
num_repeats = 10;  % Number of times to repeat CNR calculation per noise level
cnr_values = zeros(length(noise_levels), num_repeats);  % Store CNR values

% Generate noisy images, their histograms, and compute CNR
for i = 1:length(noise_levels)
    std_dev = noise_levels(i);
    noise = std_dev * randn(size(phantomImage)); % Generate noise
    noisy_image = phantomImage + noise; % Add noise to the image

    % Crop to ROI
    noisy_roi = imcrop(noisy_image, roi);

    % Display noisy images
    figure(figureNoisy); % Switch to figure for noisy images
    subplot(2, 4, i); % Arrange subplots in 2 rows and 4 columns
    imshow(noisy_roi, []);
    title(sprintf('Noise STD %.2f', std_dev));

    % Create histograms
    figure(figureHist); % Switch to figure for histograms
    subplot(2, 4, i); % Arrange subplots in 2 rows and 4 columns
    histogram(noisy_roi(:), 50); % Create histogram with 50 bins
    title(sprintf('Histogram STD %.2f', std_dev));
    xlabel('Intensity');
    ylabel('Frequency');

    % CNR calculations for the current noise level
    for j = 1:num_repeats
        noise = std_dev * randn(size(phantomImage)); % Generate noise again for consistency
        noisy_image = phantomImage + noise; % Add noise to the image
        [cnr, ~, ~] = computeAndDisplayCNR(noisy_image, roi, sprintf('Phantom Image ROI with Noise Level %.2f', noise_levels(i)), mask);
        
        cnr_values(i, j) = cnr;  % Store each CNR measurement
    end
end

% Calculate average CNR for each noise level
cnr_avg = mean(cnr_values, 2);

% Display average CNR values for each noise level
disp('Average CNR values for each noise level:');
for i = 1:length(noise_levels)
    fprintf('Noise Level %.2f: Average CNR = %.4f\n', noise_levels(i), cnr_avg(i));
end



%% PHASE 3 CNR measuremnt
path_1 = '/dtu/cfu/data/userdata18/s220464/storage/PHASE_3/stitched_sh_l/';


im_2 = imread([path_1, 'noisy_sh_l_stitched_220624.tiff']);
im_3 = imread([path_1, 'recon_block_noisy_idx_0000.tiff']);

%% crop
% Get the size of the image
[height, width, ~] = size(im_2);

% Define the size of the crop
cropSize = 512;

% Calculate the coordinates of the top-left corner of the crop
xCenter = floor(width / 2);
yCenter = floor(height / 2);
xStart = xCenter - floor(cropSize / 2);
yStart = yCenter - floor(cropSize / 2);

% Crop the image
im_2_cropped = imcrop(im_2, [xStart yStart cropSize-1 cropSize-1]);

figure;
subplot(1, 2, 1);  
imagesc(im_2_cropped, [0 1]);
%colormap("gray");
axis image;
colorbar;
title('Stitching');
%xlabel('X-axis label');
%ylabel('Y-axis label');
subplot(1, 2, 2);  
imagesc(im_3, [0 1]);
colormap("gray");
axis image;
colorbar;
title('Blocking');
%xlabel('X-axis label');
%ylabel('Y-axis label');


% Base directory for storing images
savePath = '/dtu/cfu/data/userdata18/s220464/storage/all_figures/';
%% Save the figure
figureFileName = 'stitchvsblock.png';
exportgraphics(gcf, fullfile(savePath, figureFileName), 'Resolution', 300);

%% CNR
%[cnr, ~, ~] = computeAndDisplayCNR(im_1, roi, sprintf('Phantom Image ROI'), mask);
[cnr, ~, ~] = computeAndDisplayCNR(im_2_cropped, roi, sprintf('Phantom Image ROI'), mask);
[cnr, ~, ~] = computeAndDisplayCNR(im_3, roi, sprintf('Phantom Image ROI'), mask);

