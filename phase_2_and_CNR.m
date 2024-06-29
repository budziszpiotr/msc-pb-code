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

%%
clear;
clc;
close all;

%% upload baseline recon_fbp_512 from slice 1
folderPath = '/dtu/cfu/data/userdata18/s220464/storage/PHASE_2/recons';
filePath = fullfile(folderPath, 'recon_fbp_512_idx_0000.tiff');
recon_fbp_512 = imread(filePath);

%% select a ROI in baseline image

% % ROI_1   % large
% x1 = 1000; y1 = 1070;
% x2 = 1280; y2 = 1500;
% width = x2 - x1;
% height = y2 - y1;
% roi_1 = [x1, y1, width, height];

% ROI_1    % no artifacts
x1 = 920; y1 = 1000;
x2 = 1180; y2 = 1260;
width = x2 - x1;
height = y2 - y1;
roi_1_no_artifacts = [x1, y1, width, height];

% ROI_1    
x1 = 870; y1 = 1000;
x2 = 1220; y2 = 1350;
width = x2 - x1;
height = y2 - y1;
roi_1 = [x1, y1, width, height];

% ROI_2
x1 = 1300; y1 = 1526;
x2 = 1600; y2 = 1850;
width = x2 - x1;
height = y2 - y1;
roi_2 = [x1, y1, width, height];

% ROI_3
x1 = 1000; y1 = 1070;
x2 = 1280; y2 = 1500;
width = x2 - x1;
height = y2 - y1;
roi_3 = [x1, y1, width, height];

% ROI_7
x1 = 920; y1 = 1720;
x2 = 1220; y2 = 2020;
width = x2 - x1;
height = y2 - y1;
roi_7 = [x1, y1, width, height];

% ROI_7_large
x1 = 990; y1 = 1450;
x2 = 1300; y2 = 2000;
width = x2 - x1;
height = y2 - y1;
roi_7_large = [x1, y1, width, height];

% ROI_stitched
x1 = 1250; y1 = 1510;
x2 = 1600; y2 = 1860;
width = x2 - x1;
height = y2 - y1;
roi_stitched_1_7 = [x1, y1, width, height];

%%  show ROI in the image
figure;
imshow(imadjust(recon_fbp_512));
axis on;
hold on;
rectangle('Position', roi_7, 'EdgeColor', 'r', 'LineWidth',2); % Draws a rectangle

datacursormode on

%% create baseline masks
mask_1 = createMaskFromROI(recon_fbp_512, roi_1);    
mask_2 = createMaskFromROI(recon_fbp_512, roi_2);

%% load corrected masks
%load('mask_1_corr_048.mat'); % This will load the 'mask_manual' variable into your workspace
%load('mask_1_corr_049.mat');
%load('mask_1_corr_050.mat');
%load('mask_1_corr_051.mat');
%load('mask_1_corr_052.mat');
%load('mask_1_corr_053.mat');
%load('mask_1_corr_055.mat');
%load('mask_2_corr_040.mat');
%load('mask_2_corr_040_2048.mat');
%load('mask_avg_7.mat');
%load('mask_avg_corr_7.mat');

load('mask_avg_1.mat'); 
load('mask_avg_1_no_artifacts.mat'); 
%load('mask_avg_2.mat'); 
load('mask_avg_7.mat');
load('mask_avg_1_7.mat'); 

%% display the masks

% Create a figure for displaying the masks
figure;

% Subplot 1: Display mask_avg_1
subplot(1, 3, 1);
imagesc(mask_avg_1);  % Display the mask as an image
colormap(gray);       % Set colormap to gray for better visualization of the mask
axis image;           % Ensure the axes are scaled equally
title('Mask based on ROI 1');
xlabel('Column Index');
ylabel('Row Index');

% Subplot 2: Display mask_avg_7
subplot(1, 3, 2);
imagesc(mask_avg_7);  % Display the mask as an image
colormap(gray);       % Set colormap to gray
axis image;           % Ensure the axes are scaled equally
title('Mask based on ROI 7');
xlabel('Column Index');
ylabel('Row Index');

% Subplot 3: Display mask_avg_1_7
subplot(1, 3, 3);
imagesc(mask_avg_1_7); % Display the mask as an image
colormap(gray);        % Set colormap to gray
axis image;            % Ensure the axes are scaled equally
title('Mask based on stitched ROI 1 - 7');
xlabel('Column Index');
ylabel('Row Index');

% Enhance layout
%sgtitle('Comparison of Average Masks'); % Super title for the entire figure

% Base directory for storing images
savePath = '/dtu/cfu/data/userdata18/s220464/storage/all_figures/';
%% Save the figure
figureFileName = 'average_masks.png';  % Change only the name here
exportgraphics(gcf, fullfile(savePath, figureFileName), 'Resolution', 300);

%% compare corrected mask with Otsu`s one

figure;
subplot(1, 2, 1);
imshow(mask_2); % Show the image
title('Otsu Threshold');

subplot(1, 2, 2);
imshow(mask_2_corr_040_2048); % Show the image
title('Manual Threshold');

%% upload images

clear files folderPath filePath image fieldName images


% Specify the directory containing the images
folderPath = '/dtu/cfu/data/userdata18/s220464/storage/best_measure/';

% Get a list of all files in the folder with the .tiff extension
files = dir(fullfile(folderPath, '*.tiff'));

% Initialize a structure for storing images
images = struct();

% Read and process each image file
for i = 1:length(files)
    % Full path to file
    filePath = fullfile(folderPath, files(i).name);
    
    % Read the image
    image = imread(filePath);
    
    % Store the image in the structure using the filename as a field name
    % Create a valid field name from the file name
    fieldName = matlab.lang.makeValidName(files(i).name);
    images.(fieldName) = image;

    % Display the image name
    fprintf('Stored image: %s under field %s\n', filePath, fieldName);
end

% Now, images can be accessed in the workspace by their field names
% For example, to access an image stored under the name 'image1.tiff', use:
% image1 = images.image1_tiff;

%% computeAndDisplayCNR

% Get all field names from the structure 'images'
fieldNames = fieldnames(images);

% Loop through each field name
for i = 1:length(fieldNames)
    % Get the field name
    fieldName = fieldNames{i};
    
    % Access the image from the structure using dynamic field names
    img = images.(fieldName);
    
    % Execute the function computeAndDisplayCNR with the current image
    computeAndDisplayCNR(img, roi_stitched_1_7, fieldName, mask_avg_1_7);
end

%% Presentation settings
set(groot,'defaultAxesFontName','times')
set(groot,'defaultAxesFontSize',16)
set(0,'DefaultTextFontSize',12);
set(0,'DefaultLineLineWidth',1.2);

%% line profile test


line_profile1 = images.recon_fbp_512_idx_0000_tiff(:, 1024);
line_profile2 = images.recon_fbp_4096_idx_0000_tiff(:, 1024);


figure; % Create a new figure window
plot(line_profile1, 'b'); %og line profile with clipped values
hold on;
plot(line_profile2, 'r'); %stitched, at the same level as og but without clipping 
hold off;
title('10 vs 100 CGLS iterations');
legend('CGLS 10', 'CGLS 100', 'Location', 'best'); % Add a legend
xlabel('Column Index');
ylabel('Intensity Value');

%% create one more mask for slice 7 more affected by cupping

folderPath = '/dtu/cfu/data/userdata18/s220464/storage/PHASE_1/fbp_slice7';
filePath = fullfile(folderPath, 'recon_fbp_512_idx_0000.tiff');
recon_fbp = imread(filePath);

% ROI_7
x1 = 1000; y1 = 1450;
x2 = 1300; y2 = 2000;
width = x2 - x1;
height = y2 - y1;
roi_7 = [x1, y1, width, height];

% ROI_7_5
x1 = 1000; y1 = 1800;
x2 = 1300; y2 = 2000;
width = x2 - x1;
height = y2 - y1;
roi_7_5 = [x1, y1, width, height];

% show ROI in the image
figure;
imshow(imadjust(recon_fbp));
axis on;
hold on;
rectangle('Position', roi_7_5, 'EdgeColor', 'r', 'LineWidth', 2); % Draws a rectangle

datacursormode on

mask_7_5 = createMaskFromROI(recon_fbp, roi_7_5); 

%% create mask for stitched image 3 4

folderPath = '/dtu/cfu/data/userdata18/s220464/storage/PHASE_1/fbp_stitched';
filePath = fullfile(folderPath, 'stitched_slices_3_4_pad_512.tiff');
stitched_3_4_512 = imread(filePath);

% ROI_stitched
x1 = 1460; y1 = 1820;
x2 = 1670; y2 = 1880;
width = x2 - x1;
height = y2 - y1;
roi_stitched_3_4 = [x1, y1, width, height];

% show ROI in the image
figure;
imshow(imadjust(stitched_3_4_512));
axis on;
hold on;
rectangle('Position', roi_stitched_3_4, 'EdgeColor', 'r', 'LineWidth', 2); % Draws a rectangle

datacursormode on

mask_stitched_3_4 = createMaskFromROI(stitched_3_4_512, roi_stitched_3_4); 

%% create mask for stitched image 7 1

folderPath = '/dtu/cfu/data/userdata18/s220464/storage/PHASE_1/fbp_stitched';
filePath = fullfile(folderPath, 'stitched_slices_1_7_pad_512.tiff');
stitched_1_7_512 = imread(filePath);

% ROI_stitched
x1 = 1200; y1 = 1500;
x2 = 1810; y2 = 1810;
width = x2 - x1;
height = y2 - y1;
roi_stitched_1_7 = [x1, y1, width, height];

% show ROI in the image
figure;
imshow(imadjust(stitched_1_7_512));
axis on;
hold on;
rectangle('Position', roi_stitched_1_7, 'EdgeColor', 'r', 'LineWidth', 2); % Draws a rectangle

datacursormode on

mask_stitched_1_7 = createMaskFromROI(stitched_1_7_512, roi_stitched_1_7); 

%% Entropy and STD
close all; 
%image1 = read_tif('/dtu/cfu/data/userdata18/s220464/storage/PHASE_1/pictures/rec_slice1_512_idx_0000.tiff');
%image2 = read_tif('/dtu/cfu/data/userdata18/s220464/storage/PHASE_1/pictures/rec_slice1_8192_idx_0000.tiff');
image1 = images.recon_cgls_10_idx_0000_tiff;
image2 = images.recon_cgls_100_idx_0000_tiff;
figure; 
imshow(imadjust(image1)); 
title('10 CGLS iterations');
figure; 
imshow(imadjust(image2));
title('100 CGLS iterations');
image1 = double(image1);
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

%% STD and entropy in ROI
close all; 
%image1 = read_tif('/dtu/cfu/data/userdata18/s220464/storage/PHASE_1/pictures/rec_slice1_512_idx_0000.tiff');
%image2 = read_tif('/dtu/cfu/data/userdata18/s220464/storage/PHASE_1/pictures/rec_slice1_8192_idx_0000.tiff');
image1 = images.recon_cgls_100_idx_0000_tiff;
image2 = images.recon_cgls_200_idx_0000_tiff;
image1 = double(image1);
image2 = double(image2);

x1 = 1080; y1 = 1167;
x2 = 1290; y2 = 1500;
% Calculate width and height
width = x2 - x1;
height = y2 - y1;

figure;
imshow(imadjust(image2));
axis on;
hold on;
rectangle('Position', [x1, y1, width, height], 'EdgeColor', 'r'); % Draws a rectangle

% Define the position vector for the rectangle [x, y, width, height]
rect = [x1, y1, width, height];

roiImage1 = imcrop(image1, rect);
roiImage2 = imcrop(image2, rect);
figure; imshow(imadjust(roiImage2));

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

% Calculate mean
meanRoiImage1 = mean(roiImage1(:));
meanRoiImage2 = mean(roiImage2(:));

% Calculate coefficient of variation (CV)
cvRoiImage1 = stdRoiImage1 / meanRoiImage1;
cvRoiImage2 = stdRoiImage2 / meanRoiImage2;

% Displaying the calculated values
disp(['Entropy of Image 1 ROI: ', num2str(entropyRoiImage1)]);
disp(['Entropy of Image 2 ROI: ', num2str(entropyRoiImage2)]);

disp(['Entropy of Image 1 normalized ROI: ', num2str(entropyNormRoiImage1)]);
disp(['Entropy of Image 2 normalized ROI: ', num2str(entropyNormRoiImage2)]);

disp(['Std Dev of Image 1 ROI: ', num2str(stdRoiImage1)]);
disp(['Std Dev of Image 2 ROI: ', num2str(stdRoiImage2)]);

disp(['Coefficient of Variation (CV) for Image 1 ROI: ', num2str(cvRoiImage1 * 100, '%.2f'), '%']);
disp(['Coefficient of Variation (CV) for Image 2 ROI: ', num2str(cvRoiImage2 * 100, '%.2f'), '%']);

figure; imshow(imadjust(image1)); hold on;
rectangle('Position', rect, 'EdgeColor', 'r', 'LineWidth', 2); hold off;
title('Image 1 with ROI'); 

figure; imshow(imadjust(image2)); hold on;
rectangle('Position', rect, 'EdgeColor', 'r', 'LineWidth', 2); hold off;
title('Image 2 with ROI');


%% CNR
close all; 

img = roiImage1;
figure; imshow(imadjust(img));

% scale image to [0,1] so that its fine with otsus thresholding
img_scaled = (img - min(img(:))) / (max(img(:)) - min(img(:)));
%figure; imshow(imadjust(img_scaled));

% Otsu's method
level = graythresh(img_scaled); 
%threshold = level * max(img(:)); % Scale threshold to the range of the image data
threshold = level * (max(img(:)) - min(img(:))) + min(img(:));  % Rescale back to original range
%threshold = 0.00107;
%threshold = mean(double(img(:)));

mask = img >= threshold;
mask_50 = mask; % create copy of the best mask 
%figure;
%imshow(mask);
%title('Binary Mask');

group1 = img(mask); % Bright region (above threshold)
group2 = img(~mask); % Dark region (below threshold)

mean1 = mean(double(group1));
mean2 = mean(double(group2));
std1 = std(double(group1));
std2 = std(double(group2));

cnr = abs(mean1 - mean2) / (std1 + std2);
%cnr = abs(mean1 - mean2) / sqrt(std1^2 + std2^2);
fprintf('CNR = %f\n', cnr);

% Plot histogram with the threshold
%figure;
%histogram(img, 256);
%title('Pixel Intensity Distribution with Otsu`s Threshold');
%hold on;
%xline(threshold, 'r', 'Label', sprintf('Threshold: %.5f', threshold), 'LineWidth', 2);
%hold off;
%fprintf('CNR = %f\n', cnr);



figure;

% Subplot 1: Adjusted Image
subplot(1, 3, 1);
imshow(imadjust(img));
title('Adjusted Image');

% Subplot 2: Histogram with Threshold
subplot(1, 3, 3);
histogram(img, 256);
title('Pixel Intensity Distribution');
hold on;
xline(threshold, 'r', 'Label', sprintf('Threshold: %.5f', threshold), 'LineWidth', 2);
hold off;

% Subplot 3: Binary Mask
subplot(1, 3, 2);
imshow(mask);
title('Binary Mask');


%% Use mask_50 to do some damage
figure; imshow(mask_50);
% so i guess now mask 50 will be used to calculate the cnr, so we skip the
% threshoding point a bit. 

close all; 
%image1 = read_tif('/dtu/cfu/data/userdata18/s220464/storage/PHASE_1/pictures/rec_slice1_512_idx_0000.tiff');
%image2 = read_tif('/dtu/cfu/data/userdata18/s220464/storage/PHASE_1/pictures/rec_slice1_8192_idx_0000.tiff');
image = images.recon_cgls_200_idx_0000_tiff;

image = double(image);

x1 = 1080; y1 = 1167;
x2 = 1290; y2 = 1500;
% Calculate width and height
width = x2 - x1;
height = y2 - y1;

figure;
imshow(imadjust(image));
axis on;
hold on;
rectangle('Position', [x1, y1, width, height], 'EdgeColor', 'r'); % Draws a rectangle

% Define the position vector for the rectangle [x, y, width, height]
rect = [x1, y1, width, height];

roiImage = imcrop(image, rect);
img = roiImage;

group1 = img(mask_50); % Bright region (above threshold)
group2 = img(~mask_50); % Dark region (below threshold)

mean1 = mean(double(group1));
mean2 = mean(double(group2));
std1 = std(double(group1));
std2 = std(double(group2));

cnr_1 = abs(mean1 - mean2) / (std1 + std2);
cnr_2 = abs(mean1 - mean2) / sqrt(std1^2 + std2^2);
fprintf('CNR_1 = %f\n', cnr_1);
fprintf('CNR_2 = %f\n', cnr_2);




%% compare Tikhonov alpha

%% upload images from tikhonov apha folders
%% upload images tikh_05 
% Specify the directory containing the images
folderPath = '/dtu/cfu/data/userdata18/s220464/storage/PHASE_2/cgls_tikh_05';

% Get a list of all files in the folder with the .tiff extension
files = dir(fullfile(folderPath, '*.tiff'));

% Initialize a structure for storing images
cgls_tikh_05 = struct();

% Read and process each image file
for i = 1:length(files)
    % Full path to file
    filePath = fullfile(folderPath, files(i).name);
    
    % Read the image
    image = imread(filePath);
    
    % Store the image in the structure using the filename as a field name
    % Create a valid field name from the file name
    fieldName = matlab.lang.makeValidName(files(i).name);
    cgls_tikh_05.(fieldName) = image;

    % Display the image name
    fprintf('Stored image: %s under field %s\n', filePath, fieldName);
end

% Now, images can be accessed in the workspace by their field names
% For example, to access an image stored under the name 'image1.tiff', use:
% image1 = images.image1_tiff;

%% upload images
% Specify the directory containing the images
folderPath = '/dtu/cfu/data/userdata18/s220464/storage/PHASE_2/cgls_tikh_1';

% Get a list of all files in the folder with the .tiff extension
files = dir(fullfile(folderPath, '*.tiff'));

% Initialize a structure for storing images
cgls_tikh_1 = struct();

% Read and process each image file
for i = 1:length(files)
    % Full path to file
    filePath = fullfile(folderPath, files(i).name);
    
    % Read the image
    image = imread(filePath);
    
    % Store the image in the structure using the filename as a field name
    % Create a valid field name from the file name
    fieldName = matlab.lang.makeValidName(files(i).name);
    cgls_tikh_1.(fieldName) = image;

    % Display the image name
    fprintf('Stored image: %s under field %s\n', filePath, fieldName);
end

% Now, images can be accessed in the workspace by their field names
% For example, to access an image stored under the name 'image1.tiff', use:
% image1 = images.image1_tiff;

%% upload images
% Specify the directory containing the images
folderPath = '/dtu/cfu/data/userdata18/s220464/storage/PHASE_2/cgls_tikh_5';

% Get a list of all files in the folder with the .tiff extension
files = dir(fullfile(folderPath, '*.tiff'));

% Initialize a structure for storing images
cgls_tikh_5 = struct();

% Read and process each image file
for i = 1:length(files)
    % Full path to file
    filePath = fullfile(folderPath, files(i).name);
    
    % Read the image
    image = imread(filePath);
    
    % Store the image in the structure using the filename as a field name
    % Create a valid field name from the file name
    fieldName = matlab.lang.makeValidName(files(i).name);
    cgls_tikh_5.(fieldName) = image;

    % Display the image name
    fprintf('Stored image: %s under field %s\n', filePath, fieldName);
end

% Now, images can be accessed in the workspace by their field names
% For example, to access an image stored under the name 'image1.tiff', use:
% image1 = images.image1_tiff;

%% upload images
% Specify the directory containing the images
folderPath = '/dtu/cfu/data/userdata18/s220464/storage/PHASE_2/cgls_tikh_10';

% Get a list of all files in the folder with the .tiff extension
files = dir(fullfile(folderPath, '*.tiff'));

% Initialize a structure for storing images
cgls_tikh_10 = struct();

% Read and process each image file
for i = 1:length(files)
    % Full path to file
    filePath = fullfile(folderPath, files(i).name);
    
    % Read the image
    image = imread(filePath);
    
    % Store the image in the structure using the filename as a field name
    % Create a valid field name from the file name
    fieldName = matlab.lang.makeValidName(files(i).name);
    cgls_tikh_10.(fieldName) = image;

    % Display the image name
    fprintf('Stored image: %s under field %s\n', filePath, fieldName);
end

% Now, images can be accessed in the workspace by their field names
% For example, to access an image stored under the name 'image1.tiff', use:
% image1 = images.image1_tiff;


%% upload images
% Specify the directory containing the images
folderPath = '/dtu/cfu/data/userdata18/s220464/storage/PHASE_2/cgls_tikh_20';

% Get a list of all files in the folder with the .tiff extension
files = dir(fullfile(folderPath, '*.tiff'));

% Initialize a structure for storing images
cgls_tikh_20 = struct();

% Read and process each image file
for i = 1:length(files)
    % Full path to file
    filePath = fullfile(folderPath, files(i).name);
    
    % Read the image
    image = imread(filePath);
    
    % Store the image in the structure using the filename as a field name
    % Create a valid field name from the file name
    fieldName = matlab.lang.makeValidName(files(i).name);
    cgls_tikh_20.(fieldName) = image;

    % Display the image name
    fprintf('Stored image: %s under field %s\n', filePath, fieldName);
end

% Now, images can be accessed in the workspace by their field names
% For example, to access an image stored under the name 'image1.tiff', use:
% image1 = images.image1_tiff;


%% compare different alpha for tikhonov CGLS 

line_profile0 = images.recon_fbp_idx_0000_tiff(:, 1024);
line_profile1 = cgls_tikh_05.recon_cgls_tikh_200_idx_0000_tiff(:, 1024);
line_profile2 = cgls_tikh_20.recon_cgls_tikh_200_idx_0000_tiff(:, 1024);


figure; % Create a new figure window
plot(line_profile0, 'k'); %og line profile with clipped values
hold on;
plot(line_profile1, 'b'); %og line profile with clipped values
hold on;
plot(line_profile2, 'r'); %stitched, at the same level as og but without clipping 
hold off;
title('CGLS Tikh Aplha 0.5 vs CGLS Tikh Aplha 0.5');
legend('FBP', 'CGLS Tikh 50 iter Alpha 0.5', 'CGLS Tikh 50 iter Alpha 10', 'Location', 'best'); % Add a legend
xlabel('Column Index');
ylabel('Intensity Value');


%% Use CNR to asses difference between Tikh 0.5 and 20 for 100 iterations
close all;
image1 = cgls_tikh_05.recon_cgls_tikh_25_idx_0000_tiff;
image2 = cgls_tikh_05.recon_cgls_tikh_15_idx_0000_tiff;
image1 = double(image1);
image2 = double(image2);

x1 = 1080; y1 = 1167;
x2 = 1290; y2 = 1500;
% Calculate width and height
width = x2 - x1;
height = y2 - y1;

figure;
imshow(imadjust(image2));
axis on;
hold on;
rectangle('Position', [x1, y1, width, height], 'EdgeColor', 'r'); % Draws a rectangle

% Define the position vector for the rectangle [x, y, width, height]
rect = [x1, y1, width, height];

roiImage1 = imcrop(image1, rect);
roiImage2 = imcrop(image2, rect);
figure; imshow(imadjust(roiImage2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% define rect ROI
x1 = 1000; y1 = 1070;
x2 = 1280; y2 = 1500;
% Calculate width and height
width = x2 - x1;
height = y2 - y1;

% Define the position vector for the rectangle [x, y, width, height]
rect = [x1, y1, width, height];

im = read_tif('/dtu/cfu/data/userdata18/s220464/storage/PHASE_2/recons/recon_fbp_idx_0000.tiff');

% show rectangle
figure;
imshow(imadjust(im));
axis on;
hold on;
rectangle('Position', rect, 'EdgeColor', 'r'); % Draws a rectangle

% datacursormode on

%% lets select a two ROI within one slice from fbp to later create mask

%% define rect ROI
x1 = 1010; y1 = 1070;
x2 = 1280; y2 = 1500;
% Calculate width and height
width = x2 - x1;
height = y2 - y1;

% Define the position vector for the rectangle [x, y, width, height]
rect = [x1, y1, width, height];

x1 = 1299; y1 = 1382;
x2 = 1459; y2 = 1920;
% Calculate width and height
width = x2 - x1;
height = y2 - y1;
rect_2 = [x1, y1, width, height];

im = read_tif('/dtu/cfu/data/userdata18/s220464/storage/PHASE_2/recons/recon_fbp_idx_0000.tiff');

% show rectangle
figure;
imshow(imadjust(im));
axis on;
hold on;
rectangle('Position', rect, 'EdgeColor', 'r'); % Draws a rectangle

% datacursormode on

x1 = 763; y1 = 74;
x2 = 1131; y2 = 562;
% Calculate width and height
width = x2 - x1;
height = y2 - y1;
rect_3 = [x1, y1, width, height];

im = read_tif('/dtu/cfu/data/userdata18/s220464/storage/PHASE_2/recons/recon_fbp_idx_0000.tiff');

% show rectangle
figure;
imshow(imadjust(im));
axis on;
hold on;
rectangle('Position', rect_3, 'EdgeColor', 'r'); % Draws a rectangle

%% Read images
close all;
image1 = cgls_tikh_20.recon_cgls_tikh_25_idx_0000_tiff;

%% preprocess ROI

roiImage = makeROI(image1, rect);

%% computeAndDisplayCNR

computeAndDisplayCNR(roiImage); % Call the function