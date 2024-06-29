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

clear;
clc;
close all;

%% upload baseline recons
folderPath = '/dtu/cfu/data/userdata18/s220464/storage/PHASE_1/fbp_stitched_1_7_no_clip';
filePath_1 = fullfile(folderPath, 'stitched_slices_1_7_pad_512.tiff');
filePath_2 = fullfile(folderPath, 'stitched_slices_1_7_pad_8192.tiff');
recon_fbp_512 = imread(filePath_1);
recon_fbp_8192 = imread(filePath_2);

%% select a ROI in baseline image

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

%% show ROI in the image
figure;
%imshow(imadjust(recon_fbp_512));
imagesc(recon_fbp_512, [-0.001 0.002]);
axis on;
hold on;
rectangle('Position', roi_stitched_1_7, 'EdgeColor', 'r', 'LineWidth',2); % Draws a rectangle
colormap(gray(256));
axis image;
title('ROI Stitched 1 7');

datacursormode on

%% create baseline masks

% Convert the image to double precision
imageDouble_1 = double(recon_fbp_512);
imageDouble_2 = double(recon_fbp_8192);
% Crop the image using the provided ROI
img_1 = imcrop(imageDouble_1, roi_7);
img_2 = imcrop(imageDouble_2, roi_7);

figure;
subplot(1, 2, 1)
imagesc(recon_fbp_512, [-0.001 0.002]);
axis on;
colormap(gray(256));
axis image;
hold on;
rectangle('Position', roi_7, 'EdgeColor', 'r', 'LineWidth', 2); % Draws a rectangle
hold off
title('ROI 7, 512 padding');
subplot(1, 2, 2)
imagesc(recon_fbp_8192, [-0.001 0.002]);
axis on;
colormap(gray(256));
axis image;
hold on;
rectangle('Position', roi_7, 'EdgeColor', 'r', 'LineWidth', 2); % Draws a rectangle
title('ROI 7, 8192 padding');

%% compare masks

% Scale the cropped image to [0, 1] for effective Otsu's thresholding
img_scaled_1 = (img_1 - min(img_1(:))) / (max(img_1(:)) - min(img_1(:)));
img_scaled_2 = (img_2 - min(img_2(:))) / (max(img_2(:)) - min(img_2(:)));
% Apply Otsu's method to find an optimal threshold
level_1 = graythresh(img_scaled_1); 
level_2 = graythresh(img_scaled_2);
threshold_1 = level_1 * (max(img_1(:)) - min(img_1(:))) + min(img_1(:));  % Rescale back to the original range
threshold_2 = level_2 * (max(img_2(:)) - min(img_2(:))) + min(img_2(:));  % Rescale back to the original range
% Create the mask based on the threshold
mask_7_512 = img_1 >= threshold_1;
mask_7_8192 = img_2 >= threshold_2;

%% Display the masks
figure;
subplot(1, 3, 1);
imshow(mask_7_512); % Show the first mask
title('Mask 512');

subplot(1, 3, 2);
imshow(mask_7_8192); % Show the second mask
title('Mask 8192');

% Calculate the difference between the two masks
mask_difference = xor(mask_7_512, mask_7_8192);

subplot(1, 3, 3);
imshow(mask_difference); % Show the difference mask
title('Difference');
%% Color the mask
colored_mask_1 = uint8(zeros(size(mask_7_512,1), size(mask_7_512,2), 3)); % Initialize a color mask
colored_mask_1(:,:,1) = mask_7_512 * 255; % Red channel
colored_mask_1(:,:,2) = ~mask_7_512 * 255; % Green channel, yellow = red + green

colored_mask_2 = uint8(zeros(size(mask_7_8192,1), size(mask_7_8192,2), 3)); % Initialize a color mask
colored_mask_2(:,:,1) = mask_7_8192 * 255; % Red channel
colored_mask_2(:,:,2) = ~mask_7_8192 * 255; % Green channel, yellow = red + green

%% Overlay the coloured mask on the cropped image

figure;
subplot(1, 2, 1)
imagesc(img_1, [-0.001 0.002]);
axis on;
colormap(gray(256));
%axis image;
hold on;
% Display the colored mask with transparency
h = imshow(colored_mask_1);
set(h, 'AlphaData', 0.3); % Set transparency to 30%
title('ROI 7, 512 padding');

subplot(1, 2, 2)
imagesc(img_2, [-0.001 0.002]);
%axis on;
colormap(gray(256));
axis image;
hold on;
% Display the colored mask with transparency
h = imshow(colored_mask_2);
set(h, 'AlphaData', 0.3); % Set transparency to 30%
title('ROI 7, 8192 padding');

%% average mask

% Convert masks to double for averaging
mask_7_512_double = double(mask_7_512);
mask_7_8192_double = double(mask_7_8192);

% Calculate the average of the two masks
average_mask = (mask_7_512_double + mask_7_8192_double) / 2;

% Convert the average mask back to binary
average_mask_7 = average_mask >= 0.5;  % Threshold at 0.5


%% Calculate Differences
% Difference between mask_512 and mask_8192
diff_512_8192 = xor(mask_7_512, mask_7_8192);

% Difference between average mask and mask_512
diff_avg_512 = xor(average_mask_7, mask_7_512);

% Difference between average mask and mask_512
diff_avg_8192 = xor(average_mask_7, mask_7_8192);

%% Plotting the Differences
figure;

% Subplot 1: Difference between mask_512 and mask_8192
subplot(1, 3, 1);
imagesc(diff_512_8192);  % Display the difference image
axis image; % Keep aspect ratio
colormap gray; % Use grayscale
title('Difference between Mask 512 and Mask 8192');
colorbar; % Show colorbar to indicate the presence of differences

% Subplot 2: Difference between average mask and mask_512
subplot(1, 3, 2);
imagesc(diff_avg_512);  % Display the difference image
axis image; % Keep aspect ratio
colormap gray; % Use grayscale
title('Difference between Average Mask and Mask 512');
colorbar; % Show colorbar to indicate the presence of differences

% Subplot 2: Difference between average mask and mask_512
subplot(1, 3, 3);
imagesc(diff_avg_8192);  % Display the difference image
axis image; % Keep aspect ratio
colormap gray; % Use grayscale
title('Difference between Average Mask and Mask 8192');
colorbar; % Show colorbar to indicate the presence of differences


%% Display the masks and their differences
figure;
subplot(2, 3, 2);
imshow(mask_7_512); % Show the first mask
title('Mask 512');

subplot(2, 3, 5);
imshow(mask_7_8192); % Show the second mask
title('Mask 8192');

subplot(2, 3, 3);
imshow(average_mask); % Show the averaged binary mask
title('Averaged Binary Mask');

% Calculate the difference between the two initial masks and the average mask
diff_mask_512_avg = xor(mask_7_512, average_mask);
diff_mask_8192_avg = xor(mask_7_8192, average_mask);

subplot(2, 3, 1);
imagesc(img_1, [-0.001 0.002]);
%axis on;
colormap(gray(256));
axis image;
title('512');

subplot(2, 3, 4);
imagesc(img_2, [-0.001 0.002]);
%axis on;
colormap(gray(256));
axis image;
title('8192');

% Optional: show difference between the two initial masks for reference
diff_mask_512_8192 = xor(mask_7_512, mask_7_8192);

subplot(2, 3, 6);
imshow(diff_mask_512_8192); % Show the difference between Mask 512 and 8192
title('Diff 512 & 8192');

%% export to gimp
%adjustedImage = imadjust(img);
%filename = '/dtu/cfu/data/userdata18/s220464/storage/img_roi_2.png';  % Specify the path and filename
%imwrite(adjustedImage, filename);  % Write the image file

%% make corrected mask
level_corr_1 = 0.40;
level_corr_2 = 0.40;
mask_7_512_corr = createMaskFromROI(recon_fbp_512, roi_7, level_corr_1);
mask_7_8192_corr = createMaskFromROI(recon_fbp_8192, roi_7, level_corr_2);

%% Color the corrected mask
colored_mask_corr_1 = uint8(zeros(size(mask_7_512_corr,1), size(mask_7_512_corr,2), 3)); % Initialize a color mask
colored_mask_corr_1(:,:,1) = mask_7_512_corr * 255; % Red channel
colored_mask_corr_1(:,:,2) = ~mask_7_512_corr * 255; % Green channel, yellow = red + green

colored_mask_corr_2 = uint8(zeros(size(mask_7_8192_corr,1), size(mask_7_8192_corr,2), 3)); % Initialize a color mask
colored_mask_corr_2(:,:,1) = mask_7_8192_corr * 255; % Red channel
colored_mask_corr_2(:,:,2) = ~mask_7_8192_corr * 255; % Green channel, yellow = red + green

%% Overlay corrected mask  on the cropped image
figure;
subplot(2, 2, 1)
imagesc(img_1, [-0.001 0.002]);
axis on;
colormap(gray(256));
%axis image;
hold on;
% Display the colored mask with transparency
h = imshow(colored_mask_1);
set(h, 'AlphaData', 0.3); % Set transparency to 30%
title('ROI 7, 512 padding');

subplot(2, 2, 2)
imagesc(img_2, [-0.001 0.002]);
%axis on;
colormap(gray(256));
axis image;
hold on;
% Display the colored mask with transparency
h = imshow(colored_mask_2);
set(h, 'AlphaData', 0.3); % Set transparency to 30%
title('ROI 7, 8192 padding');


subplot(2, 2, 3)
imagesc(img_1, [-0.001 0.002]);
axis on;
colormap(gray(256));
%axis image;
hold on;
% Display the colored mask with transparency
h = imshow(colored_mask_corr_1);
set(h, 'AlphaData', 0.3); % Set transparency to 30%
title('ROI 7, 512 padding, corrected');

subplot(2, 2, 4)
imagesc(img_2, [-0.001 0.002]);
%axis on;
colormap(gray(256));
axis image;
hold on;
% Display the colored mask with transparency
h = imshow(colored_mask_corr_2);
set(h, 'AlphaData', 0.3); % Set transparency to 30%
title('ROI 7, 8192 padding, corrected');


%% average mask corrected

% Convert masks to double for averaging
mask_7_512_double_corr = double(mask_7_512_corr);
mask_7_8192_double_corr = double(mask_7_8192_corr);

% Calculate the average of the two masks
average_mask_corr = (mask_7_512_double_corr + mask_7_8192_double_corr) / 2;

% Convert the average mask back to binary
average_mask_corr_7 = average_mask_corr >= 0.5;  % Threshold at 0.5


%% Display the masks and their differences
figure;
subplot(2, 3, 2);
imshow(mask_7_512_corr); % Show the first mask
title('Mask 512 corrected');

subplot(2, 3, 5);
imshow(mask_7_8192_corr); % Show the second mask
title('Mask 8192 corrected');

subplot(2, 3, 3);
imshow(average_mask_corr); % Show the averaged binary mask
title('Averaged Binary Mask corrected');

% Calculate the difference between the two initial masks and the average mask
diff_mask_512_avg_corr = xor(mask_7_512_corr, average_mask_corr);
diff_mask_8192_avg_corr = xor(mask_7_8192_corr, average_mask_corr);

subplot(2, 3, 1);
imagesc(img_1, [-0.001 0.002]);
%axis on;
colormap(gray(256));
axis image;
title('512');

subplot(2, 3, 4);
imagesc(img_2, [-0.001 0.002]);
%axis on;
colormap(gray(256));
axis image;
title('8192');

% Optional: show difference between the two initial masks for reference
diff_mask_512_8192_corr = xor(mask_7_512, mask_7_8192);

subplot(2, 3, 6);
imshow(diff_mask_512_8192_corr); % Show the difference between Mask 512 and 8192
title('Diff 512 & 8192 corrected');

%% ROI 1 function
% [mask_avg_1, mask_avg_corr_1, level_1_roi_1, level_2_roi_1] =
% createAvgCorrMask(recon_fbp_512, recon_fbp_8192, roi_1, level_corr_1, level_corr_2)
folderPath = '/dtu/cfu/data/userdata18/s220464/storage/PHASE_1/fbp_slice1';
filePath_1 = fullfile(folderPath, 'recon_fbp_512_idx_0000.tiff');
filePath_2 = fullfile(folderPath, 'recon_fbp_8192_idx_0000.tiff');
recon_fbp_512 = imread(filePath_1);
recon_fbp_8192 = imread(filePath_2);
level_corr_1 = [];
level_corr_2 = [];
[mask_avg_1, mask_avg_corr_1, mask_1, mask_2, mask_corr_1, mask_corr_2, level_1_roi_1, level_2_roi_1] = createAvgCorrMask(recon_fbp_512, recon_fbp_8192, roi_1, level_corr_1, level_corr_2);
[mask_avg_1_no_artifacts, mask_avg_corr_1, mask_1, mask_2, mask_corr_1, mask_corr_2, level_1_roi_1, level_2_roi_1] = createAvgCorrMask(recon_fbp_512, recon_fbp_8192, roi_1_no_artifacts, level_corr_1, level_corr_2);

%% ROI 2 function

% [mask_avg_1, mask_avg_corr_1, level_1_roi_1, level_2_roi_1] =
% createAvgCorrMask(recon_fbp_512, recon_fbp_8192, roi_1, level_corr_1, level_corr_2)
folderPath = '/dtu/cfu/data/userdata18/s220464/storage/PHASE_1/fbp_slice1';
filePath_1 = fullfile(folderPath, 'recon_fbp_512_idx_0000.tiff');
filePath_2 = fullfile(folderPath, 'recon_fbp_8192_idx_0000.tiff');
recon_fbp_512 = imread(filePath_1);
recon_fbp_8192 = imread(filePath_2);
level_corr_1 = [];
level_corr_2 = [];
[mask_avg_2, mask_avg_corr_2, mask_1, mask_2, mask_corr_1, mask_corr_2, level_1_roi_1, level_2_roi_1] = createAvgCorrMask(recon_fbp_512, recon_fbp_8192, roi_2, level_corr_1, level_corr_2);

%% ROI stitched function

% [mask_avg_1, mask_avg_corr_1, level_1_roi_1, level_2_roi_1] =
% createAvgCorrMask(recon_fbp_512, recon_fbp_8192, roi_1, level_corr_1, level_corr_2)
folderPath = '/dtu/cfu/data/userdata18/s220464/storage/PHASE_1/fbp_stitched_1_7_no_clip';
filePath_1 = fullfile(folderPath, 'stitched_slices_1_7_pad_512.tiff');
filePath_2 = fullfile(folderPath, 'stitched_slices_1_7_pad_8192.tiff');
stitched_fbp_512 = imread(filePath_1);
stitched_fbp_8192 = imread(filePath_2);
level_corr_1 = [];
level_corr_2 = [];
[mask_avg_1_7, mask_avg_corr_1_7, mask_1, mask_2, mask_corr_1, mask_corr_2, level_1_roi_1, level_2_roi_1] = createAvgCorrMask(stitched_fbp_512, stitched_fbp_8192, roi_stitched_1_7, level_corr_1, level_corr_2);

% Plotting the Differences
figure;

% Subplot 1: Difference between mask_512 and mask_8192
subplot(1, 3, 1);
imagesc(mask_avg_1_7);  % Display the difference image
axis image; % Keep aspect ratio
colormap gray; % Use grayscale
title('mask stitched');
colorbar; % Show colorbar to indicate the presence of differences

% Subplot 2: Difference between average mask and mask_512
subplot(1, 3, 2);
imagesc(mask_1);  % Display the difference image
axis image; % Keep aspect ratio
colormap gray; % Use grayscale
title('Mask 512');
colorbar; % Show colorbar to indicate the presence of differences

% Subplot 2: Difference between average mask and mask_512
subplot(1, 3, 3);
imagesc(mask_2);  % Display the difference image
axis image; % Keep aspect ratio
colormap gray; % Use grayscale
title('Mask 8192');
colorbar; % Show colorbar to indicate the presence of differences
%% Calculate Differences
% Difference between mask_512 and mask_8192
diff_512_8192 = xor(mask_1, mask_2);

% Difference between average mask and mask_512
diff_avg_512 = xor(mask_avg_1_7, mask_1);

% Difference between average mask and mask_512
diff_avg_8192 = xor(mask_avg_1_7, mask_2);

% Plotting the Differences
figure;

% Subplot 1: Difference between mask_512 and mask_8192
subplot(1, 3, 1);
imagesc(diff_512_8192);  % Display the difference image
axis image; % Keep aspect ratio
colormap gray; % Use grayscale
title('Difference between Mask 512 and Mask 8192');
colorbar; % Show colorbar to indicate the presence of differences

% Subplot 2: Difference between average mask and mask_512
subplot(1, 3, 2);
imagesc(diff_avg_512);  % Display the difference image
axis image; % Keep aspect ratio
colormap gray; % Use grayscale
title('Difference between Average Mask and Mask 512');
colorbar; % Show colorbar to indicate the presence of differences

% Subplot 2: Difference between average mask and mask_512
subplot(1, 3, 3);
imagesc(diff_avg_8192);  % Display the difference image
axis image; % Keep aspect ratio
colormap gray; % Use grayscale
title('Difference between Average Mask and Mask 8192');
colorbar; % Show colorbar to indicate the presence of differences

%% saving the corrected mask
% Save the manually corrected mask to a .mat file
%mask_2_corr_040_2048 = mask_manual; % Rename the variable
%save('mask_2_corr_040_2048.mat', 'mask_2_corr_040_2048'); % Save it under the new name

mask_avg_7 = average_mask_7;
mask_avg_corr_7 = average_mask_corr_7; % Rename the variable
save('mask_avg_corr_7.mat', 'mask_avg_corr_7'); % Save it under the new name
save('mask_avg_1.mat', 'mask_avg_1'); % Save it under the new name
save('mask_avg_1_no_artifacts.mat', 'mask_avg_1_no_artifacts'); % Save it under the new name
save('mask_avg_7.mat', 'mask_avg_7'); % Save it under the new name
save('mask_avg_2.mat', 'mask_avg_2'); % Save it under the new name
save('mask_avg_1_7.mat', 'mask_avg_1_7'); % Save it under the new name

%% lets just see if ROI 2 does the same result as in function

clear;
clc;
close all;

%% Presentation settings
set(groot,'defaultAxesFontName','times')
set(groot,'defaultAxesFontSize',16)
set(0,'DefaultTextFontSize',12);
set(0,'DefaultLineLineWidth',1.2);

%% upload baseline recon_fbp_512 from slice 1
%folderPath = '/dtu/cfu/data/userdata18/s220464/storage/PHASE_1/fbp_slice1';
folderPath = '/dtu/cfu/data/userdata18/s220464/storage/PHASE_1/fbp_stitched_1_7_no_clip';
%filePath_1 = fullfile(folderPath, 'recon_fbp_512_idx_0000.tiff');
filePath_1 = fullfile(folderPath, 'stitched_slices_1_7_pad_512.tiff');
%filePath_2 = fullfile(folderPath, 'recon_fbp_8192_idx_0000.tiff');
filePath_2 = fullfile(folderPath, 'stitched_slices_1_7_pad_8192.tiff');
recon_fbp_512 = imread(filePath_1);
recon_fbp_8192 = imread(filePath_2);

%% select a ROI in baseline image

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
roi_1_7 = [x1, y1, width, height];

%% show ROI in the image
figure;
%imshow(imadjust(recon_fbp_512));
imagesc(recon_fbp_512, [-0.001 0.002]);
axis on;
hold on;
rectangle('Position', roi_1, 'EdgeColor', 'r', 'LineWidth',2); % Draws a rectangle
colormap(gray(256));
axis image;
title('ROI 1');

datacursormode on

%% create baseline masks

% Convert the image to double precision
imageDouble_1 = double(recon_fbp_512);
imageDouble_2 = double(recon_fbp_8192);
% Crop the image using the provided ROI
img_1 = imcrop(imageDouble_1, roi_1_7);
img_2 = imcrop(imageDouble_2, roi_1_7);

figure;
subplot(1, 2, 1)
imagesc(recon_fbp_512, [-0.001 0.002]);
axis on;
colormap(gray(256));
axis image;
hold on;
rectangle('Position', roi_1_7, 'EdgeColor', 'r', 'LineWidth', 2); % Draws a rectangle
hold off
title('ROI 1, 512 padding');
subplot(1, 2, 2)
imagesc(recon_fbp_8192, [-0.001 0.002]);
axis on;
colormap(gray(256));
axis image;
hold on;
rectangle('Position', roi_1_7, 'EdgeColor', 'r', 'LineWidth', 2); % Draws a rectangle
title('ROI 1, 8192 padding');

%% compare masks

% Scale the cropped image to [0, 1] for effective Otsu's thresholding
img_scaled_1 = (img_1 - min(img_1(:))) / (max(img_1(:)) - min(img_1(:)));
img_scaled_2 = (img_2 - min(img_2(:))) / (max(img_2(:)) - min(img_2(:)));
% Apply Otsu's method to find an optimal threshold
level_1 = graythresh(img_scaled_1); 
level_2 = graythresh(img_scaled_2);
threshold_1 = level_1 * (max(img_1(:)) - min(img_1(:))) + min(img_1(:));  % Rescale back to the original range
threshold_2 = level_2 * (max(img_2(:)) - min(img_2(:))) + min(img_2(:));  % Rescale back to the original range
% Create the mask based on the threshold
mask_2_512 = img_1 >= threshold_1;
mask_2_8192 = img_2 >= threshold_2;

%% Display the masks
figure;
subplot(1, 3, 1);
imshow(mask_2_512); % Show the first mask
title('ROI 1-2, Mask 512');

subplot(1, 3, 2);
imshow(mask_2_8192); % Show the second mask
title('ROI 1-2, Mask 8192');

% Calculate the difference between the two masks
mask_difference = xor(mask_2_512, mask_2_8192);

subplot(1, 3, 3);
imshow(mask_difference); % Show the difference mask
title('Difference');

% Base directory for storing images
savePath = '/dtu/cfu/data/userdata18/s220464/storage/all_figures/';
%% Save the figure
figureFileName = 'roi1-2maskdifference.png';  % Change only the name here
exportgraphics(gcf, fullfile(savePath, figureFileName), 'Resolution', 300);

%% Color the masks
colored_mask_1 = uint8(zeros(size(mask_2_512,1), size(mask_2_512,2), 3)); % Initialize a color mask
colored_mask_1(:,:,1) = mask_2_512 * 255; % Red channel
colored_mask_1(:,:,2) = ~mask_2_512 * 255; % Green channel, yellow = red + green

colored_mask_2 = uint8(zeros(size(mask_2_8192,1), size(mask_2_8192,2), 3)); % Initialize a color mask
colored_mask_2(:,:,1) = mask_2_8192 * 255; % Red channel
colored_mask_2(:,:,2) = ~mask_2_8192 * 255; % Green channel, yellow = red + green

%% Overlay the coloured mask on the cropped image

figure;
subplot(1, 2, 1)
imagesc(img_1, [-0.001 0.002]);
axis on;
colormap(gray(256));
%axis image;
hold on;
% Display the colored mask with transparency
h = imshow(colored_mask_1);
set(h, 'AlphaData', 0.3); % Set transparency to 30%
title('ROI 1 - 2, 512 padding');

subplot(1, 2, 2)
imagesc(img_2, [-0.001 0.002]);
%axis on;
colormap(gray(256));
axis image;
hold on;
% Display the colored mask with transparency
h = imshow(colored_mask_2);
set(h, 'AlphaData', 0.3); % Set transparency to 30%
title('ROI 1, 8192 padding');

%% average mask

% Convert masks to double for averaging
mask_2_512_double = double(mask_2_512);
mask_2_8192_double = double(mask_2_8192);

% Calculate the average of the two masks
average_mask = (mask_2_512_double + mask_2_8192_double) / 2;

% Convert the average mask back to binary
average_mask_2 = average_mask >= 0.5;  % Threshold at 0.5

%% color the combined mask
colored_mask_com = uint8(zeros(size(average_mask_2,1), size(average_mask_2,2), 3)); % Initialize a color mask
colored_mask_com(:,:,1) = average_mask_2 * 255; % Red channel
colored_mask_com(:,:,2) = ~average_mask_2 * 255; % Green channel, yellow = red + green

%% display combined mask on to p of the image

figure;
imagesc(img_1, [-0.001 0.002]);
axis on;
colormap("gray");
%axis image;
hold on;
% Display the colored mask with transparency
h = imshow(colored_mask_com);
set(h, 'AlphaData', 0.35); % Set transparency to 30%
title('Combined mask, ROI 1 - 2');

% Base directory for storing images
savePath = '/dtu/cfu/data/userdata18/s220464/storage/all_figures/';
%% Save the figure
figureFileName = 'combined_mask_roi_1_2.png';  % Change only the name here
exportgraphics(gcf, fullfile(savePath, figureFileName), 'Resolution', 300);
%% Display the masks and their differences
figure;
subplot(2, 3, 2);
imshow(mask_2_512); % Show the first mask
title('Mask 512');

subplot(2, 3, 5);
imshow(mask_2_8192); % Show the second mask
title('Mask 8192');

subplot(2, 3, 3);
imshow(average_mask); % Show the averaged binary mask
title('Combined Mask');

% Calculate the difference between the two initial masks and the average mask
diff_mask_512_avg = xor(mask_2_512, average_mask);
diff_mask_8192_avg = xor(mask_2_8192, average_mask);

subplot(2, 3, 1);
imagesc(img_1, [-0.001 0.002]);
%axis on;
colormap(gray(256));
axis image;
title('512');

subplot(2, 3, 4);
imagesc(img_2, [-0.001 0.002]);
%axis on;
colormap(gray(256));
axis image;
title('8192');

% Optional: show difference between the two initial masks for reference
diff_mask_512_8192 = xor(mask_2_512, mask_2_8192);

subplot(2, 3, 6);
imshow(diff_mask_512_8192); % Show the difference between Mask 512 and 8192
title('Diff 512 & 8192');

% Base directory for storing images
savePath = '/dtu/cfu/data/userdata18/s220464/storage/all_figures/';
%% Save the figure
figureFileName = 'subfigure_mask_roi_1_2.png';  % Change only the name here
exportgraphics(gcf, fullfile(savePath, figureFileName), 'Resolution', 300);

%% export to gimp
%adjustedImage = imadjust(img);
%filename = '/dtu/cfu/data/userdata18/s220464/storage/img_roi_2.png';  % Specify the path and filename
%imwrite(adjustedImage, filename);  % Write the image file

%% make corrected mask
level_corr_1 = 0.40;
level_corr_2 = 0.40;
mask_2_512_corr = createMaskFromROI(recon_fbp_512, roi_2, level_corr_1);
mask_2_8192_corr = createMaskFromROI(recon_fbp_8192, roi_2, level_corr_2);

%% Color the corrected mask
colored_mask_corr_1 = uint8(zeros(size(mask_2_512_corr,1), size(mask_2_512_corr,2), 3)); % Initialize a color mask
colored_mask_corr_1(:,:,1) = mask_2_512_corr * 255; % Red channel
colored_mask_corr_1(:,:,2) = ~mask_2_512_corr * 255; % Green channel, yellow = red + green

colored_mask_corr_2 = uint8(zeros(size(mask_2_8192_corr,1), size(mask_2_8192_corr,2), 3)); % Initialize a color mask
colored_mask_corr_2(:,:,1) = mask_2_8192_corr * 255; % Red channel
colored_mask_corr_2(:,:,2) = ~mask_2_8192_corr * 255; % Green channel, yellow = red + green

%% Overlay corrected mask  on the cropped image
figure;
subplot(2, 2, 1)
imagesc(img_1, [-0.001 0.002]);
axis on;
colormap(gray(256));
%axis image;
hold on;
% Display the colored mask with transparency
h = imshow(colored_mask_1);
set(h, 'AlphaData', 0.3); % Set transparency to 30%
title('ROI 2, 512 padding');

subplot(2, 2, 2)
imagesc(img_2, [-0.001 0.002]);
%axis on;
colormap(gray(256));
axis image;
hold on;
% Display the colored mask with transparency
h = imshow(colored_mask_2);
set(h, 'AlphaData', 0.3); % Set transparency to 30%
title('ROI 2, 8192 padding');


subplot(2, 2, 3)
imagesc(img_1, [-0.001 0.002]);
axis on;
colormap(gray(256));
%axis image;
hold on;
% Display the colored mask with transparency
h = imshow(colored_mask_corr_1);
set(h, 'AlphaData', 0.3); % Set transparency to 30%
title('ROI 2, 512 padding, corrected');

subplot(2, 2, 4)
imagesc(img_2, [-0.001 0.002]);
%axis on;
colormap(gray(256));
axis image;
hold on;
% Display the colored mask with transparency
h = imshow(colored_mask_corr_2);
set(h, 'AlphaData', 0.3); % Set transparency to 30%
title('ROI 2, 8192 padding, corrected');


%% average mask corrected

% Convert masks to double for averaging
mask_2_512_double_corr = double(mask_2_512_corr);
mask_2_8192_double_corr = double(mask_2_8192_corr);

% Calculate the average of the two masks
average_mask_corr = (mask_2_512_double_corr + mask_2_8192_double_corr) / 2;

% Convert the average mask back to binary
average_mask_corr_2 = average_mask_corr >= 0.5;  % Threshold at 0.5


%% Display the masks and their differences
figure;
subplot(2, 3, 2);
imshow(mask_2_512_corr); % Show the first mask
title('Mask 512 corrected');

subplot(2, 3, 5);
imshow(mask_2_8192_corr); % Show the second mask
title('Mask 8192 corrected');

subplot(2, 3, 3);
imshow(average_mask_corr); % Show the averaged binary mask
title('Averaged Binary Mask corrected');

% Calculate the difference between the two initial masks and the average mask
diff_mask_512_avg_corr = xor(mask_2_512_corr, average_mask_corr);
diff_mask_8192_avg_corr = xor(mask_2_8192_corr, average_mask_corr);

subplot(2, 3, 1);
imagesc(img_1, [-0.001 0.002]);
%axis on;
colormap(gray(256));
axis image;
title('512');

subplot(2, 3, 4);
imagesc(img_2, [-0.001 0.002]);
%axis on;
colormap(gray(256));
axis image;
title('8192');

% Optional: show difference between the two initial masks for reference
diff_mask_512_8192_corr = xor(mask_2_512, mask_2_8192);

subplot(2, 3, 6);
imshow(diff_mask_512_8192_corr); % Show the difference between Mask 512 and 8192
title('Diff 512 & 8192 corrected');
%% save
mask_avg_2 = average_mask_2; % Rename the variable
save('mask_avg_2.mat', 'mask_avg_2'); % Save it under the new name

mask_avg_corr_2 = average_mask_corr_2; % Rename the variable
save('mask_avg_corr_2.mat', 'mask_avg_corr_2'); % Save it under the new name



%% display all the ROI

folder_1 = '/dtu/cfu/data/userdata18/s220464/storage/PHASE_1/fbp_slice1';
folder_7 = '/dtu/cfu/data/userdata18/s220464/storage//PHASE_1/fbp_slice7';
folder_1_7 = '/dtu/cfu/data/userdata18/s220464/storage/PHASE_1/fbp_stitched_1_7_no_clip';
filePath_1 = fullfile(folder_1, 'recon_fbp_512_idx_0000.tiff');
filePath_7 = fullfile(folder_7, 'recon_fbp_512_idx_0000.tiff');
filePath_1_7 = fullfile(folder_1_7, 'stitched_slices_1_7_pad_512.tiff');
recon_fbp_1 = imread(filePath_1);
recon_fbp_7 = imread(filePath_7);
recon_fbp_1_7 = imread(filePath_1_7);

% Create figure
figure;

% Display Image 1 with ROI 1
subplot(1, 3, 1);
imagesc(recon_fbp_1, [-0.001 0.002]);
colormap(gray);
hold on;
rectangle('Position', roi_1, 'EdgeColor', 'r', 'LineWidth', 2);
hold off;
title('ROI 1');
xlabel('Column index');
ylabel('Row index');
axis on;
colorbar; 
axis image;

% Display Image 7 with ROI 7
subplot(1, 3, 2);
imagesc(recon_fbp_7, [-0.001 0.002]);
colormap(gray);
hold on;
rectangle('Position', roi_7, 'EdgeColor', 'r', 'LineWidth', 2);
hold off;
title('ROI 2');
xlabel('Column index');
ylabel('Row index');
colorbar;
axis on;
axis image;

% Display Stitched Image 1_7 with ROI stitched
subplot(1, 3, 3);
imagesc(recon_fbp_1_7, [-0.001 0.002]);
colormap(gray);
hold on;
rectangle('Position', roi_stitched_1_7, 'EdgeColor', 'r', 'LineWidth', 2);
hold off;
title('ROI 1-2');
xlabel('Column index');
ylabel('Row index');
colorbar;
axis on;
axis image;

% Adjust subplot spacing
%set(gcf, 'Position', [100, 100, 1500, 500]);  % Adjust figure size as needed

% Base directory for storing images
savePath = '/dtu/cfu/data/userdata18/s220464/storage/all_figures/';
%% Save the figure
figureFileName = 'ROI_marked.png';  % Change only the name here
exportgraphics(gcf, fullfile(savePath, figureFileName), 'Resolution', 300);

%% Presentation settings
set(groot,'defaultAxesFontName','times')
set(groot,'defaultAxesFontSize',16)
set(0,'DefaultTextFontSize',12);
set(0,'DefaultLineLineWidth',1.2);