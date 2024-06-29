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

%% Display sorted recons and ROI
clear;
clc;
close all;

%% Presentation settings
set(groot,'defaultAxesFontName','times')
set(groot,'defaultAxesFontSize',16)
set(0,'DefaultTextFontSize',12);
set(0,'DefaultLineLineWidth',1.2);

%% Display sorted recons
% Define the directory containing the images
folderPath = '/dtu/cfu/data/userdata18/s220464/storage/PHASE_2/cgls_stitched_1_7_no_clip';

% Get a list of all .tiff files in the folder
files = dir(fullfile(folderPath, '*.tiff')); % Adjust the extension if needed

% Check if files were found
if isempty(files)
    error('No files found. Check the directory and file type.');
end

% Names of files to exclude
%excludedFiles = {'recon_fista_tv_0_idx_0000.tiff', 'recon_fista_tv_10_idx_0000.tiff'};

% Filter out excluded files
%files = files(~ismember({files.name}, excludedFiles));

% Extract numbers from filenames to sort by
numbers = zeros(length(files), 1);
for i = 1:length(files)
    filename = files(i).name;
    numStr = regexp(filename, '(?<=recon_fbp_)\d+', 'match', 'once');  
    %numStr = regexp(filename, '(?<=data_pad_)\d+', 'match', 'once');  
    %numStr = regexp(filename, '(?<=stitched_slices_1_7_pad_)\d+', 'match', 'once'); 
    %numStr = regexp(filename, '(?<=recon_fista_tv_)\d*\.?\d*(?=_idx)', 'match', 'once');
    %numStr = regexp(filename, '(?<=recon_cgls_)\d*\.?\d*(?=_idx)', 'match', 'once');
    %numStr = regexp(filename, '(?<=recon_sirt_)\d*\.?\d*(?=_idx)', 'match', 'once');
    
    numbers(i) = str2double(numStr);
end

% Sort files based on the extracted numbers
[~, idx] = sort(numbers);
files = files(idx);

% Initialize a figure for displaying images
figure;
colormap gray; % Set the colormap to gray for all subplots

% Super title for the entire figure
%sgtitle('FBP Stitched Slice 1 and 7');
%sgtitle('FISTA TV Slice 1 500 Iterations');
figureFileName = 'empty_padding.png';  % Change only the name here
% Number of rows and columns in subplot grid
numRows = 2;
numCols = 3;

% Prepare to adjust subplot sizes
left = 0.04; % Left margin
bottom = 0.04; % Bottom margin
inter_row_spacing = 0.1; % Space between rows
% Adjust width and height
width = (1 - left * 2) / numCols; % Width of each subplot
height = (1 - 2 * bottom - (numRows - 1) * inter_row_spacing) / numRows; % Adjust height to account for space between rows

% Loop through each file and load/display the image
for i = 1:length(files)
    % Full path to file
    filePath = fullfile(folderPath, files(i).name);
    
    % Read the image
    img = imread(filePath);
    
    % Calculate position of the subplot
    col = mod(i-1, numCols); % 0-indexed column
    row = floor((i-1) / numCols); % 0-indexed row

    % Calculate exact position [left, bottom, width, height]
    % Adjust 'bottom' to move each row up or down
    bottom_pos = 1 - bottom - (row * (height + inter_row_spacing));
    
    position = [left + col * width, bottom_pos, width, height];
    
    % Create axes and show image
    %axes('Position', position);
    % Display the image in a subplot
    subplot(numRows, numCols, i);
    %imagesc(img);
    imagesc(img, [-0.001 0.002]); % uncomment for intensity range
    colorbar; % Add a colorbar to each subplot
    %title(sprintf('alpha = %.3f', numbers(idx(i)))); % Use the sorted number as the title with prefix "alpha ="
    title(sprintf('padding = %.0f', numbers(idx(i)))); % Use the sorted number as the title with prefix "padding ="
    %title(sprintf('iterations = %.0f', numbers(idx(i)))); % Use the sorted number as the title with prefix "padding ="
    
    axis image; % Keep aspect ratio of the image
    %axis off;  % Optionally turn off axis for a cleaner look
end

% Save the figure
% Specify the figure file name

% Base directory for storing images
savePath = '/dtu/cfu/data/userdata18/s220464/storage/all_figures/';
%% Save the figure
exportgraphics(gcf, fullfile(savePath, figureFileName), 'Resolution', 300);

%% Presentation settings
set(groot,'defaultAxesFontName','times')
set(groot,'defaultAxesFontSize',16)
set(0,'DefaultTextFontSize',12);
set(0,'DefaultLineLineWidth',1.2);

%% Display sorted ROI


clear;
clc;
close all;

% ROI_1    
x1 = 870; y1 = 1000;
x2 = 1220; y2 = 1350;
width = x2 - x1;
height = y2 - y1;
roi_1 = [x1, y1, width, height];

% ROI_7
x1 = 920; y1 = 1720;
x2 = 1220; y2 = 2020;
width = x2 - x1;
height = y2 - y1;
roi_7 = [x1, y1, width, height];

% ROI_stitched
x1 = 1250; y1 = 1510;
x2 = 1600; y2 = 1860;
width = x2 - x1;
height = y2 - y1;
roi_1_7 = [x1, y1, width, height];


% Define the directory containing the images
folderPath = '/dtu/cfu/data/userdata18/s220464/storage/PHASE_2/cgls_tikh_stitched_50_iter/';
roi = roi_1_7;
% Get a list of all .tiff files in the folder
files = dir(fullfile(folderPath, '*.tiff')); % Adjust the extension if needed

% Check if files were found
if isempty(files)
    error('No files found. Check the directory and file type.');
end

% Names of files to exclude
%excludedFiles = {'recon_fista_tv_0_idx_0000.tiff', 'recon_fista_tv_10_idx_0000.tiff'};
%excludedFiles = {'stitched_slices_1_7_pad_512_iter_5.tiff'};

% Filter out excluded files
%files = files(~ismember({files.name}, excludedFiles));

% Extract numbers from filenames to sort by
numbers = zeros(length(files), 1);
for i = 1:length(files)
    filename = files(i).name;
    %numStr = regexp(filename, '(?<=recon_fbp_)\d+', 'match', 'once');  
    %numStr = regexp(filename, '(?<=data_pad_)\d+', 'match', 'once');  
    %numStr = regexp(filename, '(?<=stitched_slices_1_7_pad_)\d+', 'match', 'once'); 
    %numStr = regexp(filename, '(?<=stitched_slices_1_7_iter_)\d+', 'match', 'once'); 
    %numStr = regexp(filename, '(?<=stitched_slices_1_7_pad_512_iter_)\d+', 'match', 'once'); 
    numStr = regexp(filename, '(?<=stitched_slices_1_7_alpha_)\d+\.\d+', 'match', 'once');
    %numStr = regexp(filename, '(?<=recon_fista_tv_)\d*\.?\d*(?=_idx)', 'match', 'once');
    %numStr = regexp(filename, '(?<=recon_cgls_)\d*\.?\d*(?=_idx)', 'match', 'once');
    %numStr = regexp(filename, '(?<=recon_sirt_)\d*\.?\d*(?=_idx)', 'match', 'once');
    %numStr = regexp(filename, '(?<=recon_cgls_tikh_)\d*\.?\d*(?=_idx)', 'match', 'once');
    %numStr = regexp(filename, '(?<=recon_sirt_)\d*\.?\d*(?=_idx)', 'match', 'once');
    numbers(i) = str2double(numStr);
end

% Sort files based on the extracted numbers
[~, idx] = sort(numbers);
files = files(idx);

% Initialize a figure for displaying images
figure;
colormap gray; % Set the colormap to gray for all subplots

% Super title for the entire figure
%sgtitle('Slice 1, 200 Iterations');

% Number of rows and columns in subplot grid
numRows = 2;
numCols = 4;

% Loop through each file and load/display the image
for i = 1:length(files)
    % Full path to file
    filePath = fullfile(folderPath, files(i).name);
    
    % Read the image
    img = imread(filePath);
    
    % Crop image to the defined ROI
    croppedImg = imcrop(img, roi);
    
    % Create axes and show image
    subplot(numRows, numCols, i);
    imagesc(croppedImg, [-0.001 0.002]);
    colorbar; % Add a colorbar to each subplot
    title(sprintf('alpha = %.1f', numbers(idx(i)))); % Use the sorted number as the title with prefix "alpha ="
    %title(sprintf('padding = %.0f', numbers(idx(i)))); % Use the sorted number as the title with prefix "padding ="
    %title(sprintf('iterations = %.0f', numbers(idx(i)))); % Use the sorted number as the title with prefix "iterations ="
    axis image; % Keep aspect ratio of the image
    axis off;  
end

%Automatic Scaling by imagesc
%The imagesc function in MATLAB scales the color data in the image to the full range 
% of the current colormap and adjusts the CLim property of the axes accordingly. 
% This means that for each image displayed with imagesc, the pixel intensity values 
% are mapped to span the entire colormap range based on the minimum and maximum values 
% within that specific image or ROI.

%When you crop an image to a specific ROI, if this ROI contains a narrower 
% range of pixel values than the full image, the contrast appears to be enhanced 
% because the colormap is now spread across a smaller range of intensity values. 
% Essentially, the dynamic range of the smaller image (ROI) is different, potentially 
% making the image appear brighter or more contrasted.

% Base directory for storing images
savePath = '/dtu/cfu/data/userdata18/s220464/storage/all_figures/';
%% Save the figure
figureFileName = 'cgls_tikh_1_7_roi.png';
exportgraphics(gcf, fullfile(savePath, figureFileName), 'Resolution', 300);


%% code with caxis range based on minimum and maximum pixel values

clear;
clc;
close all;

% Define the directory containing the images
folderPath = '/dtu/cfu/data/userdata18/s220464/storage/PHASE_2/fista_tv_slice1_200_iter';

% Get a list of all .tiff files in the folder
files = dir(fullfile(folderPath, '*.tiff')); % Adjust the extension if needed

% Check if files were found
if isempty(files)
    error('No files found. Check the directory and file type.');
end

% Names of files to exclude
excludedFiles = {'recon_fista_tv_0_idx_0000.tiff', 'recon_fista_tv_10_idx_0000.tiff'};

% Filter out excluded files
files = files(~ismember({files.name}, excludedFiles));

% Extract numbers from filenames to sort by
numbers = zeros(length(files), 1);
for i = 1:length(files)
    filename = files(i).name;
    numStr = regexp(filename, '(?<=recon_fista_tv_)\d*\.?\d*(?=_idx)', 'match', 'once');
    numbers(i) = str2double(numStr);
end


% Sort files based on the extracted numbers
[~, idx] = sort(numbers);
files = files(idx);

% Initialize variables to find global min and max
globalMin = inf;
globalMax = -inf;

% Loop through each file to find the global min and max
for i = 1:length(files)
    % Full path to file
    filePath = fullfile(folderPath, files(i).name);
    
    % Read the image
    img = imread(filePath);
    
    % Update global min and max
    globalMin = min(globalMin, min(img(:)));
    globalMax = max(globalMax, max(img(:)));
end

% Define the ROI
x1 = 1000; y1 = 1070;
x2 = 1280; y2 = 1500;
roi = [x1, y1, x2 - x1, y2 - y1];  % [X, Y, Width, Height]

% Initialize a figure for displaying images
figure;
colormap gray; % Set the colormap to gray for all subplots

% Super title for the entire figure
sgtitle('Slice 1, 200 Iterations');

% Number of rows and columns in subplot grid
numRows = 2;
numCols = 6;

% Loop through each file and load/display the image
for i = 1:length(files)
    % Full path to file
    filePath = fullfile(folderPath, files(i).name);
    
    % Read the image
    img = imread(filePath);
    
    % Crop image to the defined ROI
    croppedImg = imcrop(img, roi);

    % Display the image in a subplot
    subplot(numRows, numCols, i);
    imagesc(croppedImg);
    caxis([globalMin, globalMax]);  % Apply the global min and max
    colorbar;
    title(sprintf('alpha = %.3f', numbers(idx(i))));  % Title with alpha value
    
    axis image;  % Keep aspect ratio of the image
    axis off;    % Optionally turn off axis for a cleaner look
end