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

%% line profile through empty slice
% Define the path and filename
folderPath = '/dtu/cfu/data/userdata18/s220464/storage/PHASE_1/fbp_slice3/'; % Update this to the path where your file is located
filename = 'recon_fista_tv_05_idx_0000.tiff';

% Full path to the file
filePath = fullfile(folderPath, filename);

% Check if the file exists
if ~isfile(filePath)
    error('File does not exist. Check the directory and file name.');
end

% Read the image
img = imread(filePath);

% Extract the line profile at row 1024
profile = double(img(1024, :));

% Create a figure for the subplot arrangement
figure;

% Subplot 1: Show the image with the highlighted row
subplot(1, 2, 1);
imagesc(img);  % Display the image
colormap(gray);  % Use gray colormap
hold on;
plot([1, size(img, 2)], [1024, 1024], 'r', 'LineWidth', 2);  % Highlight row 1024
hold off;
title('Slice with Highlighted Row 1024');
xlabel('Column Index');
ylabel('Row Index');
colorbar;  % Optional: include a colorbar
axis on;  % Turn on axis ticks

% Subplot 2: Show the line profile
subplot(1, 2, 2);
plot(profile, 'b-', 'LineWidth', 2);  % Plot the profile in blue
title('Line Profile Through Row 1024');
xlabel('Column Index');
ylabel('Intensity Value');
grid on;  % Enable grid for better visualization

% Finalize the plot with labels and legend
title('Line Profiles at Row 1024');
xlabel('Column Index');
ylabel('Intensity Value');
xlim([1 2048]);
legend(legends, 'Location', 'best', 'Interpreter', 'none');
grid on;  
hold off; 

% Base directory for storing images
savePath = '/dtu/cfu/data/userdata18/s220464/storage/all_figures/';
%% Save the figure
figureFileName = 'empty_padding_profiles.png';  % Change only the name here
exportgraphics(gcf, fullfile(savePath, figureFileName), 'Resolution', 300);



%% multiple line profiles

clear;
clc;
close all;

% Define the directory containing the TIFF images
folderPath = '/dtu/cfu/data/userdata18/s220464/storage/PHASE_2/sirt_stitched_no_clip/';

% Get a list of all TIFF files in the folder
files = dir(fullfile(folderPath, '*.tiff'));  % Adjust the extension if needed

if isempty(files)
    error('No TIFF files found. Check the directory and file type.');
end

% Names of files to exclude
%excludedFiles = {'recon_fista_tv_0_idx_0000.tiff', 'recon_fista_tv_10_idx_0000.tiff'};
excludedFiles = {'stitched_slices_1_7_alpha_0.000_500_iter.tiff', 'stitched_slices_1_7_alpha_0.000_2000_iter.tiff','stitched_slices_1_7_alpha_0.000_4000_iter.tiff','stitched_slices_1_7_alpha_0.000_5000_iter.tiff','stitched_slices_1_7_alpha_0.000_3000_iter.tiff''stitched_slices_1_7_alpha_0.000_6000_iter.tiff','stitched_slices_1_7_alpha_0.000_8000_iter.tiff', 'stitched_slices_1_7_alpha_0.000_9000_iter.tiff'};
% Filter out excluded files
files = files(~ismember({files.name}, excludedFiles));


% Extract padding numbers from filenames
numbers = zeros(length(files), 1);
legends = cell(length(files), 1);
for i = 1:length(files)
    filename = files(i).name;
    %numStr = regexp(filename, '(?<=pad_)\d+', 'match', 'once');
    %numStr = regexp(filename, '(?<=recon_fista_tv_)\d*\.?\d*(?=_idx)', 'match', 'once');
    %numStr = regexp(filename, '(?<=recon_fbp_)\d+', 'match', 'once');  
    %numStr = regexp(filename, '(?<=data_pad_)\d+', 'match', 'once');  
    %numStr = regexp(filename, '(?<=stitched_slices_1_7_pad_)\d+', 'match', 'once'); 
    %numStr = regexp(filename, '(?<=stitched_slices_3_4_pad_)\d+', 'match', 'once'); 
    %numStr = regexp(filename, '(?<=stitched_slices_1_7_iter_)\d+', 'match', 'once'); 
    %numStr = regexp(filename, '(?<=stitched_slices_1_7_pad_512_iter_)\d+', 'match', 'once');
    numStr = regexp(filename, '(?<=stitched_slices_1_7_alpha_\d+\.\d+_)\d+(?=_iter)', 'match', 'once');
    %numStr = regexp(filename, '(?<=stitched_slices_1_7_alpha_)\d+', 'match', 'once'); 
    %numStr = regexp(filename, '(?<=stitched_slices_1_7_alpha_)\d+\.\d+', 'match', 'once');
    %numStr = regexp(filename, '(?<=recon_cgls_)\d*\.?\d*(?=_idx)', 'match', 'once');
    %numStr = regexp(filename, '(?<=recon_sirt_)\d*\.?\d*(?=_idx)', 'match', 'once');
    %numStr = regexp(filename, '(?<=recon_cgls_tikh_)\d*\.?\d*(?=_idx)', 'match', 'once');
    %numStr = regexp(filename, '(?<=recon_sirt_)\d*\.?\d*(?=_idx)', 'match', 'once');
    numbers(i) = str2double(numStr);
    %legends{i} = sprintf('%d padding', numbers(i));
    %legends{i} = sprintf('%.3f alpha', numbers(i));
    legends{i} = sprintf('%.0f iterations', numbers(i));
end

% Sort files based on the extracted numbers
[~, idx] = sort(numbers);

% Create a figure for plotting line profiles
figure;
hold on;  % Hold on for multiple plots

% Define colors for plotting, optional
colors = lines(length(files));  % Generates a colormap with distinct colors

% Loop through each sorted file to process and plot line profiles
for i = 1:length(files)
    % Full path to file using sorted index
    filePath = fullfile(folderPath, files(idx(i)).name);
    
    % Read the image
    img = imread(filePath);
    
    % Extract the line profile at row 1024
    profile = double(img(1024, :));
    
    % Plot the line profile using the color cycle
    plot(profile, 'LineWidth', 2, 'Color', colors(i, :));
end

% Finalize the plot with labels and legend
title('Line profiles extracted at row 1024');
xlabel('Column Index');
ylabel('Intensity Value');
legend(legends(idx), 'Location', 'best', 'Interpreter', 'none');
grid on;  % Enable grid for better visualization
xlim([0 2000]);  % Set x-axis limits based on the width of the last image processed
hold off;  % Release the hold

% Base directory for storing images
savePath = '/dtu/cfu/data/userdata18/s220464/storage/all_figures/';
%% Save the figure
figureFileName = 'sirt_profiles.png';  % Change only the name here
exportgraphics(gcf, fullfile(savePath, figureFileName), 'Resolution', 300);


%% Presentation settings
set(groot,'defaultAxesFontName','times')
set(groot,'defaultAxesFontSize',16)
set(0,'DefaultTextFontSize',12);
set(0,'DefaultLineLineWidth',1.2);


%%

% Create a figure for plotting line profiles
figure;
% Subplot for the 512 padding image
subplot(1, 3, 1);
img_512 = imread(fullfile(folderPath, 'recon_fbp_512_idx_0000.tiff'));
imagesc(img_512, [-0.001 0.002]);
colormap(gray);
hold on;
plot([1, size(img_512, 2)], [1024, 1024], 'r', 'LineWidth', 2); % Mark row 1024
hold off;
title('Single image (512 padding)');
xlabel('Column Index');
ylabel('Row Index');
%colorbar;
axis image

% Subplot for line profiles
subplot(1, 3, [2 3]);  % Use remaining space for line profiles
hold on;
colors = lines(length(files));  % Generates a colormap with distinct colors
for i = 1:length(files)
    % Full path to file using sorted index
    filePath = fullfile(folderPath, files(idx(i)).name);
    img = imread(filePath);
    profile = double(img(1024, :));
    plot(profile, 'LineWidth', 2, 'Color', colors(i, :));
end

% Finalize the plot with labels and legend
title('Line profiles extracted at row 1024');
xlabel('Column Index');
ylabel('Intensity Value');
legend(legends(idx), 'Location', 'best', 'Interpreter', 'none');
grid on;
xlim([1 2048]);  % stitched
%xlim([1000 1350]);  % stitched
hold off;

% Base directory for storing images
savePath = '/dtu/cfu/data/userdata18/s220464/storage/all_figures/';
%% Save the figure
figureFileName = 'empty_profiles_subfigure.png';  % Change only the name here
exportgraphics(gcf, fullfile(savePath, figureFileName), 'Resolution', 300);
