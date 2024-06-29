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

%% Comparison of STP and CIL reconstructions on selected data, 
% based on line profiles and their difference comparison

%% import images
clear;
close all;
clc;

path_1 = '/dtu/cfu/data/userdata18/s220464/storage/PHASE_3/';
path_2 = '/dtu/cfu/data/userdata18/s220464/storage/PHASE_3/recons_rat_77_full_Z2.35_X-1.17_Y1.83_2/';
slice_1 = imread([path_1, 'rec_2_pad_mask_idx_0000.tiff']);    % offset - 2 done in matlab
slice_2 = imread([path_2, 'slice_1022.tif']);                  % my elettra recon with offset to 2
cil = slice_1;
stp = slice_2;

% Check if images are RGB, convert to grayscale
if size(cil, 3) == 3
    cil = rgb2gray(cil);
end
if size(stp, 3) == 3
    stp = rgb2gray(stp);
end

profile1 = stp(1023, :);
profile2 = cil(1023, :);

% Create a figure to plot the line profiles
figure;
plot(profile1, 'b-', 'LineWidth', 2); % Plot the first profile in blue
hold on; % Hold the plot to overlay the second profile
plot(profile2, 'r-', 'LineWidth', 2); % Plot the second profile in red
% Calculate and plot the difference between the two profiles
difference = profile1 - profile2;
plot(difference, 'k-', 'LineWidth', 2); % Plot the difference in black

hold off;
% Set x-axis limits
xlim([1 length(profile1)]);  % Assuming you want to show all columns
% Add labels and legend
xlabel('Column Index');
ylabel('Intensity Value');
title('Horizontal line profiles and their difference');
legend('STP reconstruction', 'CIL reconstruction', 'Intensity difference');

% Ensure the plot displays grid lines for better visualization
grid on;

% Base directory for storing images
savePath = '/dtu/cfu/data/userdata18/s220464/storage/all_figures/';
%% Save the figure
figureFileName = 'stp_cil_profiles.png';  % Change only the name here
exportgraphics(gcf, fullfile(savePath, figureFileName), 'Resolution', 300);


%% line profiles difference
% Trim 20 elements from the start and end of each profile
profile1 = profile1(501:end-500);
profile2 = profile2(501:end-500);
% Calculate the absolute differences
abs_diff = abs(profile1 - profile2);

% Choose profile1 as the baseline for percentage calculation
% Avoid division by zero by adding a small epsilon where profile1 is zero
epsilon = 1e-5;  % Small constant to avoid division by zero
baseline = profile1 + (profile1 == 0) * epsilon;

% Calculate percentage differences
percentage_diff = (abs_diff ./ baseline) * 100;

% Calculate the average percentage difference across the profile
average_percentage_diff = mean(percentage_diff);

% Display the result
fprintf('The average intensity difference between the two profiles is %.2f%%.\n', average_percentage_diff);

%% plot the two images side by side

figure;
subplot(1, 2, 1);
imagesc(slice_2, [-0.001 0.002]);
colormap(gray(256));
colorbar;
axis image;
title('Reconstructed with STP');
subplot(1, 2, 2);
imagesc(slice_1, [-0.001 0.002]);
colormap(gray(256));
colorbar;
axis image;
title('Reconstructed with CIL');

% Base directory for storing images
savePath = '/dtu/cfu/data/userdata18/s220464/storage/all_figures/';
%% Save the figure
figureFileName = 'stp_cil_recons.png';  % Change only the name here
exportgraphics(gcf, fullfile(savePath, figureFileName), 'Resolution', 300);

%% Presentation settings
set(groot,'defaultAxesFontName','times')
set(groot,'defaultAxesFontSize',16)
set(0,'DefaultTextFontSize',12);
set(0,'DefaultLineLineWidth',1.2);