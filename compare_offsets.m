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


%% compare different offset techniques:

close all;

% Load slices
path_1 = '/dtu/cfu/data/userdata18/s220464/storage/PHASE_3/';

path_2 = '/dtu/cfu/data/userdata15/panum_rats/rat_77/elettra_ct_data/recons/rat_77_full_Z2.35_X-1.17_Y1.83/slices/';
path_3 = '/dtu/cfu/data/userdata18/s220464/storage/PHASE_3/recons_rat_77_full_Z2.35_X-1.17_Y1.83/';
slice_1 = imread([path_1, 'rec_1_pad_mask_idx_0000.tiff']);
slice_2 = imread([path_1, 'rec_2_pad_mask_idx_0000.tiff']);
slice_3 = imread([path_1, 'rec_cil_mask_idx_0000.tiff']);
slice_4 = imread([path_2, 'slice_1022.tif']);
slice_5 = imread([path_3, 'slice_1022.tif']);

%% Plots

% plot slices
figure;imagesc(slice_1);colormap(gray(256));axis image;
figure;imagesc(slice_2);colormap(gray(256));axis image;
figure;imagesc(slice_3);colormap(gray(256));axis image;
figure;imagesc(slice_4);colormap(gray(256));axis image;
figure;imagesc(slice_5);colormap(gray(256));axis image;
%%%%%%%%%%%%%%%%%%%%

% Define the paths and filenames
path_1 = '/dtu/cfu/data/userdata18/s220464/storage/PHASE_3/';
path_2 = '/dtu/cfu/data/userdata15/panum_rats/rat_77/elettra_ct_data/recons/rat_77_full_Z2.35_X-1.17_Y1.83/slices/';
path_3 = '/dtu/cfu/data/userdata18/s220464/storage/PHASE_3/recons_rat_77_full_Z2.35_X-1.17_Y1.83/';
filenames = {'rec_1_pad_mask_idx_0000.tiff', 'rec_2_pad_mask_idx_0000.tiff', 'rec_cil_mask_idx_0000.tiff', 'slice_1024.tif', 'rec_cil_mask_2_idx_0000.tiff', 'slice_1024.tif'};
paths = {path_1, path_1, path_1, path_2, path_1, path_3};

% Load slices
slices = cell(1, numel(filenames));
for i = 1:numel(filenames)
    slices{i} = imread([paths{i}, filenames{i}]);
end

% Choose the row for line profile analysis
chosen_row = 1024;  % Example: you can change this as needed

% Create a new figure for the line plots
figure;
hold on;  % Hold on to plot all lines in the same graph
colors = lines(numel(slices));  % Get some distinct colors for each line

% Loop through the slices and plot the line profiles
for i = 1:numel(slices)
    slice = slices{i};
    line_profile = slice(chosen_row, :);  % Extract the chosen row from each slice
    plot(line_profile, 'Color', colors(i, :), 'LineWidth', 2);  % Plot with a unique color and thicker line
end

% Add legend and labels
legend({'0 offset', '-2 offset', 'cil offset -2', 'elettra 1022', 'cil offset +2' 'elletra mine 1022'}, 'Location', 'northeast');
title('Line Profiles Through Chosen Row Across Four Slices');
xlabel('Column Index');
ylabel('Intensity');
hold off;  % Release the hold on the current figure

% Customize the colormap and axis properties if necessary
colormap(gray(256));
axis tight;  % Tighten the axis to fit the line profiles


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% other section
close all;
path_1 = '/dtu/cfu/data/userdata18/s220464/storage/PHASE_3/';

path_2 = '/dtu/cfu/data/userdata15/panum_rats/rat_77/elettra_ct_data/recons/rat_77_full_Z2.35_X-1.17_Y1.83/slices/';
path_3 = '/dtu/cfu/data/userdata18/s220464/storage/PHASE_3/recons_rat_77_full_Z2.35_X-1.17_Y1.83/';
path_4 = '/dtu/cfu/data/userdata18/s220464/storage/PHASE_3/recons_rat_77_full_Z2.35_X-1.17_Y1.83_0/';
path_5 = '/dtu/cfu/data/userdata18/s220464/storage/PHASE_3/recons_rat_77_full_Z2.35_X-1.17_Y1.83_2/';
slice_1 = imread([path_1, 'rec_1_pad_mask_idx_0000.tiff']);    % normal fbp cil reconstruction
slice_2 = imread([path_1, 'rec_2_pad_mask_idx_0000.tiff']);    % offset - 2 done in matlab
slice_3 = imread([path_1, 'rec_cil_mask_idx_0000.tiff']);      % offset -2 done in cil
slice_4 = imread([path_1, 'rec_cil_mask_2_idx_0000.tiff']);    % cil offset but set to +2
slice_5 = imread([path_2, 'slice_1022.tif']);                  % original elettra recon from their folder
slice_6 = imread([path_3, 'slice_1022.tif']);                  % my elettra recon with offset to -2
slice_7 = imread([path_4, 'slice_1022.tif']);                  % my elettra recon with offset to 0
slice_8 = imread([path_5, 'slice_1022.tif']);                  % my elettra recon with offset to 2
% add matlab +2 to compare it with original slice
cil = slice_2;
stp = slice_8;

% conclusions:
% % offset - 2 done in matlab is corresponding to % my elettra recon with
% offset to 2 slice 1022 but cil is darker
% slice 2 correesponds to 8. 
% slice 4 corresponds to 6. 
% slice 1 corresponds to 7
% slice 3 correspodns to 8
% slice 2 is the same as 3. 
% slice 5 corresponds to 2
% slice 5 correspodns to 3
% slice 8 IS THE SAME AS 5     so my +2 reconsturciton is good because
% original offset was also set to 2
%set FILENAMES[6]=rat_77_full_Z2.35_X-1.17_Y1.83
%set OFFSETS[6]=2

% 

% Check if images are RGB, convert to grayscale
if size(cil, 3) == 3
    cil = rgb2gray(cil);
end
if size(stp, 3) == 3
    stp = rgb2gray(stp);
end

% Extract line profiles at row 1023 from each image
%row_index = 1023;
%profile1 = double(stp(1023, :));
%profile2 = double(cil(1023, :));

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
%title('Line profiles through row 1023 of STP and CIL reconstructions, slice Z2.35 X-1.17 Y1.83');
legend('STP reconstruction', 'CIL reconstruction', 'Intensity difference');

% Ensure the plot displays grid lines for better visualization
grid on;


%% line profiles difference\

% Assuming 'profile1' and 'profile2' have been previously defined and extracted
% Assuming 'profile1' and 'profile2' have been previously defined and extracted

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
axis image;
title('Reconstructed with STP');
subplot(1, 2, 2);
imagesc(slice_8);
colormap(gray(256));
axis image;
title('Reconstructed with CIL');



figure;
subplot(1, 2, 1);
imagesc(slice_8, [-0.001 0.002]);
colormap(gray(256));
colorbar;
axis image;
title('Reconstructed with STP');
subplot(1, 2, 2);
imagesc(slice_2, [-0.001 0.002]);
colormap(gray(256));
colorbar;
axis image;
title('Reconstructed with CIL');



%% Presentation settings
set(groot,'defaultAxesFontName','times')
set(groot,'defaultAxesFontSize',16)
set(0,'DefaultTextFontSize',12);
set(0,'DefaultLineLineWidth',1.2);