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


% Define the folder containing the data
folderPath = '/dtu/cfu//data/userdata18/s220464/storage/PHASE_3/projections_rat_77_full_Z2.35_X-1.17_Y1.83/';

% Initialize matrices to hold the extracted data
X_proj = [];
X_dark = [];
X_flat = [];

% Define the offset
offset = -2; % Example: shift 10 columns to the left

% Load images if they exist using the updated loadImages function
X_proj = loadImages(folderPath, 'tomo', 0:1799, 1023, offset);

% Save the extracted slices into a .mat file
saveFilePath = '/dtu/cfu/data/userdata18/s220464/storage/PHASE_3/slice5_offset_-2_1023_data.mat';
save(saveFilePath, 'X_proj', 'X_dark', 'X_flat');

% Function to check and load images with offset application
function X = loadImages(folder, prefix, range, row, offset)
    X = [];
    for i = range
        fileName = fullfile(folder, sprintf('%s_%04d.tif', prefix, i));
        if isfile(fileName)
            fprintf('Processing %s\n', fileName);  % Print the name of the file being processed
            img = imread(fileName);
            img = shiftImage(img, offset);  % Apply the offset before extracting the row
            X(:, end+1) = img(row, :);  % Extract the specified row
        end
    end
end

% Function to apply columnwise shift to an image
function shiftedImg = shiftImage(img, offset)
    if offset > 0
        shiftedImg = [img(:, offset+1:end), repmat(img(:, end), 1, offset)];
    elseif offset < 0
        shiftedImg = [repmat(img(:, 1), 1, -offset), img(:, 1:end+offset)];
    else
        shiftedImg = img;  % No offset applied
    end
end