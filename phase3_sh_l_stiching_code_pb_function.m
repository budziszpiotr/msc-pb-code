%  
%
%
%  THE STITCHING PART WITHIN THE CODE IS BORROWED FROM IMAN TAGHAVI, THE ORIGINAL AUTHOR. 
%  THE CODE HAS BEEN ADJUSTED FOR THE PROJECT'S NEEDS.
%
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
%% start
% Directories
%tif_directory = '/dtu/cfu/data/userdata18/s220464/storage/PHASE_1/pictures/';
tif_directory = '/dtu/cfu/data/userdata18/s220464/storage/PHASE_3/';
save_dir = '/dtu/cfu/data/userdata18/s220464/storage/PHASE_3/stitched_sh_l/';

% Define the slices to be processed
% use sl_200_pad_no_mask to gather the tforms and use them to stitch the
% other data. 

%slices = [ "sl_4_pad_idx_0000", "sl_3_pad_idx_0000"];   %- this works!
%slices = ["sl_4_200_pad_mask_idx_0000", "sl_2_200_pad_mask_idx_0000", "sl_3_200_pad_mask_idx_0000", "sl_1_200_pad_mask_idx_0000"];
%slices = ["sl_4_200_pad_no_mask_idx_0000", "sl_2_200_pad_no_mask_idx_0000", "sl_3_200_pad_no_mask_idx_0000", "sl_1_200_pad_no_mask_idx_0000"];  % this works!
%slices = ["noisy_sl_1_200_idx_0000", "noisy_sl_3_200_idx_0000", "noisy_sl_2_200_idx_0000", "noisy_sl_4_200_idx_0000"]; 
%slices = ["sl_4_200_idx_0000", "sl_2_200_idx_0000", "sl_3_200_idx_0000", "sl_1_200_idx_0000"]; 
%slices = ["sl_4_400_pad_no_mask_idx_0000", "sl_2_400_pad_no_mask_idx_0000", "sl_3_400_pad_no_mask_idx_0000", "sl_1_400_pad_no_mask_idx_0000"];
slices = ["sl_2_400_pad_mask_idx_0000", "sl_4_400_pad_mask_idx_0000", "sl_3_400_pad_mask_idx_0000", "sl_1_400_pad_mask_idx_0000"]; %this works
%slices = ["noisy_sl_2_400_idx_0000", "noisy_sl_4_400_idx_0000", "noisy_sl_3_400_idx_0000", "noisy_sl_1_400_idx_0000"];
% Initialize an array to store file names for stitching
fileset = strings(size(slices));

% Generate file paths for each slice
for i = 1:numel(slices)
    filename = sprintf('%s%s.tiff', tif_directory, slices(i));
    fileset(i) = string(filename);
end

numImages = numel(fileset);

% Initialize variable to hold image sizes.
imageSize = zeros(numImages,2);
%% transform
trans = [];
if isempty(trans)

    for n = 1:numImages
        if n == 1
            %             fixed = readimage(imds,n);  % not supported for older than matlab2020
            fixed = read_tif(fileset(n));
            
            continue;
        end
        
        moving = read_tif(fileset(n));
        figure; imshow(fixed, []); title(sprintf('Moving, Iteration %d', n));
        figure; imshow(moving, []); title(sprintf('Moving, Iteration %d', n));
        Reg = registerImages_surf(moving,fixed);
        tforms(n) = Reg.Transformation;
        % Compute T(n) * T(n-1) * ... * T(1)
        tforms(n).T = tforms(n).T * tforms(n-1).T;

        % Save image size.
        imageSize(n,:) = size(moving);

        fixed = moving;
    end

    % Compute the output limits for each transform.
    for i = 1:numel(tforms)
        [xlim(i,:), ylim(i,:)] = outputLimits(tforms(i), [1 imageSize(i,2)], [1 imageSize(i,1)]);
    end

    avgXLim = mean(xlim, 2);
    [~,idx] = sort(avgXLim);
    centerIdx = floor((numel(tforms)+1)/2);
    centerImageIdx = idx(centerIdx);

    %%
    Tinv = invert(tforms(centerImageIdx));
    for i = 1:numel(tforms)
        tforms(i).T = tforms(i).T * Tinv.T;
    end
    trans = tforms;
else
    %tforms = trans;
    %     fixed = readimage(imds,1);  % not supported for older than matlab2020
    %fixed = read_tif(fileset(1));
    %imageSize = repmat(size(fixed),numImages,1);
end

%%
for i = 1:numel(tforms)
    [xlim(i,:), ylim(i,:)] = outputLimits(tforms(i), [1 imageSize(i,2)], [1 imageSize(i,1)]);
end

maxImageSize = max(imageSize);

% Find the minimum and maximum output limits.
xMin = min([1; xlim(:)]);
xMax = max([maxImageSize(2); xlim(:)]);

yMin = min([1; ylim(:)]);
yMax = max([maxImageSize(1); ylim(:)]);

% Width and height of panorama.
width  = round(xMax - xMin);
height = round(yMax - yMin);
%%
% Initialize the "empty" panorama.
stitched_image = zeros([height width], 'like', fixed);

%%
blender = vision.AlphaBlender('Operation', 'Binary mask', ...
    'MaskSource', 'Input port');

% Create a 2-D spatial reference object defining the size of the panorama.
xLimits = [xMin xMax];
yLimits = [yMin yMax];
stitched_image_view = imref2d([height width], xLimits, yLimits);

%% Create the panorama.
for i = 1:numImages
    %     I = readimage(imds, i);  % not supported for older than matlab2020
    I = read_tif(fileset(i));

     % shepp logan:
    center_x = 85;
    center_y = 85;
    radius = 85;

    % set c
    %1e-3 clips
    c = 1e-5;

    % Create a grid of x, y coordinates 
    % for slices 1:2048 1:2048
    [x, y] = meshgrid(1:170, 1:170);

    % Calculate the distance of each point from the center
    distances = sqrt((x - center_x).^2 + (y - center_y).^2);

    % Create a mask for points within the circle
    circle_mask = distances <= radius;

    % Add the constant c to every value inside the circle
    modified_image = I; % Create a copy of the image to modify
    modified_image(circle_mask) = modified_image(circle_mask) + c;

    % At this point, modified_image contains the original image values
    % with c added to the values inside the circle

    % Transform I into the panorama.
    warpedImage = imwarp(I, tforms(i), 'OutputView', stitched_image_view);
   
    mask = imwarp(true(size(I,1),size(I,2)).*modified_image>0.00005, tforms(i), 'OutputView', stitched_image_view);
    figure(i); imshow(mask, []); title(sprintf('Mask for Image %d', i));
    % Debug: Visualize the mask
    %figure; imshow(mask, []); title(sprintf('Mask for Image %d', i));
    % Overlay the warpedImage onto the panorama.
    stitched_image = step(blender, stitched_image, warpedImage, mask);
end

figure; imshow(imadjust(stitched_image));

%% save
output_file_name = 'sh_l_stitched_220624.tiff';
output_filepath = fullfile(save_dir, output_file_name);
save_tif(output_filepath, stitched_image);