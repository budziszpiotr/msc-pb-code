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


% function
clc; clear; close all;

%% load slices
tif_directory_1 = '/dtu/cfu/data/userdata18/s220464/storage/PHASE_2/sirt_slice1/';
tif_directory_2 = '/dtu/cfu/data/userdata18/s220464/storage/PHASE_2/sirt_slice7/';
save_dir = '/dtu/cfu/data/userdata18/s220464/storage/PHASE_2/sirt_stitched_no_clip';

slices = ["recon_sirt_5_idx_0000", "recon_sirt_5_idx_0000"];
alpha = 0;
iterations = 5;
pad_value = 512;


% Initialize an array to store file names for stitching
fileset(1) = fullfile(tif_directory_1, slices(1) + ".tiff");
fileset(2) = fullfile(tif_directory_2, slices(2) + ".tiff");

numImages = numel(fileset);

% Initialize variable to hold image sizes.
imageSize = zeros(numImages,2);
%% main loop
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

    %
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

%
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
%
% Initialize the "empty" panorama.
stitched_image = zeros([height width], 'like', fixed);

%
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
    %I_scaled = (I - min(I(:))) / (max(I(:)) - min(I(:)));

    % Define the circle's center and radius
    
    center_x = 1025;
    center_y = 1025;
    radius = 1023;

    c = 8e-2;

    % Create a grid of x, y coordinates 
    % for slices 1:2048 1:2048
    [x, y] = meshgrid(1:2048, 1:2048);

    % Calculate the distance of each point from the center
    distances = sqrt((x - center_x).^2 + (y - center_y).^2);

    % Create a mask for points within the circle
    circle_mask = distances <= radius;

    % Add the constant c to every value inside the circle
    modified_image = I; % Create a copy of the image to modify
    modified_image(circle_mask) = modified_image(circle_mask) + c;

    warpedImage = imwarp(I, tforms(i), 'OutputView', stitched_image_view);
    
    mask = imwarp(true(size(I,1),size(I,2)).*modified_image>0, tforms(i), 'OutputView', stitched_image_view);
    %mask = imwarp(true(size(I,1),size(I,2)).*(I > 0 & I < 0.2), tforms(i), 'OutputView', stitched_image_view);
    
    % Debug: Visualize the mask
    figure(i); imshow(mask, []); title(sprintf('Mask for Image %d', i));
    % Overlay the warpedImage onto the panorama.
    stitched_image = step(blender, stitched_image, warpedImage, mask);
end

% line profile check for clipping

figure; 
plot(stitched_image(1024,:), 'b'); % Plot with blue color

figure; imshow(imadjust(stitched_image));

%% Save the stitched image
%output_file_name = sprintf('%srat_77_3_4_2048_stitched_Z%2.2f.tif', save_dir, z);
slice_set = [1 2 3 4 5 6 7]';
Z = 2.35;
z = 2.35;
%output_file_name = sprintf('stitched_slices_%d_%d_pad_%d_iter_%d_05_1.tiff', slice_set(1), slice_set(7), pad_value, iterations);
%output_file_name = sprintf('stitched_slices_%d_%d_pad_%d.tiff', slice_set(3), slice_set(4), pad_value);
%output_file_name = sprintf('stitched_slices_%d_%d_alpha_%.3f.tiff', slice_set(1), slice_set(7), alpha);
output_file_name = sprintf('stitched_slices_%d_%d_alpha_%.3f_%d_iter.tiff', slice_set(1), slice_set(7), alpha, iterations);
output_filepath = fullfile(save_dir, output_file_name);
save_tif(output_filepath, stitched_image);
%print('-f1', [output_file_name '.png'], '-dpng', '-r100');