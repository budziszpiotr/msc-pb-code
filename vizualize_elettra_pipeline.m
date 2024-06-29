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

%% demonstration of elettra_pipeline output
% display images from each stage: raw, preprocessed, phase retrieved and
% reconstructed.
clear;
close all;
clc;

%% settings
set(groot,'defaultAxesFontName','times')
set(groot,'defaultAxesFontSize',16)
set(0,'DefaultTextFontSize',12);
set(0,'DefaultLineLineWidth',1.2);

%% Define common axis limits
xLimits = [1 2048]; % Assuming maximum width of 2048
yLimits = [1 2048]; % Set maximum height to 2048 for uniform display

%% raw projections
path = '/dtu/cfu/data/userdata18/s220464/storage/elettra_pipeline/raw_projections/';
dark = imread([path, 'dark_0005.tif']);
flat = imread([path, 'flat_0005.tif']);
proj = imread([path, 'tomo_0005.tif']);

figure;

subplot(1, 2, 1);   
imagesc(dark);
colormap(gray(256));
colorbar;
axis image;
title('Dark image');

subplot(1, 2, 2); 
imagesc(flat);
colormap(gray(256));
colorbar;
axis image;
title('Flat image');

%subplot(1, 3, 3);  
%imagesc(proj);
%colormap(gray(256));
%colorbar;
%axis image;
%title('tomo 1023');
% Base directory for storing images
savePath = '/dtu/cfu/data/userdata18/s220464/storage/all_figures/';
%% Save the figure
figureFileName = 'stp_dark_flat.png';
exportgraphics(gcf, fullfile(savePath, figureFileName), 'Resolution', 300);

%% raw sinogram

clear;
clc; 
close all;

path = '/dtu/cfu/data/userdata18/s220464/storage/elettra_pipeline/raw_sinograms/';
sino = imread([path, 'tomo_1023.tif']);

figure;
imagesc(sino);
colormap(gray(256));
axis image;
title('sino 1023');

%% preprocessing effect

clear;
clc; 
close all;

path_1 = '/dtu/cfu/data/userdata18/s220464/storage/elettra_pipeline/raw_projections/';
path_2 = '/dtu/cfu/data/userdata18/s220464/storage/elettra_pipeline/preprocessed_projections/';
proj_1 = imread([path_1, 'tomo_1023.tif']);
proj_2 = imread([path_2, 'tomo_1023.tif']);

figure;
subplot(1, 2, 1);  
imagesc(proj_1);
colormap(gray(256));
axis image;
title('Projection before pre-processing');
subplot(1, 2, 2);  
imagesc(proj_2);
colormap(gray(256));
axis image;
title('Projection after pre-processing');

% Base directory for storing images
savePath = '/dtu/cfu/data/userdata18/s220464/storage/all_figures/';
%% Save the figure
figureFileName = 'stp_preprocessing.png';
exportgraphics(gcf, fullfile(savePath, figureFileName), 'Resolution', 300);
%% phase retrieval

path_1 = '/dtu/cfu/data/userdata18/s220464/storage/elettra_pipeline/phaseretrieved_projections/';
path_2 = '/dtu/cfu/data/userdata18/s220464/storage/elettra_pipeline/phaseretrieved_sinograms/';
proj = imread([path_1, 'tomo_1023.tif']);
sino = imread([path_2, 'tomo_1023.tif']);

figure;
subplot(1, 2, 1);  
imagesc(proj);
colormap(gray(256));
axis image;
title('projection 1023 after phase retrieval');
subplot(1, 2, 2);  
imagesc(sino);
colormap(gray(256));
axis image;
set(gca, 'XLim', xLimits, 'YLim', yLimits); % Apply defined limits
title('sinogram 1023 after phase retrieval');

%% reconstruction

path = '/dtu/cfu/data/userdata18/s220464/storage/elettra_pipeline/reconstructions/';
recon = imread([path, 'slice_1023.tif']);

figure;
subplot(1, 2, 1);
imagesc(recon);
colormap(gray(256));
axis image;
title('recon 1023');
subplot(1, 2, 2);
imagesc(recon, [-0.001 0.002]);
colormap(gray(256));
axis image;
title('recon 1023 with [-0.001 0.002] range');


%% reconstruction after phase retrieval vs raw reconstruction

path_1 = '/dtu/cfu/data/userdata18/s220464/storage/elettra_pipeline/reconstructions/';
path_2 = '/dtu/cfu/data/userdata18/s220464/storage/elettra_pipeline/raw_reconstructions/';
recon_1 = imread([path_1, 'slice_1023.tif']);
recon_2 = imread([path_2, 'slice_1023.tif']);

figure;
subplot(1, 2, 1);
imagesc(recon_2);
colormap(gray(256));
axis image;
colorbar;
title('Reconstruction without phase retrieval');
subplot(1, 2, 2);
imagesc(recon_1);
colormap(gray(256));
axis image;
colorbar;
title('Reconstruction with phase retrieval');

% Base directory for storing images
savePath = '/dtu/cfu/data/userdata18/s220464/storage/all_figures/';
%% Save the figure
figureFileName = 'stp_phase_recon.png';
exportgraphics(gcf, fullfile(savePath, figureFileName), 'Resolution', 300);



%% CIL SHOW & TELL GRAPH
path = '/dtu/cfu/data/userdata18/s220464/storage/elettra_pipeline/raw_projections/';
proj = imread([path, 'tomo_1023.tif']);
path = '/dtu/cfu/data/userdata18/s220464/storage/elettra_pipeline/raw_sinograms/';
sino = imread([path, 'tomo_1023.tif']);

path = '/dtu/cfu/data/userdata18/s220464/storage/elettra_pipeline/preprocessed_projections/';
proj_1 = imread([path, 'tomo_1023.tif']);
path = '/dtu/cfu/data/userdata18/s220464/storage/elettra_pipeline/preprocessed_sinograms/';
sino_1 = imread([path, 'tomo_1023.tif']);

path = '/dtu/cfu/data/userdata18/s220464/storage/elettra_pipeline/phaseretrieved_projections/';
proj_2 = imread([path, 'tomo_1023.tif']);
path = '/dtu/cfu/data/userdata18/s220464/storage/elettra_pipeline/phaseretrieved_sinograms/';
sino_2 = imread([path, 'tomo_1023.tif']);


% Define common axis limits
xLimits = [1 2048]; % Assuming maximum width of 2048
yLimits = [1 2048]; % Set maximum height to 2048 for uniform display

figure;

subplot(2, 3, 1);   
imagesc(proj);
colormap(gray(256));
colorbar;
axis image;
set(gca, 'XLim', xLimits, 'YLim', yLimits); % Apply defined limits
title('Raw projection');

subplot(2, 3, 2); 
imagesc(proj_1);
colormap(gray(256));
colorbar;
axis image;
set(gca, 'XLim', xLimits, 'YLim', yLimits); % Apply defined limits
title('Preprocessed projection');

subplot(2, 3, 3);  
imagesc(proj_2);
colormap(gray(256));
colorbar;
axis image;
set(gca, 'XLim', xLimits, 'YLim', yLimits); % Apply defined limits
title('Phase retrieved projection');

subplot(2, 3, 4);   
imagesc(sino);
colormap(gray(256));
colorbar;
axis image;
set(gca, 'XLim', xLimits, 'YLim', yLimits); % Apply defined limits
title('Raw sinogram');

subplot(2, 3, 5); 
imagesc(sino_1);
colormap(gray(256));
colorbar;
axis image;
set(gca, 'XLim', xLimits, 'YLim', yLimits); % Apply defined limits
title('Preprocessed sinogram');

subplot(2, 3, 6);  
imagesc(sino_2);
colormap(gray(256));
colorbar;
axis image;
set(gca, 'XLim', xLimits, 'YLim', yLimits); % Apply defined limits
title('Phase retrieved sinogram');

% Base directory for storing images
savePath = '/dtu/cfu/data/userdata18/s220464/storage/all_figures/';
%% Save the figure
figureFileName = 'elettra_pipeline.png';
exportgraphics(gcf, fullfile(savePath, figureFileName), 'Resolution', 300);

%% figurecm

clear;
clc;
close all;

% Load the image
path = '/dtu/cfu/data/userdata18/s220464/storage/elettra_pipeline/raw_sinograms/';
sino = imread([path, 'tomo_1023.tif']);

% Image dimensions
[imgHeight, imgWidth] = size(sino);

% Desired width based on document layout
desiredWidthCm = 15.99773; 

% Calculate the desired height to maintain the image's aspect ratio
desiredHeightCm = desiredWidthCm * (imgHeight / imgWidth);

% Create the figure using the figurecm function
figure_dimensions = [desiredWidthCm, desiredHeightCm];
fig = figurecm(figure_dimensions);
imagesc(sino);
colormap(gray(256));
axis image;
title('Sinogram 1023');

% Saving the figure as PNG with a specific resolution
% Define the file path and name
filename = '/dtu/cfu/data/userdata18/s220464/storage/elettra_pipeline/FormattedSinogram.png';

% Use exportgraphics for better control over the output
% Specify resolution in dots per inch (DPI)
exportgraphics(fig, filename, 'Resolution', 300);