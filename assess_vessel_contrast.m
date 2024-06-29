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



function assess_vessel_contrast(line_profile1, line_profile2)
    % Find peaks and troughs in the line profiles
    % Peaks represent vessel walls, troughs represent inside of vessels
    [pks1, locs_pks1] = findpeaks(line_profile1);
    [pks2, locs_pks2] = findpeaks(line_profile2);
    
    [troughs1, locs_troughs1] = findpeaks(-line_profile1);
    troughs1 = -troughs1;
    [troughs2, locs_troughs2] = findpeaks(-line_profile2);
    troughs2 = -troughs2;
    
    % Calculate average intensity of peaks and troughs
    avg_peak_intensity1 = mean(pks1);
    avg_peak_intensity2 = mean(pks2);
    avg_trough_intensity1 = mean(troughs1);
    avg_trough_intensity2 = mean(troughs2);
    
    % Calculate contrast
    contrast1 = avg_peak_intensity1 - avg_trough_intensity1;
    contrast2 = avg_peak_intensity2 - avg_trough_intensity2;
    
    % Display results
    fprintf('Image 1 - Average Peak Intensity: %.2f, Average Trough Intensity: %.2f, Contrast: %.2f\n', ...
        avg_peak_intensity1, avg_trough_intensity1, contrast1);
    fprintf('Image 2 - Average Peak Intensity: %.2f, Average Trough Intensity: %.2f, Contrast: %.2f\n', ...
        avg_peak_intensity2, avg_trough_intensity2, contrast2);
    
    % Assess which image has better contrast
    if contrast1 > contrast2
        fprintf('Image 1 has better contrast of the walls and inside of the blood vessels.\n');
    elseif contrast2 > contrast1
        fprintf('Image 2 has better contrast of the walls and inside of the blood vessels.\n');
    else
        fprintf('Both images have equal contrast of the walls and inside of the blood vessels.\n');
    end
end