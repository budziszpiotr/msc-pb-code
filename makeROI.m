function roiImage = makeROI(image, rect)
    % Convert the image to double precision format
    imageDouble = double(image);
    
    % Crop the image to the specified rectangle
    % rect is expected to be in the format [x y width height]
    roiImage = imcrop(imageDouble, rect);
end