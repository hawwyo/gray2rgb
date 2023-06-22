clc, clear;

grayImage = imread("starfield.jpg");
refImage = imread("starfield_ref.jpg");


h = 512;
w = 911;
K = 0.0005;

grayImage = imresize(grayImage, [h w]);
grayImage = im2double(grayImage);
grayImage = im2gray(grayImage);


refImage = imresize(refImage, [h w]);
refImage = im2double(refImage);
refGrayImage = im2gray(refImage);

imshow(grayImage);
return;

size(grayImage), size(refImage)


% ----------------------MAKE FILTERS-------------------------------------

imageSize = size(refImage);
numRows = imageSize(1);
numCols = imageSize(2);

wavelengthMin = 4 / sqrt(2);
wavelengthMax = hypot(numRows,numCols);
n = floor(log2(wavelengthMax/wavelengthMin));
wavelength = 2.^(0:(n-2)) * wavelengthMin;

deltaTheta = 30;
orientation = 0 : deltaTheta : (180 - deltaTheta);

g = gabor(wavelength, orientation);


% --------------------FIT--------------------------------------------

refFeatures = imgaborfilt(refGrayImage, g);
% for i = 1:length(g)
%     sigma = 0.5*g(i).Wavelength;
%     refFeatures(:,:,i) = imgaussfilt(refFeatures(:,:,i),K*sigma); 
% end

R = refImage(:,:,1);
G = refImage(:,:,2);
B = refImage(:,:,3);


refFeatures = reshape(refFeatures, [], length(g));
refFeatures = cat(2, refFeatures, reshape(refGrayImage(:,:,1), [], 1));
refFeatures = cat(2, refFeatures, ones(length(refFeatures), 1));

R_coeff = mldivide(refFeatures, reshape(R, [], 1));
G_coeff = mldivide(refFeatures, reshape(G, [], 1));
B_coeff = mldivide(refFeatures, reshape(B, [], 1));

% -------------------------PREDICT---------------------------------------

grayImageFeatures = imgaborfilt(grayImage, g);
% for i = 1:length(g)
%     sigma = 0.5*g(i).Wavelength;
%     grayImageFeatures(:,:,i) = imgaussfilt(grayImageFeatures(:,:,i),K*sigma); 
% end

grayImageFeatures = reshape(grayImageFeatures, [], length(g));
grayImageFeatures = cat(2, grayImageFeatures, reshape(grayImage(:,:,1), [], 1));
grayImageFeatures = cat(2, grayImageFeatures, ones(length(grayImageFeatures), 1));

% grayImage(:,:,1) = reshape(grayImageFeatures * R_coeff, h, []);
grayImage(:,:,2) = reshape(grayImageFeatures * G_coeff, h, []);
grayImage(:,:,3) = reshape(grayImageFeatures * B_coeff, h, []);

imshow(grayImage);
colorbar;
