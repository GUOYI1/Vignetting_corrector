function demo_bias()
% demo for bias correction

% path and name of input image file
fn='..\data\bias\IMG_8282.jpg';

% setting
addpath(genpath('.'));
checkMex;% compile poly2d.cpp if neccesary. You need to set up your compiler configuration with "mex -setup" beforehand.
addpath(genpath('.'));

% read and process the image
im=imread(fn);

tic
[im_bias_corrected,im_bias]=biasCorrection_bipoly(im);%you may need to input the second parameter and adapt it to the bias effect in the input image
toc

% show results
figure; 
subplot(1,3,1); imshow(im); title('Given Image');
subplot(1,3,2); imshow(im_bias_corrected); title('Bias Corrected Image')
subplot(1,3,3); imshow(im_bias); title('Estimated Bias')
