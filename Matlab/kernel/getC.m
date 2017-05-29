function C=getC(sz)
% C=getC(sz); Get the correspondence matrix for which each row has one
% nonzero value (value is 1) showing which radius this pixel corresponds to.
% Currently, it can work only for 2D image. 
% Input: 
% sz:   size of the image, which can be got using size()
% Output:
% C:    N_pxN_b sparse matrix showing what radius each pixel corresponds
% 
% Note:
% If the maximum value of radius is maxR, then possible radii are from 0 to
% maxR. Therefore, N_b=maxR+1
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Designed and implemented by Yuanjie Zheng at picsl of UPenn on 07/25/2009
% contact: zheng.vision@gmail.com
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% get radius value of each pixel
centx=round((sz(2)+1)/2); 
centy=round((sz(1)+1)/2);
[X,Y]=meshgrid(1-centx:sz(2)-centx,1-centy:sz(1)-centy);
R=round(sqrt(X.^2+Y.^2));

% max R and number of R
maxR=max(R(:));
numR=maxR+1;

% ind of R value for each pixel
R_ind=R+1;

% number of pixels
numPixels=sz(1)*sz(2);

% nonzero values of the correspondence matrix are all 1
allones=ones(numPixels,1);

% number of the nonzero elements in the sparse matrix of C
nnz=numPixels;

% form the sparse matrix
C = sparse(1:numPixels,R_ind(:),allones,numPixels,numR,nnz);