% function Z=poly2d(sz,D,a)
% 
% Z=poly2d(sz,D,a) does 2-variat polynomial interpolation given the
% parameters of the polynmomial, with equation (14) in the journal version
% of NICESIGN
% 
% Input arguments:
% 
% sz: 2-component vector specifying the size of the image, e.g. [254 360]
% D: scalar value shwoing the degree of the polynormial
% a: vector containing the polynomial parameters in the order determined by
%    equation (14) in the journal version of NICESIGN  
% 
% Ouput argument:
% 
% Z: the interpolation image.
% 
% Note:
% 1. x,y value of the polynormial are 1,2,....,sz(2), 1,2,...,sz(1), repectively
% 2. this function will call the mex file in name poly2d.dll (windows) etc.
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Designed and implemented by Yuanjie Zheng at picsl of UPenn
% Sep. 07, 2009. contact: zheng.vision@gmail.com
% All rights reserved
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
