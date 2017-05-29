function [szCP]=sizeCP(sz,res)
% [szCP]=sizeCP(sz,res) gets size of the grid of control points
% 
% Input arguments:
% sz:      2-component vector specifying the size of the image: row and
%          columns respectively
% res:     2-component vector specifying the resolution of the grid of
%          control points. rows and colums, respectively. It's it the
%          number of pixels of sz(1) and sz(2) respectively.
% 
% Output:
% szCP:    2-component vector specifying the size of the grid of CPs
% 
% Note
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Designed and implemented by Yuanjie Zheng on May 03, 2010 @ PICSL@UPenn
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

szCP=floor((sz(:)-1)./res(:))+1+3;