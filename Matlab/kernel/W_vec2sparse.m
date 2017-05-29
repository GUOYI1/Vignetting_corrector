function sW=W_vec2sparse(W)
% Inputing vector or matrix representing weight of each pixel, and ouput
% the sparse matrix for which the diagonal elements are each pixel's weight
% Input:
% W: a MxN weight matrix of MNx1 weight vector
% Output:
% sW: a MNxMN sparse matrix for which the diagnal elements are weight
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Designed and implemented by Yuanjie Zheng at picsl of UPenn on 07/25/2009
% contact: zheng.vision@gmail.com
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

numPixels=length(W(:)); nnz=numPixels;
ii=1:numPixels; jj=1:numPixels; ss=W(:);
sW=sparse(ii,jj,ss,numPixels,numPixels,nnz);