function Lxx=getLxx(sz,dim)
% get Lxx (or Lyy, Lzz) according to the input size. dim=1 means Lxx is the one when
% horizontal direvative is considered, dim=2 means when the vertical
% direvative is considered,dim=3 means the derivative on the third axis is
% considered
% Input:
% sz:   size of the image, which can be got using size()
% dim:  on which axis the first derivative will be computed
% Output:
% Lxx:   the sparse matrix of Lv in eq. (18) of the miccai09 paper. No
%        resolution will be considered
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Designed and implemented by Yuanjie Zheng at picsl of UPenn on 07/25/2009
% contact: zheng.vision@gmail.com
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

if ~((dim>=1)&(dim<=length(sz)))
    error('The input value of dim is not correct!');
end

if sz(dim)<2
    Lxx=0;
    return;
end

% number of pixels
numPixels=1;
for i=1:length(sz)
    numPixels=numPixels*sz(i);
end
numPixels_used=1;
for i=1:length(sz)
    if dim==i
        thissz=sz(i)-2;
    else
        thissz=sz(i);
    end
    numPixels_used=numPixels_used*thissz;
end

nnz=numPixels_used*3;%number of the nonzero elements in the sparse matrix Lxx

% inds
if length(sz)>1
    inds=reshape(1:numPixels,sz);
else
    inds=(1:numPixels);
end

% compute inds_middle, inds_former, inds_latter
% inds_middle
func_l='inds_middle';

func_r='inds(';
for i=1:dim-1
    func_r=[func_r ':,'];
end
func_r=[func_r '2:end-1'];
for i=dim+1:length(sz)
    func_r=[func_r ',:'];
end
func_r=[func_r ')'];

func=[func_l '=' func_r ';'];
eval(func);

% inds_former
func_l='inds_former';

func_r='inds(';
for i=1:dim-1
    func_r=[func_r ':,'];
end
func_r=[func_r '3:end'];
for i=dim+1:length(sz)
    func_r=[func_r ',:'];
end
func_r=[func_r ')'];

func=[func_l '=' func_r ';'];
eval(func);

% inds_latter
func_l='inds_latter';

func_r='inds(';
for i=1:dim-1
    func_r=[func_r ':,'];
end
func_r=[func_r '1:end-2'];
for i=dim+1:length(sz)
    func_r=[func_r ',:'];
end
func_r=[func_r ')'];

func=[func_l '=' func_r ';'];
eval(func);

% % % creat the sparse matrix
ii=zeros(nnz,1); jj=zeros(nnz,1); ss=zeros(nnz,1);
% ii(1:numPixels_used)=inds_middle(:); jj(1:numPixels_used)=inds_middle(:); ss(1:numPixels_used)=2;
% ii(1+numPixels_used:numPixels_used+numPixels_used)=inds_middle(:); jj(1+numPixels_used:numPixels_used+numPixels_used)=inds_latter(:); ss(1+numPixels_used:numPixels_used+numPixels_used)=-1;
% ii(1+numPixels_used*2:numPixels_used+numPixels_used*2)=inds_middle(:); jj(1+numPixels_used*2:numPixels_used+numPixels_used*2)=inds_former(:); ss(1+numPixels_used*2:numPixels_used+numPixels_used*2)=-1;
ii(1:numPixels_used)=inds_middle(:); jj(1:numPixels_used)=inds_middle(:); ss(1:numPixels_used)=-2;
ii(1+numPixels_used:numPixels_used+numPixels_used)=inds_middle(:); jj(1+numPixels_used:numPixels_used+numPixels_used)=inds_latter(:); ss(1+numPixels_used:numPixels_used+numPixels_used)=1;
ii(1+numPixels_used*2:numPixels_used+numPixels_used*2)=inds_middle(:); jj(1+numPixels_used*2:numPixels_used+numPixels_used*2)=inds_former(:); ss(1+numPixels_used*2:numPixels_used+numPixels_used*2)=1;
Lxx = sparse(ii,jj,ss,numPixels,numPixels,nnz);
    