function Lv=getLv(sz,dim,inds2process)
% get Lv according to the input size. dim=1 means Lv is the one when
% horizontal direvative is considered, dim=2 means when the vertical
% direvative is considered,dim=3 means the derivative on the third axis is
% considered
% Input:
% sz:   size of the image, which can be obtained using size()
% dim:  on which axis the first derivative will be computed
% Output:
% Lv:   the sparse matrix of Lv in eq. (17) of the miccai09 paper. No weight
%       will be considered
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Designed and implemented by Yuanjie Zheng at picsl of UPenn on 07/25/2009
% contact: zheng.vision@gmail.com
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

if ~((dim>=1)&(dim<=length(sz)))
    error('The input value of dim is not correct!');
end

if sz(dim)<2
    Lv=0;
    return;
end

% number of pixels
numPixels=1;
for i=1:length(sz)
    numPixels=numPixels*sz(i);
end

if nargin==2
    inds2process=1:numPixels;
end
% number of the nonzero elements in the sparse matrix Lv
% for the pixels on the first line, there is no derivative, we can simply
% set the row for them to be all zero
% nnz=numPixels*2;

% inds
if length(sz)>1
    inds=reshape(1:numPixels,sz);
else
    inds=(1:numPixels);
end

% % % % creat inds_latter
if length(sz)==1
    inds_latter=zeros(sz,1);
else
    inds_latter=zeros(sz);
end
% inds_latter(:,2:end,:)=inds(:,1:end-1,:)
func_r='inds(';
for i=1:dim-1
    func_r=[func_r ':,'];
end
func_r=[func_r '1:end-1'];
for i=dim+1:length(sz)
    func_r=[func_r ',:'];
end
func_r=[func_r ')'];

func_l='inds_latter(';
for i=1:dim-1
    func_l=[func_l ':,'];
end
func_l=[func_l '2:end'];
for i=dim+1:length(sz)
    func_l=[func_l ',:'];
end
func_l=[func_l ')'];

func=[func_l '=' func_r ';'];
eval(func);

% inds_latter(:,1,:)=inds(:,1,:)
func_r='inds(';
for i=1:dim-1
    func_r=[func_r ':,'];
end
func_r=[func_r '1'];
for i=dim+1:length(sz)
    func_r=[func_r ',:'];
end
func_r=[func_r ')'];

func_l='inds_latter(';
for i=1:dim-1
    func_l=[func_l ':,'];
end
func_l=[func_l '1'];
for i=dim+1:length(sz)
    func_l=[func_l ',:'];
end
func_l=[func_l ')'];

func=[func_l '=' func_r ';'];
eval(func);

% % % % % only process the pixels indexed with inds2process
numPixels2Process=length(inds2process);
nnz=numPixels2Process*2;
inds4Pixel2Proc=1:numPixels2Process;
inds_latter=inds_latter(inds2process);

% % % sparse matrix
ii=zeros(nnz,1); jj=zeros(nnz,1); ss=zeros(nnz,1);

ii(1:numPixels2Process)=inds4Pixel2Proc(:); 
jj(1:numPixels2Process)=inds(inds2process); 
ss(1:numPixels2Process)=1;

ii(1+numPixels2Process:numPixels2Process+numPixels2Process)=inds4Pixel2Proc(:); 
jj(1+numPixels2Process:numPixels2Process+numPixels2Process)=inds_latter(:); 
ss(1+numPixels2Process:numPixels2Process+numPixels2Process)=-1;

Lv = sparse(ii,jj,ss,numPixels2Process,numPixels,nnz);
    
    