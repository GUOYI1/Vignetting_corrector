function [im_corrected,im_bias]=biasCorrection_bipoly(im_given,varargin)
% % 
% % biasCorrection_bipoly can esitmate and remove bias field in gray or
% % color images. Bias field is represented paramatrically with the bipolynomial
% % model as shown in the reference paper.
% % 
% % Examples:
% % biasCorrection_nonPara(im_given)
% % biasCorrection_nonPara(im_given,D)
% % biasCorrection_nonPara(im_given,D,itrNum)
% % biasCorrection_nonPara(im_given,D,itrNum,alpha)
% % 
% % Input Arguments:
% % im_given: Matrix of image data, obtained with imread or other ways.
% % D:        Degree of polynomial. Default is 5.
% % itrNum:   Number of iterations for the IRLS. Default is 4.
% % alpha:    Exponent value of gradient distribution
% % 
% % Ouput Arguments:
% % im_corrected:       image with vignetting corrected
% % im_bias:            image of the estimated bias field
% % 
% % Reference:
% % Y. Zheng, M. Grossman, S. Awate, J. Gee “Automatic Correction of
% % Intensity Nonuniformity From Sparseness of Gradient Distribution in
% % Medical Images,” in MICCAI 2009: the 12th International Conference on
% % Medical Image Computing and Computer Assisted Intervention, London, UK,
% % September 20-24, 2009    
% % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % Re-implemented by Yuanjie Zheng on April 16, 2010 @PICSL@UPenn.
% % Contact: zheng.vision@gmail.com
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % Modified on May 01, 2010 to concatenate x, y gradients instead of
% % adding them 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% % % % % % % % % % % % % % % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% check input
% % % % % % % % % % % % % % % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if nargin==0
    error('You need to provide image data for biasCorrection_bipoly.m. Try "help biasCorrection_bipoly" to get more information.');
end

% % % % % % % % % % % % % % % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% set up parameters
% % % % % % % % % % % % % % % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% variable paras
[D,itrNum,alpha]=parseInputs(nargin-1,varargin);

epsilon=0.0001;% perturbation on B
shift=1; % shift to aoivd computation of log(0)

% % % % % % % % % % % % % % % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% downsample to inrease speed 
% % % % % % % % % % % % % % % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dsfact=0.25;
im_given_sampled=imresize(im_given,dsfact);
% % % % % % % % % % % % % % % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% convert color to gray
% % % % % % % % % % % % % % % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
numsz=length(size(im_given_sampled));
if numsz==3
    DIM=3;
    im_gray=rgb2gray(uint8(im_given_sampled));
else
    DIM=1;
    im_gray=im_given_sampled;
end
im_data=double(im_gray);

sz=size(im_data);
numPixels=sz(1)*sz(2);

% % % % % % % % % % % % % % % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% preparitions
% % % % % % % % % % % % % % % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
disp('Preparing ...')
Z_shift=log(im_data+shift); %%%%%%%%%%important

% compute Lvi and Lvii
% 1st derivative
Lvxny=[getLv(sz,2);getLv(sz,1)];


numCoeff=(D+1)*(D+2)/2;
myI=speye(numCoeff,numCoeff);
% initial W for which all pixels are 1
vector_W=ones(numPixels*2,1);
W=W_vec2sparse(vector_W); 

% get C
[Cx,Cy]=biPoly_getC(sz(1:2),D);
C=[Cx';Cy'];

Gamma1=epsilon*myI;


% % % % % % % % % % % % % % % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% IRLS
% % % % % % % % % % % % % % % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
colNum=itrNum+1;
ims{1}=im_given_sampled; tts{1}='Given Image';
ims{colNum+1}=uint8(ones(sz(1),sz(2))*255); tts{colNum+1}='Intial Bias';
ims{2*colNum+1}=uint8(reshape(vector_W(1:numPixels),sz)*255); tts{2*colNum+1}='Intial Weight';
for i=1:4
    disp(['IRLS: iteration ' num2str(i) ' ... ...']);
    
    %solve B
    right=W*Lvxny*Z_shift(:);
    A=W*C;
    warning off all;
    a=((A'*A+Gamma1))\(A'*right);
    B=poly2d(sz(1:2),D,a);
    
    %update W
    S1=abs(Lvxny*B(:)-Lvxny*Z_shift(:)); S2=alpha*(S1).^(alpha-1);
    vector_W=exp(-S1).*(1-exp(-S2));   
    %average w for x and y
    vector_W_temp=[vector_W(numPixels+1:end);vector_W(1:numPixels)];
    vector_W=(vector_W+vector_W_temp)/2;
    W=W_vec2sparse(vector_W);
    
    image_W=uint8(reshape(vector_W(1:numPixels),sz)*255);
    
    B=B-mean(B(:));
    
    b=exp(B-(max(B(:))));% such that b is within [0 1], for display
    image_b=reshape(b,sz);
    image_b=uint8(image_b*255);
    
    X=log(double(im_given_sampled)+shift)-repmat(reshape(B,sz(1),sz(2)),[1 1 DIM]);
    x=exp(X)-shift;
    x=uint8(x);
    
    % put all images in a cell array in order to display
    ims{1+i}=x; tts{1+i}=['Corrected: Iter ' num2str(i)];
    ims{colNum+1+i}=image_b; tts{colNum+1+i}=['Bias: Iter ' num2str(i)];
    ims{2*colNum+1+i}=image_W; tts{2*colNum+1+i}=['Weight: Iter ' num2str(i)];
end
figure,imdisp(ims,'Size',[3 1+itrNum],'Title',tts);

% % % % % % % % % % % % % % % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% upsample
% % % % % % % % % % % % % % % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% bias
im_bias=imresize(image_b,[size(im_given,1) size(im_given,2)],'bicubic');

% bias corrected image
im_bias_4compute=reshape(exp(B),sz)*255;
im_bias_4compute=imresize(im_bias_4compute,[size(im_given,1) size(im_given,2)],'bicubic');

data_bias=double(im_bias_4compute)/255; %data_bias=data_bias-mean(data_bias(:));
X=log(double(im_given)+shift)-repmat(log(data_bias),[1 1 DIM]);
x=exp(X)-shift;
im_corrected=uint8(x);

function [D,itrNum,alpha]=parseInputs(numVar,vars)
% default
D=5;% smoothness
itrNum=4;
alpha=0.1;%diff 2 weight %smaller alpha means more contingent or sharper

% read input
if numVar>=1
    D=vars{1};
end
if numVar>=2
    itrNum=vars{2};
end
if numVar>=3
    alpha=vars{3};
end
