function [im_corrected,im_vignetting]=vignCorrection_nonPara(im_given,varargin)
% % 
% % vignCorrection_nonPara can esitmate and remove vignetting effect in gray
% % or color images. Vignetting image is represented nonparamatrically as
% % shown in the reference paper.
% % 
% % Examples:
% % vignCorrection_nonPara(im_given)
% % vignCorrection_nonPara(im_given,labmda)
% % vignCorrection_nonPara(im_given,labmda,itrNum)
% % vignCorrection_nonPara(im_given,labmda,itrNum,alpha)
% % 
% % Input Arguments:
% % im_given: Matrix of image data, obtained with imread or other ways.
% % labmda:   Parameter adjusting the smoothness of the resulting vignetting
% %           image. Default is 0.3. See Eq. (16)
% % itrNum:   Number of iterations for the IRLS. Default is 4.
% % alpha:    Exponent value of gradient distribution, in Eq. (14)  
% % 
% % Ouput Arguments:
% % im_corrected:       image after vignetting correction
% % im_vignetting:      image of the estimated vignetting
% % 
% % Reference:
% % Y. Zheng, J. Yu, S. B. Kang, S. Lin, and C. Kambhamettu, “Single-image vignetting
% % correction using radial gradient symmetry,?in IEEE Conference on
% % Computer Vision and Pattern Recognition, June 2008, pp. 1?.
% % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % Re-implemented by Yuanjie Zheng on April 16, 2010 @PICSL@UPenn,
% % based on the original work for the conference paper.
% % Contact: zheng.vision@gmail.com
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% % % % % % % % % % % % % % % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% check input
% % % % % % % % % % % % % % % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if nargin==0
    error('You need to provide image data for vignCorrection_nonPara.m. Try "help vignCorrection_nonPara" to get more information.');
end

% % % % % % % % % % % % % % % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% set up parameters
% % % % % % % % % % % % % % % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% variable paras
[lambda,itrNum,alpha]=parseInputs(nargin-1,varargin);
% fixed paras
epsilon=0.000001;% perturbation on B
shift=1; % shift to aoivd computation of log(0)

% % % % % % % % % % % % % % % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% downsample to increase speed 
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
%Z_shift=log(im_data+shift+eps); %%%%%%%%%%important
Z_shift=log(double(im_data)+shift); %%%%%%%%%%important
% Z_orn=log(im_data+eps); 

% initial W to 1 for each pixel
vector_W=ones(numPixels,1);
W=W_vec2sparse(vector_W); 

rg=zeros(sz(1),sz(2));
for j=2:1:sz(1)
	for i=2:1:sz(2)
		cx = int32(i - sz(1)*0.5);
		cy = int32(j - sz(2)*0.5);
         %calculate the gradient
        dx = Z_shift(j,i)-Z_shift(j,i-1);
		dy = Z_shift(j,i)-Z_shift(j-1,i);
		%calculate the radial gradient
        cx=double(cx);
        cy=double(cy);
		rg_value = (cx*dx+cy*dy) / sqrt( cx*cx+cy*cy + epsilon );
		rg(j,i) = rg_value;
    end
end
%rg=reshape(rg,sz(1)*sz(2),1);

% get C
C=getC(sz);
numR=size(C,2);

% compute Lvi
R=zeros(sz(1),sz(2));
for j=1:1:sz(1)
	for i=1:1:sz(2)
		cx = int32(i - sz(1)*0.5);
		cy = int32(j - sz(2)*0.5);
        R(j,i) = int32(sqrt(double(cx*cx + cy*cy)));
      end
end
R=reshape(R,sz(1)*sz(2),1);
A=zeros(sz(1)*sz(2),numR);
for i=1:1:sz(1)*sz(2)
    if(R(i,1)<numR+1&&R(i,1)>1)
        A(i,R(i,1))=1;
        A(i,R(i,1)-1)=-1;
    end
  
end
%Lvxny=getLv(sz,1)+getLv(sz,2);
%Lvxny=getLv(sz,1);


% compute Lvii
Lxxnyy=getLxx(numR,1);

%  identiy matrix
myI=speye(numR,numR);

epsilon=0.03;
Gamma1=epsilon*myI;
%Gamma2=lambda*(2*numPixels/numR)*Lxxnyy;
Gamma2=lambda*(2*numPixels/numR)*Lxxnyy;
% % % % % % % % % % % % % % % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Iterative reweighted least squares
% % % % % % % % % % % % % % % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
colNum=itrNum+1;
ims{1}=im_given; tts{1}='Given Image';
ims{colNum+1}=uint8(ones(sz(1),sz(2))*255); tts{colNum+1}='Intial Vignetting';
ims{2*colNum+1}=uint8(reshape(vector_W,sz)*255); tts{2*colNum+1}='Intial Weight';
for k=1:itrNum
    disp(['IRLS: iteration ' num2str(k) ' ... ...'])
    
    %solve B
    %right=W*Lvxny*Z_shift(:);
    %L=Lvxny*Z_shift(:);
    %A=W*Lvxny;
    %B_r=((C'*A'*A*C+Gamma1'*Gamma1+Gamma2'*Gamma2))\(C'*A'*right); 
    
    
    right=W*rg(:); 
    %A=Lvxny*C;
    B_r=((A'*W*W*A+Gamma1'*Gamma1+Gamma2'*Gamma2))\(A'*W*right);
    B=C*B_r;
    
    bg=zeros(sz(1),sz(2));
    for j=2:1:sz(1)
        for i=2:1:sz(2)
            cx = int32(i - sz(1)*0.5);
            cy = int32(j - sz(2)*0.5);
            radius = int32(sqrt( double(cx*cx + cy*cy ) ));
            radius=max(2, min(numR, radius) );
            bg(j,i)=B_r(radius,1)-B_r(radius-1,1);
            
        end
    end
    %bg=reshape(bg,sz(1)*sz(2),1);
    % update W
    %S1=abs(Lvxny*B(:)-Lvxny*Z_shift(:)); S2=alpha*(S1).^(alpha-1);
    S1=abs(bg-rg); S2=alpha*(S1).^(alpha-1);
    vector_W=exp(-S1).*(1-exp(-S2));    
    W=W_vec2sparse(vector_W);
    
    image_W=uint8(reshape(vector_W,sz)*255);
    

    
    % 
    %B=B-mean(B(:));%this can make the mean value of the corrected image equals to the one of the original image, this also makes the estimated bias's values around 1   
    %get image for showing estimated vignetting
      
    b=exp(B-(max(B(:))));% such that b is within [0 1], for display   

    image_b=reshape(b,sz);   
    image_b=uint8(image_b*255);
    
    
    % get corrected image
    X=log(double(im_given_sampled))-repmat(log(reshape(b,sz(1),sz(2))),[1 1 DIM]);
    X_sample=Z_shift./reshape(b,sz(1),sz(2));
    %X=log(double(im_given_sampled)+shift)-repmat(reshape(B,sz(1),sz(2)),[1 1 DIM]);
    x=exp(X);
    %x=exp(X)-shift;
    x=uint8(x);
    
    % put all images in a cell array in order to display them together
    ims{1+k}=x; tts{1+k}=['Corrected: Iter ' num2str(k)];
    ims{colNum+1+k}=image_b; tts{colNum+1+k}=['Vignetting: Iter ' num2str(k)];
    ims{2*colNum+1+k}=image_W; tts{2*colNum+1+k}=['Weight: Iter ' num2str(k)];
end

figure,imdisp(ims,'Size',[3 1+itrNum],'Title',tts);
image_W=imresize(image_W,[size(im_given,1) size(im_given,2)],'bicubic');

% % % % % % % % % % % % % % % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% upsample
% % % % % % % % % % % % % % % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% bias
im_vignetting=imresize(image_b,[size(im_given,1) size(im_given,2)],'bicubic');

% bias corrected image

%im_bias_4compute=reshape(exp(B),sz)*255;
im_bias_4compute=reshape(b,sz)*255;
im_bias_4compute=imresize(im_bias_4compute,[size(im_given,1) size(im_given,2)],'bicubic');

data_bias=double(im_bias_4compute)/255;

X=log(double(im_given))-repmat(log(data_bias),[1 1 DIM]);
x=exp(X);

%X=log(double(im_given)+shift)-repmat(log(data_bias),[1 1 DIM]);
%x=exp(X)-shift;


im_corrected=uint8(x);

%gradient of im_corrected
rg_corrected=zeros(sz(1),sz(2));
rg_corrected=rg-bg;

%Show histogram of gradient
rg_whole=imresize(rg,[size(im_given,1) size(im_given,2)],'bicubic');
rg_whole=roundn(rg_whole,-2);
range_z=int32(max(abs(max(max(rg_whole))*100),abs(min(min(rg_whole))*100)));
range_z_1=range_z*2+1;
Z_R=zeros(range_z_1,1);
for j=1:1:size(im_given,1)
	for i=1:1:size(im_given,2)
        Z_R(int32(rg_whole(j,i)*100+range_z+1))=Z_R(int32(rg_whole(j,i)*100+range_z+1))+1;
    end
end
Z_R=Z_R/(size(im_given,1)*size(im_given,2));

corrected_whole=imresize(rg_corrected,[size(im_given,1) size(im_given,2)],'bicubic');
corrected_whole=roundn(corrected_whole,-2);
range_I=int32(max(abs(max(max(corrected_whole))*100),abs(min(min(corrected_whole))*100)));
range_I_1=range_I*2+1;
I_R=zeros(range_I_1,1);
for j=1:1:size(im_given,1)
	for i=1:1:size(im_given,2)
        I_R(int32(corrected_whole(j,i)*100+range_I+1))=I_R(int32(corrected_whole(j,i)*100+range_I+1))+1;
    end
end
I_R=I_R/(size(im_given,1)*size(im_given,2));

figure;
x=1:1:range_z;
plot(double(x-range_z)/10,Z_R(x+1),'-r',double(x-1)/10,Z_R(x+range_z),'-b');
plot(log(1+abs(double(x-range_z)/10)),Z_R(x+1),'-r',log(1+abs(double(x-1)/10)),Z_R(x+range_z),'-b');
set (gcf,'Position',[0,0,300,200]);
axis normal;

% figure;
% x=1:1:range_I;
% plot(double(x-range_I)/10,I_R(x+1),'r',double(x-1)/10,I_R(x+range_I),'b');
% plot(log(1+abs(double(x-range_I)/10)),I_R(x+1),'r',log(1+abs(double(x-1)/10)),I_R(x+range_I),'b');
% set (gcf,'Position',[0,0,300,200]);
% axis normal;




function [lambda,itrNum,alpha]=parseInputs(numVar,vars)
% default
lambda=0.3;% smoothness
itrNum=4;
alpha=0.6;%diff 2 weight %smaller alpha means more contingent or sharper

% read input
if numVar>=1
    lambda=vars{1};
end
if numVar>=2
    itrNum=vars{2};
end
if numVar>=3
    alpha=vars{3};
end

