function [Cx,Cy]=FFD_getC(sz,res,inds2process)
% [Cx,Cy]=FFD_getC(sz,D,inds2process) get bi-polynomial coefficients for
% 2D image
% 
% Input:
% sz: size of image. vector in length 2
% res:  2-components vector: number pixels for resolution of control points of ffd b-spline. Default is [15 15].
% inds2process: optional; inds of pixels to processed. It can be used to add mask
% 
% Output:
% Cx,Cy: FFD coefficients along x, and y directions.
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Designed and implemented by Yuanjie Zheng at picsl of UPenn on 04/28/2010
% contact: zheng.vision@gmail.com
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% parse inputs
if nargin<2
    res=[15 15];
end
numpixel=sz(1)*sz(2);
if nargin<3
    inds2process=1:numpixel;
end
numpixel2process=length(inds2process);

% prepare
deltx=1; delty=1;

[szCP]=sizeCP(sz,res);
% numCtps=szCP(1)*szCP(2);

% 
% % % % % compuate image of B0, B1, B2, B3 for both x and y
% % compute b spline--v-->y-->sz(1), u-->x-->sz(2)
v_cand=((0:res(1)-1)/res(1))'; B_1=cell(4,1);

temp=repmat(B0(v_cand),ceil(sz(1)/res(1)),sz(2));
B_1{1}=temp(1:sz(1),:); 

temp=repmat(B1(v_cand),ceil(sz(1)/res(1)),sz(2));
B_1{2}=temp(1:sz(1),:); 

temp=repmat(B2(v_cand),ceil(sz(1)/res(1)),sz(2));
B_1{3}=temp(1:sz(1),:); 

temp=repmat(B3(v_cand),ceil(sz(1)/res(1)),sz(2));
B_1{4}=temp(1:sz(1),:); 

% 
u_cand=((0:res(2)-1)/res(2)); B_2=cell(4,1);

temp=repmat(B0(u_cand),sz(1),ceil(sz(2)/res(2)));
B_2{1}=temp(:,1:sz(2)); 

temp=repmat(B1(u_cand),sz(1),ceil(sz(2)/res(2)));
B_2{2}=temp(:,1:sz(2));

temp=repmat(B2(u_cand),sz(1),ceil(sz(2)/res(2)));
B_2{3}=temp(:,1:sz(2));

temp=repmat(B3(u_cand),sz(1),ceil(sz(2)/res(2)));
B_2{4}=temp(:,1:sz(2));

% % % % % compuate image of B0_grad, B1_grad, B2_grad, B3_grad for both x and y
% % compute b spline--v-->y-->sz(1), u-->x-->sz(2)
v_cand=((0:res(1)-1)/res(1))'; B_grad_1=cell(4,1);

temp=repmat(B0_grad(v_cand),ceil(sz(1)/res(1)),sz(2));
B_grad_1{1}=temp(1:sz(1),:); 

temp=repmat(B1_grad(v_cand),ceil(sz(1)/res(1)),sz(2));
B_grad_1{2}=temp(1:sz(1),:); 

temp=repmat(B2_grad(v_cand),ceil(sz(1)/res(1)),sz(2));
B_grad_1{3}=temp(1:sz(1),:); 

temp=repmat(B3_grad(v_cand),ceil(sz(1)/res(1)),sz(2));
B_grad_1{4}=temp(1:sz(1),:); 

% 
u_cand=((0:res(2)-1)/res(2)); B_grad_2=cell(4,1);

temp=repmat(B0_grad(u_cand),sz(1),ceil(sz(2)/res(2)));
B_grad_2{1}=temp(:,1:sz(2)); 

temp=repmat(B1_grad(u_cand),sz(1),ceil(sz(2)/res(2)));
B_grad_2{2}=temp(:,1:sz(2));

temp=repmat(B2_grad(u_cand),sz(1),ceil(sz(2)/res(2)));
B_grad_2{3}=temp(:,1:sz(2));

temp=repmat(B3_grad(u_cand),sz(1),ceil(sz(2)/res(2)));
B_grad_2{4}=temp(:,1:sz(2));


[X,Y]=meshgrid(1:sz(2),1:sz(1));

% % % % % specify image control points's indices
% from up to down
indCpt1{1}=floor((Y(:)-1)/res(1))-1+2; indCpt1{2}=indCpt1{1}+1;
indCpt1{3}=indCpt1{1}+2; indCpt1{4}=indCpt1{1}+3;

% from left to right
indCpt2{1}=floor((X(:)-1)/res(2))-1+2; indCpt2{2}=indCpt2{1}+1;
indCpt2{3}=indCpt2{1}+2; indCpt2{4}=indCpt2{1}+3;

indCpt=cell(4,4);
for i=1:4
    for j=1:4
        indCpt{i,j}=indCpt1{i}+(indCpt2{j}-1)*szCP(1);
    end
end

% % do interpolation by representing linear combination with matrix
% % operation
mm=sz(1)*sz(2); nn=szCP(1)*szCP(2); nnz=sz(1)*sz(2)*4*4;

indpix=(1:sz(1)*sz(2))';

ii=[]; jj=[]; ss=[];
for i=1:4
    for j=1:4
        ii=[ii;indpix];
        jj=[jj;indCpt{i,j}(:)];
        ss=[ss;B_1{i}(:).*B_grad_2{j}(:)];
    end
end
Cx=sparse(ii,jj,ss,mm,nn,nnz);

ii=[]; jj=[]; ss=[];
for i=1:4
    for j=1:4
        ii=[ii;indpix];
        jj=[jj;indCpt{i,j}(:)];
        ss=[ss;B_grad_1{i}(:).*B_2{j}(:)];
    end
end
Cy=sparse(ii,jj,ss,mm,nn,nnz);

Cx=Cx./(deltx*res(2));
Cy=Cy./(delty*res(1)); %inds2process
Cx=Cx(inds2process,:);
Cy=Cy(inds2process,:);

% ffd b-spline functions
function val=B0(u)
val=(1-u).^3/6;

function val=B1(u)
val=(3*u.^3-6*u.^2+4)/6;

function val=B2(u)
val=(-3*u.^3+3*u.^2+3*u+1)/6;

function val=B3(u)
val=u.^3/6;

% ffd b-spline functions
function val=B0_grad(u)
val=-(1-u).^2/2;

function val=B1_grad(u)
val=(3*u.^2-4*u)/2;

function val=B2_grad(u)
val=(-3*u.^2+2*u+1)/2;

function val=B3_grad(u)
val=u.^2/2;