function Z=interp2D_FFD_mat(sz,res,phi)
% Z=interp2D_FFD_mat(sz,res,phi), interpolate an image with values of control
% points based on matrix operations.
% 
% Input arguments:
% sz:      2-component vector specifying the size of the image: row and
%          columns respectively
% res:     2-component vector specifying the resolution of the grid of
%          control points. rows and colums, respectively. It's it the
%          number of pixels of sz(1) and sz(2) respectively.
% phi:     mxn matrix denoting value of the grid of control points.
% 
% this version is different from interp2D_FFD in the sense that this one is
% computed with Z=A*phi, representing a linear combination of control
% points' values. I found that matrix operation can increase speed by 10
% times


% % % % % check the inputs of phi
[szCP]=sizeCP(sz,res);
szphi=size(phi);
if szCP(:)~=szphi(:)
    error('size of the input phi matrix looks not correct!');
end

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


[X,Y]=meshgrid(1:sz(2),1:sz(1));
% % % specify u and v for all pixels
% V=mod((Y(:)-1),res(1))./res(1);
% U=mod((X(:)-1),res(2))./res(2);

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
        ss=[ss;B_1{i}(:).*B_2{j}(:)];
    end
end

A=sparse(ii,jj,ss,mm,nn,nnz);

Z=reshape(A*phi(:),sz(1),sz(2));

% ffd b-spline functions
function val=B0(u)
val=(1-u).^3/6;

function val=B1(u)
val=(3*u.^3-6*u.^2+4)/6;

function val=B2(u)
val=(-3*u.^3+3*u.^2+3*u+1)/6;

function val=B3(u)
val=u.^3/6;