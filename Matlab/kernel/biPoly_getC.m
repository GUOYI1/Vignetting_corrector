function [Cx,Cy]=biPoly_getC(sz,D)

deltx=1; delty=1;

numcoeff=(D+1)*(D+2)/2;
numpixel=sz(1)*sz(2);

% get coordinates of X and Y
[X,Y]=meshgrid(1:sz(2),1:sz(1));

% compute C for gradient X
X_left=X; X_left(:,2:end)=X(:,1:end-1);%make the left most column of pixels to be 0 for all coeffs
Cx=zeros(numcoeff,numpixel);
num=0;
for t=0:D
    for l=0:D-t
        num=num+1;
        Cx(num,:)=((X(:).^t-X_left(:).^t).*(Y(:).^l))'/deltx;
    end
end

% compute C for gradient Y
Y_up=Y; Y_up(2:end,:)=Y(1:end-1,:);%make the up most column of pixels to be 0 for all coeffs
Cy=zeros(numcoeff,numpixel);
num=0;
for t=0:D
    for l=0:D-t
        num=num+1;
        Cy(num,:)=((X(:).^t).*(Y(:).^l-Y_up(:).^l))'/delty;
    end
end
