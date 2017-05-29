function Estimate()
%im=imread('..\data\vignetting\gallery3_Est.jpg');
im=imread('..\data\Examples\M\white2_Est.png');
im2=imread('..\data\Examples\M\white2_Est_truth.png');
%im3=imread('..\data\vignetting\M\white_Est2.jpg');
%im3=imread('..\data\vignetting\M\Estimate.png');
%im4=imread('..\data\vignetting\M\sea_Est3.jpg');

[length,width]=size(im);
s_im=min(length,width);
Reshape1=im((length-s_im)/2+1:(length+s_im)/2,(width-s_im)/2+1:(width+s_im)/2);
Reshape2=im2((length-s_im)/2+1:(length+s_im)/2,(width-s_im)/2+1:(width+s_im)/2);
%Reshape3=im3((length-s_im)/2+1:(length+s_im)/2,(width-s_im)/2+1:(width+s_im)/2);
%Reshape4=im4((length-s_im)/2+1:(length+s_im)/2,(width-s_im)/2+1:(width+s_im)/2);

x=1:1:s_im;
scale=1/255;
%temp=(x-1)*length+x
y1=double(Reshape1((x-1)*s_im+x))*scale;
y2=double(Reshape2((x-1)*s_im+x))*scale;
%y3=double(Reshape3((x-1)*s_im+x))*scale;
%y4=double(Reshape4((x-1)*s_im+x))*scale;

figure;
plot(1.414*(x-s_im/2),y1,1.414*(x-s_im/2),y2);

axis normal;
%plot(x,y1)



