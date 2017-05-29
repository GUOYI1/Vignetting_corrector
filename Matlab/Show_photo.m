function Show_photo()
im=imread('..\data\screen cut\25.jpg');
figure; 
imshow(im,'border','tight','initialmagnification','fit');
set (gcf,'Position',[0,0,240,185]);
axis normal;
