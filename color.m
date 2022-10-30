function [ output]= color( input, lambuta)

 img3= min( 255, max( 0, input));
 img3= img3./ 255;
 img2= 0* img3;

for  i= 1: 3

 mean1= mean( mean( img3(:,:, i)));
 var1= std2( img3(:,:, i));
 min1= mean1- lambuta* var1;
 max1= mean1+ lambuta* var1;
 x=( img3(:,:, i)- min1)./( max1- min1);
 x= min( 1, max( x, 0));
 img2(:,:, i)= x/ 1.3;
end

 output= uint8( img2* 255);

