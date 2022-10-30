%% Underwater Image Enhancement Using Hyper-Laplacian Reflection Priors 
%% Peixian Zhuang, Jiamin Wu, Fatih Porikli, and Chongyi Li
%% Accepted by IEEE Transactions on Image Processing, 2022
%% The parameter ecoff should be adjusted for different underwater images.


%% input raw image
% img1 = double(imread('Coral.png'));
 img1 = double(imread('Y5.jpg')); 
 figure,imshow(uint8(img1));
 
 
 %% Enhance result
 ecoff = 2.0;     
 enhanced4 = HLRP(img1, ecoff);
 figure,imshow(uint8(enhanced4));
 
