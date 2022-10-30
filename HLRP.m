function  enhanced4= HLRP( img1, ecoff)

 up_scale= 10;
 img1= modcrop( img1, up_scale);

 img= double( color( img1, 2.5));%颜色校正
imtool(uint8(img))
 
 cform= makecform( 'srgb2lab');  %RGB转lab
 lab= applycform( uint8( img), cform);
 LL= lab(:,:, 1);
 A= lab(:,:, 2);
 B= lab(:,:, 3);

 img= double( img);
 V_S= double( LL);  % 0`255 double型

 fL=[ 0, 1, 0; 1,- 4, 1; 0, 1, 0];  %二阶梯度
[ row, col, dim]= size( img);
 eigsDtD= getC( V_S);         
 eigsLtL= getL( V_S);         

 beta= 0.01;
 beta2= 0.001;

 lambda= 0.0001;
 lambda2= 0.01;
 lambda3= 0.1;
 lambda4= 0.001;

 gamma= 1;

 alpha1= 0.00001;
 alpha2= 0.001;

 ffL= psf2otf( fL,[ row, col]);
 H= fspecial( 'gaussian', 51, 2);
 L1= filter2( H, V_S);  

 r= 8;
 eps= 0.01;
 L1= GF( L1, V_S, r, eps);   

 x= L1;
 R1= 0* V_S;

 alpha= 1/ 2;
 pbeta= 1;
 pbeta_rate= sqrt( 2);
 pbeta_max= 2^ 8;

 iter= 1;

 tic
for  i= 1: iter


 t1=[ R1(:,end ,:)- R1(:, 1,:),- diff( R1, 1, 2)]; %反射一阶水平差分
 t2=[ R1(end ,:,:)- R1( 1,:,:);- diff( R1, 1, 1)]; %反射一阶垂直差分

 rho2= imfilter( R1, fL, 'circular');   
while  pbeta< pbeta_max
 t11= Solver( t1, pbeta, alpha);    
 t22= Solver( t2, pbeta, alpha);    
 t33= Solver( rho2, pbeta, alpha);  
 pbeta= pbeta* pbeta_rate;
end


 Normin2=[ t22(:,end ,:)- t22(:, 1,:),- diff( t22, 1, 2)];
 Normin2= Normin2+[ t11(end ,:,:)- t11( 1,:,:);- diff( t11, 1, 1)]; %一阶水平垂直

 S1= fft2( V_S./ max( L1, 5)+ beta* lambda*( Normin2)+ beta2* lambda4* conj( ffL).* fft2( t33));
 R1= ifft2( S1./( 1+ beta* lambda.* eigsDtD+ beta2* lambda4.* eigsLtL));
 R1= real( R1);
 R1= max( 0.001, min( R1, 1));


 d1=[ L1(:,end ,:)- L1(:, 1,:),- diff( L1, 1, 2)];    % 水平差分
 d2=[ L1(end ,:,:)- L1( 1,:,:);- diff( L1, 1, 1)];    % 垂直差分

 d11= sign( d1).* max( 0, abs( d1)- 1/( 2* lambda2));
 d22= sign( d2).* max( 0, abs( d2)- 1/( 2* lambda2));

 Nom2=[ d22(:,end ,:)- d22(:, 1,:),- diff( d22, 1, 2)];
 Nom2= Nom2+[ d11(end ,:,:)- d11( 1,:,:);- diff( d11, 1, 1)];

 rho= imfilter( L1, fL, 'circular');
 d33= sign( rho).* max( 0, abs( rho)- 1/( 2* lambda3));


 S11= fft2( V_S./ R1+( 1* x)* gamma+ alpha1.* fft2( Nom2)+ alpha2.* conj( ffL).* fft2( d33));
 L1= ifft2( S11./(( 1+ gamma)+ alpha1.* eigsDtD+ alpha2.* eigsLtL));
 L1= real( L1);
 L1= min( 255, max( L1, V_S));

end
 toc

 R= real( R1);
 L= real( L1);

[ L11]= IA( L, 15, 230);
 R11= adapthisteq( R);

 lab(:,:, 1)= L.* R;
 cform= makecform( 'lab2srgb');
 enhanced1= applycform( lab, cform);
 enhanced1= cast( enhanced1, 'uint8');

 lab(:,:, 1)= L11.* R;
 cform= makecform( 'lab2srgb');
 enhanced2= applycform( lab, cform);
 enhanced2= cast( enhanced2, 'uint8');

 lab(:,:, 1)= L11.* R11;
 cform= makecform( 'lab2srgb');
 enhanced3= applycform( lab, cform);
 enhanced4= cast( enhanced3, 'uint8');

 coff= 0.5;
 enhanced3= enhanced3* coff+ uint8( img)*( 1- coff);


 enhanced4= nature( uint8( enhanced3), ecoff);



 img1= shave( uint8( img1),[ up_scale, up_scale]);
 enhanced4= shave( uint8( enhanced4),[ up_scale, up_scale]);


end

