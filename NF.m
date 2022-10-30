






function  Nf= NF( im, theta)




if  nargin< 2
 theta= 0.5;
end



[ m, n, t]= size( im);
 PDF= zeros( 2041, t);
 offset= 1021;
 eps= 0.0001;

 im= double( im);
for  k= 1: t
 index= del2( im(:,:, k))* 4+ offset;
 index= reshape( index, m* n, 1);
[ h, h_x]= hist( index, 1: 2041);
 PDF( h_x, k)= h/ sum( h);
end



[ mm, nn]= size( PDF);
 CDF= PDF;
for  i= 2: mm
 CDF( i,:)= CDF( i- 1,:)+ PDF( i,:);
end


 T_d= zeros( t, 1);

for  k= 1: t
 T_d( k)= SearchT2( CDF(:, k), eps);
end


 T_priorD= 0.1446;
 Nd= T_d/ T_priorD;




 PDF= zeros( 512, 512, t);
[ m, n, t]= size( im);

 im= double( im);

for  k= 1: t
 PDF_tmp= zeros( 512);
 r= conv2( im(:,:, k),[ 0,- 1, 1], 'same')+ 256;
 r= reshape( r, m* n, 1);
 c= conv2( im(:,:, k),[ 0,- 1, 1]', 'same')+ 256;
 c= reshape( c, m* n, 1);

 index= sub2ind([ 512, 512], r, c);
[ h, h_x]= hist( index, 1:( 512^ 2));

 PDF_tmp( h_x)= h/ sum( h);
 PDF(:,:, k)= PDF_tmp;
end


[ mm, nn, kk]= size( PDF);
 CDF= zeros( mm, nn, kk);
for  k= 1: kk
 CDF_tmp= cumsum( PDF(:,:, k));
 CDF(:,:, k)= cumsum( CDF_tmp, 2);
end





 T_f= zeros( t, 1);

for  k= 1: t
 T_f( k)= SearchT1_project( CDF(:,:, k), eps);
end

 T_priorG= 0.3754;
 Ng= T_f/ T_priorG;

 Nf= zeros( size( theta, 1), t);
for  i= 1: size( theta, 1)
 Nf( i,:)=( 1- theta( i))* Ng+( theta( i))* Nd;
end

function  T1= SearchT1_project( CDF, eps)
 x=- 256: 255;
 y= x';

 CDFx= pi*( sum( CDF)/ 255- 0.5);
 CDFy= pi*( sum( CDF, 2)/ 255- 0.5);

 right= 1;
 left= 0;
while ( right- left> eps)
 midleft= left+( right- left)/ 3;
 midright= right-( right- left)/ 3;
 tt= CDFx- atan( midleft* x);
 tt= tt* tt';
 tt2= CDFx- atan( midright* x);
 tt2= tt2* tt2';
if  tt<= tt2
 right= midright;
else 
 left= midleft;
end
end

 Tx=( left+ right)/ 2;

 right= 1;
 left= 0;
while ( right- left> eps)
 midleft= left+( right- left)/ 3;
 midright= right-( right- left)/ 3;
 tt= CDFy- atan( midleft* y);
 tt= tt'* tt;
 tt2= CDFy- atan( midright* y);
 tt2= tt2'* tt2;
if  tt<= tt2
 right= midright;
else 
 left= midleft;
end
end

 Ty=( left+ right)/ 2;

 T1=( Tx+ Ty)/ 2;


function  T2= SearchT2( CDF, eps)
 left= 0;
 right= 1;

 x=- 1020: 1020;
 x= x';

 CDF= pi*( CDF- 0.5);

while ( right- left> eps)
 midleft= left+( right- left)/ 3;
 midright= right-( right- left)/ 3;

if  CDF1D_error( CDF, x, midleft)<= CDF1D_error( CDF, x, midright)
 right= midright;
else 
 left= midleft;
end
end

 T2=( left+ right)/ 2;

function  v= CDF1D_error( CDF, x, T)
 model= atan( T* x);

 tmp= CDF- model;
 v= tmp'* tmp;
