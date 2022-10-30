
function  eigsDtD= getC( F)
 fx=[ 1,- 1];
 fy=[ 1;- 1];
[ N, M, D]= size( F);
 sizeI2D=[ N, M];
 otfFx= psf2otf( fx, sizeI2D);
 otfFy= psf2otf( fy, sizeI2D);
 eigsDtD= abs( otfFx).^ 2+ abs( otfFy).^ 2;
