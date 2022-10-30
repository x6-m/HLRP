function  eigsLtL= getL( F)
 f=[ 0, 1, 0; 1,- 4, 1; 0, 1, 0];
[ N, M, D]= size( F);
 sizeI2D=[ N, M];
 otfF= psf2otf( f, sizeI2D);
 eigsLtL= abs( otfF).^ 2;