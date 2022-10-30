function [ w]= Solver( v, beta, alpha)














 range= 10;
 step= 0.0001;

persistent  lookup_v known_beta xx known_alpha
 ind= find( known_beta== beta& known_alpha== alpha);
if  isempty( known_beta| known_alpha)
 xx=[- range: step: range];
end
if  any( ind)


if ( exist( 'pointOp')== 3)

 w= pointOp( double( v), lookup_v( ind,:),- range, step, 0);
else 
 w= interp1( xx', lookup_v( ind,:)', v(:), 'linear', 'extrap');
 w= reshape( w, size( v, 1), size( v, 2));
end;
else 

 tmp= compute_w( xx, beta, alpha);
 lookup_v=[ lookup_v; tmp(:)'];
 known_beta=[ known_beta, beta];
 known_alpha=[ known_alpha, alpha];


if ( exist( 'pointOp')== 3)

 w= pointOp( double( v), lookup_v(end ,:),- range, step, 0);
else 
 w= interp1( xx', lookup_v(end ,:)', v(:), 'linear', 'extrap');
 w= reshape( w, size( v, 1), size( v, 2));
end;


end






function  w= compute_w( v, beta, alpha)

if ( abs( alpha- 1)< 1e-9)

 w= compute_w1( v, beta);
return ;
end;

if ( abs( alpha- 2/ 3)< 1e-9)

 w= compute_w23( v, beta);
return ;
end;

if ( abs( alpha- 1/ 2)< 1e-9)

 w= compute_w12( v, beta);
return ;
end;



 w= newton_w( v, beta, alpha);


function  w= compute_w23( v, beta)



 epsilon= 1e-6;

 k= 8/( 27* beta^ 3);
 m= ones( size( v))* k;







 v2= v.* v;
 v3= v2.* v;
 v4= v3.* v;
 m2= m.* m;
 m3= m2.* m;


 alpha=- 1.125* v2;
 beta2= 0.25* v3;


 q=- 0.125*( m.* v2);
 r1=- q/ 2+ sqrt(- m3/ 27+( m2.* v4)/ 256);

 u= exp( log( r1)/ 3);
 y= 2*(- 5/ 18* alpha+ u+( m./( 3* u)));

 W= sqrt( alpha./ 3+ y);


 root= zeros( size( v, 1), size( v, 2), 4);
 root(:,:, 1)= 0.75.* v+ 0.5.*( W+ sqrt(-( alpha+ y+ beta2./ W)));
 root(:,:, 2)= 0.75.* v+ 0.5.*( W- sqrt(-( alpha+ y+ beta2./ W)));
 root(:,:, 3)= 0.75.* v+ 0.5.*(- W+ sqrt(-( alpha+ y- beta2./ W)));
 root(:,:, 4)= 0.75.* v+ 0.5.*(- W- sqrt(-( alpha+ y- beta2./ W)));





 v2= repmat( v,[ 1, 1, 4]);
 sv2= sign( v2);
 rsv2= real( root).* sv2;



 root_flag3= sort((( abs( imag( root))< epsilon)&(( rsv2)>( abs( v2)/ 2))&(( rsv2)<( abs( v2)))).* rsv2, 3, 'descend').* sv2;

 w= root_flag3(:,:, 1);


function  w= compute_w12( v, beta)



 epsilon= 1e-6;

 k=- 0.25/ beta^ 2;
 m= ones( size( v))* k.* sign( v);


 t1=( 2/ 3)* v;

 v2= v.* v;
 v3= v2.* v;


 t2= exp( log(- 27* m- 2* v3+( 3* sqrt( 3))* sqrt( 27* m.^ 2+ 4* m.* v3))/ 3);

 t3= v2./ t2;


 root= zeros( size( v, 1), size( v, 2), 3);
 root(:,:, 1)= t1+( 2^( 1/ 3))/ 3* t3+( t2/( 3* 2^( 1/ 3)));
 root(:,:, 2)= t1-(( 1+ i* sqrt( 3))/( 3* 2^( 2/ 3)))* t3-(( 1- i* sqrt( 3))/( 6* 2^( 1/ 3)))* t2;
 root(:,:, 3)= t1-(( 1- i* sqrt( 3))/( 3* 2^( 2/ 3)))* t3-(( 1+ i* sqrt( 3))/( 6* 2^( 1/ 3)))* t2;

 root( find( isnan( root)| isinf( root)))= 0;



 v2= repmat( v,[ 1, 1, 3]);
 sv2= sign( v2);
 rsv2= real( root).* sv2;
 root_flag3= sort((( abs( imag( root))< epsilon)&(( rsv2)>( 2* abs( v2)/ 3))&(( rsv2)<( abs( v2)))).* rsv2, 3, 'descend').* sv2;

 w= root_flag3(:,:, 1);


function  w= compute_w1( v, beta)


 w= max( abs( v)- 1/ beta, 0).* sign( v);

function  w= newton_w( v, beta, alpha)





 iterations= 4;

 x= v;

for  a= 1: iterations
 fd=( alpha)* sign( x).* abs( x).^( alpha- 1)+ beta*( x- v);
 fdd= alpha*( alpha- 1)* abs( x).^( alpha- 2)+ beta;

 x= x- fd./ fdd;
end;

 q= find( isnan( x));
 x( q)= 0;


 z= beta/ 2* v.^ 2;
 f= abs( x).^ alpha+ beta/ 2*( x- v).^ 2;
 w=( f< z).* x;
