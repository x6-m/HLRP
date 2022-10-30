function [ Lm]= IA( Y, min1, max1)

 Lr= round( Y);

 Lg= atan( double( Lr));

 L= max( Lr(:))+ 1;
 mp= zeros( 1, L);

 LI= LC( double( Lr), 5);
 sum1= sum( sum( Lg.* LI));

for  k= 1: L

 t=( Lr== k- 1);
 mp( k)= sum( sum( Lg( t).* LI( t)))/ sum1;
 clear t;
end

 cL= cumsum( mp);
 z= 0: 1: max1;
 z= double( z);
 s= atan(( z- min1));

 sum2= sum( s);
 cf= cumsum( s)/ sum2;

 tmp= cell( 1, L);
 index= zeros( 1, L);
 a= zeros( 1, L);
for  j= 1: L
 tmp{ j}= abs( cf- cL( j));
[ a( j), index( j)]= min( tmp{ j});
end

 Lm= zeros( size( Lr));
for  i= 1: L
 Lm( Lr== i- 1)=( index( i)- 1);
end


end

