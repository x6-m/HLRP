function  q= GF( I, p, r, eps)







[ hei, wid]= size( I);
 N= BF( ones( hei, wid), r);

 mean_I= BF( I, r)./ N;
 mean_p= BF( p, r)./ N;
 mean_Ip= BF( I.* p, r)./ N;
 cov_Ip= mean_Ip- mean_I.* mean_p;

 mean_II= BF( I.* I, r)./ N;
 var_I= mean_II- mean_I.* mean_I;

 a= cov_Ip./( var_I+ eps);
 b= mean_p- a.* mean_I;

 mean_a= BF( a, r)./ N;
 mean_b= BF( b, r)./ N;

 q= mean_a.* I+ mean_b;
end