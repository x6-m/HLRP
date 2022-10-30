function  imDst= BF( imSrc, r)








[ hei, wid]= size( imSrc);
 imDst= zeros( size( imSrc));


 imCum= cumsum( imSrc, 1);

 imDst( 1: r+ 1,:)= imCum( 1+ r: 2* r+ 1,:);
 imDst( r+ 2: hei- r,:)= imCum( 2* r+ 2: hei,:)- imCum( 1: hei- 2* r- 1,:);
 imDst( hei- r+ 1: hei,:)= repmat( imCum( hei,:),[ r, 1])- imCum( hei- 2* r: hei- r- 1,:);


 imCum= cumsum( imDst, 2);

 imDst(:, 1: r+ 1)= imCum(:, 1+ r: 2* r+ 1);
 imDst(:, r+ 2: wid- r)= imCum(:, 2* r+ 2: wid)- imCum(:, 1: wid- 2* r- 1);
 imDst(:, wid- r+ 1: wid)= repmat( imCum(:, wid),[ 1, r])- imCum(:, wid- 2* r: wid- r- 1);
end

