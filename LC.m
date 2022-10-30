function  RawDarkChannel= LC( OriImage, sizeOfBlock)



[ mImage, nImage]= size( OriImage);
 PatchSize= sizeOfBlock;
 PadWidth=( PatchSize- 1)/ 2;
 RawDPaddedImage= padarray( OriImage,[ PadWidth, PadWidth], 'symmetric');

 RawDarkChannel= zeros( mImage, nImage);


for  i= 1: mImage




for  j= 1: nImage
 TmpPatch= RawDPaddedImage( i:( i+ PatchSize- 1), j:( j+ PatchSize- 1),:);
 RawDarkChannel( i, j)=( RawDPaddedImage( i, j)- min( TmpPatch(:)))./ max( 0.00001,( max( TmpPatch(:))- min( TmpPatch(:))));
end
end



end

