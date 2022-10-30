



function [ enhanced]= nature( enhance, ecoff)

 im= enhance;


 tic
 Nf= NF( im);


 LinearScaling= double( im);
 tic
for  i= 1: size( im, 3)
 tmp_c= LinearScaling(:,:, i);

 LinearScaling(:,:, i)= ecoff*( tmp_c- mean( tmp_c(:)))* Nf( i)+ mean( tmp_c(:));
end

 LinearScaling= uint8( LinearScaling);

 enhanced= LinearScaling;

end

