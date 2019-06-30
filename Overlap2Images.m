function [ h ] = Overlap2Images( Ibottom, Iup )
%Overlap2Images shows the image of two double image matrix
%   Input: two images
%   output: figure handle h
Ibottom=NormalizeImage(Ibottom)*2^16;
Iup=NormalizeImage(Iup)*2^16;

figure();
imshow(uint16(Ibottom));
%generate a green mask with the same size as up image. 
green = cat(3, zeros(size(Ibottom)),ones(size(Ibottom)), zeros(size(Ibottom))); 
hold on 
h = imagesc(green); 
hold off 
set(h, 'AlphaData', uint16(Iup)); %shows as green image
end

