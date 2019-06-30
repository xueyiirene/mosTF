clear all;close all;clc

sigma=2;
[x,y]=meshgrid(linspace(-5,5,100));
Pex=exp(-x.^2);
Pem=exp(-x.^2/sigma^2-y.^2/sigma^2);
Pem1=sum(Pem,2);
[x1,y1]=meshgrid(Pem1/max(Pem1(:)));
Peff=Pex.*y1;

figure(1);subplot(1,4,1);imagesc(Pex);title('P_{ex}: 1D Gaussian');axis image;
subplot(1,4,2);imagesc(Pem);title('P_{em}: 2D Gaussian');axis image;
subplot(1,4,3);imagesc(y1);title('integral P_{em} in 1D');axis image;
subplot(1,4,4);imagesc(Peff);title('reconstructed I_{psf}');
axis image;
% figure(4);imagesc(Peff/max(Peff(:))-Pem/max(Pem(:)));