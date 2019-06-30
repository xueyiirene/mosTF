clear all;close all;clc
%% overlapping vertical and horizontal reconstruction images
path='D:\1035-multilineTF\Prime0.5umperpixel_version';
A = imread([path '\redchannel_beads.tif'], 1) ; 
B = imread([path '\yellowchannel_beads.tif'], 1) ; 
load([path '\BG_PCO.mat']);
load([path '\BG_prime.mat']);
A1=double(A)-BG_prime;
B1=double(B)-BG_PCO;
A1=fliplr(A1);

A2=(A1-min(min(A1)))/max(max((A1-min(min(A1)))));
B2=(B1-min(min(B1)))/max(max((B1-min(min(B1)))));
[x,y]=meshgrid(1:500);
[X,Y]=meshgrid(linspace(1,500,800));
A3=interp2(x,y,A2,X,Y);

%binary imag    e for easy registration
Vthres=0.2*max(A3(:));
Hthres=0.05*max(B2(:));
A3(A3>Vthres)=max(A3(:));
A3(A3<=Vthres)=min((A3(:)));
B2(B2>Hthres)=max(B2(:));
B2(B2<=Hthres)=min(B2(:));
% % show overlapping two images
% figure();
% imshow(uint16(IrecV));
% green = cat(3, zeros(size(IrecV)),ones(size(IrecV)), zeros(size(IrecV))); 
% hold on 
% h = imagesc(green); 
% hold off 
% set(h, 'AlphaData', uint16(IrecH)); 
% registration of two images
% registrationEstimator;%automatically registration
cpselect(A3,B2);%manually select control points
%%
mytform = fitgeotrans(movingPoints, fixedPoints, 'projective');%change points names
A4= imwarp(A3, mytform,'OutputView',imref2d(size(A3)));%keep the same size after transform
%check registration effect
figure();
imshow(A4);
green = cat(3, zeros(size(B2)),ones(size(B2)), zeros(size(B2))); 
hold on 
h = imagesc(green); 
hold off 
set(h, 'AlphaData', A4);
%% If the registration looks good, transform the original image 
IrecVr=IrecV(1:512,:);
IrecHr= imwarp(IrecH(:,1:512), mytform,'OutputView',imref2d(size(IrecV1)));
%% Fourier domain

