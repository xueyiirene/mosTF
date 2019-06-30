% reconstruction full frame from individual line scan images
clear all;close all;clc
%% import images
ScanDirection='V';
WrapUpSlide=50;%V: 18; H: 50
% path='D:\1035-multilineTF\4line-mosTF\2direction_2\YFPScarletGephyrinMouse_dendrite2_H';
path='D:\1035-multilineTF\4line-mosTF\2direction_2\LaserOn_63fps_BG';
for ii = 1 : 128
        if ii<11
            Iname=['LaserOn_400ms_BG00000' num2str(ii-1) '.tif'];
        elseif ii>=11 && ii<101
            Iname=['LaserOn_400ms_BG0000' num2str(ii-1) '.tif'];
        else
            Iname=['LaserOn_400ms_BG000' num2str(ii-1) '.tif'];
        end
    
temp_tiff  = imread([path '\' Iname]) ; % read in first image
%concatenate each successive tiff to tiff_stack
if ii==1
    I0=temp_tiff;
else
    I0 = cat(3 , I0, temp_tiff);
end
end

% load('IbgGain5000.mat');%double precision, background matrix

%% rotate, remove background
I1=double(I0)-Ibg;
switch ScanDirection
    case 'V'
% angleV=2;
for ii=1:size(I1,3)
    I2(:,:,ii)=imrotate(I1(:,:,ii),angleV,'bilinear','crop');
end
    case 'H'
% angleH=-97;
for ii=1:size(I1,3)
    I2(:,:,ii)=imrotate(I1(:,:,ii),angleH,'bilinear','crop');
end 
end

%% Wrap the image stack
I3=cat(3,I2(:,:,WrapUpSlide:128),I2(:,:,1:(WrapUpSlide-1)));
I=padarray(I3,[0,5,0],0,'both');
%% Find the position of xstart, xspace
slideNum1=1;
figure();imagesc(I(:,:,slideNum1));
figure();plot(sum(I(:,:,slideNum1),1));
LinePos1=[6,130,256,377];
slideNum2=61;
figure();imagesc(I(:,:,slideNum2));
figure();plot(sum(I(:,:,slideNum2),1));
LinePos64=[74,200,325,446];
slideNum3=116;
figure();imagesc(I(:,:,slideNum3));
figure();plot(sum(I(:,:,slideNum3),1));
LinePos128=[138,263,388,510];
%% medfilt saturation pixels to adjacent area
% I=Iraw;
% Iraw=I;
% for kk=123:127%vertical
for kk=93:104
%     I(270:320,110:155,kk)=medfilt2(I(270:320,110:155,kk),[20,1]);%vertical
    I(93:164,225:253,kk)=medfilt2(I(93:164,225:253,kk),[1,20]);%horizon
    figure();imagesc(I(:,:,kk));
end

% for kk=95
%     figure();imagesc(I(:,:,kk));axis image;colorbar;
% %     saveas(gcf,['D:\Dropbox (MIT)\research\4linesTF\experimentData3\MouseH\' num2str(kk) '.tif']);
% end

%%  Reconstruction
xspace=mean([diff(LinePos1,1,2),diff(LinePos64,1,2),diff(LinePos128,1,2)]);% pixels in between two lines
xstep1=mean((LinePos64-LinePos1)/(slideNum2-slideNum1));
xstep2=mean((LinePos128-LinePos64)/(slideNum3-slideNum2));

% w=3;%the bandwidth or window of the filter, roughly FWHM
xstart=6;
%%
load('Iin.mat');
load('wLUT.mat');
Irec=zeros(size(I(:,:,1)));
w=zeros(size(I(:,:,1)));
for ii=1:size(I,3)
    if ii<size(I,3)/2
        xstep=xstep1;
    else
        xstep=xstep2;
    end
    
    x0=round((ii-1)*xstep+xstart);%center position of the filter
    x1=round((ii-1)*xstep+xstart+xspace);
    x2=round((ii-1)*xstep+xstart+2*xspace); 
    x3=round((ii-1)*xstep+xstart+3*xspace);  
    
    for jj=1:size(Irec,1)
        [temp,Indmin]=min(abs(Iin-I(jj,x0,ii)));
        w(jj,x0)=wLUT(Indmin);
      if (x0-w(jj,x0))<1
         Irec(jj,x0)=sum(I(jj,1:(x0+w(jj,x0)),ii),2);
      else
         Irec(jj,x0)=sum(I(jj,(x0-w(jj,x0)):(x0+w(jj,x0)),ii),2);          
      end
    
      [temp,Indmin]=min(abs(Iin-I(jj,x1,ii)));
      w(jj,x1)=wLUT(Indmin);    
      Irec(jj,x1)=sum(I(jj,(x1-w(jj,x1)):(x1+w(jj,x1)),ii),2);   
    
      [temp,Indmin]=min(abs(Iin-I(jj,x2,ii)));
      w(jj,x2)=wLUT(Indmin);  
      Irec(jj,x2)=sum(I(jj,(x2-w(jj,x2)):(x2+w(jj,x2)),ii),2);
    
      [temp,Indmin]=min(abs(Iin-I(jj,x3,ii)));
      w(jj,x3)=wLUT(Indmin);  
      if (x3+w(jj,x3))>size(I,2)
         Irec(jj,x3)=sum(I(jj,(x3-w(jj,x3)):end,ii),2);
      else
        Irec(jj,x3)=sum(I(jj,(x3-w(jj,x3)):(x3+w(jj,x3)),ii),2);
      end
    end
    
    %remove crosstalk scattering
    switch ScanDirection
    case 'V'
  if (ii>13 && ii<34) || (ii>119 && ii<128) %vertical
      Irec(:,x0)=Irec(:,x0)-Irec(:,x3)*1.5;%-x3*1.5 in vertical,-x0 in horizon
      Irec(:,x2)=Irec(:,x2)-Irec(:,x3)*1.5;
      Irec(:,x1)=Irec(:,x1)-Irec(:,x3)*1.5;
      Irec(:,x3)=Irec(:,x3-20);
  end
      case 'H'
   if ii>93 && ii<104   %horizon       
      Irec(:,x3)=Irec(:,x3)-Irec(:,x0)*0.05;%-x3*1.5 in vertical,-x0 in horizon
      Irec(:,x2)=Irec(:,x2)-Irec(:,x0)*0.05;
      Irec(:,x1)=Irec(:,x1)-Irec(:,x0)*0.5;
      Irec(:,x0)=Irec(:,x0-20);
  end
   end
% figure(1);imagesc(Irec);caxis([0 5*10^4]);title(num2str(ii));pause(1);
% axis image;
% saveas(gcf,['D:\Dropbox (MIT)\research\4linesTF\presentation\process water beads\' num2str(ii) '.tif']);
disp(ii);
end
 figure();imagesc(w);
%%
% remove artificial skips because round xstep
ind_0=find(sum(Irec,1)==0);
for kk=1:length(ind_0)
    ind=ind_0(kk);
    if ind>1 && ind<size(Irec,2)
    Irec(:,ind)=(Irec(:,ind-1)+Irec(:,ind+1))/2;
    end
end
Irec1=Irec(:,6:517);
% Irec1=medfilt2(Irec1,[1 3]);
figure();imagesc(Irec1);
% caxis([0 10^5]);
%% simply sum
Isum=sum(I(:,:,1:128),3);
Isum1=Isum(:,6:517);
%% Plot results
figure();
subplot(1,2,1);
imagesc(Irec1);caxis([0 4.8e5]);
colorbar;title('sum back to lines');
axis image;
subplot(1,2,2);
imagesc(Isum1);caxis([0 4.8e5]);
colorbar;title('directly sum');
axis image;
%% Rotate H and V back to their original direction for overlapping
switch ScanDirection
    case 'H'
IrecH=imrotate(Irec1,-angleH,'bilinear','crop');
IsumH=imrotate(Isum1,-angleH,'bilinear','crop');
    case 'V'
        IrecV=imrotate(Irec1,-angleV,'bilinear','crop');
        IsumV=imrotate(Isum1,-angleV,'bilinear','crop');
end
%% Crop the rotation edge zeros
EdgeW=30;
IrecH1=IrecH(EdgeW:(513-EdgeW),EdgeW:(512-EdgeW));
IrecV1=IrecV(EdgeW:(513-EdgeW),EdgeW:(512-EdgeW));
IsumH1=IsumH(EdgeW:(513-EdgeW),EdgeW:(512-EdgeW));
IsumV1=IsumV(EdgeW:(513-EdgeW),EdgeW:(512-EdgeW));
%% Fourier domain analysis 
cutoff=8;
[fw,fhp,flp]=CreateFilter(IrecV1,cutoff);
border=cutoff*6; 
IrecH1 = padarray(IrecH1,[border,border],'replicate','both');
IrecV1 = padarray(IrecV1,[border,border],'replicate','both');
IsumH1 = padarray(IsumH1,[border,border],'replicate','both');
IsumV1 = padarray(IsumV1,[border,border],'replicate','both');

FrecV=fft2(IrecV1);
FrecH=fft2(IrecH1);
FsumV=fft2(IsumV1);
FsumH=fft2(IsumH1);
%% Image registration
figure();imagesc(IrecV1);
figure();imagesc(IrecH1);
Overlap2Images(IrecV1,IrecH1);

ROI1=IrecV1(150:360,60:480);
ROI2=IrecH1(150:360,60:480);
% use FastFindPeaks2D.m GUI to find the peaks
r=1;
load('FindPeaks_ROI1.mat');
ROIxy1=h;
load('FindPeaks_ROI2.mat');
ROIxy2=h;

%generate binary images for registration
ROI1_1= TargetMap( ROIxy1,r,size(ROI1));
ROI2_1 = TargetMap(ROIxy2,r,size(ROI2));
figure();subplot(1,2,1);imagesc(ROI1);axis image
subplot(1,2,2);imagesc(ROI2);axis image
registrationEstimator;

ROI1reg = registerImages(ROI1,ROI2);
ROI1r=ROI1reg.RegisteredImage;

Overlap2Images(ROI1,ROI2);
Overlap2Images(ROI1r,ROI2);

FrecV=fft2(ROI1r);
FrecH=fft2(ROI2);

%% Filter and combine in Fourier domain
cutoff=8;
[fw,fhp,flp]=CreateFilter(ROI1r,cutoff);
border=cutoff*6; 
ROI1r = padarray(ROI1r,[border,border],'replicate','both');
ROI2 = padarray(ROI2,[border,border],'replicate','both');

FrecV=fft2(ROI1r);
FrecH=fft2(ROI2);

FrecV=FrecV/max(abs(FrecV(:)));
FrecH=FrecH/max(abs(FrecH(:)));

Fw=0.5;
Frec_lp=(FrecV*(1-Fw)+FrecH*Fw).*flp;
Frec_hp=(FrecV*(1-Fw)+FrecH*Fw).*fhp*2;
Frec=Frec_lp+Frec_hp;
Irec2D=ifft2(Frec);
% Fsum=(FsumV*(1-Fw)+FsumH*Fw).*flp+(FsumV*(1-Fw)+FsumH*Fw).*fhp*2;
% Isum2D=ifft2(Fsum);

Irec2D=Irec2D((1+border):border+453,(1+border):border+452);
Isum2D=Isum2D((1+border):border+453,(1+border):border+452);

figure();
% subplot(1,2,1);
imagesc(abs(Irec2D));
colorbar;title('sum back to lines');
axis image;
subplot(1,2,2);
imagesc(abs(Isum2D));
colorbar;title('directly sum');
axis square;

imwrite(uint16(Irec2D/400), ['experimentData3\200nmBeads_Water_reconstruct.tif']);
imwrite(uint16(Isum2D/400), ['experimentData3\200nmBeads_Water_sum.tif']);


%% noise analysis
Ibg_rec=Irec1(286:305,376:395);
Ibg_sum=Isum1(286:305,376:395);

mean_rec=mean(Ibg_rec(:));
std_rec=std(Ibg_rec(:));
mean_sum=mean(Ibg_sum(:));
std_sum=std(Ibg_sum(:));
%% Plot histogram 
figure();
subplot(1,2,1);
histogram(Ibg_rec(:),20);title('background:sum back to lines');
str1={['mean=' num2str(mean_rec)],['std=' num2str(std_rec)]};
xlim([-1 1]*10^4);
annotation('textbox',[0.3,0.8,0.1,0.1],'String',str1,'FitBoxToText','on');
subplot(1,2,2);
histogram(Ibg_sum(:),20);title('background:directly sum');
str2={['mean=' num2str(mean_sum)],['std=' num2str(std_sum)]};
annotation('textbox',[0.8,0.8,0.1,0.1],'String',str2,'FitBoxToText','on');
xlim([-1 1]*10^4);

%% PSF compare
x=(-8:8)*0.4;
Pos=[249,342];
yr=Irec1(Pos(1)-8:Pos(1)+8,Pos(2));
xr=Irec1(Pos(1),Pos(2)-8:Pos(2)+8);

ys=Isum1(Pos(1)-8:Pos(1)+8,Pos(2));
xs=Isum1(Pos(1),Pos(2)-8:Pos(2)+8);

figure();
plot(x,xr);hold on;
plot(x,yr,'r');hold on;
plot(x,xs,'b--');hold on
plot(x,ys,'r--');
grid on;
xlabel('\mum');
ylabel('camera counts');
axis tight
legend('sum back to lines: x','sum back to lines: y','directly sum: x','directly sum: y');
title('PSF comparison');
% imwrite(uint16(Irec_interp), [path '\beads_reconstruct.tif']);



