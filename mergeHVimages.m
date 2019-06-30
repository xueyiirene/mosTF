%reconstruction of mosTF after first step finished using
%reconstruction_4lines.m
%input IrecV, IrecH
%% initialize
need2reg=0;
Norm=0;
flatten=1;%flatten=0 don't flatten; flatten=1 flatten using function flattenfield.m; 
%flatten=2 flatten image using prior skyline
filter=1;
%% 1. flat field

figure();imagesc(IrecV);axis image;colorbar;title('IrecV');
figure();imagesc(IrecH);axis image;colorbar;title('IrecH');
figure();imagesc(IsumV);axis image;colorbar;title('IsumV');
figure();imagesc(IsumH);axis image;colorbar;title('IsumH');
% 

if flatten==1
IV=flattenfield(IrecV,1,50,100);
IH=flattenfield(IrecH,2,50,100);
figure();imagesc(IV);axis image;colorbar;
figure();imagesc(IH);axis image;colorbar;
elseif flatten==0
     IV=IrecV;
     IH=IrecH;
%     IV=IsumV;
%     IH=IsumH;
elseif flatten==2
    IV=IrecV./ repmat(skylineV, [size(IrecV,1), 1]);
    IH = IrecH ./ repmat(skylineH, [1, size(IrecH,2)]);
    IV=max(IrecV(:))/max(IV(:))*IV;
    IH=max(IrecH(:))/max(IH(:))*IH;
end

% figure();plot(sum(IV,2));hold on;plot(sum(IH,2),'r');title('row boundry');
% figure();plot(sum(IV,1));hold on;plot(sum(IH,1),'r');title('column boundry');

figure();imagesc(IV);axis image;colorbar;title('IV');
figure();imagesc(IH);axis image;colorbar;title('IH');
%% 2. find overlapping area

rowedge=[185,315];
coledge=[120,390];
IVo=IV(rowedge(1):rowedge(2),coledge(1):coledge(2));
IHo=IH(rowedge(1):rowedge(2),coledge(1):coledge(2));
figure();imagesc(IVo);axis image;
figure();imagesc(IHo);axis image;
Overlap2Images(IVo,IHo);
%% generate peak map
if need2reg==1%use IVo and IHo instead of ROI1_1 &ROI2_1 in registrationEstimator
save('IVo.mat','IVo');
save('IHo.mat','IHo');
%run FastFindPeaks2D.m
%%
r=1;
load('FindPeaks_IVo_M.mat');
hv=h;
load('FindPeaks_IHo_M.mat');
hh=h;

%generate binary images for registration
ROI1_1= TargetMap( hv,r,size(IVo));
ROI2_1 = TargetMap(hh,r,size(IHo));
figure();subplot(1,2,1);imagesc(ROI1_1);axis image
subplot(1,2,2);imagesc(ROI2_1);axis image
registrationEstimator;
end
%save the new function registerImages.m
%% register original images using the registerImages.m
IHo_result = registerImages_Lr(IHo,IVo);
IHor=IHo_result.RegisteredImage;
% IHor=imwarp(IHo,movingReg.Transformation,'OutputView',imref2d(movingReg.SpatialRefObj.ImageSize));
% if necessary, iterate the registration process using IVor
% registrationEstimator;
% IVo_result = registerImages2(IVor,IHo);
% IVor=IVo_result.RegisteredImage;
% Norm=1;%normalized image mark
%  cpselect(NormalizeImage(IHo),NormalizeImage(IVo));%control points for mouse data
% tform = fitgeotrans(movingPoints3,fixedPoints3,'affine');
% IHor=imwarp(IHo,movingReg.Transformation,'OutputView',imref2d(movingReg.SpatialRefObj.ImageSize));

figure();subplot(1,2,1);imagesc(IHor);axis image;colorbar;
subplot(1,2,2);imagesc(IVo);axis image;colorbar;
Overlap2Images(IHor,IVo);


%% stitching image

%beads
IHr=IH;
IVr=IV;
if Norm==1
IHr(rowedge(1):rowedge(2),coledge(1):coledge(2))=IHor*max(IHr(:))/max(IHor(:));
elseif Norm==0
IHr(rowedge(1):rowedge(2),coledge(1):coledge(2))=IHor;
end
% IVr=max(max(IV))/max(max(IVr))*IVr;
% IHr=max(max(IH))/max(max(IHr))*IHr;

%mouse
% IHr=IHor;
% IVr=IVo;
figure();imagesc(IVr);axis image;colorbar;title('IVr');
figure();imagesc(IHr);axis image;colorbar;title('IHr');
Overlap2Images(IVr,IHr);
%% overlap images in Fourier plane
% IVr=IV;
% IHr=IHor*(max(max(IH-min(IH(:)))))+min(IH(:));
cutoff=8;%cutoff frequency of the filter, cutoff is smaller, lp band is broader.
border=cutoff*6;
if filter==1
[fw,fhp,flp]=CreateFilter(IVr,cutoff);
IVr = padarray(IVr,[border,border],'replicate','both');
IHr = padarray(IHr,[border,border],'replicate','both');
end
FrecV=fft2(IVr);
FrecH=fft2(IHr);

% FrecV=fftshift(FrecV);
% FrecH=fftshift(FrecH);

% FrecV=FrecV/max(abs(FrecV(:)))*10^5;
% FrecH=FrecH/max(abs(FrecH(:)))*10^5;

%% combine in Fourier domain

Fw=0.50;%Fw=0.50 Wr;
if filter==1
Frec_lp=(FrecV*(1-Fw)+FrecH*Fw).*flp;
Frec_hp=(FrecV*(1-Fw)+FrecH*Fw).*fhp;
else
    Frec_lp=FrecV*(1-Fw)+FrecH*Fw;
Frec_hp=FrecV*(1-Fw)+FrecH*Fw;
end

Frec=Frec_lp+Frec_hp;
Irec2D=ifft2(Frec);

if flatten==1
Irec2D1=flattenfield(Irec2D,1,128,400);%cutedge=128, window=400, Wr;
Irec2D1=flattenfield(Irec2D1,2,128,150);%x direction,1. cutedge=128, window=150, Wr;
elseif flatten==0
    Irec2D1=Irec2D;
else
    Irec2D1=Irec2D./repmat(skyline, [size(Irec2D,1), 1]);
end

% figure();
% imagesc(Irec2D1);
% colorbar;title('sum back to lines');
% axis image;
% figure();
% imagesc(IrecV);
% colorbar;title('IrecV');
% axis image;
%% Plot Fourier domain of H, V, and merged image
kx=linspace(-256/204,256/204,512);
ky=linspace(-256/204,256/204,513);
figure();
subplot(1,3,1);
imagesc(kx,ky,abs(fftshift(fft2(IH))));axis image;
caxis([0 2*10^7]);
xlabel('k_x');
ylabel('k_y');
title('F(H)')
subplot(1,3,2);
imagesc(kx,ky,abs(fftshift(fft2(IV))));axis image;
caxis([0 2*10^7]);
xlabel('k_x');
ylabel('k_y');
title('F(V)')
subplot(1,3,3);
imagesc(kx,ky,abs(fftshift(Frec)));axis image;
caxis([0 2*10^7]);
xlabel('k_x');
ylabel('k_y');
title('F(H+\alphaV)');

saveas(gcf,'manuscript figures/figure4_FourierPlane_Waters_fhp1.fig','fig');
saveas(gcf,'manuscript figures/figure4_FourierPlane_Waters_fhp1.pdf','pdf');
%% save results
Irec2D2=Irec2D1((border+1):(end-border),(border+1):(end-border));
%remove bg from Isum
% Irec2D2=Irec2D2-mean(mean(Irec2D2(420:470,50:100)));
Irec2D3=imrotate(Irec2D2,0,'crop');
Irec2D3(Irec2D3<0)=0;

figure();
imagesc(Irec2D3);
colorbar;title('sum back to lines');
axis image;
%% save tiff image
    imdata = uint16(Irec2D3_Mw/404);%convert to photon count
        t = Tiff('manuscript figures\YFPScarletMouse_directsum_photoncount_Fw5.tif','w');
        t.setTag('ImageLength',size(imdata,1));
        t.setTag('ImageWidth', size(imdata,2));
        t.setTag('Photometric', Tiff.Photometric.MinIsBlack);
        t.setTag('BitsPerSample', 16);
        t.setTag('SamplesPerPixel', size(imdata,3));
%         t.setTag('TileWidth', 128);
%         t.setTag('TileLength', 128);
        t.setTag('Compression', Tiff.Compression.None);
        t.setTag('PlanarConfiguration', Tiff.PlanarConfiguration.Chunky);
        t.setTag('Software', 'MATLAB');
        t.write(imdata);
        t.close();

%% extract signal
% h=h-1;
Iw=Irec2D3_Lw;
Ir=Irec2D3_Lr;
hw=h;
hr=h;

figure();imagesc(Ir);axis image;
hold on;
scatter(hw(:,2),hw(:,1),'ro');
% h=cat(1,[[285:295]',[147:157]'],h);
hr=hr(1:min(size(hr,1),size(hw,1)),:);
hw=hw(1:min(size(hr,1),size(hw,1)),:);

Ir_signal=zeros(1,size(hr,1));
Iw_signal=zeros(1,size(hw,1));

for ii=1:size(hr,1)
Iw_signal(ii)=Iw(hr(ii,1),hr(ii,2));
Ir_signal(ii)=Ir(hr(ii,1),hr(ii,2));
end
 
% ILr(ILr<10)=[];
% ILw(ILw<10)=[];
% ILr(ILr>10^5)=[];

figure();
G=[zeros(size(Ir_signal)),ones(size(Iw_signal))];
boxplot([sqrt(Ir_signal/404),sqrt(Iw_signal/404)],G);
title('SNR comparison Mouse dendrites(photon count)');
grid on;
set(gca,'FontSize',16);
set(gcf,'Position',[100 100 300 500]);
saveas(gcf,'manuscript figures/figure4_SNRcompareMousedendrite.fig','fig');
saveas(gcf,'manuscript figures/figure4_SNRcompareMousedendrite.pdf','pdf');
%% histogram
edges=0:0.5:11;
figure();
h1 = histogram(Iw_signal/404,edges);
h1.FaceColor = [0 0 1];
hold on
h2 = histogram(Ir_signal/404,edges);
h2.FaceColor = [1 0 0];
grid on;
set(gca,'FontSize',16);
legend('IntensityW','IntensityR');
% xlim([-4 12]);
saveas(gcf,'manuscript figures/figure4_hist_LwLr.fig','fig');
saveas(gcf,'manuscript figures/figure4_hist_LwLr.pdf','pdf');

%% hot plot
figure();
imagesc(Ir(24:24+465-1,24:24+465-1)/404);
axis image;
set(gcf,'Position',[100 100 600 600]);
colorbar;
caxis([0 300]);%[0 250]Mouse
colormap('hot');
title('Wr');
saveas(gcf,'manuscript figures/figure3_Wr_fhp1.fig','fig');
saveas(gcf,'manuscript figures/figure3_Wr_fhp1.pdf','pdf');

figure();
imagesc(Iw(24:24+465-1,24:24+465-1)/404);
axis image;
set(gcf,'Position',[100 100 600 600]);
colorbar;
caxis([0 300]);
colormap('hot');
title('Ww');
saveas(gcf,'manuscript figures/figure3_Ww_fhp1.fig','fig');
saveas(gcf,'manuscript figures/figure3_Ww_fhp1.pdf','pdf');

ps=imread('manuscript figures/mouse_pointscan.tif','tif');
ps1=ps(118:(118+560),160:(560+160));
figure();
imagesc(ps1*2);
axis image;
set(gcf,'Position',[100 100 600 600]);
colorbar;
caxis([0 250]);
colormap('hot');
title('Mps');
saveas(gcf,'manuscript figures/figure6_Mps.fig','fig');
saveas(gcf,'manuscript figures/figure6_Mps.pdf','pdf');
%% calculate MSE to evaluate how much information reconstructed
%registrate reconstruction image to w/o image
%registrationEstimator
Lr_reg=imwarp(Lr_norm,movingReg.Transformation,'OutputView',movingReg.SpatialRefObj);
Lw_reg=imwarp(Lw_norm,movingReg.Transformation,'OutputView',movingReg.SpatialRefObj);

Lr_reg=imwarp(Lr_reg,movingReg.Transformation,'OutputView',movingReg.SpatialRefObj);
Lw_reg=imwarp(Lw_reg,movingReg.Transformation,'OutputView',movingReg.SpatialRefObj);

% Overlap2Images(Ir_reg1,Iw);

%normalize image
% Ir_reg_norm=NormalizeImage(Ir_reg1);
% Iw_norm=NormalizeImage(Iw);
% Iw_norm=Iw_norm*mean(mean(Ir_reg_norm))/mean(mean(Iw_norm));
% Iw_norm=Iw*mean(mean(Ir))/mean(mean(Iw));
Lw_norm=Iw/mean(mean(Iw));
Lr_norm=Ir/mean(mean(Ir));
%MSE between Ir_reg1 and Iw_norm
LrMSE = immse(medfilt2(Wr_norm(200:300,100:200)), medfilt2(Lr_reg(200:300,100:200)));
LwMSE = immse(medfilt2(Wr_norm(200:300,100:200)), medfilt2(Lw_reg(200:300,100:200)));
WwMSE = immse(medfilt2(Wr_norm(200:300,200:300)), medfilt2(Ww_norm(200:300,200:300)));
% WwMSE=sqrt(sum(sum((NormalizeImage(Ww_norm(150:350,150:350))-NormalizeImage(Wr_norm(150:350,150:350))).^2,1),2));
save('experimentData3/MSE_beads.mat','LrMSE','LwMSE','WwMSE','WrMSE','Lr_reg','Lw_reg','Ww_norm','Wr_norm');
%% Binary MSE
% figure();imagesc(Lr_norm);
% % p=98;
% % Y = prctile(Lr_norm(:),p);
% Wr_binary=im2bw(Wr_norm,0.9999999999);
% figure();imagesc(Wr_binary);
WwMSEbinary=sqrt(sum(sum((Ww_binary(200:300,200:300)-Wr_binary(200:300,200:300)).^2,1),2));
LrMSEbinary=sqrt(sum(sum((Lr_binary(200:300,100:200)-Wr_binary(200:300,100:200)).^2,1),2));
LwMSEbinary=sqrt(sum(sum((Lw_binary(200:300,100:200)-Wr_binary(200:300,100:200)).^2,1),2));

%% plot MSE
figure();
c = categorical({'Water_{rec}','Water_{sum}'});
bar(c,[WrMSE WwMSE],0.5);
ylim([0 10]);
title('MSE comparison');
grid on;
set(gca,'FontSize',16);
set(gcf,'Position',[100 100 300 500]);
saveas(gcf,'manuscript figures\figure3_MSEcompareWwWr.fig','fig');
saveas(gcf,'manuscript figures\figure3_MSEcompareWwWr.pdf','pdf');

figure();
c = categorical({'Lipid_{rec}','Lipid_{sum}'});
bar(c,[LrMSE LwMSE],0.5);
ylim([0 10]);
title('MSE comparison');
grid on;
set(gca,'FontSize',16);
set(gcf,'Position',[100 100 300 500]);
saveas(gcf,'manuscript figures\figure4_MSEcompareLwLr.fig','fig');
saveas(gcf,'manuscript figures\figure4_MSEcompareLwLr.pdf','pdf');

figure();
% c = categorical({'Lipid_{rec}','Lipid_{sum}'});
bar([WrMSE WwMSE LrMSE LwMSE]);
ylim([0 10]);
% title('MSE comparison');
grid on;
set(gca,'FontSize',16);
set(gcf,'Position',[100 100 200 300]);
saveas(gcf,'manuscript figures\figure4_MSEcompare.fig','fig');
saveas(gcf,'manuscript figures\figure4_MSEcompare.pdf','pdf');
%% calculate SSIM 
[LrSSIM,Lrmap]=ssim(Lr_reg(200:300,200:300),Wr_norm(200:300,200:300));
[LwSSIM,Lwmap]=ssim(Lw_reg(200:300,200:300),Wr_norm(200:300,200:300));
[WwSSIM,Wwmap]=ssim(Ww_norm(200:300,200:300),Wr_norm(200:300,200:300));
[WrSSIM,Wrmap]=ssim(Wr_norm(200:300,200:300),Wr_norm(200:300,200:300));

figure();subplot(1,3,1);imshow(Wwmap);subplot(1,3,2);imshow(Lrmap);
subplot(1,3,3);imshow(Lwmap);

figure();
% c = categorical({'Lipid_{rec}','Lipid_{sum}'});
bar([WrSSIM WwSSIM LrSSIM LwSSIM]);
ylim([0 1]);
% title('SSIM comparison');
grid on;
set(gca,'FontSize',16);
set(gcf,'Position',[100 100 200 300]);
saveas(gcf,'manuscript figures\figure4_SSIMcompare.fig','fig');
saveas(gcf,'manuscript figures\figure4_SSIMcompare.pdf','pdf');

save('experimentData3/SSIM_beads.mat','LrSSIM','LwSSIM','WwSSIM','WrSSIM','Lrmap','Lwmap','Wwmap','Wrmap');

%% binary SSIM
[LrSSIM,Lrmap]=ssim(double(im2bw(Lr_reg(200:300,200:300),0.9999)),double(im2bw(Wr_norm(200:300,200:300),0.9999)));
[LwSSIM,Lwmap]=ssim(double(im2bw(Lw_reg(200:300,200:300),0.9999)),double(im2bw(Wr_norm(200:300,200:300),0.9999)));
[WwSSIM,Wwmap]=ssim(double(im2bw(Ww_norm(200:300,200:300),0.9999)),double(im2bw(Wr_norm(200:300,200:300),0.9999)));
[WrSSIM,Wrmap]=ssim(double(im2bw(Wr_norm(200:300,200:300),0.9999)),double(im2bw(Wr_norm(200:300,200:300),0.9999)));

%% Contrast
%background mask
%add soma coordinate
% h=cat(1,[[280:300]',[142:162]'],h);
r=7;
maskflag='w';
switch maskflag
    case 'w'
    h=hw;
    case 'r'
        h=hr;
end

mask=ones(size(Iw));
for ii=1:size(h,1)
    if h(ii,1)-r >0 && h(ii,2)-r>0 && h(ii,1)+r<size(Ir,1) && h(ii,2)+r<size(Ir,2)
    mask(h(ii,1)-r:h(ii,1)+r,h(ii,2)-r:h(ii,2)+r)=0;
    end
end
mask([1:150,350:end],:)=0;
mask(:,[1:150,350:end])=0;
figure();imagesc(mask);axis image;
[i,j,s]=find(mask);
Ind=sub2ind(size(mask),i,j);
switch maskflag
    case 'w'
    Iw_bg=Iw(Ind); 
    Iw_bg(Iw_bg>10^3)=[];
    BG_Lw=mean(Iw_bg);
    case 'r'
       Ir_bg=Ir(Ind);
       Ir_bg(Ir_bg>10^3)=[];
       BG_Lr=mean(Ir_bg);
end
%%
figure();
boxplot([Ir_signal/BG_Lr,Iw_signal/BG_Lw],[zeros(size(Ir_signal)),ones(size(Iw_signal))]);
title('Contrast comparison');
grid on;
set(gca,'FontSize',16);
set(gcf,'Position',[100 100 300 500]);
saveas(gcf,'manuscript figures\figure4_SBRcompareLwLr.fig','fig');
saveas(gcf,'manuscript figures\figure4_SBRcompareLwLr.pdf','pdf');

%% plot PSF
 Itemp=Irec2D3;
% figure();imagesc(Itemp1_Lr);axis image;
% hold on;
% scatter(h(:,2),h(:,1),'ro');
% % 
[xq,yq]=meshgrid(linspace(1,512,512),linspace(1,513,513));
[Xq,Yq]=meshgrid(linspace(1,512,5120),linspace(1,513,5130));
Itemp1_Lr=interp2(xq,yq,Itemp,Xq,Yq);
save('Itemp1_Lr.mat','Itemp1_Lr');

% for ii=1:size(h,1)
%     row=h(ii,1);col=h(ii,2);
%     psfrow(ii,:)=interp(Itemp(row-4:row+4,col),10);
%     psfcol(ii,:)=interp(Itemp(row,col-4:col+4),10);
%     psfdiag(ii,:)=interp(diag(Itemp(row-4:row+4,col-4:col+4)'),14);
% end
%%
load('experimentData3\FindPeaks_Lr.mat');
Itemp1=Itemp1_Lr;
for ii=1:size(h,1)
    row=h(ii,1);col=h(ii,2);
    psfrow(ii,:)=Itemp1(row-40:row+40,col);
    psfcol(ii,:)=Itemp1(row,col-40:col+40);
    psfdiag(ii,:)=diag(Itemp1(row-40:row+40,col-40:col+40)');
end

[peak, ind]=max(mean(psfrow,1));
psfrow=psfrow(:,1:2*ind-1)/404;
[peak, ind]=max(mean(psfcol,1));
psfcol=psfcol(:,1:2*ind-1)/404;
[peak, ind]=max(mean(psfdiag,1));
psfdiag=psfdiag(:,1:2*ind-1)/404;

figure();
% e1=errorbar(linspace(-size(psfrow,2)*0.02,size(psfrow,2)*0.02,size(psfrow,2)),mean(psfrow,1),std(psfrow,0,1));
x=linspace(-size(psfcol,2)*0.02,size(psfcol,2)*0.02,size(psfcol,2));
yrowneg=mean(psfrow,1)-std(psfrow,0,1);
yrowpos=mean(psfrow,1)+std(psfrow,0,1);
figure();
h1=area(x,[yrowneg',yrowpos'],-50);
h1(1).FaceColor = [1 1 1];
h1(2).FaceColor = [0 0 1];
xlim([-size(psfrow,2)*0.02 size(psfrow,2)*0.02]);
ylim([-2 10]);
saveas(gcf,'manuscript figures\figure4_PSFrowLr_symm.fig','fig');
saveas(gcf,'manuscript figures\figure4_PSFrowLr_symm.pdf','pdf');

figure();
e2=errorbar(linspace(-size(psfcol,2)*0.02,size(psfcol,2)*0.02,size(psfcol,2)),mean(psfcol,1),std(psfcol,0,1),'r');
xlim([-size(psfrow,2)*0.02 size(psfrow,2)*0.02]);
ylim([-2 10]);
grid on;
saveas(gcf,'manuscript figures\figure4_PSFcolLr_symm.fig','fig');
saveas(gcf,'manuscript figures\figure4_PSFcolLr_symm.pdf','pdf');

figure();
errorbar(linspace(-size(psfdiag,2)*0.02*sqrt(2),size(psfdiag,2)*0.02*sqrt(2),size(psfdiag,2)),mean(psfdiag,1),std(psfdiag,0,1),'g');
xlim([-size(psfrow,2)*0.02 size(psfrow,2)*0.02]);
ylim([-2 10]);
saveas(gcf,'manuscript figures\figure4_PSFdiagLr_symm.fig','fig');
saveas(gcf,'manuscript figures\figure4_PSFdiagLr_symm.pdf','pdf');

%% calculate FWHM
FWHMrow = FWHMcalculator( mean(psfrow,1),0.04,1);
FWHMcol = FWHMcalculator( mean(psfcol,1),0.04,1);
FWHMdiag = FWHMcalculator( mean(psfdiag,1),0.04*sqrt(2),1);

save('PSF_Lr_symm.mat','FWHMcol','FWHMrow','FWHMdiag','psfrow','psfcol','psfdiag','Itemp1','xq','yq','Xq','Yq','h');

%% zoom-in version
center=[235,326];ROI=25;
Irzoom=Ir(center(1)-ROI:center(1)+ROI,center(2)-ROI:center(2)+ROI);
Iwzoom=Iw(center(1)-ROI:center(1)+ROI,center(2)-ROI:center(2)+ROI);

figure();
imagesc(Irzoom/404);
axis image;
set(gcf,'Position',[100 100 600 600]);
colorbar;
caxis([0 600]);
colormap('hot');
set(gca,'FontSize',16);
title('WHzoom');
saveas(gcf,'manuscript figures/figureS3_W_V_zoomin.fig','fig');
saveas(gcf,'manuscript figures/figureS3_W_V_zoomin.pdf','pdf');

figure();
imagesc(Iwzoom/404);
axis image;
set(gcf,'Position',[100 100 600 600]);
colorbar;
caxis([0 850]);
colormap('hot');
title('WVzoom');
set(gca,'FontSize',16);
saveas(gcf,'manuscript figures/figureS3_WV_zoomin.fig','fig');
saveas(gcf,'manuscript figures/figureS3_WV_zoomin.pdf','pdf');
%% crosssection
Mrcross=csvread('experimentData3/cross_section_Mps_verticaldendrite.csv',1,1);
Mrcross=Mrcross(1:end);
Mwcross=csvread('experimentData3/cross_section_Mps_horizontaldendrite.csv',1,1);
Mwcross=Mwcross(1:end);
x=linspace(-length(Mrcross)*5/3,length(Mrcross)*5/3,length(Mrcross));
figure();
plot(x,Mwcross*1.2^2*2,'k--','LineWidth',2);
set(gca,'FontSize',16);
grid on;
axis tight;
% hold on;
figure();
plot(x,Mrcross*1.2^2*2,'k--','LineWidth',2);
set(gca,'FontSize',16);
grid on;
axis tight;
% legend('Mw','Mr');
% saveas(gcf,'manuscript figures/figure6_Mps_zoomincross_horizontalDendrite.fig','fig');
% saveas(gcf,'manuscript figures/figure6_Mps_zoomincross_horizontalDendrite.pdf','pdf');
%% background image

figure();
imagesc(Isum/404);
axis image;
set(gcf,'Position',[100 100 600 600]);
colorbar;
caxis([-1 9]);
colormap('hot');
title('BGsum');
% saveas(gcf,'manuscript figures/figureS1_BGsum.fig','fig');
% saveas(gcf,'manuscript figures/figureS1_BGsum.pdf','pdf');

figure();
imagesc(Irec1/404);
axis image;
set(gcf,'Position',[100 100 600 600]);
colorbar;
caxis([-1 9]);
colormap('hot');
title('BGrec');
% saveas(gcf,'manuscript figures/figureS1_BGrec.fig','fig');
% saveas(gcf,'manuscript figures/figureS1_BGrec.pdf','pdf');