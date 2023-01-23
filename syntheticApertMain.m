%% Synthetic Aperture - Main File 

%This code is a demo. It works under the assumption you've provided your
% own data 

clear all; close all; clc;

%% Load in Data - Input variables

load("data.mat");           %raw data file where dimensions are [fast time x number of channels x number of beams]  
t0 = veraStrct.timeZero;    %start time 
data = data(t0:end,:,:);    
fs = 20*10^6;               %[Hz] - sampling rate
numChannels = 128;          % number of active channels
numElements = 256;          % total number of elements                      
pitch = 0.20;               %[mm] - element spacing  
c = 1540;                   %[m/s] - sound speed 

%% Run synthAperture

Xf = (0:numChannels-1) .* pitch;                        % transmit beam location
Xf = Xf - mean(Xf);

for beamLoc = 1:numChannels
   
    oneTx = data(:,:,beamLoc); 
    Xfocus = Xf(beamLoc); 
    
    %calls synthAperture function
    elementSum = synthAperture(oneTx, Xfocus, numChannels, numElements, pitch, fs, c);
    %sum synApt across the pnts 
    %dataAllTx = [depths, rxLoc, beamLoc] 
    dataAllTx(:,:,beamLoc) = elementSum;
end 

%% Make mask from a single transmit 

%make axis'
xAxis = (0:numElements-1) .* pitch; 
xAxis = xAxis - mean(xAxis);
zAxis = ((0:size(data,1)-1) .* 1/fs .* c/2);           

bModeImage = 20*log10(abs(hilbert(dataAllTx(:,:,64))));

figure(1); imagesc(xAxis , zAxis, bModeImage);
title("Beamformed Image from all Transmits without Mask"); 
xlabel("Distance from Center (mm) "); ylabel("Depth in Tissue (mm)"); 
colormap("gray"); axis image;

beamMask = roipoly; 

figure; 
imagesc(xAxis, zAxis, beamMask); axis image;
title("Beam Mask (Drawn)"); colormap("gray");

%% Get Masks for All Beams

[x1, y1] = find(beamMask(round(1:end/2),:),1,'first');
[x2, y2] = find(beamMask(round(1:end/2),:),1,'last');

%eliminating 0's at the top 
beamMaskMod = beamMask;
beamMaskMod(1:x1,y1:y2) = 1;

for i = 1:numChannels
    beamMaskPad = zeros(size(beamMaskMod));
   
    if i <=65 %shift left
        shift = 64-(i-1); 
        beamMaskPad(:,1:end-shift) = beamMaskMod(:,shift+1:end); 
        shift = 1; 
    elseif i>=66 %shift right 
        beamMaskPad(:,shift+1:end) = beamMaskMod(:,1:end-shift);
        shift = shift + 1;
    end
    allMasks(:,:,i) = beamMaskPad; 
end

%% graphing for check 

bModeImage = 20*log10(abs(hilbert(dataAllTx(:,:,128))));

figure; subplot(1,2,1); imagesc(xAxis, zAxis, bModeImage);
title("Beamformed Image from Tx=128"); 
xlabel("Distance from Center (mm)"); ylabel("Depth in Tissue (mm)"); 
colormap("gray"); axis image; 

subplot(1,2,2);
imagesc(xAxis, zAxis, allMasks(:,:,128)); axis image;
title("Beam Mask for Tx=128"); colormap("gray");

%% multiply each beam by correct mask 

for i = 1:numChannels 
    maskedAllTx(:,:,i) = dataAllTx(:,:,i).*allMasks(:,:,i); 
end 
sumAllTransmits = sum(maskedAllTx,3);
bModeImage = 20*log10(abs(hilbert(sumAllTransmits)));

figure; imagesc(xAxis , zAxis, bModeImage, [15 120]);
title("Synthetic Aperture from all Transmits with Mask"); 
xlabel("Distance from Center (mm) "); ylabel("Depth in Tissue (mm)"); 
colormap("gray"); axis image;



