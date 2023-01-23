function [delayedData, delayTimes] = dynamicReceive(data, numChannels, pitch, fs, c)

% Dynamic Receive/Continuous delay and sum beamforming for linear array 
% INPUTS: 
% - data: 3D RF data where the dimensions are = 
% [fast time x number of transmit channels x number of receive beams] 
% - numChannels: number of Transmit
% channels - pitch: distance between each channel (mm) 
% - fs: sampling frequency (Hz) 
% - c: sound speed (m/s) 
% OUTPUTS: 
% - delayedData: original with delayTimes applied 
% - delayTimes: the delay times for each channel for the focusDist 

c = c * 1000;                                 % mm/s - speed of sound in tissue                                 % 
Zf = ((0:size(data,1)-1) .* 1/fs .* c/2);     % mm - depth in tissue 
Xf =  (0:numChannels-1) .* pitch; 
Xf = Xf - mean(Xf);                           % mm - lateral channel distance from center

% mm - distance calculation for each channel and depth 
for i = 1:length(Xf)
    for j = 1:length(Zf)
        dist(j,i) = sqrt(Zf(j).^2 + Xf(i).^2); 
    end 
end 

timeofflight = dist./c;                        % seconds 

for i = 1:numChannels
    for j = 1:length(Zf)
        delayTimes(j,i) = timeofflight(j,i) + Zf(j)/c;  % [depth , channels] = delayDR
    end
end

time = ((0:size(data,1)-1) ./ fs)';                % seconds 

for i = 1:numChannels 
    for j = 1:size(data,3)
       
        delayUnsum(:,i,j) = interp1(time, data(:,i,j), delayTimes(:,i),"spline",0);
     
    end
end
 
% sum across the channels  
delayedData = sum(delayUnsum,2);

% Plotting
hilbertTrans = 20*log10(abs(hilbert(delayedData)));
figure; imagesc(Xf, Zf, squeeze(hilbertTrans));
title("Continuous Delay and Sum"); axis image;
xlabel("Distance (mm)"); ylabel("Depth in Tissue (mm)"); colormap("gray");
