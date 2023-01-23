function [delayedData, delayTimes] = parallelReceive(data, focusDist, numChannels, pitch, fs, c, N)

% Parallel Receive Delay and Sum Beamforming for Linear Array 
% INPUTS: 
% - data: 3D RF data where the dimensions are = 
% [fast time x number of transmit channels x number of receive beams] 
% - focusDist: focusing depth (mm) 
% - numChannels: number of Transmit
% channels - pitch: distance between each channel (mm) 
% - fs: sampling frequency (Hz) 
% - c: sound speed (m/s) 
% - N: number of receive beams per transmit beam (even number) 
% OUTPUTS: 
% - delayedData: original with delayTimes applied 
% - delayTimes: the delay times for each channel for the focusDist 

c = c * 1000;                                 % mm/s - speed of sound in tissue                                 % 
Zf = ((0:size(data,1)-1) .* 1/fs .* c/2);     % mm - depth in tissue 
Xf =  (0:numChannels-1) .* pitch; 
Xf = Xf - mean(Xf);                           % mm - lateral channel distance from center

time = ((0:size(data,1)-1) ./ fs)';           % seconds 

dataN = data(:,:,1:N:end);

% delay = [depth x channels x N]
for i = 1:numChannels           % delay for each depth, channel for each transmit beam location 
    for j = 1:length(Zf)
        for k = 1:N
            delay(j,i,k) = (1/c) * sqrt((Xf(i) - pitch*-(2*k-1-N))^2 + Zf(j)^2) + Zf(j)/c;        
        end
    end 
end

delayTimes = delay;

% apply the delay 
for i = 1:numChannels
    for n = 1:N
        flippedN = N:-1:1; 
        delayUnsum(:,i,n,:) = interp1(time, dataN(:,i,:), delay(:,i,flippedN(n)),"spline",0);
   
    end 
end 
% restructure to 3-D
delayUnsum = delayUnsum(:,:,:);
 
% sum across the channels  
delayedData = sum(delayUnsum,2);

% Plotting
hilbertTrans = 20*log10(abs(hilbert(delayedData)));
figure; imagesc(Xf, Zf, squeeze(hilbertTrans));
title("Parallel Receive"); axis image;
xlabel("Distance (mm)"); ylabel("Depth in Tissue (mm)"); colormap("gray");
