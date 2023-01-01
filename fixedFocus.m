function [delayedData, delayTimes] = fixedFocus(data, focusDist, numChannels, pitch, fs, c)

% Fixed focus delay and sum beamforming 
% Inputs: - data: 3D RF data where the dimensions are = 
% [fast time x number of transmit channels x number of receive beams] 
% - focusDist: focusing depth (mm) - numChannels: number of Transmit
% channels - pitch: distance between each channel (mm) - fs: sampling
% frequency (Hz) - c: sound speed (m/s) 
% Outputs: - delayedData: original with delayTimes applied - delayTimes:
% the delay times for each channel for the focusDist 

c = c * 1000;                                 % mm/s - speed of sound in tissue                                 % 
Zf = ((0:size(data,1)-1) .* 1/fs .* c/2)      % mm - depth in tissue 
Xf =  (0:numChannels-1) .* pitch; 
Xf = Xf - mean(Xf)                            % mm - lateral channel distance from center

for i = 1:length(Xf)
    dist(i) = sqrt(focusDist^2 + Xf(i).^2);   % mm - distance equation
end  

timeofflight = dist./c;                       % sec - time of flight 

for i = 1:length(timeofflight)
    delayTimes(i) = timeofflight(i) + focusDist/c;   % sec -  delay profile
end

delayTimes = delayTimes - min(delayTimes);    

% Before applying delays you have to upsample some in order to make the
% data dimensions the same once the delays are applied 

upSamp = 8;                                     % up sampling factor 
z1 = linspace(1,length(Zf),length(Zf));
z2 = linspace(1,length(Zf), length(Zf)*upSamp); 
dataUp = interp1(z1,data,z2,'linear');          % up samples data 
 
% convert delayTimes to pixels 
pixelDel = round(delayTimes .* fs .* upSamp) + 1;   
delayUnsummed = zeros(size(dataUp));            % allocate data space

% apply delays to data 
for i = 1:numChannels
  delayUnsummed(1:end-pixelDel(i)+1,i,:) = dataUp(pixelDel(i):end,i,:);
end 

% sum across the channels 
delayedData = sum(delayUnsummed,2);            % sum across each channel

% Plotting
hilbertTrans = 20*log10(abs(hilbert(delayedData)));
figure; imagesc(Xf, Zf, squeeze(hilbertTrans),[-100 100]);
title("Fixed Focus Delay and Sum"); axis image;
xlabel("Lateral Distance (mm)"); ylabel("Depth in Tissue (mm)"); colormap("gray");

