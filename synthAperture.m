function elementSum = synthAperture(oneTx, Xfocus, numChannels, numElements, pitch, fs, c)

Xf = (0:numChannels-1) .* pitch;                    % mm, 128 transmit beam location
Xf = Xf - mean(Xf);
Xe = (0:numElements-1) .* pitch;                    % mm, 256 receive beam locations 
Xe = Xe - mean(Xe); 
Zf = 40;                                            % mm, focus depth 
Zp = ((0:size(oneTx,1)-1) .* 1/fs .* c/2);          % mm, grid point locations
time = ((0:size(oneTx,1)-1) ./ fs)';                % sec, original time vector 

Xp = Xe;                                            % mm, point locations

Xp = repmat(Xp, [length(Zp) 1]);                    %2353x256 matrix of x locations
Zp = repmat(Zp',[1, numElements]);                  %2353x256 matrix of z locations

%Transmit Delay = [depths x pntlocs] 
Txdist = sqrt((Xp-Xfocus).^2 + (Zp-Zf).^2);         % mm , distance between Tx element and point locations

Zswitch = ones(size(Txdist)); 
Zswitch(Zp<Zf)=-1;                                  % depth correction factor   
Txdelay = (Zswitch.*Txdist)/c + Zf/c; 

%Receive delay = [depths x pntlocs x rxlocs]
for i = 1:numElements  
    Rxdelay(:,:,i) = (1/c).*sqrt((Xe(i)-Xp).^2 + Zp.^2);
end

% %Full Delay at each rx loc need to make sure the points match from txdelay
for i = 1:numElements 
    fullDelay(:,:,i) = repmat(Txdelay(:,i), [1 numElements]) + Rxdelay(:,:,i); 
end 

%add zeros according to where tx is
ind = find(Xfocus == Xf);
if ind == 1 
    rightPad = zeros(size(Zp,1),numChannels);
    dataPad = cat(2,oneTx, rightPad);
elseif ind == 128
    leftPad = zeros(size(Zp,1),numChannels);
    dataPad = cat(2,leftPad,oneTx);
else  
    leftPad = zeros(size(Zp,1), ind-1);
    Ncols = size(leftPad,2);
    rightPad = zeros(size(Zp,1),numChannels);
    rightPad(:,1:Ncols)=[]; 
    dataPad = cat(2,leftPad, oneTx, rightPad);
end 

%apply the delay 
%delay = [depths x pntlocs x rxLoc
for pntLoc = 1:numElements
    for rxLoc = 1:numElements 
        delayUnsum(:,pntLoc,rxLoc) = interp1(time,dataPad(:,pntLoc),fullDelay(:,pntLoc,rxLoc),"spline",0);
    end 
end    

%elementSum is [2353 x 256] - sums across all 256 receive locations , so you should have [Zp x Xp] 
elementSum = sum(delayUnsum,2); 


 