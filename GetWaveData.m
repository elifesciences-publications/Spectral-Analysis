function [waveData] = GetWaveData(Wxy,freq)
%GetWaveData Extract some pertinent information from wavelet coefficient
%   matrix and append as fields to the structure variable 'waveData'
% Inputs:
% Wxy - 2D matrix of wavelet coefficients. Each row corresponds to a 
%   particular frequency and each column to a particular time.
% freq - Vector of freq values corresponding to each row of Wxy
% Outputs:
% waveData - Structure variable where each field stores some useful
%   information extracted from Wxy.


if any(isreal(Wxy))
    errordlg('Entries in Wxy must be complex!')
elseif isempty(Wxy)
    disp('Empty coefficient matrix')
    return
end

%% Extracting frequency info
W = abs(Wxy);
ftmat = repmat(freq(:), 1, size(Wxy,2));
% coimat = repmat(1./Wxy.coiLine(:)',length(Wxy.freq), 1);
sumFP = sum(ftmat(:).*W(:));
sumP = sum(W(:));

waveData.freq_mean = round((sum(ftmat(:).*W(:))/sumP)*100)/100;
waveData.freq_std = circ_std(ftmat(W~=0),W(W~=0));

try
    waveData.freq_mostPow = round(freq(find(sum(W,2)==max(sum(W,2))))*100)/100;
catch
    waveData.freq_mostPow = nan;
end

%% Computing pow spectrum
waveData.pow_spec = {mean(abs(W),2)};
blah = waveData.pow_spec{:};
blah(blah==0) = nan;
waveData.pow_spec_log = {log2(blah)};
% [maxtab,~] = peakdet(waveData.powSpec,max(pkDetThr*max(Wxy.powSpec{file,chNum}),0.1));
% 
% if ~isempty(maxtab)
%     pv = find(maxtab(:,2)== max(maxtab(:,2)));
%     pf = round(Wxy.freq(maxtab(pv,1))*100)/100;
% else
%     [maxtab(:,1),maxtab(:,2), pv, pf] = deal(nan);
% end

%% Estimating phase
W = Wxy;
W(abs(W)==0)=[];
nPhaseBins = min([numel(angle(W)), 90]);
[waveData.phase_hist,theta] = hist(angle(W(:)),nPhaseBins); % Unweighted phase histogram
[ph_dist3,vals] = hist3([angle(W(:)) abs(W(:))],[nPhaseBins,nPhaseBins]);
powmat = repmat(vals{2},size(ph_dist3,1),1);
waveData.phase_hist_wt = ph_dist3.*powmat;
waveData.phase_hist_wt = sum(waveData.phase_hist_wt,2)'; % Power-weighted phase histogram
mphase = angle(sum(W));
mphase(mphase<0) = mphase(mphase<0)+ 2*pi; % 0 to 360 instead of -180 to +180 deg
waveData.phase_mean = round(mphase*(180/pi)*100)/100;
waveData.phase_std = round(100*circ_std(angle(W(:)),abs(W(:)))*(180/pi))/100;

if ~isempty(waveData.phase_hist)
    waveData.phase_hist = [waveData.phase_hist(:); waveData.phase_hist(1)]; % Ligates the plot
    waveData.phase_hist = {waveData.phase_hist/max(waveData.phase_hist)};
    
    waveData.phase_hist_wt = [waveData.phase_hist_wt(:); waveData.phase_hist_wt(1)];
    waveData.phase_hist_wt = {waveData.phase_hist_wt/max(waveData.phase_hist_wt)};
    
    waveData.phase_angles = {[theta(:); theta(1)]};
    maxInd = find(waveData.phase_hist_wt{:} == max(waveData.phase_hist_wt{:}));
    
    if ~isempty(maxInd)
        maxInd = maxInd(1);
        waveData.phase_mostPow = round(waveData.phase_angles{:}(maxInd)*180/pi);
    else
        waveData.phase_mostPow = nan;
    end
else
    waveData.phase_mostPow = nan;
end

if waveData.phase_mostPow < 0
    waveData.phase_mostPow = waveData.phase_mostPow + 360;
end

theta(theta<0) = theta(theta<0) + 360;

%% Estimating power
waveData.pow_tot  = round(sum(abs(W)));
waveData.pow_tot_alt = round(sum(abs(W(angle(W) > pi/2 | angle(W) < -pi/2))));
waveData.pow_tot_synch = round(sum(abs(W(angle(W) < pi/2 & angle(W) > -pi/2))));
waveData.pow_mean_alt = round(mean(abs(W(angle(W) > pi/2 | angle(W) < -pi/2))));
waveData.pow_mean_synch = round(mean(abs(W(angle(W) < pi/2 & angle(W) > -pi/2))));
waveData.pow_std_alt = round(std(abs(W(angle(W) > pi/2 | angle(W) < -pi/2))));
waveData.pow_std_synch = round(std(abs(W(angle(W) < pi/2 & angle(W) > -pi/2))));

end

