
function varargout = GetWaveFreqTimeseries(Wxy,freq)
% GetWaveFreqTimeseries Given Wxy, the matrix of cross wavelet
%    coefficients, returns the timeseries of mean frequencies, max-power
%    frequencies and the power and phase at these values
% freq_mean = GetWaveFreqTimeseries(Wxy,freq);
% [, freq_peakpower] = GetWaveFreqTimeseries(Wxy,freq)
% Inputs:
% Wxy - 2D matrix of complex crosswavelet coefficients.
% freq - Freq vector with values corresponding to each of the rows of Wxy
% Outputs:
% freq_mean - Structure variable with the following fields
%   .freq -  A length T vector where T is the total # of time points. Each
%       element corresponds to mean frequency (weighted by XW power) at a
%       given time point
%   .pow - A vector of XW power values at each of the frequencies in the
%       .freq vector
%   .phase - A vector of phase values
%   .std - Standard deviation of frequency at each time point
% freq_power - Structure variable with same fields as freq_mean, but .freq
%   contains freq with most power at each time point
% 
% Avinash Pujala, JRC/HHMI, 2016

if nargin < 2
    errordlg('At least 2 input variables required')
    return
end

if isreal(Wxy)
%         errordlg('First input variable must be a matrix of complex wavelet coefficients!') 
elseif numel(freq)==1
    errordlg('2nd input variable must be a vector!')    
elseif (size(freq,1)>1 && size(freq,2)>1)
    errordlg('2nd input variable must be a vector!')    
end


freq = flipud(sort(freq(:))); % Ensures that 'freq' is a col vec with values in descending order
fmat  = repmat(freq,1,size(Wxy,2)); % Matrix where each col is the 
        % freq vector and # of cols = length(time)

Wxy_abs = abs(Wxy);
fmat2 = fmat;
fmat(Wxy_abs==0)=0;
tvpow =sum(Wxy_abs,1); % Vector of length = length(time), where each
% entry is the total power at all frequencies at each time point
Wxy_tvpow = repmat(tvpow,size(Wxy,1),1);
Wxy_tvpow_prob = Wxy_abs./Wxy_tvpow;

mfvec = sum((fmat.*Wxy_tvpow_prob),1);
mfvec(isnan(mfvec))=0;

[~,freq_inst_std] = WeightedStats(fmat,Wxy_tvpow_prob);

maxmat = max(Wxy_abs);
maxmat = repmat(maxmat,size(Wxy,1),1);
diffmat = Wxy_abs-maxmat;
diffmat(diffmat==0)=1;
diffmat(diffmat>1)=0;
diffmat(diffmat<0)=0;
pfvec = sum(diffmat.*fmat);
pfvec(isnan(pfvec))=0;

% [~,maxInds] = max(Wxy_abs,[],1);
% pfvec2 = SparseIndex(fmat2,maxInds,1:size(fmat,2));

freq_mean = struct;
freq_pow = struct;
for tt = 1:length(mfvec)
    [~, ind] = min(abs(freq - mfvec(tt)));
    freq_mean.freq(tt) = mfvec(tt);
    freq_mean.pow(tt) = Wxy_abs(ind,tt);
    freq_mean.phase(tt) = angle(Wxy(ind,tt))*180/pi;
    freq_mean.std(tt) = freq_inst_std(tt);
    
    [~, ind] = min(abs(freq - pfvec(tt)));
    freq_pow.freq(tt) = pfvec(tt);
    freq_pow.pow(tt) = Wxy_abs(ind,tt);
    freq_pow.phase(tt) = angle(Wxy(ind,tt))*180/pi;
    freq_pow.std(tt) = freq_inst_std(tt);
end




varargout{1} = freq_mean;
varargout{2} = freq_pow;


if nargout ==3
    maxF = zeros(1,size(Wxy,2));
    for t = 1:size(Wxy,2)
        gps = sqrt(abs(Wxy(:,t)));
        pks = GetPks(gps);
        pks(gps(pks)<20) = [];
        if ~isempty(pks)         
            maxFInd = find(freq == max(freq(pks)));
            maxF(t) = freq(maxFInd);
        else
            maxF(t) = 0;
        end
    end
varargout{3} = maxF;    
end


