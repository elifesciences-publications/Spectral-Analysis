
function varargout = instantaneouswavefreq(Wxy,freq)
% INSTANTANEOUSWAVEFREQ Computes instantaneous freq from a matrix of
%   wavelet coefficients
% mean_instantaneous_freq = instantaneouswavefreq(Wxy,freq);
% [...,peakpower_instantaneous_frequency] = instantaneouswavefreq(Wxy,freq)
% [~,~,std_instantaneous_freq,maxPow] = instantaneouswavefreq(Wxy,freq)
% Inputs:
% Wxy - Matrix of wavelet coefficients of size [F,T], where F is number of
%   frequency scales and T is # of time points
% freq - Freq vector of length equal to the number of rows of Wxy
% Outputs:
% mean_instantaneous_freq - Vector of length T, where each entry is the
%   power weighted frequency at a given time
% peakpower_instantaneous_freq - Vector of length T, where each element is
%   the freq at max power at an instant in time
% 
% maxPow - Max power value at each time point

if nargin < 2
    errordlg('At least 2 input variables required')
    return
end

if isreal(Wxy)
    %     errordlg('First input variable must be a matrix of complex wavelet coefficients!')
    %     return
elseif numel(freq)==1
    errordlg('2nd input variable must be a vector!')
    return
elseif (size(freq,1)>1 && size(freq,2)>1)
    errordlg('2nd input variable must be a vector!')
    return
end


freq = flipud(sort(freq(:))); % Ensures that 'freq' is a col vec with values in descending order
fmat  = repmat(freq,1,size(Wxy,2)); % Matrix where each col is the 
        % freq vector and # of cols = length(time)

Wxy_abs = abs(Wxy);
fmat(Wxy_abs==0)=0;
tvpow =sum(Wxy_abs); % Vector of length = length(time), where each
% entry is the total power at all frequencies at each time point
Wxy_tvpow = repmat(tvpow,size(Wxy,1),1);
Wxy_tvpow_prob = Wxy_abs./Wxy_tvpow;

mfvec = sum(fmat.*Wxy_tvpow_prob);
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


[maxPow,maxPowInds] = max(abs(Wxy),[],1);
maxPowFreq = SparseIndex(fmat, maxPowInds,1:size(fmat,2));

varargout{1} = mfvec;
varargout{2} = maxPowFreq;
varargout{3} = freq_inst_std;
varargout{4} = maxPow;


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


% varargout{3} = freq_inst_std;
