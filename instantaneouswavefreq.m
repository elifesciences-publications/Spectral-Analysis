
function varargout = instantaneouswavefreq(Wxy,freq)
% INSTANTANEOUSWAVEFREQ Computes instantaneous freq from a matrix of
%   wavelet coefficients
% mean_instantaneous_freq = instantaneouswavefreq(Wxy,freq);
% [...,peakpower_instantaneous_frequency] = instantaneouswavefreq(Wxy,freq)
% [~,~,std_instantaneous_freq] = instantaneouswavefreq(Wxy,freq)

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
<<<<<<< HEAD
fmat  = repmat(freq,1,size(Wxy,2)); % Matrix where each col is the 
        % freq vector and # of cols = length(time)
=======

fmat  = repmat(freq,1,size(Wxy,2)); % Matrix where each col is the
% freq vector and # of cols = length(time)
>>>>>>> c66f2a4db027781d7893de3571bf637aec587bff
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

varargout{1} = mfvec;
varargout{2} = pfvec;
<<<<<<< HEAD
varargout{3} = freq_inst_std;
=======

if nargout ==3
    maxF = zeros(1,size(Wxy,2));
    for t = 1:size(Wxy,2)
        gps = sqrt(abs(Wxy(:,t)));
        pks = findpeaks_hht(gps);
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

>>>>>>> c66f2a4db027781d7893de3571bf637aec587bff
