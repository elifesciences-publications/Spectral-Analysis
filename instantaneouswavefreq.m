
function varargout = instantaneouswavefreq(Wxy,freq)
% INSTANTANEOUSWAVEFREQ Computes instantaneous freq from a matrix of
%   wavelet coefficients
% mean_instantaneous_freq = instantaneouswavefreq(Wxy,freq);
% [...,peakpower_instantaneous_frequency] = instantaneouswavefreq(Wxy,freq)
% [~,~,std_instantaneous_freq,maxPow, freqPks] = instantaneouswavefreq(Wxy,freq)
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
% freqPks - A matrix of size(Wxy). Each non-zero entry for a given time
%   point (# of cols of Wxy) corresponds to the peaks of the local power
%   spectrum at that time point.
% 
% Avinash Pujala, HHMI, 2016

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



freq = sort(freq(:),'descend'); % Ensures that 'freq' is a col vec with values in descending order
fmat  = repmat(freq,1,size(Wxy,2)); % Matrix where each col is the
% freq vector and # of cols = length(time)

Wxy_abs = abs(Wxy);
% fmat(Wxy_abs==0)=0;
% entry is the total power at all frequencies at each time point
Wxy_tvpow = repmat(sum(Wxy_abs,1),size(Wxy,1),1);
Wxy_tvpow_prob = Wxy_abs./Wxy_tvpow;

mfvec = sum(fmat.*Wxy_tvpow_prob,1);
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
% varargout{2} = pfvec;
varargout{3} = freq_inst_std;
varargout{4} = maxPow;

if nargout > 4
    disp('Getting frequency peaks...')
    freqPks = zeros(size(Wxy));
    meanFreqPks = zeros(1,size(Wxy,2));
    dispChunk = round(size(Wxy,2)*0.2);
    for t = 1:size(Wxy,2)
        if mod(t,dispChunk)==0
            disp([num2str(round(t*100/size(Wxy,2))) ' %'])
        end
        gps = abs(Wxy(:,t));
        %         gps = sqrt(abs(Wxy(:,t)));
%         pks = GetPks(gps, 'polarity',1, 'peakThr',mean(gps));
        pks = GetPks(gps, 'polarity',1);
        %         pks(gps(pks)<20) = [];
        if ~isempty(pks)
            %             maxFInd = find(freq == max(freq(pks)));
            %             maxF(pks,t) = freq(maxFInd);
            freqPks(pks,t) = gps(pks);        
        else
            freqPks(:,t) = 0;  
        end
    end
    disp('Done!')
    varargout{5} = freqPks;
end

end

