
function [varargout] = ComputeMultibandXWT(varargin)
% ComputeMultibandXWT - For a given time seris, computes & displays the XWT and a variant of it
% where only frequency peaks at each time point are kept.
% 
% [Wxy,Wxy_fPeaks,freq,coi,sig95] = ComputeMultibandXWT(signals,timeVec,freqRange,dj,stringency,phaseType);
% [Wxy,Wxy_fPeaks,freq,coi,sig95] = ComputeMultibandXWT(signals,samplingInt,freqRange,dj,stringency,phaseType);

if nargin < 2
    errordlg('At least a minimum of 2 inputs reqd'); 
elseif nargin < 3
     if numel(varargin{2} < 2)
        dt = varargin{2};
        timeVec = 0:dt:length(varargin{1})-1;
    else
        dt = mode(diff(varargin{2}));
        timeVec = varargin{2};
     end
 freqRange = [];
 freqRange(1) = 1/length(varargin{1});
 freqRange(2) = 0.5*dt;
 dj = 1/2^5;
 stringency = 1;
 phaseType = 'all';
elseif nargin < 4
    dj = 1/2^5;
    stringency = 1;
    phaseType = 'all';
elseif nargin < 5
      stringency = 1;
    phaseType = 'all';
elseif nargin < 6   
    phaseType = 'all';
elseif nargin > 6
    errordlg('Too many inputs!');
else
    freqRange = varargin{3};
    if numel(freqRange)~=2
        errordlg('Frequency Range input must contain 2 values!')
    end
    dj = varargin{4};
    stringency = varargin{5};
    phaseType = varargin{6};
end
signal = varargin{1};

if ndims(signal)> 2
    errordlg('Input signal size cannot exceed 2 dimensions!');
elseif size(signal,2) > 2
    errordlg('Signal input must be a matrix with no more than 2 cols, with each col being a different timeseries!')
elseif any(size(signal)==1)
   [Wxy,freq,coi, sig95]  = ComputeXWT(signal(:),signal(:),timeVec,freqRange,dj,stringency,phaseType);
 else
    [Wxy,freq,coi, sig95]  = ComputeXWT(signal(:,1),signal(:,2),timeVec,freqRange,dj,stringency,phaseType);
end

R = zeros(size(Wxy));
B = abs(Wxy);
for tt = 1:size(B,2)
 blah = B(:,tt);
 dBlah = diff(blah);
 peaks = dBlah(1:end-1)>0 & dBlah(2:end)<=0;
  peakInds = find(peaks);
 R(peakInds,tt) = B(peakInds,tt);
end

varargout{1} = Wxy;
varargout{2} = R;
varargout{3} = freq;
varargout{4} = coi;
varargout{5} = sig95;
