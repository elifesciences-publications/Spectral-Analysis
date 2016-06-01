function pks = GetPks(x,varargin)
% GetPks - Given a signal x, returns all peaks in the signal
% pks = GetPks(x)
% pks = GetPks(x,'polarity',polarity,'peakThr',peakThr,'thrType',
%   thrType, 'minPkDist',minPkDist);
% Inputs:
% x - Signal in which to find peaks
% 'polarity' - 0, 1, or -1;
%       0 - Returns maxima and minima
%       1 - Returns maxima only
%      -1 - Returns minima only
% 'peakThr' - Amplitude threshold over which to return peaks
% 'thrType' - 'abs' or 'rel'; 'abs' results in returning of peaks over an
%   absolute value, 'rel' returns in returning of peaks whose amplitude is
%   determined relative to maxima or minima flanking the peak.
% 'minPkDist' - Starting from the largest peak, returns only those peaks
%   that are separated by a minimum distance.
%
% Avinash Pujala, Koyama lab/HHMI, 2016

peakThr = 0;
minPkDist = 0;
thrType = 'abs';
pol = 1;

%# Parse inputs
for jj = 1:numel(varargin)
    switch lower(varargin{jj})
        case 'polarity'
            pol = varargin{jj+1};
        case 'peakthr'
            peakThr = varargin{jj+1};
        case 'minpkdist'
            minPkDist = varargin{jj+1};
        case 'thrtype'
            thrType = varargin{jj+1};
    end
end

%# Pks by polarity
if pol == 1
    pks = GetCrests(x);
elseif pol == -1
    pks = GetCrests(-x);
elseif pol ==0
    crests = GetCrests(x);
    troughs = GetCrests(-x);
    pks = union(crests,troughs);
end

%# Pks by threshold
if peakThr == 0 && strcmpi(thrType,'abs')
    pks(abs(x(pks)) < peakThr) =[];
elseif strcmpi(thrType,'abs')
     pks(abs(x(pks)) < peakThr) =[];
elseif peakThr ~=0 && strcmpi(thrType,'rel')   
    x = x(:) -((GetMaximalEnvelope(x)- GetMaximalEnvelope(-x))/2)';   % Not time-efficient, peak computation redundant, must change later
    pks(abs(x(pks)) < peakThr)=[];
end

%# Pks by min distance
[~, inds] = sort(x(pks),'descend');
pks_sort = pks(inds);
dPks = diff(pks_sort);
pks_sort(find(abs(dPks) < minPkDist)+1)=[];
pks = unique(pks_sort);


end

function crests = GetCrests(x)
dx = diff(x);
crests = find(dx(1:end-1)>0 & dx(2:end)<=0)+1;
% crests    = find(diff(diff(x) > 0) < 0);
% u    = find(x(crests+1) > x(crests));
% crests(u) = crests(u)+1;

end
