
function [varargout] = ComputeMultiBandXWT(varargin)
% ComputeMultiBandXWT - For a given time seris, computes & displays the XWT and a variant of it
% where only frequency peaks at each time point are kept.
%
% Wxy = ComputeMultiBandXWT(signals,timeVec);
% [Wxy,Wxy_fPeaks,freq,coi,sig95] = ComputeMultiBandXWT(signals,dt);
% [....] = ComputeMultiBandXWT(signals,timeVec,'freqRange', freqRange, 'dj', dj, 
%           'timeRange', timeRange, 'phaseType', phaseType,'stringency', stringency,
%           'plotsSwitch', plotSwitch)
% Inputs:
% plotSwitch = 1; will result in plotting of the Wxy
% plotSwitch = 2; will result in plotting of Wxy_fPeaks,
% plotSwitch = 3; will result in plotting both
% 
% Avinash Pujala, JRC/HHMI, 2016

freqRange = [];
timeRange = [];
dj = 1/24;
stringency  = 1;
phaseType = 'all';
plotSwitch = 3;

if numel(varargin)<2
    errordlg('Minimum 2 inputs reqd');
else
    signal = varargin{1};
    if any(size(signal)==1)
        signal = signal(:);
    elseif size(signal,1)==2
        signal = signal'; % So as to make arrange the different signals along columns
    end
    if numel(varargin{2})==2
        dt = varargin{2};
        timeVec = (0:size(signal,1)-1)*dt;        
        
    else
        timeVec = varargin{2};       
        dt = timeVec(2)-timeVec(1);
    end
end


for jj = 1:numel(varargin)
    if ischar(varargin{jj})
        switch lower(varargin{jj})
            case 'freqrange'
                freqRange = varargin{jj+1};
            case 'timerange'
                timeRange = varargin{jj+1};
            case 'dj'
                dj = varargin{jj+1};
            case 'stringency'
                stringency = varargin{jj+1};
            case 'phasetype'
                phaseType = varargin{jj+1};
                if sum(strcmpi(phaseType,{'alt','sync','all'}))==0
                    errordlg('phaseType input must be "alt", "sync", or "both"!')
                end
            case 'plotswitch'
                plotSwitch = varargin{jj+1};
                if isempty(intersect(plotSwitch,[1 2 3]))
                    errordlg('plotSwitch input can only be 1, 2, or 3');
                end
        end
    end
end
if ~isempty(timeRange)
    inds = find(timeVec>= timeRange(1) & timeVec<=timeRange(2));
    timeVec = timeVec(inds);
    signal = signal(inds,:);
end

maxF = 1/(0.5*dt);
minF = 1/size(signal,1);
if isempty(freqRange)
    freqRange(1) = minF;
    freqRange(2) = maxF;
end


if ndims(signal)> 2
    errordlg('Input signal size cannot exceed 2 dimensions!');
elseif size(signal,2) > 2
    errordlg('Signal input must be a matrix with no more than 2 cols, with each col being a different timeseries!')
elseif any(size(signal)==1)
    [Wxy,freq,coi, sig95]  = ComputeXWT(signal,signal,timeVec,freqRange,dj,stringency,phaseType);
else
    [Wxy,freq,coi, sig95]  = ComputeXWT(signal(:,1),signal(:,2),timeVec,freqRange,dj,stringency,phaseType);
end

% s_norm = signal/max(signal);
% s_norm = log2(signal);
% sigMat = repmat(signal(:)',size(Wxy,1),1);
R = zeros(size(Wxy));
% B = log2(abs(Wxy)).*sigMat;
% B = log2(abs(Wxy));
B = abs(Wxy);
Wxy = abs(Wxy).^2;
for tt = 1:size(B,2)
    blah = B(:,tt);
    blah = SubtractMinimalEnvelope(blah);
    dBlah = diff(blah);
    peaks = dBlah(1:end-1)>0 & dBlah(2:end)<=0;
    peakInds = find(peaks);
    blah = blah./max(blah);
    peakInds(blah(peakInds)<0.1)=[];
    R(peakInds,tt) = B(peakInds,tt);
end

varargout{1} = Wxy;
varargout{2} = R;
varargout{3} = freq;
varargout{4} = coi;
varargout{5} = sig95;

switch plotSwitch
    case 1
        figure('Name','Wavelet Spectrogram')
        subaxis(2,1,1),colormap(hot), hold on
        imagesc(timeVec,log2(freq),log2(abs(Wxy)))
        cLim = get(gca,'clim');
        dC = diff(cLim);
        set(gca,'clim',[cLim(1)+ 0*dC cLim(1) + 1*dC], 'xtick',[],'ydir','normal');
        xlim([timeVec(1) timeVec(end)])
        ylim(log2([freqRange]))
        set(gca,'yticklabel',num2str(round(2.^(str2num(get(gca,'yticklabel'))))))
        box off
        for ch = 1:size(signal,2)
            if ch == 2
                col = 'r';
            else
                col = 'b';
            end
            subaxis(2,1,2), plot(timeVec,signal(:,1),'color',col)
            hold on
            set(gca,'color','k','tickdir','out'), box off
        end
        xlim([timeVec(1) timeVec(end)])
        ylim([-inf inf])
    case 2
        figure('Name','Wavelet Spectrogram_peaks only')
        subaxis(2,1,1),colormap(hot), hold on
        imagesc(timeVec,log2(freq),R)
        cLim = get(gca,'clim');
        dC = diff(cLim);
        set(gca,'clim',[cLim(1)+ 0*dC cLim(1) + 1*dC], 'xtick',[],'ydir','normal');
        xlim([timeVec(1) timeVec(end)])
        ylim(log2([freqRange]))
        set(gca,'yticklabel',num2str(round(2.^(str2num(get(gca,'yticklabel'))))))
        box off
        for ch = 1:size(signal,2)
            if ch == 2
                col = 'r';
            else
                col = 'b';
            end
            subaxis(2,1,2), plot(timeVec,signal(:,1),'color',col)
            hold on
            set(gca,'color','k','tickdir','out'), box off
        end
        xlim([timeVec(1) timeVec(end)])
        ylim([-inf inf])
    case 3
      figure('Name','Wavelet Spectrogram')
        subaxis(2,1,1),colormap(hot), hold on
%         imagesc(timeVec,log2(freq),log2(abs(Wxy)))
        imagesc(timeVec,log2(freq),abs(Wxy.^0.5))
        cLim = get(gca,'clim');
        dC = diff(cLim);
        set(gca,'clim',[cLim(1)+ 0*dC cLim(1) + 1*dC], 'xtick',[],'ydir','normal');
        xlim([timeVec(1) timeVec(end)])
        ylim(log2([freqRange]))
        yTick = get(gca,'ytick');
        set(gca,'ytick',yTick(1):0.5:yTick(end));
        set(gca,'yticklabel',num2str(round(2.^(str2num(get(gca,'yticklabel'))))))
        box off
        for ch = 1:size(signal,2)
            if ch == 2
                col = 'r';
            else
                col = 'b';
            end
            subaxis(2,1,2), plot(timeVec,signal(:,1),'color',col)
            hold on
            set(gca,'color','k','tickdir','out'), box off
        end
        xlim([timeVec(1) timeVec(end)])
        ylim([-inf inf])  
        
        figure('Name','Wavelet Spectrogram_peaks only')
        subaxis(2,1,1),colormap(hot), hold on
        imagesc(timeVec,log2(freq), abs(R.^0.5))
        cLim = get(gca,'clim');
        dC = diff(cLim);
        set(gca,'clim',[cLim(1)+ 0*dC cLim(1) + 1*dC], 'xtick',[],'ydir','normal');
        xlim([timeVec(1) timeVec(end)])
        ylim(log2([freqRange]))
        yTick = get(gca,'ytick');
        set(gca,'ytick',yTick(1):0.5:yTick(end));
        set(gca,'yticklabel',num2str(round(2.^(str2num(get(gca,'yticklabel'))))))
        box off
        for ch = 1:size(signal,2)
            if ch == 2
                col = 'r';
            else
                col = 'b';
            end
            subaxis(2,1,2), plot(timeVec,signal(:,1),'color',col)
            hold on
            set(gca,'color','k','tickdir','out'), box off
        end
        xlim([timeVec(1) timeVec(end)])
        ylim([-inf inf])
end


function s = getspline(x)
N = length(x);
p = findpeaks_hht(x);
s = spline([0 p N+1],[0 x(p) 0],1:N);