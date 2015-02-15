
function [varargout] = ComputeMultiBandXWT(varargin);
% ComputeMultiBandXWT - For a given time seris, computes & displays the XWT and a variant of it
% where only frequency peaks at each time point are kept.
%
% [Wxy,Wxy_fPeaks,freq,coi,sig95] = ComputeMultiBandXWT(signals,timeVec,freqRange,dj,stringency,phaseType,plotSwitch);
% [Wxy,Wxy_fPeaks,freq,coi,sig95] = ComputeMultiBandXWT(signals,samplingInt,freqRange,dj,stringency,phaseType,plotSwitch);
% Inputs:
% plotSwitch = 1; will result in plotting of the Wxy
% plotSwitch = 2; will result in plotting of Wxy_fPeaks,
% plotSwitch = 3; will resulti in plotting both


if nargin < 2
    errordlg('At least a minimum of 2 inputs reqd');
    return;
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
    plotSwitch = 0;
elseif nargin < 4
    dj = 1/2^5;
    stringency = 1;
    phaseType = 'all';
    plotSwitch = 0;
elseif nargin < 5
    stringency = 1;
    phaseType = 'all';
    plotSwitch = 0;
elseif nargin < 6
    phaseType = 'all';
    plotSwitch = 0;
elseif nargin > 7
    errordlg('Too many inputs!');
else
    freqRange = varargin{3};
    if numel(freqRange)~=2
        errordlg('Frequency Range input must contain 2 values!')
    end
    dj = varargin{4};
    stringency = varargin{5};
    phaseType = varargin{6};
    plotSwitch = varargin{7};
end
signal = varargin{1};
if any(size(signal))==1
    signal = signal(:);
elseif size(signal,1)==2
    signal = signal'; % So as to make arrange the different signals along columns
end
% s1 = getspline(signal);
% s2 = getspline(-signal);
% s = ((s1+s2)/2)';

timeVec = varargin{2};
% f_max = findpeaks_hht(signal);
% f_min = findpeaks_hht(signal);
% f_all = union(f_min(:), f_max(:));
% s1 = spline(timeVec(f_max),signal(f_max),timeVec);
% s2 = spline(timeVec(f_min),signal(f_min),timeVec);  % Spline-method does
% some weird stuff at high frequencies
% s = (s1+s2)/2;
% s = interp1(timeVec(f_all),signal(f_all), timeVec);

% s = zeros(size(signal));
% for cc = 1:size(signal,2);
%     [~,s_min,~] = SubtractMinimalEnvelope(signal(:,cc));
%     [~,s_max,~] = SubtractMinimalEnvelope(-signal(:,cc));
%     blah = (s_min(:) - s_max(:))/2;
%     [~,s_min,~] = SubtractMinimalEnvelope(blah);
%     [~,s_max,~] = SubtractMinimalEnvelope(-blah);
%     s(:,cc)  = signal(:,cc) - (s_min(:) - s_max(:))/2;
% end

s = signal;

if ndims(signal)> 2
    errordlg('Input signal size cannot exceed 2 dimensions!');
elseif size(signal,2) > 2
    errordlg('Signal input must be a matrix with no more than 2 cols, with each col being a different timeseries!')
elseif any(size(signal)==1)
    [Wxy,freq,coi, sig95]  = ComputeXWT(s,s,timeVec,freqRange,dj,stringency,phaseType);
else
    [Wxy,freq,coi, sig95]  = ComputeXWT(s(:,1),s(:,2),timeVec,freqRange,dj,stringency,phaseType);
end

% s_norm = s/max(s);
% s_norm = log2(s);
% sigMat = repmat(s(:)',size(Wxy,1),1);
R = zeros(size(Wxy));
% B = log2(abs(Wxy)).*sigMat;
% B = log2(abs(Wxy));
B = abs(Wxy).^1;

for tt = 1:size(B,2)
    blah = B(:,tt);
    blah = SubtractMinimalEnvelope(blah);
    dBlah = diff(blah);
    peaks = dBlah(1:end-1)>0 & dBlah(2:end)<=0;
    peakInds = find(peaks);
    blah = blah./max(blah);
    peakInds(blah(peakInds)<0.75)=[];
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