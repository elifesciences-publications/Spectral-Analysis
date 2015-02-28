function varargout = ComputeSparseSpectrogram(varargin)
%ComputeSparseSpectrogram - In addition to plotting, returns a
% sparse spectrogram that is computed using minimal EMD (Empirical Mode Decomposition).
% Frequencies are computed by detecting peaks of minimal (1 iteration, 3 decompositions)
% IMFs (Intrinsic Mode Functions).
%
% S = ComputeSparseSpectrogram(signal,timeVec,peakDetThresh,freqRange,colorMap,plotOrNot)
% S = ComputeSparseSpectrogram(signal,samplingInt,peakDetThresh,freqRange,colorMap,plotOrNot)

x = varargin{1};
if all(size(x)> 1)
    errordlg('First input is a matrix: please enter a vector!')
    return;
end
x = x(:);
N = length(x);
switch nargin
    case 1
        errordlg('At least 2 inputs reqd!')
        return;
    case 2
        if numel(varargin{2})==1; % Assuming this is sampling interval
            samplingInt = varargin{2};
            time = (0:length(x)-1)*samplingInt;
        else
            time = varargin{2};
            samplingInt = mode(diff(time));
        end
        thresh = 5;
        freqRange(1) = 0.5/(length(x)*samplingInt);
        freqRange(2) = 1./(2*samplingInt);
        colorMap = hot; 
        plotOrNot = 'y';
    case 3
          if numel(varargin{2})==1; % Assuming this is sampling interval
            samplingInt = varargin{2};
            time = (0:length(x)-1)*samplingInt;
        else
            time = varargin{2};
            samplingInt = mode(diff(time));
        end
        thresh = varargin{3};
        freqRange(1) = 0.5/(length(x)*samplingInt);
        freqRange(2) = 1/(2*samplingInt);
        colorMap = hot; 
        plotOrNot = 'y';
    case 4
          if numel(varargin{2})==1; % Assuming this is sampling interval
            samplingInt = varargin{2};
            time = (0:length(x)-1)*samplingInt;
        else
            time = varargin{2};
            samplingInt = mode(diff(time));
        end
        thresh = varargin{3};
        freqRange = varargin{4};
        if numel(freqRange) < 2
            errordlg('Freq Range input must contain at least 2 values!')
            return;
        end
         colorMap = hot;
         plotOrNot = 'y';
    case 5
        if numel(varargin{2})==1; % Assuming this is sampling interval
            samplingInt = varargin{2};
            time = (0:length(x)-1)*samplingInt;
        else
            time = varargin{2};
            samplingInt = mode(diff(time));
        end
        thresh = varargin{3};
        freqRange = varargin{4};
        if numel(freqRange) < 2
            errordlg('Freq Range input must contain at least 2 values!')
            return;
        end
        colorMap = varargin{5};
        plotOrNot = 'y';
    case 6
        if numel(varargin{2})==1; % Assuming this is sampling interval
            samplingInt = varargin{2};
            time = (0:length(x)-1)*samplingInt;
        else
            time = varargin{2};
            samplingInt = mode(diff(time));
        end
        thresh = varargin{3};
        freqRange = varargin{4};
        if numel(freqRange) < 2
            errordlg('Freq Range input must contain at least 2 values!')
            return;
        end
        colorMap = varargin{5}; 
        plotOrNot = varargin{6};
end

if length(time)~= length(x)
    errordlg('Length of time vector does not match length of signal, please check inputs!')
    return;
end

% thresh = 40;
% newSamplingRate = round(3*freqRange(2));
% newSamplingInt = 1./newSamplingRate;
% dt = ceil(samplingRate/newSamplingRate);
% x = x(1:dt:end);
% time = time(1:dt:end);

% x = x(:)';

% specMat = zeros(length(freqRange(1):freqRange(2)),N);


%% Getting 3 IMFs using cubic interpolation

[imf1,r1,pks1] = GetCubic(x);
amp1 = abs(imf1(pks1));

[imf2,r2,pks2] = GetCubic(r1);
amp2 = abs(imf2(pks2));

[imf3,r3,pks3] = GetCubic(r2);
amp3 = abs(imf3(pks3));


maxFrac = 0.2;
%% Frequency from IMF1
thresh = mean(amp1) + thresh*std(amp1);
pks1(amp1 < thresh)= [];
amp1(amp1 < thresh) = [];
t1 = time(pks1);
f1 = 1./diff(t1);
t_adj = (t1(1:end-1) + t1(2:end))/2;
t1 = t_adj;
pks1 = round((pks1(1:end-1) + pks1(2:end))/2);
outInds = find(t1 > time(end));
t1(outInds)=[];

amp_adj = abs(sqrt(amp1(1:end-1).*amp1(2:end))); % Geometric mean of a point and the next
amp1 = amp_adj;
amp1(outInds)=[];
f1(outInds)=[];
pks1(outInds)=[];
outInds = find(f1 > freqRange(2) | f1 < freqRange(1));
f1(outInds) = [];
t1(outInds) =[];
amp1(outInds) = [];
pks1(outInds)=[];

weakInds = find(amp1 < (maxFrac*max(amp1)));
amp1(weakInds)=[];
pks1(weakInds) = [];
f1(weakInds) = [];
t1(weakInds) = [];

%% Frequency from IMF2
thresh = mean(amp2) + 0.9*thresh*std(amp2);
pks2(amp2 < thresh/20)= [];
amp2(amp2 < thresh/20) = [];
t2 = time(pks2);
f2 = 1./diff(t2);
t_adj = (t2(1:end-1) + t2(2:end))/2;
t2 = t_adj;
pks2 = round((pks2(1:end-1) + pks2(2:end))/2);
outInds = find(t2 > time(end));
t2(outInds)=[];
amp_adj = abs(sqrt(amp2(1:end-1).*amp2(2:end))); % Geometric mean of a point and the next
amp2 = amp_adj;
amp2(outInds)=[];
f2(outInds)=[];
pks2(outInds)=[];
outInds = find(f2 > freqRange(2) | f2 < freqRange(1));
t2(outInds) =[];
f2(outInds) = [];
amp2(outInds) = [];
pks2(outInds) =[];

weakInds = find(amp2 < (maxFrac*max(amp2)));
amp2(weakInds)=[];
pks2(weakInds) = [];
f2(weakInds) = [];
t2(weakInds) = [];

%% Frequency from IMF3
thresh = mean(amp3) + 0.8*thresh*std(amp3);
pks3(amp3 < thresh/40)= [];
amp3(amp3 < thresh/40) = [];
t3 = time(pks3);
f3 = 1./diff(t3);
t_adj = (t3(1:end-1) + t3(2:end))/2;
t3 = t_adj;
pks3 = round((pks3(1:end-1) + pks3(2:end))/2);
outInds = find(t3 > time(end));
t3(outInds)=[];
amp_adj = abs(sqrt(amp3(1:end-1).*amp3(2:end))); % Geometric mean of a point and the next
amp3 = amp_adj;
amp3(outInds)=[];
f3(outInds)=[];
pks3(outInds) = [];
outInds = find(f3 > freqRange(2) | f3 < freqRange(1));
t3(outInds) =[];
f3(outInds) = [];
amp3(outInds) = [];
pks3(outInds) = [];

weakInds = find(amp3 < (maxFrac*max(amp3)));
amp3(weakInds)=[];
pks3(weakInds) = [];
f3(weakInds) = [];
t3(weakInds) = [];


%% Gathering all points together
t_all = [t1(:); t2(:); t3(:)];
f_all = [f1(:); f2(:); f3(:)];
amp_all = [amp1(:); amp2(:); amp3(:)];
pks_all = [pks1(:); pks2(:); pks3(:)];

%% Arrange all values in chronological order.
[t_all,inds] = sort(t_all,'ascend');
f_all = f_all(inds);
amp_all = amp_all(inds);
pks_all = pks_all(inds);

inds = zeros(size(t_all));
for jj = 2:length(t_all)
    if (t_all(jj) - t_all(jj-1)) <= (2*samplingInt) && (f_all(jj) - f_all(jj-1) <=2)
        t_all(jj-1) = median([t_all(jj), t_all(jj-1)]);
        f_all(jj-1) = mean([f_all(jj), f_all(jj-1)]);
        amp_all(jj-1) = sqrt(amp_all(jj)* amp_all(jj-1));
        pks_all(jj-1) = round(mean([pks_all(jj), pks_all(jj-1)]));
        inds(jj) = jj;
    end
end
inds(inds ==0)=[];
t_all(inds) = [];
f_all(inds) = [];
amp_all(inds) = [];
pks_all(inds) = [];
S = [f_all amp_all, pks_all,t_all];

varargout{1} = S;

% weakInds = find(amp_all < (0.08*max(amp_all)));
% S(weakInds,:)=[];


%  [colVals,LUT,cmap_new] = MapValsToColors(S(:,3),colorMap);

colVals = MapValsToColorsVer4(S(:,2),colorMap);
%  cmap_new = cmap;


if strcmpi(plotOrNot,'y')
    figure('Name','Mimal Spectrogram')
    subaxis(2,1,1,'mb',0.1,'mt', 0.1)
    hold on
    set(gca,'color','k','tickdir','out', 'xtick',[]), box off
    plot(time,x,'color','g')
    title('Minimal Spectrogram')
    ylabel('Amp (z-scores)')
    ylim([-inf inf]), xlim([time(1) time(end)])
    hold off
    subaxis(2,1,2,'mt',0.001,'mb',0.1)
    hold on
    set(gca,'color','k','tickdir','out'), box off
    ylabel('Frequency (Hz)'), xlabel('Time (sec)')
    axis([time(1) time(end) min(f_all) max(f_all)])
    % ylim(freqRange), xlim([time(1) time(end)])
    for jj = 1:size(S,1)
        plot(S(jj,4),S(jj,1),'color',colVals(jj,:),'marker','o','markersize',5,'markerfacecolor',colVals(jj,:))
    end
    colormap(colorMap)
    ch = colorbar;
    set(ch,'location','NorthOutside')
end


if nargout == 2
    freqVec  = zeros(size(time));
    ampVec = zeros(size(time));
    freqVec(pks_all) = f_all;
    ampVec(pks_all) = amp_all;
    SS =[freqVec(:), ampVec(:)]; 
    varargout{2} = SS;
end





