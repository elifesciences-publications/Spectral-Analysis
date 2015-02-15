function S = ComputeSparseSpectrogram(varargin)
%ComputeSparseSpectrogram - In addition to plotting, returns a
% sparse spectrogram that is computed using minimal EMD (Empirical Mode Decomposition).
% Frequencies are computed by detecting peaks of minimal (1 iteration, 3 decompositions)
% IMFs (Intrinsic Mode Functions).
%
% S = ComputeSparseSpectrogram(signal,timeVec,peakDetThresh,freqRange,colorMap)
% S = ComputeSparseSpectrogram(signal,samplingInt,peakDetThresh,freqRange,colorMap)

x = varargin{1};
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


%% Frequency from IMF1
pks1(amp1 < thresh)= [];
amp1(amp1 < thresh) = [];
t1 = time(pks1);
f1 = 1./diff(t1);
t_adj = (t1(1:end-1) + t1(2:end))/2;
t1 = t_adj;
outInds = find(t1 > time(end));
t1(outInds)=[];
amp_adj = abs(sqrt(amp1(1:end-1).*amp1(2:end))); % Geometric mean of a point and the next
amp1 = amp_adj;
amp1(outInds)=[];
f1(outInds)=[];
outInds = find(f1 > freqRange(2) | f1 < freqRange(1));
f1(outInds) = [];
t1(outInds) =[];
amp1(outInds) = [];


%% Frequency from IMF2
pks2(amp2 < thresh/2)= [];
amp2(amp2 < thresh/2) = [];
t2 = time(pks2);
f2 = 1./diff(t2);
t_adj = (t2(1:end-1) + t2(2:end))/2;
t2 = t_adj;
outInds = find(t2 > time(end));
t2(outInds)=[];
amp_adj = abs(sqrt(amp2(1:end-1).*amp2(2:end))); % Geometric mean of a point and the next
amp2 = amp_adj;
amp2(outInds)=[];
f2(outInds)=[];
outInds = find(f2 > freqRange(2) | f2 < freqRange(1));
t2(outInds) =[];
f2(outInds) = [];
amp2(outInds) = [];

%% Frequency from IMF3
pks3(amp3 < thresh/4)= [];
amp3(amp3 < thresh/4) = [];
t3 = time(pks3);
f3 = 1./diff(t3);
t_adj = (t3(1:end-1) + t3(2:end))/2;
t3 = t_adj;
outInds = find(t3 > time(end));
t3(outInds)=[];
amp_adj = abs(sqrt(amp1(1:end-1).*amp1(2:end))); % Geometric mean of a point and the next
amp1 = amp_adj;
amp1(outInds)=[];
f1(outInds)=[];

outInds = find(f3 > freqRange(2) | f3 < freqRange(1));
t3(outInds) =[];
f3(outInds) = [];
amp3(outInds) = [];

%% Gathering all points together
t_all = [t1(:); t2(:); t3(:)];
f_all = [f1(:); f2(:); f3(:)];
amp_all = [amp1(:); amp2(:); amp3(:)];
S = [t_all f_all amp_all];

weakInds = find(amp_all < (0.1*max(amp_all)));
S(weakInds,:)=[];


%  [colVals,LUT,cmap_new] = MapValsToColors(S(:,3),colorMap);

colVals = MapValsToColorsVer4(S(:,3),colorMap);
%  cmap_new = cmap;



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
ylim(freqRange), xlim([time(1) time(end)])
for jj = 1:size(S,1)
    plot(S(jj,1),S(jj,2),'color',colVals(jj,:),'marker','o','markersize',5,'markerfacecolor',colVals(jj,:))
end
% colormap(cmap_new)
ch = colorbar;
set(ch,'location','NorthOutside')





