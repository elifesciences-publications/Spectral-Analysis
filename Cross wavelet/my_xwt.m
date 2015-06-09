function [varargout] = my_xwt(x,y,t,varargin)
% My custom written XWT based on xwt.m by Grinsted et al
% [Wxy,period,scale,coi,sig95] = my_xwt(x,y,time);
% [Wxy,period,scale,coi,sig95] = my_xwt(x,y,time,freqRange,stringency,threshold,plotOrNot,sigmaxy);

% plotOrNot = 0 ==> Don't plot; 1 ==> Plot.

%% Reading inputs
if nargin < 3
    errordlg('Minimum 3 inputs reqd');
    return;
elseif nargin < 4
    freqRange = [10 150];
    stringency = 1;
    threshold = 3;
    plotOrNot = 1;
elseif nargin < 5
    freqRange = varargin{1};
    stringecy = 1;
    threshold =3;
    plotOrNot = 1;
elseif nargin < 6
    freqRange = varargin{1};
    stringency = varargin{2};
    threshold = 3;
    plotOrNot = 1;
elseif nargin < 7
    freqRange = varargin{1};
    stringency = varargin{2};
    threshold = varargin{3};
    plotOrNot = 1;
elseif nargin < 8
     freqRange = varargin{1};
    stringency = varargin{2};
    threshold = varargin{3};
    plotOrNot = varargin{4};
else
    freqRange = varargin{1};
    stringency = varargin{2};
    threshold = varargin{3};
    plotOrNot = varargin{4};
    sigmaxy = varargin{5};
end

noiseType = 'white'; %('red' or 'white')

%% Fixed parameters & adjustable parameters
fourier_factor = 1.0330; % Conversion factor for changing wavelet scales to periods.

peakDetectionThreshold = 0.01; % Determines the peak detection for global wavelet spectrum plotted to the right of XWT.
% Lower value results in detection of smaller peaks.

scaleRange = 1./(freqRange*fourier_factor);
S0 = min(scaleRange);
MaxScale = max(scaleRange);
Args=struct('Pad',1,...      % pad the time series with zeroes (recommended)
    'Dj',1/32, ...    % this will do 48 sub-octaves per octave
    'S0',S0,...    % this says start at a scale of 2*dt
    'J1',[],...
    'Mother','Morlet', ...
    'MaxScale',MaxScale,...   % a more simple way to specify J1
    'MakeFigure',(nargout==0),...
    'BlackandWhite',0,...
    'AR1','auto',...
    'ArrowDensity',[30 30],...
    'ArrowSize',1,...
    'ArrowHeadSize',1);


x = [t(:) x(:)];
y = [t(:) y(:)];

% ------validate and reformat timeseries.
[x,dt]=formatts(x);
[y,dty]=formatts(y);
if dt~=dty
    errordlg('Timestep must be equal for time series')
end

t=(max(x(1,1),y(1,1)):dt:min(x(end,1),y(end,1)))'; % Common time period
if length(t)<4
    errordlg('The two time series must overlap.')
end

n=length(t);
dt = t(2)-t(1);
dty= dt;

%% Verifying arguments for wavelet transform
if isempty(Args.J1)
    if isempty(Args.MaxScale)
        Args.MaxScale=(n*.17)*2*dt; %auto maxscale
    end
    Args.J1=round(log2(Args.MaxScale/Args.S0)/Args.Dj);
end

ad=mean(Args.ArrowDensity);
Args.ArrowSize=Args.ArrowSize*30*.03/ad;
Args.ArrowHeadSize=Args.ArrowHeadSize*Args.ArrowSize*220;


if strcmpi(Args.AR1,'auto')
    Args.AR1=[ar1nv(x(:,2)) ar1nv(y(:,2))];
    if any(isnan(Args.AR1))
        error('Automatic AR1 estimation failed. Specify them manually (use the arcov or arburg estimators).')
    end
end


if ~ exist('sigmaxy')
% sigmax = std(x(:,2));
% sigmay = std(y(:,2));
[~,muX,sigmax] = ZscoreByHist(x(:,2)); 
[~,muY,sigmay] = ZscoreByHist(y(:,2));
else
    sigmax = sigmaxy;
    sigmay = sigmaxy;
end


% ----- Wavelet Transforms
[X,period,scale,coix] = wavelet(x(:,2),dt,Args.Pad,Args.Dj,Args.S0,Args.J1,Args.Mother);
[Y,period,scale,coiy] = wavelet(y(:,2),dt,Args.Pad,Args.Dj,Args.S0,Args.J1,Args.Mother);
freq = 1./period;

% ----- Truncate X,Y to common time interval (this is first done here so that the coi is minimized)
dte=dt*.01; %to cricumvent round off errors with fractional timesteps
idx=find((x(:,1)>=(t(1)-dte))&(x(:,1)<=(t(end)+dte)));
X=X(:,idx);
coix=coix(idx);

idx=find((y(:,1)>=(t(1)-dte))&(y(:,1)<=(t(end)+dte)));
Y=Y(:,idx);
coiy=coiy(idx);


coi=min(coix,coiy);

% -------- Cross
Wxy=X.*conj(Y);

[~,muWxy,sigWxy] = ZscoreByHist(abs(Wxy(:)));

% muWxy = mean(Wxy(:));
% sigWxy = std(Wxy(:));
threshWxy = muWxy + threshold*sigWxy;


%---- Significance levels
if strcmpi(noiseType,'red')
    Pkx=ar1spectrum(Args.AR1(1),period./dt);
    Pky=ar1spectrum(Args.AR1(2),period./dt);
elseif strcmpi(noiseType,'white')
    Pkx=ar1spectrum(0,period./dt);
    Pky=ar1spectrum(0,period./dt);
end

V=2; %(default: V = 2)
Zv = 3.9999; %(default: Zv = 3.9999)
signif=sigmax*sigmay*sqrt(Pkx.*Pky)*Zv/V;
sig95 = (signif')*(ones(1,n));  % expand signif --> (J+1)x(N) array
sig95 = abs(Wxy) ./ sig95;
if ~strcmpi(Args.Mother,'morlet')
    sig95(:)=nan;
end
Wxy(sig95 < stringency)=0;
% Wxy(abs(Wxy) < threshWxy) = 0;
varargout={Wxy,period,scale,coi,sig95};
varargout=varargout(1:nargout);



if plotOrNot       
    %% Calculating XW power spectrum
    powerSpectrum = sum(abs(Wxy),2);
    normPowerSpectrum = powerSpectrum/max(powerSpectrum);
    logPowerSpectrum = log2(powerSpectrum);
    normLogPowerSpectrum = logPowerSpectrum/max(logPowerSpectrum);
    [maxtab,~] = peakdet(normLogPowerSpectrum,peakDetectionThreshold);
    if numel(maxtab) == 0
        maxtab(:,1)= find(normPowerSpectrum==max(normPowerSpectrum));
        maxtab(:,2)= normPowerSpectrum(normPowerSpectrum == max(normPowerSpectrum));
    end
        
    
    %% Plotting Figures
    figure('Name','XWT','color','w')
    figPos = get(gcf,'position');
    set(gcf,'position',[figPos(1) figPos(2)-figPos(4)/2 figPos(3)* 1.5...
        figPos(4)+figPos(4)/2]);
    
    %% XWT Axes
    cData = Wxy;
    cData(cData==0) = NaN;
    ax1 = axes; box off
    aPos = get(ax1,'position');
    aPos = [aPos(1)-0.02 aPos(2)+ aPos(4)*(1/3) aPos(3)*(0.8) aPos(4)*(0.75)];
    set(ax1,'position', aPos)
    plotwavelin(cData,t,period,coi,sig95,sigmax, sigmay,stringency)
    set(ax1,'color','k','xtick',[], 'xcolor','w','ycolor','k')
    xlim([t(1) t(end)]) %%%% This line is NECESSARY to ensure that x-axis is aligned with traces below
    xlabel('')
    ylims1 = get(ax1,'ylim');
    yticklabels_ax1 = str2num(get(ax1,'yticklabel'));
    dyticklabels_ax1 = diff(yticklabels_ax1);
    if dyticklabels_ax1(1)== dyticklabels_ax1, yScale = 'linear'; else yScale = 'log'; end
    yl = ylabel('Frequency (Hz)','fontsize',14);
    ylpos = get(yl,'pos');
    aPos = get(ax1,'position');
    colormap('jet')
    %% XW Power Spectrum Axes
    ax2 = axes; hold on, box off
    aPos2 = get(ax2,'pos');
    aPos2 = [aPos(1) + aPos(3) aPos(2) aPos2(3)*(0.3) aPos(4)];
    set(ax2,'pos', aPos2, 'color','w','tickdir','out','fontsize',14);
    xlabel([{'Normalized'}; {'Power'}])
    peakFreqs = freq(maxtab(:,1));
    if strcmpi(yScale,'linear'), plot(normPowerSpectrum, freq,'k','linewidth',2);
        hold on
        plot(normLogPowerSpectrum, freq,'k:','linewidth',2)
        logPeakFreqs = peakFreqs;
    else
        plot(normPowerSpectrum,log2(freq),'k','linewidth',2)
        hold on
        plot(normLogPowerSpectrum,log2(freq),'k:','linewidth',2)
        logPeakFreqs = log2(peakFreqs);
    end
    set(ax2,'ylim',ylims1,'xlim',[0 1],'xtick',[0.5 1],'ytick',[])
    ylims2 = get(ax2,'ylim');
    
    xvals = 0.2*ones(size(peakFreqs)); % X coordinates for text displaying peak frequencies
    
    for pf = 1:length(peakFreqs)
        txt{pf} = [num2str(round(peakFreqs(pf)*100)/100) ' Hz'];
    end
    
    text(xvals,logPeakFreqs,txt,'fontsize',14,'color','r');
    
    
    ax2b = axes('position',aPos2,'xaxislocation','top',...
        'yaxislocation','right','color','none','xscale','log','ytick',[],'ycolor','w','fontsize',14);
    xt = [0.25 0.5 1];
    set(ax2b,'xtick',xt,'xlim',[0 1]);
    leg = {'Linear';'Log'};
    legend(ax2,leg)
    
    
    
    
    %% Time Series Axes
    ax3=axes; hold on, box off
    aPos3 = get(ax3,'position');
    aPos3 = [aPos(1) aPos3(2) aPos(3) aPos3(4)*(1/3)];
    set(ax3,'position',aPos3,'tickdir','out','color','w','ycolor','w');
    
    combAmp = std(x(:,2)+ y(:,2));
    
    plot(x(:,1),x(:,2) + 1.5*combAmp ,'k','linewidth',1.5)
    
    plot(y(:,1),y(:,2) - 1.5*combAmp,'k','linewidth',1.5)
    
    yl2 = ylabel([{'Normalized'};{'Amplitude'}],'fontsize',14,'color','k');
    % ylpos2 = get(yl2,'pos');
    % ylpos2_mod = ylpos2;
    % ylpos2_mod(1) = t(1)-abs(ylpos(1))+ 0.65; ylpos2_mod(2)=0;
    % set(yl2,'pos',ylpos2_mod);
    
    set(ax3,'ytick',[])
    axis([t(1) t(end) -inf inf])
    xlabel('Time (s)','fontsize',14)
    set(gca,'fontsize',14)
    set(ax3,'ytick',[])
    hold off;
    
    
end

%% Pending fixes
% 1)It would be nice to not display arrows in insignificant regions
% 2)Ticks on colorbar seem a little screwy