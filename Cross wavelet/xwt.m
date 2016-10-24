function varargout=xwt(x,y,varargin)
%% Cross wavelet transform
% Creates a figure of cross wavelet power in units of
% normalized variance.
%
% USAGE: [Wxy,period,scale,coi,sig95]=xwt(x,y,[,settings])
%
% x & y: two time series
% Wxy: the cross wavelet transform of x against y
% period: a vector of "Fourier" periods associated with Wxy
% scale: a vector of wavelet scales associated with Wxy
% coi: the cone of influence
%
% Settings: Pad: pad the time series with zeros?
% .         Dj: Octaves per scale (default: '1/12')
% .         S0: Minimum scale
% .         J1: Total number of scales
% .         Mother: Mother wavelet (default 'morlet')
% .         MaxScale: An easier way of specifying J1
% .         MakeFigure: Make a figure or simply return the output.
% .         BlackandWhite: Create black and white figures
% .         AR1: the ar1 coefficients of the series
% .              (default='auto' using a naive ar1 estimator. See ar1nv.m)
% .         ArrowDensity (default: [30 30])
% .         ArrowSize (default: 1)
% .         ArrowHeadSize (default: 1)
%
% Settings can also be specified using abbreviations. e.g. ms=MaxScale.
% For detailed help on some parameters type help wavelet.
%
% Example:
%    t=1:200;
%    xwt(sin(t),sin(t.*cos(t*.01)),'ms',16)
%
% Phase arrows indicate the relative phase relationship between the series
% (pointing right: in-phase; left: anti-phase; down: series1 leading
% series2 by 90°)
%
% Please acknowledge the use of this software in any publications:
%   "Crosswavelet and wavelet coherence software were provided by
%   A. Grinsted."
%
% (C) Aslak Grinsted 2002-2004
%
% http://www.pol.ac.uk/home/research/waveletcoherence/

% -------------------------------------------------------------------------
%   Copyright (C) 2002-2004, Aslak Grinsted
%   This software may be used, copied, or redistributed as long as it is not
%   sold and this copyright notice is reproduced on each copy made.  This
%   routine is provided as is without any express or implied warranties
%   whatsoever.

% Custom modifications by Avinash Pujala, Janelia Research Campus, 2015

%% Choose background noise type (for statistical significance testing)
noise_type = 'red'; % ('red' or 'white'; default: 'white')

%% Validate and reformat timeseries
[x,dt]=formatts(x);
[y,dty]=formatts(y);
if dt~=dty
    error('timestep must be equal between time series')
end
% t=(max(x(1,1),y(1,1)):dt:min(x(end,1),y(end,1)))'; 
t = linspace(max(x(1,1),y(1,1)),min(x(end,1),y(end,1)),size(x,1)); %common time period
t = t(:);
if length(t)<4
    error('The two time series must overlap.')
end

n = length(t);
if n < size(x,1)
    t = linspace(t(1),t(end), size(x,1));
end
dt = t(2)-t(1);
dty = dt;
n = length(t);

dj = varargin{2};
%----------default arguments for the wavelet transform-----------
Args=struct('Pad',1,...      % pad the time seri                                                                                                                                                                                                                                   es with zeroes (recommended)                                                                                                                                                                                                                                                                                                                                                  
    'Dj',dj, ...    % this will do 64 sub-octaves per octave
    'S0',2*dt,...    % this says start at a scale of 2 years                                                                                                                                                                                                                                                                                                                                   
    'J1',[],...
    'Mother','DOG', ...
    'MaxScale',[],...   %a more simple way to specify J1
    'MakeFigure',(nargout==0),...
    'BlackandWhite',0,...
    'AR1','auto',...
    'ArrowDensity',[30 30],...
    'ArrowSize',1,...
    'ArrowHeadSize',1);

Args=parseArgs(varargin,Args,{'BlackandWhite'});
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

% [~,~,sigmax] = ZscoreByHist(x(:,2));
% [~,~,sigmay] = ZscoreByHist(y(:,2));
%# Although, ZscoreByHist() will likely provide a better estimate of the true std of a timeseries, I will
%#  switch back to std() so as to make it consistent with std estimates in other dependent scripts.
sigmax = std(x(:,2));
sigmay = std(y(:,2));


%% Compute crosswavelet
[X,period,scale,coix] = wavelet(x(:,2),dt,Args.Pad,Args.Dj,Args.S0,Args.J1,Args.Mother);
[Y,period,scale,coiy] = wavelet(y(:,2),dt,Args.Pad,Args.Dj,Args.S0,Args.J1,Args.Mother);

%# Truncate X,Y to common time interval (this is first done here so that the coi is minimized)
dte=dt*.01; %to cricumvent round off errors with fractional timesteps
idx=find((x(:,1)>=(t(1)-dte))&(x(:,1)<=(t(end)+dte)));
X=X(:,idx);
coix=coix(idx);

idx=find((y(:,1)>=(t(1)-dte))&(y(:,1)<=(t(end)+dte)));
Y=Y(:,idx);
coiy=coiy(idx);

coi=min(coix,coiy);

Wxy=(X.*conj(Y)); % Not being normalized by division by sigmax*sigmay here.

switch lower(noise_type)
    case 'red'
        %       Pkx=ar1spectrum_ap(Args.AR1(1),period); % This version uses eqn 16
        %           from Torrence & Compo, 1998, which should be equivalent to
        %           Grinsted's version, but seems not to be (????)
        Pkx =ar1spectrum(Args.AR1(1),period./dt);
        Pky=ar1spectrum(Args.AR1(2),period./dt);
    case 'white'
        %       Pkx=ar1spectrum_ap(0,period);
        Pkx = ar1spectrum(0,period./dt);
        Pky=ar1spectrum(0,period./dt);
end

V=2;
Zv=3.9999; %(default:Zv = 3.999; Grinsted et al., 2004, eqn (5))
signif=sigmax*sigmay*sqrt(Pkx.*Pky)*Zv/V; % Eqn (5)
sig95 = (signif')*(ones(1,n));  % expand signif --> (J+1)x(N) array
sig95 = abs(Wxy)./sig95;
if ~strcmpi(Args.Mother,'morlet')
    sig95(:)=nan;
end

if 0 %Args.MakeFigure
    Yticks = 2.^(fix(log2(min(period))):fix(log2(max(period))));
    if Args.BlackandWhite
        levels = [0.25,0.5,1,2,4,8,16];
        [cout,H]=safecontourf(t,log2(period),log2(abs(Wxy/(sigmax*sigmay))),log2(levels));%,log2(levels));  %*** or use 'contourf3ill'
        cout(1,:)=2.^cout(1,:);
        HCB=colorbarf(cout,H);
        barylbls=rats([0 levels 0]');
        barylbls([1 end],:)=' ';
        barylbls(:,all(barylbls==' ',1))=[];
        set(HCB,'yticklabel',barylbls);
        cmap=(1:-.01:.5)'*.9;
        cmap(:,2:3)=cmap(:,[1 1]);
        %cmap(:,1:2)=cmap(:,1:2)*.8;
        colormap(cmap);
        set(gca,'YLim',log2([min(period),max(period)]), ...
            'YDir','reverse', ...
            'YTick',log2(Yticks(:)), ...
            'YTickLabel',num2str(Yticks'), ...
            'layer','top')
        %xlabel('Time')
        ylabel('Period')
        hold on
        
        aWxy=angle(Wxy);
        
        phs_dt=round(length(t)/Args.ArrowDensity(1)); tidx=max(floor(phs_dt/2),1):phs_dt:length(t);
        phs_dp=round(length(period)/Args.ArrowDensity(2)); pidx=max(floor(phs_dp/2),1):phs_dp:length(period);
        phaseplot(t(tidx),log2(period(pidx)),aWxy(pidx,tidx),Args.ArrowSize,Args.ArrowHeadSize);
        
        if strcmpi(Args.Mother,'morlet')
            [c,h] = contour(t,log2(period),sig95,[1 1],'k');%#ok
            set(h,'linewidth',3)
        else
            warning('XWT:sigLevelNotValid','XWT Significance level calculation is only valid for morlet wavelet.')
            %TODO: alternatively load from same file as wtc (needs to be coded!)
        end
        
        %tt=[t([1 1])-dt*.5;t;t([end end])+dt*.5];
        %hcoi=patch(tt,log2([period([end 1]) coi period([1 end])]),ones(size(tt))*0,'w');
        %set(hcoi,'alphadatamapping','direct','facealpha',.8)
        
        plot(t,log2(coi),'k','linewidth',3)
        %hcoi=fill([t([1 1:end end])],log2([period(end) coi period(end)]),'r')
        %set(hcoi,'alphadatamapping','direct','facealpha',.3)
        hold off
    else
        H=imagesc(t,log2(period),log2(abs(Wxy/(sigmax*sigmay))));%#ok
        %logpow=log2(abs(Wxy/(sigmax*sigmay)));
        %[c,H]=safecontourf(t,log2(period),logpow,[min(logpow(:)):.25:max(logpow(:))]);
        %set(H,'linestyle','none')
        
        clim=get(gca,'clim'); %center color limits around log2(1)=0
        clim=[-1 1]*max(clim(2),3);
        set(gca,'clim',clim)
        
        HCB=safecolorbar;
        set(HCB,'ytick',-7:7);
        barylbls=rats(2.^(get(HCB,'ytick')'));
        barylbls([1 end],:)=' ';
        barylbls(:,all(barylbls==' ',1))=[];
        set(HCB,'yticklabel',barylbls);
        
        set(gca,'YLim',log2([min(period),max(period)]), ...
            'YDir','reverse', ...
            'YTick',log2(Yticks(:)), ...
            'YTickLabel',num2str(Yticks'), ...
            'layer','top')
        %xlabel('Time')
        ylabel('Period')
        hold on
        
        aWxy=angle(Wxy);
        
        phs_dt=round(length(t)/Args.ArrowDensity(1)); tidx=max(floor(phs_dt/2),1):phs_dt:length(t);
        phs_dp=round(length(period)/Args.ArrowDensity(2)); pidx=max(floor(phs_dp/2),1):phs_dp:length(period);
        phaseplot(t(tidx),log2(period(pidx)),aWxy(pidx,tidx),Args.ArrowSize,Args.ArrowHeadSize);
        
        if strcmpi(Args.Mother,'morlet')
            [c,h] = contour(t,log2(period),sig95,[1 1],'k');%#ok
            set(h,'linewidth',2)
        else
            warning('XWT:sigLevelNotValid','XWT Significance level calculation is only valid for morlet wavelet.')
            %TODO: alternatively load from same file as wtc (needs to be coded!)
        end
        tt=[t([1 1])-dt*.5;t;t([end end])+dt*.5];
        hcoi=fill(tt,log2([period([end 1]) coi period([1 end])]),'w');
        set(hcoi,'alphadatamapping','direct','facealpha',.5)
        hold off
    end
end

varargout={Wxy,period,scale,coi,sig95};
varargout=varargout(1:nargout);





function [cout,H]=safecontourf(varargin)
vv=sscanf(version,'%i.');
if (version('-release')<14)|(vv(1)<7)
    [cout,H]=contourf(varargin{:});
else
    [cout,H]=contourf('v6',varargin{:});
end

function hcb=safecolorbar(varargin)
vv=sscanf(version,'%i.');

if (version('-release')<14)|(vv(1)<7)
    hcb=colorbar(varargin{:});
else
    hcb=colorbar('v6',varargin{:});
end