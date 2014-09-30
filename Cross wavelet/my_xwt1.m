function [varargout] = my_xwt1(x,y,t,varargin)
% My custom written XWT based on xwt.m by Grinsted et al
% [Wxy,period,scale,coi,sig95] = my_xwt1(x,y,time);

%% Fixed parameters
fourier_factor = 1.0330; % Conversion factor for changing wavelet scales to periods.

%% Adjustable parameter
peakDetectionThreshold = 0.2; % Determines the peak detection for global wavelet spectrum plotted to the right of XWT.
                     % Lower value results in detection of smaller peaks.
noise_type = 'white'; % ('red' or 'white'; default: 'red')                 
freqRange = [0.2 10]; % Power for frequencies only within this range will be calculated and displayed
powspec_coi = 'y'; %('y' = only computes power spectrum for regions within the COI).
ROI_time = [5 15]; % Time range within which to compute power. If [0 0], then computes power for whole timeseries
ROI_freq  = [0.2 8]; % Freq range within which to computer power.

scaleRange = 1./(freqRange*fourier_factor);
S0 = min(scaleRange);
MaxScale = max(scaleRange);
Args=struct('Pad',1,...      % pad the time series with zeroes (recommended)
    'Dj',1/48, ...    % this will do 48 sub-octaves per octave
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
    error('timestep must be equal between time series')
end % AP: Requires the two inputs to have time cols with identical time steps
t=(max(x(1,1),y(1,1)):dt:min(x(end,1),y(end,1)))'; %common time period
if length(t)<4
    error('The two time series must overlap.')% AP: Overlap by at least 4 time pts
end

n=length(t);
dt = t(2)-t(1);
% dty= dt;

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

sigmax = std(x(:,2));
sigmay = std(y(:,2));

%% ----- Wavelet Transforms
[X,period,scale,coix] = wavelet(x(:,2),dt,Args.Pad,Args.Dj,Args.S0,Args.J1,Args.Mother);
[Y,period,scale,coiy] = wavelet(y(:,2),dt,Args.Pad,Args.Dj,Args.S0,Args.J1,Args.Mother);
freq = 1./period;

%% ----- Truncate X,Y to common time interval (this is first done here so that the coi is minimized)
dte=dt*.01; %to cricumvent round off errors with fractional timesteps
idx=find((x(:,1)>=(t(1)-dte))&(x(:,1)<=(t(end)+dte)));
X=X(:,idx);
coix=coix(idx);

idx=find((y(:,1)>=(t(1)-dte))&(y(:,1)<=(t(end)+dte)));
Y=Y(:,idx);
coiy=coiy(idx);


coi=min(coix,coiy);

%% -------- Cross-wavelet
Wxy=X.*conj(Y);

%% ---- Significance levels
%Pk1=fft_theor(freq,lag1_1);
%Pk2=fft_theor(freq,lag1_2);
switch lower(noise_type)
    case 'red'
        %         Pkx=ar1spectrum(Args.AR1(1),period./dt);
        %         Pky=ar1spectrum(Args.AR1(2),period./dt); % This is from
        %         Grinsted et al (2004), but looks fishy and in disagreement with
        %         Torrence & Compo(1998). So, I modified it (see below)
        
        Pkx=ar1spectrum_ap(Args.AR1(1),period);
        Pky=ar1spectrum_ap(Args.AR1(2),period);
    case 'white'
        %         Pkx=ar1spectrum(0,period./dt);
        %         Pky=ar1spectrum(0,period./dt);
        
        Pkx=ar1spectrum_ap(0,period);
        Pky=ar1spectrum_ap(0,period);
end


V=2; %(default: V = 2)
Zv = 3.9999; %(default: Zv = 3.9999)
signif=sigmax*sigmay*sqrt(Pkx.*Pky)*Zv/V;
sig95 = (signif')*(ones(1,n));  % expand signif --> (J+1)x(N) array
sig95 = abs(Wxy) ./ sig95;
if ~strcmpi(Args.Mother,'morlet')
    
    sig95(:)=nan;
end
varargout={Wxy,period,scale,coi,sig95};
varargout=varargout(1:nargout);


%% Calculating XW power spectrum
freq =1./period;
ftmat = repmat(freq(:), 1, n); % n = length(time)
coimat = repmat(1./coi(:)',length(freq), 1);
Wxy_full = Wxy;
Wxy_coi = Wxy;
Wxy_coi(ftmat<coimat) = 0; % Setting power for regions outside coi to zero
Wxy_sig = Wxy_coi;
Wxy_sig(sig95<1)=0;



%% Frequencies, Phases, XW Powers, Etc
Wxy = Wxy_sig;
timeRange = [t(2) t(end)];
temp = [];
if ROI_time(1)==0, temp(1) = timeRange(1); else temp(1)= ROI_time(1); end
if ROI_time(2)==0, temp(2) = timeRange(2); else temp(2)= ROI_time(2); end
temp = ceil(temp/dt);

blah = [];
if ROI_freq(1)==0 || ROI_freq (1) < min(freq), blah(1) = min(freq); else blah(1)= ROI_freq(1); end
if ROI_freq(2)==0 || ROI_freq(2)> max(freq), blah(2) = max(freq); else blah(2)= ROI_freq(2); end

pts =[];
pts(1) = max(find(freq>=blah(2))); pts(2) = max(find(freq>=blah(1)));
zmat = zeros(size(Wxy));
zmat(pts(1):pts(2),temp(1):temp(2))=1;
Wxy_ROI = Wxy.*zmat;

if strcmpi(powspec_coi,'y')
% powerSpectrum = sum(abs(Wxy_sig),2);
  powerSpectrum = sum(abs(Wxy_ROI),2);
else
    powerSpectrum = sum(abs(Wxy_full),2);
end

normPowerSpectrum = powerSpectrum/max(powerSpectrum);
logPowerSpectrum = log2(powerSpectrum);
normLogPowerSpectrum = logPowerSpectrum/max(logPowerSpectrum);

[maxtab,mintab] = peakdet(normLogPowerSpectrum,peakDetectionThreshold);
if numel(maxtab) == 0
    maxtab(:,1)= find(normPowerSpectrum==max(normPowerSpectrum));
    maxtab(:,2)= normPowerSpectrum(normPowerSpectrum == max(normPowerSpectrum));
end

Pxy_ROI = sum(abs(Wxy_ROI(:))) % Total XW power within the ROI
Mxy_ROI = mean(abs(Wxy_ROI(:))) % Mean XW power within the ROI
% Wxy_prob = abs(Wxy./sum(Wxy(:))); % Each entry in Wxy is expressed as the probability of occurrence of that value within the distribution of all values of Wxy to create Wxy_prob

% Calculating mean & std of frequencies
%         Wxy_s = sum(Wxy_prob,2); 
%         meanFreqs = round((freq(:)'*Wxy_s)*100)/100; % Mean frequency is obtained by weighting the frequency coordinate of each value in Wxy with the corresponding value in Wxy_prob and taking the sum of all these values
%         Fxy = ftmat; Fxy(Wxy==0)=0;
%         Fxy_pow_weighted = Fxy.*Wxy_prob;
%         Fxy_pow_weighted(Fxy_pow_weighted==0)=[];
%         mf = mean(Fxy_pow_weighted);
%         sf = std(Fxy_pow_weighted);
%         stdFreqs = sf * meanFreqs/mf;
%         stdFreqs = round(stdFreqs*100)/100;
%         
%         % Calculating XW power spectrum & Peak Frequency(frequency at which coherent power is highest)
%         power_distribution_along_freq_axis = sum(abs(Wxy),2);
%         norm_power_distribution_along_freq_axis =...
%             power_distribution_along_freq_axis/...
%             max(power_distribution_along_freq_axis); % Power units indicate cumulative probability
%         normlog_power_distribution_along_freq_axis = ...
%             log2(power_distribution_along_freq_axis);
%         normlog_power_distribution_along_freq_axis = ...
%             normlog_power_distribution_along_freq_axis/max(normlog_power_distribution_along_freq_axis);
%         [maxtab,mintab] = peakdet(norm_power_distribution_along_freq_axis,peakDetectionThresh);
%         
%         if numel(maxtab) == 0
%             maxtab(:,1)= find(norm_power_distribution_along_freq_axis == max(norm_power_distribution_along_freq_axis));
%             maxtab(:,2)= norm_power_distribution_along_freq_axis(maxtab(:,1));
%         end
%         
%         pv = find(maxtab(:,2)== max(maxtab(:,2)));
%         pf = round(freq(maxtab(pv,1))*100)/100;
%         
%         if isempty(pv) | isempty(pf)
%             maxtab =[1 1+i];
%             pv = 1
%             pf= 1+i;
%         end
%         % Calculating mean and std of phases within significant regions
%         Wxy_nonzero_lin = Wxy; Wxy_nonzero_lin(Wxy==0)=[]; % This removes all zero elements from the matrix and vectorizes it.
%         if isempty(Wxy_nonzero_lin)
%             Wxy_nonzero_lin = rand(10,1)*(1+i);
%         end
%                 Axy = angle(Wxy_nonzero_lin); % Axy will also be a vector now       
%         nPhaseBins = min([numel(Axy), 90]);
%         mphase = circ_mean(Axy(:),abs(Wxy_nonzero_lin(:)));
%         %%%%% Alternatively, mphase = angle(sum(Wxy_nonzero_lin));
%         mphase(mphase<0) = mphase(mphase<0)+ 2*pi; % Addition of 2*pi to -ve values converts angle range from 0 to 360 rather than -180 to +180
%         meanPhases = mphase*180/pi; % Converts radians to degrees.
%         meanPhases = round(meanPhases*100)/100;
%         sphase = circ_std(Axy(:),abs(Wxy_nonzero_lin(:)));
%         stdPhases = sphase*180/pi;
%         stdPhases = round(stdPhases*100)/100;
%         [ph_dist,th] = hist(Axy(:),nPhaseBins);
%         pow_dist = hist(abs(Wxy_nonzero_lin(:)),nPhaseBins);
%         
%         %         [sortph,ind] = sort(Axy(:));
%         %         Z = abs(Wxy_nonzero_lin); Z = Z(:);
%         %         histWxy = Z(ind); % This essentially amounts to sorting by the ascending order of values in Axy
%         %         num_phases_per_bin = floor(numel(histWxy)/nPhaseBins);
%         %         histWxy = decimate(histWxy(:),num_phases_per_bin);
%         %         histWxy = histWxy(1:nPhaseBins);
%         ph_dist_wt = ph_dist(:).*pow_dist(:);
%         if isempty(ph_dist)
%             theta = zeros(size(ph_dist_wt));
%         else
%             ph_dist = [ph_dist(:); ph_dist(1)]; % This will close the loop in polar plot by circularizing the vector.
%             ph_dist_wt = [ ph_dist_wt(:); ph_dist_wt(1)];
%             phase_dist = ph_dist./max(ph_dist);
%             phase_dist_weight = ph_dist_wt./max(ph_dist_wt);
%             theta = [th(:);th(1)];
%             phf = find(phase_dist_weight == max(phase_dist_weight));
%             if numel(phf)~=0
%                 phf = phf(1);
%             else phf = 1;
%             end           
%         end        
%          peakPhase = round(theta(phf)*180/pi);
%         if peakPhase < 0, peakPhase = peakPhase + 360; end
%            
%         % Calculating mean and std of xw power in significant regions
%         altPhases = find(angle(Wxy_nonzero_lin)>pi/2 | angle(Wxy_nonzero_lin)< -pi/2);
%         synchPhases = find(angle(Wxy_nonzero_lin)<= pi/2 & angle(Wxy_nonzero_lin)>= -pi/2);
%         sblah = abs(Wxy_nonzero_lin(synchPhases));
%         synchPowers = round(sum(sblah(:))*100)/100;
%         ablah = abs(Wxy_nonzero_lin(altPhases));
%         altPowers = round(sum(ablah(:))*100)/100;
%         meanPowers = round(mean(abs(Wxy(:)))*100)/100;
%         stdPowers = round(std(abs(Wxy(:)))*100)/100;
%         totalPowers = round(sum(abs(Wxy(:))));
%         altSynchPowRatio = round((altPowers/synchPowers)*100)/100;
%         if round(totalPowers -(synchPowers + altPowers))>10
%             errordlg('Synch Pow + Alt Pow ~= Tot Pow');
%         end
% %         if cellNum==2
% %             firstPow = totalPowers;
% %         end
% %         totPowRatio = round(100*totalPowers/firstPow)/100;
% %         statMat(:,cellNum)= deal({pf, meanFreqs,...
% %             stdFreqs,peakPhase,meanPhases,...
% %             stdPhases,meanPowers,...
% %             stdPowers,synchPowers,...
% %             altPowers,altSynchPowRatio,totalPowers,...
% %             totPowRatio,num2str(timeRange),num2str(freqRange)});
% %         chLabelMat{cellNum} = ['f' num2str(fileNum) ' ch' num2str(ch(chNum))...
% %             num2str(ch(chNum+1))];
% %         fNamesMat{cellNum} = fNames(fileNum,:);
% %         
        
%% Plotting Figures
figure('Name','XWT','color','w')
figPos = get(gcf,'position');
set(gcf,'position',[figPos(1) figPos(2)-figPos(4)/2 figPos(3)* 1.5...
    figPos(4)+figPos(4)/2]);

%% XWT Axes
ax1 = axes; box off
aPos = get(ax1,'position');
aPos = [aPos(1)-0.02 aPos(2)+ aPos(4)*(1/3) aPos(3)*(0.8) aPos(4)*(0.75)];
set(ax1,'position', aPos)
% plotwave3(Wxy_coi,t,period,coi,sig95,sigmax, sigmay); % Does not plot regions outside COI
% plotwave3(Wxy_sig,t,period,coi,sig95,sigmax, sigmay); % Does not plot
% regions outide COI and below significance
% plotwave3(Wxy_full,t,period,coi,sig95,sigmax, sigmay);
plotwave3(Wxy_ROI,t,period,coi,sig95,sigmax, sigmay);
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
figure('Name','XW power spectrum - linear frequency scale')
plot(freq,normPowerSpectrum)
hold on
 
freq = flipud(freq(:)); % Freq col vector with ascending vals
normPowerSpectrum = flipud(normPowerSpectrum(:));% XW power col vec with ascending vals
df = freq(2)-freq(1); % Smallest freq sampling interval
% df = df/100; % I tried interpolating with finer resolution, but makes
               % little difference
freq2 = freq(1):df:freq(end); % New freq vector with uniform sampling interval = df
nps = spline(freq,normPowerSpectrum,freq2); % New XW power col vec with ascending vals
% plot(freq2,nps,'r:')
hold off
% figure
% % plot(freq,cumsum(normPowerSpectrum))
ld = log2(diff(freq));
% ld = diff(freq);
ld = ld - min(ld);
aps = normPowerSpectrum(2:end).* ld;
% % plot(ld)
%   plot(freq(2:end),aps)
figure ('Name','Cumulative power over uneven frequency intervals')
plot(freq,cumsum(normPowerSpectrum))
figure('Name','Cumulative power over even frequency intervals')
plot(freq2,cumsum(nps))

%% Pending fixes
% 1)It would be nice to not display arrows in insignificant regions
% 2)Ticks on colorbar seem a little screwy