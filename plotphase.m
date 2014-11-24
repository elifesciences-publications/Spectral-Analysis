% PLOTPHASE Plots rose diagrams for phases obtained from XWT in
% XWFREQSPECPLOT. Displays both weighted(weighted by log2(XW power)) and
% unweighted rose diagrams
% XWFREQSPECPLOT must always be run before to load the requisite variable
% into the workspace
% Author: AP
% Last modified 22-Mar-2011 18:52:43



%% Display parameters
border_toggle = 'on';
ttdelta = 90; % Spacing in degrees between angular ticks
tl = [0:ttdelta:360];
ttl= [];
for jj = 1:length(tl)
    ttl{jj} = [num2str(tl(jj)) '\circ'];
    ttl = ttl(:);
end


Wxy_nonzero_lin = Wxy_avg_files(:,:,1,cp);
Wxy_nonzero_lin(Wxy_nonzero_lin==0)=[]; % This removes all zero elements from the matrix and vectorizes it.
Axy = angle(Wxy_nonzero_lin);
nPhaseBins = min([numel(Axy), 90]); % Number of bins circular for phase histograms
[ph_dist,th] = hist(Axy(:),nPhaseBins); % Unweighted phase histogram
mag = abs(Wxy_nonzero_lin);
[ph_dist3,vals] = hist3([Axy(:) mag(:)],[nPhaseBins,nPhaseBins]); % 3D bivariate (phases and power) histogram
powmat = repmat(vals{2},size(ph_dist3,1),1);
ph_dist_wt = ph_dist3.*powmat; % Scaling the number of elements in each bin by power (power-weighting)
ph_dist_wt = sum(ph_dist_wt,2)'; % Power-weighted phase histogram

mphase = angle(sum(Wxy_nonzero_lin)); % Mean phase is the angle of the resultant vector obtained by summing all the wavelet coefficients
mphase(mphase<0) = mphase(mphase<0)+ 2*pi; % Addition of 2*pi to -ve values converts angle range from 0 to 360 rather than -180 to +180
meanPhases.avg = mphase*180/pi; % Converts radians to degrees.
meanPhases.avg = round(meanPhases.avg*100)/100;
sphase = circ_std(Axy(:),abs(Wxy_nonzero_lin(:)));
stdPhases.avg = sphase*180/pi;
stdPhases.avg = round(stdPhases.avg*100)/100;

if isempty(ph_dist)
    theta.avg = zeros(size(ph_dist_wt));
    phf = 1;
else
    ph_dist = [ph_dist(:); ph_dist(1)]; % This will close the loop in polar plot by circularizing the vector.
    ph_dist_wt = [ ph_dist_wt(:); ph_dist_wt(1)];
    phase_dist.avg = ph_dist./max(ph_dist);
    phase_dist_weight.avg = ph_dist_wt./max(ph_dist_wt);
    theta.avg = [th(:);th(1)];
    phf = find(phase_dist_weight.avg == max(phase_dist_weight.avg));
    if numel(phf)~=0
        phf = phf(1);
    else phf = 1;
    end
end
peakPhase.avg = round(theta.avg(phf)*180/pi);
if peakPhase.avg <0, peakPhase.avg = peakPhase.avg + 360; end

rtl = [max(phase_dist_weight.avg)/2 max(phase_dist_weight.avg)];

figure('Name', ['Rose Diagram, AVERAGED, Channels '...
    num2str(ch(cp)) ' vs ' num2str(ch(cp+1))],'color','w')

handle = mmpolar(theta.avg,phase_dist_weight.avg,...
    'b-','border',border_toggle,'ttickdelta',ttdelta,...
    'rtickvalue',rtl,'border','on','tticklabel',ttl);
set(gcf,'color','w'), set(handle,'linewidth',2)
hold on
handle = mmpolar(theta.avg,phase_dist.avg,...
    'color',[0 0.5 0],'linestyle',':','border',border_toggle,'ttickdelta',ttdelta,...
    'rtickvalue',rtl,'border','on','tticklabel',ttl);
set(gcf,'color','w'), set(handle,'linewidth',2,'markersize',15)

mp = meanPhases.avg*pi/180;
pol = polar(mp*ones(1,2),[0 1]);
set(pol,'color','r','linewidth',2,'linestyle','--','linewidth',2)
legend('Weighted', 'Unweighted','Mean Phase')
hold off



