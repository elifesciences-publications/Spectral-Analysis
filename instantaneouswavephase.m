
function [mean_instantaneous_phase, peakpow_instantaneous_phase,varargout] = instantaneouswavephase(Wxy)
% [mean_instantaneous_phase, peakpow_instantaneous_phase, phaseStd] = instantaneouswavephase(Wxy)

if isreal(Wxy)
    %     errordlg('First input variable must be a matrix of complex wavelet coefficients!')
    %     return
    
end

Wxy_abs = abs(Wxy);

mpv = angle(sum(Wxy)/size(Wxy,1));
mpv(mpv<0) = 2*pi+(mpv(mpv<0));
mean_instantaneous_phase = (360/(2*pi))*mpv;

ppv = angle(max(Wxy));
ppv(ppv<0) = 2*pi+(ppv(ppv<0));
peakpow_instantaneous_phase = (360/(2*pi))*(ppv);

nouts = nargout;
if nouts > 2
    stdPhases = zeros(size(Wxy,1),1);
    for tt = 1:size(Wxy,2)
        dezero = Wxy(:,tt);
        dezero(dezero==0)=[];
        dezero = angle(dezero);
        if isempty(dezero)
            stdPhases(tt) = 0;
        else
            stdPhases(tt) = std(dezero);
        end
    end
    stdPhases(stdPhases<0)= 2*pi + (stdPhases(stdPhases<0));
    stdPhases = (180/pi)*stdPhases;
    varargout{1} = stdPhases';
end
