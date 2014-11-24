% % 
% % 
% % signal = signal1;
% % close all
% % winSize = 30; % In seconds
% % nPts = floor(winSize/samplingInt);
% % winShift = 1; % In seconds
% % wsPts = round(1/samplingInt);
% % ls = length(signal);
% % for jj = 1:wsPts:length(signal)-nPts
% %     plot(time(jj:jj+nPts),signal(jj:jj+nPts,1),'b',time(jj:jj+nPts),signal(jj:jj+nPts,2),'r','linewidth',2)
% %     axis([-inf inf -15 15])
% %     drawnow
% %     pause
% end
%    
for fn = 1:nFiles
  plotxwt(W_gm(:,:,fn),time_reduced,freq,coi),
 pause
end
for fn = 1:nFiles
 plotxwt(Wxy_avg_channels(:,:,fn),time_reduced,freq,coi)
 set(gcf,'Name',['Arithmetic mean ' num2str(fn)])
 pause
end