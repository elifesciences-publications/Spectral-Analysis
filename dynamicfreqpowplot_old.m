
if avgCheck == 1
    figure('Name', ['Dynamic Freq and Pow, AVERAGED , Channels '...
                num2str(ch(chNum)) ' vs ' num2str(ch(chNum+1))],'color','w')
else
figure('Name', ['Dynamic Freq and Pow, File ' num2str(fileNum) ', Channels '...
                num2str(ch(chNum)) ' vs ' num2str(ch(chNum+1))],'color','w')
end
title('Dynamic Freq and Pow')
plot(time_reduced(1:10:end),mfvec(1:10:end),'-k.','markersize',10),box off,
% plot(tempTime(1:10:end),mfvec(1:10:end),'-k.','markersize',20),box off,
ylabel('Freq'),xlabel('Time')
hold on

plot(time_reduced(1:10:end),pfvec(1:10:end),'.r-','markersize',10)
% plot(tempTime(1:10:end),pfvec(1:10:end),'.r-','markersize',20)
set(gca,'tickdir','out'),shg
title('Dynamic Freq and Pow')


plot(time_reduced(1:10:end),tvpower(1:10:end),'.b-','markersize',2)
% plot(tempTime(1:10:end),powvec(1:10:end),'ob-','markersize',2)
legend('Mean freq','Max-pow freq','Power','Location','Best')
xlim([timeRange])
ylim([min(freq) max(freq)])


%% ####### PLOTTING DYNAMIC MEAN AND MAX-POWER PHASES #######

[meanPhaseVec,maxPhaseVec,stdPhaseVec] = instantaneouswavephase(Wxy);

if avgCheck == 1
    figure('Name', ['Dynamic Mean and Max-Pow Phases, AVERAGED, Channels '...
                num2str(ch(chNum)) ' vs ' num2str(ch(chNum+1))],'color','w')
else
figure('Name', ['Dynamic Mean and Max-Pow Phases, File ' num2str(fileNum) ', Channels '...
                num2str(ch(chNum)) ' vs ' num2str(ch(chNum+1))],'color','w')
end

plot(time_reduced, meanPhaseVec,'k.','linewidth',2)
hold on,
plot(time_reduced,maxPhaseVec,'r.','linewidth',2)
legend('Mean Phases','Max-Pow Phases','location','best')
plot(time_reduced, 180*ones(size(meanPhaseVec)),'k:')
plot(time_reduced, 180*ones(size(meanPhaseVec))+60,'r:')
plot(time_reduced, 180*ones(size(meanPhaseVec))-60,'r:')
title('Dynamic Mean and Max-Pow phases')
box off, axis([-inf inf -inf inf]), set(gca,'tickdir','out','ytick',[0:60:360],'ydir','reverse')
xlabel('Time (sec)'), ylabel('Phase (Deg)')
hold off
xlim(timeRange)
ylim([0 360])

% arrowGap = round(length(uu)/40);
%  figure('Name', ['Dynamic Mean Phases, File ' num2str(fileNum) ', Channels '...
%                 num2str(ch(chNum)) ' vs ' num2str(ch(chNum+1))],'color','w')
% 

% quiver(time_reduced(1:arrowGap:end),mfvec(1:arrowGap:end),...
%     u(1:arrowGap:end),v(1:arrowGap:end),0.25,'color',[0 0 0])
% hold on,
% quiver(time_reduced(1:arrowGap:end),pfvec(1:arrowGap:end),...
%     uu(1:arrowGap:end),vv(1:arrowGap:end),0.25,'color',[1 0 0])
% title('Dynamic Mean Phases in Frequency-Time Space')
% legend('Mean Phases','Max-Pow Phases')
% set(gca,'tickdir','out'), box off
% xlim([timeRange])
% ylim([-inf inf])
% xlabel('Time (sec)')
% ylabel('Frequency (Hz)')
% hold off
% 
