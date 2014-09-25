%% DYNAMICFREQPOWPLOT - Plots dynamic frequencies, power, and phase

 figure('Name', ['Dynamic Freq and Pow, AVERAGED , Channels '...
                num2str(ch(cp)) ' vs ' num2str(ch(cp+1))],'color','w')

title('Dynamic Frequency & Power')
plot(time(1:10:end),tvmf(cp,1:10:end),'-k.','markersize',10),box off,
% plot(tempTime(1:10:end),mfvec(1:10:end),'-k.','markersize',20),box off,
xlabel('Time (sec)')
hold on

plot(time(1:10:end),tvpf(cp,1:10:end),'r.','markersize',10)
% plot(tempTime(1:10:end),pfvec(1:10:end),'.r-','markersize',20)
set(gca,'tickdir','out'),shg
title('Dynamic Frequency and Power')

blah = tvpow(cp,:);
scalingFactor = max(blah(:));
plot(time(1:10:end),tvpow(cp,1:10:end)/scalingFactor,'.b-','markersize',2,'linewidth',1.5)
legend('Mean freq','Max-pow freq','Power','Location','Best')
axis([time(1) time(end) 0 max(freq)])
ylabel('Frequency (Hz)/\color{blue}Power (Max-Normalized)')
%% ####### PLOTTING DYNAMIC MEAN AND MAX-POWER PHASES #######

figure('Name', ['Dynamic Mean and Max-Pow Phases, AVERAGED, Channels '...
                num2str(ch(cp)) ' vs ' num2str(ch(cp+1))],'color','w')

plot(time, tvmph(cp,:),'k.','linewidth',2)
hold on,
plot(time,tvpph(cp,:),'r.','linewidth',2)
plot(time, 180*ones(size(tvmph(cp,:))),'k:')
plot(time, 180*ones(size(tvmph(cp,:)))+60,'r:')
plot(time, 180*ones(size(tvmph(cp,:)))-60,'r:')
title('Dynamic Mean and Max-Pow phases')
box off, axis([-inf inf -inf inf]), set(gca,'tickdir','out','ytick',[0:60:360],'ydir','reverse')
xlabel('Time (sec)'), ylabel('Phase (Deg)')
hold off
xlim([time(1) time(end)])
ylim([0 360])
legend('Mean Phases','Max-Pow Phases','location','best')

% arrowGap = round(length(uu)/40);
%  figure('Name', ['Dynamic Mean Phases, File ' num2str(fileNum) ', Channels '...
%                 num2str(ch(cp)) ' vs ' num2str(ch(cp+1))],'color','w')
% 

% quiver(time(1:arrowGap:end),mfvec(1:arrowGap:end),...
%     u(1:arrowGap:end),v(1:arrowGap:end),0.25,'color',[0 0 0])
% hold on,
% quiver(time(1:arrowGap:end),pfvec(1:arrowGap:end),...
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
