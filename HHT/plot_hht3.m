function plot_hht3(x,Ts)
% Plot the HHT.
% plot_hht(x,Ts)
% 
% :: Syntax
%    The array x is the input signal and Ts is the sampling period.
%    Example on use: [x,Fs] = wavread('Hum.wav');
%                    plot_hht(x(1:6000),1/Fs);
% Func : emd

% Get HHT.

freqRange = [10 100];
dj = 1;

imf = emd(x);
I = cell2mat(imf);
for k = 1:length(imf)
   b(k) = sum(imf{k}.*imf{k});
   th   = angle(hilbert(imf{k}));
   d{k} = diff(th)/Ts/(2*pi);
end
D = cell2mat(d);
D(D <= 0)=[];
D = round(D/dj)*dj;
D = sort(unique(D),'ascend');
D_ind = round(D/dj)+1;
% D_ind = round(D/dj);
freqVec  = linspace(min(D), max(D), numel(D));
freqInd = 1:length(freqVec);
timeInd = 1:length(imf{k})-1;

hsMat = zeros(D_ind(end), length(timeInd),length(imf)-3);
for kk = 1:length(imf)-3;
    blah = d{kk};
    temp = imf{kk};
    temp = temp(1:end-1);
    blah(blah <= 0)= nan;
%     blah(blah < 0)=[];
    blah = round(blah/dj)*dj;
%       blah(blah < freqRange(1)) = nan;
%     blah(blah > freqRange(2)) = nan;
    blah_ind = round(blah/dj) +1;
    blah_ind(isnan(blah_ind)) = D_ind(end);
    tempMat = sparse(blah_ind,timeInd,temp.^2);
    hsMat(:,:,kk) = full(tempMat);        
end 
hES = sum(hsMat,3);
N = length(x);
time = linspace(0,(N-1)*Ts,N);
freqVec  = D(1):dj:D(end);

figure('Name', 'Hilbert Energy Spectrum')
imagesc(time(1:end-1),freqVec,hES)
set(gca,'clim',[0 0.5*max(hES(:))],'ydir','normal')

% [u,v] = sort(-b);
% b     = 1-b/max(b);
% 
% % Set time-frequency plots.
% N = length(x);
% c = linspace(0,(N-2)*Ts,N-1);
% for k = v(1:2)
%    figure, plot(c,d{k},'k.','Color',b([k k k]),'MarkerSize',3);
%    set(gca,'FontSize',8,'XLim',[0 c(end)],'YLim',[0 1/2/Ts]); xlabel('Time'), ylabel('Frequency');
% end
% 
% % Set IMF plots.
% M = length(imf);
% N = length(x);
% c = linspace(0,(N-1)*Ts,N);
% for k1 = 0:4:M-1
%    figure
%    for k2 = 1:min(4,M-k1), subplot(4,1,k2), plot(c,imf{k1+k2}); set(gca,'FontSize',8,'XLim',[0 c(end)]); end
%    xlabel('Time');
% end
