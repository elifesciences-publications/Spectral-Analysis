function plot_hht2(x,Ts)
% Plot the HHT.
% plot_hht(x,Ts)
% 
% :: Syntax
%    The array x is the input signal and Ts is the sampling period.
%    Example on use: [x,Fs] = wavread('Hum.wav');
%                    plot_hht(x(1:6000),1/Fs);
% Func : emd

% Get HHT.

B = sum(x.*x);
imf = emd2(x);

% Eliminate imfs with low energy (i.e. noise contribution)
inds = [];
for k = 1:length(imf)
    b(k) = sum(imf{k}.*imf{k});  
    if b(k)/B <0.3
      inds = [inds; k];
   end
end
imf(inds)=[];
b = [];
sig = zeros(size(imf{1}));
for k = 1:length(imf)
   b(k) = sum(imf{k}.*imf{k});  
   th   = angle(hilbert(imf{k}));
   d{k} = diff(th)/Ts/(2*pi);
   sig = sig + imf{k};
end
sig = sig/length(k);
imf{k+1} = sig
b(k+1) = sum(sig.*sig);
th = angle(hilbert(sig));
d{k+1} = diff(th)/Ts/(2*pi);

[u,v] = sort(-b);
b     = 1-b/max(b);

% Set time-frequency plots.
N = length(x);
c = linspace(0,(N-2)*Ts,N-1);
for k = v(1:max(2,length(imf)))
   figure, plot(c,d{k},'k.','Color',b([k k k]),'MarkerSize',3);
   set(gca,'FontSize',8,'XLim',[0 c(end)],'YLim',[0 1/2/Ts]); xlabel('Time'), ylabel('Frequency');
end

% Set IMF plots.
M = length(imf);
N = length(x);
c = linspace(0,(N-1)*Ts,N);
for k1 = 0:4:M-1
   figure
   for k2 = 1:min(4,M-k1), subplot(4,1,k2), plot(c,imf{k1+k2}); set(gca,'FontSize',8,'XLim',[0 c(end)]); end
   xlabel('Time');
end
