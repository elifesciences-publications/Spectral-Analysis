function imf = emd2(x)
% Empiricial Mode Decomposition (Hilbert-Huang Transform)
% imf = emd(x)
% Func : findpeaks_hht2

peakThresh = mean(x) + 1*std(x);
x   = transpose(x(:));
imf = [];
while ~ismonotonic(x,peakThresh)
   x1 = x;
   sd = Inf;
   while (sd > 0.1) | ~isimf(x1,peakThresh)
      s1 = getspline(x1,peakThresh);
      s2 = -getspline(-x1,peakThresh);
      x2 = x1-(s1+s2)/2;      
      sd = sum((x1-x2).^2)/sum(x1.^2);
      x1 = x2;
   end
   
   imf{end+1} = x1;
   x          = x-x1;
end
imf{end+1} = x;

% FUNCTIONS

function u = ismonotonic(x,peakThresh)

u1 = length(findpeaks_hht2(x,peakThresh))*length(findpeaks_hht2(-x,peakThresh));
if u1 > 0, u = 0;
else,      u = 1; end

function u = isimf(x,peakThresh)

N  = length(x,peakThresh);
u1 = sum(x(1:N-1).*x(2:N) < 0);
u2 = length(findpeaks_hht2(x,peakThresh))+length(findpeaks_hht2(-x,peakThresh));
if abs(u1-u2) > 1, u = 0;
else,              u = 1; end

function s = getspline(x,peakThresh)

N = length(x);
p = findpeaks_hht2(x,peakThresh);
s = spline([0 p N+1],[0 x(p) 0],1:N);
