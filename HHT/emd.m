function imf = emd(x,varargin)
% Empiricial Mode Decomposition (Hilbert-Huang Transform)
% imf = emd(x);
% imf = emd(x,'interp',interp);
% Inputs:
% x - Timeseries on which to perform EMD
% interp - 'spline' or 'cubic' [default]
% nComps - Manually specifying number of components
% Outputs:
% imf - Cell array where each cell holds an imf.
% Required functions: GetPks
% 
% Modified by Avinash Pujala, JRC/HHMI, 2017

nComps = Inf;

for jj = 1:numel(varargin)
    if ischar(varargin{jj})
        switch varargin{jj}          
            case 'ncomps'
                nComps = varargin{jj+1};
        end
    end
end

x   = transpose(x(:));
imf = [];
c = Inf;
x_orig = x;
count_outer = 0;
while ~ismonotonic(x) 
    x1 = x;
    sd = Inf;
    while (sd > 0.1) & ~isimf(x1)        
        %             s1 = getspline(x1);
        %             s2 = -getspline(-x1);
        %             x2 = x1-(s1+s2)/2;        
        x2 = GetCubic(x1);        
        sd = sum((x1-x2).^2)/sum(x1.^2);
        x1 = x2;
    end    
    imf{end+1} = x1;
    x  = x-x1;
    count_outer = count_outer + 1;
end
imf{end+1} = x;

% FUNCTIONS

function u = ismonotonic(x)

u1 = length(GetPks(x))*length(GetPks(-x));
if u1 > 0, u = 0;
else
    u = 1;
end

function u = isimf(x)

N  = length(x);
u1 = sum(x(1:N-1).*x(2:N) < 0);
u2 = length(GetPks(x))+length(GetPks(-x));
if abs(u1-u2) > 1, u = 0;
else
    u = 1;
end

function s = getspline(x)
N = length(x);
p = GetPks(x);
s = spline([0 p N+1],[0 x(p) 0],1:N);
