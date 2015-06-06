function P=ar1spectrum_ap(ar1,period)
% AR1 power spectrum... 
%
% power=ar1spectrum_ap(ar1,period)
%
%
% (c) Aslak Grinsted 2002-2004
% Modified by Avinash Pujala to make eqns look more like those seen in 
% Torrence & Compo (1998), although they are equivalent


% -------------------------------------------------------------------------
%   Copyright (C) 2002-2004, Aslak Grinsted
%   This software may be used, copied, or redistributed as long as it is not
%   sold and this copyright notice is reproduced on each copy made.  This
%   routine is provided as is without any express or implied warranties
%   whatsoever.


freq=1./period;

% P0=(1-ar1.^2)./(abs(1-ar1.*exp(-2*pi*i*freq))).^2; % Eqn 3 from Grinsted et al, 2004 
%http://www.madsci.org/posts/archives/may97/864012045.Eg.r.html
%fixed typo in numerical recipes

N = round(max(freq)*2);
P=(1-ar1.^2)./(1+ar1.^2-2*ar1*cos(2*pi*freq/N)); % Eqn 16 from Torrence & Compo, 1998
% The two above eqns should be equivalent, but I prefer the latter for
% reasons of simplicity (-Avinash Pujala)