
function totpow = instantaneouswavepow(Wxy)

% 
% if isreal(Wxy)
% %     errordlg('First input variable must be a matrix of complex wavelet coefficients!')
% %     return
% 
% end

Wxy_abs = abs(Wxy);

totalpow_norm = sum(Wxy_abs,1);
totalpow_norm = totalpow_norm/max(totalpow_norm);

totpow = totalpow_norm;
