function S = ComputeLocomotorStrength(varargin)
% ComputeLocomotorStrength - Computes locomotor strengths
% strengths = ComputeLocomotorStrength(allData);
% strengths = ComputeLocomotorStrength(allData,phaseRange)
% Inputs:
% allData - Variable containing all the relevant information created by ComputeAndPlotXWTs
% phaseRange - [a,b]; where a and b must take some value between 0 to 2*pi;
% Outputs:
% strengths - A vector of strength values computed for each file data stored in allData

if nargin == 1
    phaseRange = [pi-pi/3 pi+pi/3];
else
    phaseRange = sort(varargin{2});
end
allData = varargin{1};
phaseRange(phaseRange < 0) = phaseRange(phaseRange<0)+ 2*pi;

data = allData.data;
Wxy = allData.Wxy;

S = nan(size(data,1),1);
for fn = 1:size(data,1)
    W = mean(Wxy.sig(:,:,fn,:),4);
    W(abs(W)==0) =[];
    A = angle(W(:));
    A(A<0) = A(A<0)+ 2*pi;
    W(A < min(phaseRange) | A > max(phaseRange)) = [];
    S(fn) = round(sqrt(abs(sum(W))));
end
