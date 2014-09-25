function np = pnbt(vals)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
ncols = size(vals,2);
% % vals(:,ncols+1) = log2(vals(:,ncols));
% vals(:,ncols+1) = 1*(vals(:,ncols));
n1 = find(vals(:,1)==1);
if isempty(n1);
    thr = min(vals(:,1));
    if ~exist('IS_mat','var')
    global IS_mat
%     errordlg('Please load IS_mat');
    % return
%     load LUT % loads LUT.mat which contains the variable IS_mat;
    end
    ind = find(IS_mat(:,1)==thr);
    scalingFactor  = mean(IS_mat(ind,2:3));
else thr = 1;
    scalingFactor = 1;
end

indices = find(vals(:,1)==thr & (vals(:,2)>=4 & vals(:,2)<=5));
muvals = mean(vals(indices,ncols));


np = vals(:,ncols)/muvals;
np = np*scalingFactor;
np(abs(np)==inf) = 0;
np = (round(np*100))/100;
np = [vals np(:)];

end


