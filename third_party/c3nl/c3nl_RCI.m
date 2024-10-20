function CI = c3nl_RCI(R,N,oneTailed)
%    ___________ _   ____   
%   / ____|__  // | / / /   
%  / /     /_ </  |/ / /    
% / /___ ___/ / /|  / /___  
% \____//____/_/ |_/_____/  
% Richard E. Daws NOV 2019
% 
% Calculate 95% confidence interval of a person correlation coefficient
%

if exist('oneTailed','var') && oneTailed
    thr = 1.65;
else
    thr = 1.96;
end

Z=c3nl_fisher(R, 'r2z');% Convert r to z
SE=1/sqrt(N-3);% pull Z standard error

% Compute a confidence interval in terms of z & convert to r
CI=c3nl_fisher([Z - (thr * SE) Z + (thr * SE)], 'z2r');

end
