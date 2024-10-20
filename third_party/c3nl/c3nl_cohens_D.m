function d = c3nl_cohens_D(X,I)
%    ___________ _   ____   
%   / ____|__  // | / / /   
%  / /     /_ </  |/ / /    
% / /___ ___/ / /|  / /___  
% \____//____/_/ |_/_____/  
% Richard E. Daws APR19
%
% Cohen's d = (M2 - M1) / SDpooled
%
% d = c3nl_cohens_D(X)
%
% Cohen's d = (M2 - M1) ? SDpooled
% SDpooled = ?((SD12 + SD22) ? 2)
%
% Paired Sample
% X - An nx2 matrix
% One Sample
% X - An nx1 matrix
%
% I - Grouping index for 2 samples
%%
if ~exist('I','var')
    if size(X,2)==2
        d = diff(mean(X)) / std(diff(X,[],2)); % Difference in mean divided by std of the change
    end

    if size(X,2)==1
        d = mean(X) / sqrt(mean(std(X).^2));
    end
else
    G = unique(I);
    d = diff([mean(X(I==G(1))) mean(X(I==G(2)))]) / sqrt(mean(std(X).^2));
end
