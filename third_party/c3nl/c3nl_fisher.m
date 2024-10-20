function y = c3nl_fisher(x,type)
%    ___________ _   ____   
%   / ____|__  // | / / /   
%  / /     /_ </  |/ / /    
% / /___ ___/ / /|  / /___  
% \____//____/_/ |_/_____/  
% Richard Daws Imperial
% Fisher transform
%
% z = .5[ln(1+r)-ln(1-r)];
%
% e.g. y = c3nl_fisher(x, 'r2z')
%   
%
y=zeros(size(x));
if strcmp(type,'r2z')
    for ii = 1:numel(x)
        y(ii) = .5*[log(1+x(ii))-log(1-x(ii))];
    end
elseif strcmp(type,'z2r')
    for ii = 1:numel(x)
        y(ii) = ((exp(2 * x(ii))) - 1) / ((exp(2 * x(ii))) + 1);
    end
end
