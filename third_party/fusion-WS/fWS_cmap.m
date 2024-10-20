function map = fWS_cmap(color,n,w)
%% GEN.CMAP: creates a gradient colormap
%
%   __           _             
%  / _|         (_)            
% | |_ _   _ ___ _  ___  _ __    
% |  _| | | / __| |/ _ \| `_ \    :- Functional and Structural 
% | | | |_| \__ \ | (_) | | | |      Integration of Neuroimages
% |_|  \__,_|___/_|\___/|_| |_|
%
%
%% AUTHOR:  Eyal Soreq
%  EMAIL:  e.soreq14@imperial.ac.uk
%  AFFILIATION:  Imperial College London
%  VERSION:  0.0 CREATED: 10-Jun-2017 09:28:18
%
%% INPUTS:
%    color - at least two colours are necessery 
%    n -  number of colours needed
%    weights - interpulation weights
%
%
%% OUTPUT:
%
%% EXAMPLES:
%{
map = gen.cmap([244,235,37;251,198,85;0,0,0;30,152,213;70,64,153])
figure;imagesc(-6:0.1:6);
colormap(map)
%}
%
%% DEPENDENCIES:
%
% This file is part of Fusion Pipeline
% Fusion Pipeline is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% Fusion Pipeline is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with Fusion Pipeline.If not, see <http://www.gnu.org/licenses/>.
%------------- BEGIN CODE --------------
%

if nargin < 2; n = 64;end

if ischar(color) && strcmp(color, 'soreq')
    map = get.cmap(flipud([255 205 155;255 180 104;255 155 54;255 130 4;191 98 3;51 51 51;14 83 164;49 92 143;62 117 182;110 152 200;207 221 237]./255),n);
else
    if nargin < 3; w = linspace(0,1,size(color,1));end
    if (range(min(color(:),max(color(:))))>1||any(color(:))>1)
        color = c3nl.scale(color);
    end
    map = interp1(w,color,linspace(0,1,n));
end



end
%------------- END OF CODE --------------
