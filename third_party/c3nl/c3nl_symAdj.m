function Adj = c3nl_symAdj(A,type)
%    ___________ _   ____   
%   / ____|__  // | / / /   
%  / /     /_ </  |/ / /    
% / /___ ___/ / /|  / /___  
% \____//____/_/ |_/_____/  
% Eyal Soreq Imperial College London
%
% symmetrize a square matrix using different ways
%
% e.g. Adj = c3nl_symAdj(mat,'mean')
%

sz=size(A);
Adj = zeros(sz);
if numel(sz)==2
    AA =   fliplr(rot90(A,-1));
    ix_l = tril(ones(size(A)),-1)>0;
    ix_u = triu(ones(size(A)),1)>0;
    
    switch type
        case 'lower'; Adj = ix_l.*A+ix_u.*AA;
        case 'upper'; Adj = ix_u.*A+ix_l.*AA;
        case 'mean';  Adj = mean(cat(3,c3nl_symAdj(A,'lower'),c3nl_symAdj(A,'upper')),3);
        otherwise; disp('type is not recognized');return;
    end
elseif numel(sz)==3
    AA =   fliplr(rot90(A,-1));
    ix_l = repmat(tril(ones(sz(1),sz(2),1),-1)>0,1,1,sz(3));
    eye3 = repmat(logical(eye(sz(1),sz(2))),1,1,sz(3));
    %ix_u = repmat(triu(ones(sz(1),sz(2),1),1)>0,1,1,sz(3));
    
    switch type
        case 'mean'  
            Adj(ix_l) = mean([A(ix_l) AA(ix_l)],2);
            Adj = Adj + flipud(rot90(Adj));
            Adj(eye3) = A(eye3);
        otherwise; disp('type is not recognized');return;
    end    
end

end


%{
TEST
tmp = rand(5,5,3);
A=tmp;
Adj=zeros(size(A))

%}



