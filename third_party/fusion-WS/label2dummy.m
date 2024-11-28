function [D,nc] = label2Dummy(L)
dim = size(L);
L = double(L(:));
L((isnan(L)))=0;
nc = unique(L(L~=0));
nc(nc==0) = [];
D = false([size(L(:),1),numel(nc)]);
for ii=1:numel(nc)
    D(:,ii) = L==nc(ii);
end
D = reshape(D,[dim,numel(nc)]);

end
