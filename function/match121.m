function [match, iA, uMatch] = match121(ssmatch,idx0,idx1)
%   Detailed explanation goes here
% Find neighbors of unique match in original image idx
if ~isempty(ssmatch)
    
    [u,ia,~] = unique(ssmatch);
    n = histc(ssmatch,u);
    n = n==1;
    ia = nonzeros(ia.*n);
    iA = idx0(ia);
    
    match = ssmatch(ia);
else
    match = [];
    iA = [];
end

% match = zeros(size(ssmatch));
% match(ia) = ssmatch(ia);
uMatch = 0;

end