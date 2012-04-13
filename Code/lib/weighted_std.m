function [wstd] = weighted_std(x, w)

wmean = nansum(w.*x)/nansum(w);

N = length(find(w~=0)); %number of non zero weights
wstd = sqrt( nansum(w.*((x-wmean).^2))/( (N-1)*nansum(w)/N) );
end