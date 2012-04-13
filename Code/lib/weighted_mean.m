function [wmean] = weighted_mean(x, w)

wmean = nansum(w.*x)/nansum(w);

end