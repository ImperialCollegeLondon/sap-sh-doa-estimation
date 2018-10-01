function[smoothed_hist] = smooth_histogram_pd(counts,sigma,azg,incg)
% sigma: standard deviation of normal distribution filter kernel [degrees]

persistent weightmat
persistent pd

if ~isempty(pd)
    if ~isequal(pd.sigma,sigma) || ~isequal(pd.azg,azg) || ~isequal(pd.incg,incg)
        pd = [];
    end
end
if isempty(pd) || isempty(weightmat)
    distmat = az_inc_bin_distances(azg,incg);
    weightmat = normpdf(real(distmat)*180, 0,sigma);
    weightmat = bsxfun(@rdivide,weightmat,sum(weightmat,2));
    weightmat(weightmat<0.001*max(weightmat(:))) = 0;
    weightmat = sparse(weightmat);
    
    pd.azg = azg;
    pd.incg = incg;
    pd.sigma = sigma;
end


[nInc,nAz] =  size(azg);
smoothed_hist = zeros(nInc,nAz);
for i_az = 1:nAz
    tmp = circshift(weightmat,[0,(i_az-1)*nInc]);
    smoothed_hist(:,i_az) = tmp * counts(:);
end