function[est_doa_vec,nFoundSrcs,found_peak_vals] = sph_local_peaks_to_doa_vecs_with_peak_vals(X,azg,incg,nSrcs,min_separation_deg)

[iinc,iaz] = spherical_find_local_maxima(X);
est_doa_vec = NaN(nSrcs,3);
found_peak_vals = NaN(nSrcs,1);
nFoundSrcs = min(length(iinc),nSrcs);
peak_vals = X(sub2ind(size(X),iinc,iaz));
[~,sort_order] = sort(peak_vals,'descend');
idc = sub2ind(size(X),iinc(sort_order(1:nFoundSrcs)),iaz(sort_order(1:nFoundSrcs)));
%radtodeg([azg(idc),incg(idc)])
[tmpx,tmpy,tmpz] = mysph2cart(azg(idc),incg(idc),ones(nFoundSrcs,1));
est_doa_vec(1:nFoundSrcs,:) = [tmpx,tmpy,tmpz];
found_peak_vals(1:nFoundSrcs) = peak_vals(sort_order(1:nFoundSrcs));


% [ninc,naz] = size(X);
%
% buf = zeros(ninc+2,naz+2);
% buf(2:end-1,2:end-1) = X;
% buf(2:end-1,1) = X(:,end);
% buf(2:end-1,end) = X(:,1);
% buf(1,2:end-1) = fliplr(X(1,:));
% buf(end,2:end-1) = fliplr(X(end,:));
%
% [iinc,iaz] = find( X >= buf(1:end-2,2:end-1) & ... %compare to the top
%                    X >= buf(3:end,2:end-1) & ...  %...bottom
%                    X >= buf(2:end-1,1:end-2) & ... %...left
%                    X >= buf(2:end-1,3:end) ); %...right
%
% %check for adjacent points forming a plateau
% keyboard
%
%
% peak_vals = X(sub2ind(size(X),iinc,iaz));
% [sorted_peaks,sort_order] = sort(peak_vals,'descend');
% found_peaks = 0
% idc = 1;
% peak_diff = diff(full(sorted_peaks));
% i_equal = find(peak_diff==0);
% while found_peaks<max_peaks & idc<length(sort_order)
%    this_peak = sorted_peaks(idc);
%    equal_peaks = pea
%
%
%
% end

