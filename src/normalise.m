function[out,gain] = normalise(in,nbits)
%NORMALISE  normalises each column of a matrix to +/- 1
%
%[out] = normalise(in)
%
%[out] = normalise(in,nbits)
%
%Alastair Moore, September 2006

if nargin > 1 && ((nbits==16) || (nbits==24))
    m = 2^(nbits-1); %gives 2^15 or 2^23
    limit_max = (m - 1)/m;
else
    limit_max = 1;
end
limit_min = -1;

%normalise
max_val = max(in(:));
min_val = min(in(:));
scale_pos = limit_max/max(max_val,eps);
scale_neg = limit_min/min(min_val,-eps);
gain = min(scale_pos,scale_neg);
out = gain * in;
