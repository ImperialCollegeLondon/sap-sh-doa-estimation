function[hout] = plot_az_inc(az,inc,varargin)
h = plot(radtodeg(wrapTo2Pi(az)),radtodeg(inc),varargin{:});
if nargout>0
    hout = h;
end