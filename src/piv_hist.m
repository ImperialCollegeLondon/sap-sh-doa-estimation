function[counts,az_grid,inc_grid,extras] = piv_hist(pivs,az_grid,inc_grid)



%% deal with inputs
if nargin<2 || isempty(az_grid)
    %use default grid
    az_grid=2;
    inc_grid=2;
end
if nargin<3 || isempty(inc_grid)
    %use default grid
    inc_grid=2;
end

if numel(az_grid)==1
    az_grid = pi/180 .* [0:az_grid:360-az_grid];%wrapTo360([-180:1:179]);
end
if numel(inc_grid)==1
    inc_grid = pi/180 .* [0:inc_grid:180];

end
if isvector(az_grid) && isvector(inc_grid)
    [az_grid,inc_grid] = meshgrid(az_grid,inc_grid);
end

%%
nAz = size(az_grid,2);
nInc = size(az_grid,1);

[az_final,inc_final,r] = mycart2sph(pivs);
%az_final(toremove) = [];
%inc_final(toremove) = [];

az_bin = round(wrapTo2Pi(az_final)./(2*pi) * nAz)+1;
az_bin(az_bin>nAz) = az_bin(az_bin>nAz) - nAz;
inc_bin = round(inc_final./pi * (nInc-1))+1;

toremove = find(r==0 | isnan(r));
az_bin_safe = az_bin; az_bin_safe(toremove) = [];
inc_bin_safe = inc_bin; inc_bin_safe(toremove) = [];

nPoints = size(az_bin_safe,1);

%counts = sparse(az_bin,inc_bin,ones(nPoints,1),nAz,nInc);
counts = sparse(inc_bin_safe,az_bin_safe,ones(nPoints,1),nInc,nAz);

if nargout<1 
imagesc(180/pi*az_grid(1,:),180/pi*inc_grid(:,1),counts)
set(gca,'ydir','normal')
colorbar
    xlabel('Azimuth [degrees]');
    ylabel('Inclination [degrees]');
end

if nargout>3
    % return the raw pivs with their associated bin assignments
    extras.pivs = pivs;   extras.pivs(r==0,:)=NaN;
    extras.az_bin = az_bin; extras.az_bin(r==0)=NaN;
    extras.inc_bin = inc_bin; extras.inc_bin(r==0)=NaN;
end