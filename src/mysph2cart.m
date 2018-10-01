function [x,y,z] = mysph2cart(az,inc,r)

%   Spherical Coordinates to 3D Cartesian
%
%   [x,y,z] = mysph2cart(az,inc,r)
%
%   Inputs:
%       az      Azimuth angle in radians, counterclockwise from xy plane 
%               from the positive x axis (otherwise referred to as phi)
%       inc     Inclination angle in radians, from positive z axis 
%               (otherwise referred to as theta) 
%       r       Radius
%
%   Outputs:
%       x       x-coordinate
%       y       y-coordinate
%       z       z-coordinate
%
%   Notes:
%       The MATLAB function cart2sph reverses phi and theta.
%
%**************************************************************************
% Author:           E. A. P. Habets, M. R. P. Thomas and D. P. Jarrett
% Date:             27 July 2010
% Version: $Id: mysph2cart.m 351 2011-07-21 15:10:09Z dpj05 $
%**************************************************************************

z = r .* cos(inc);
rcosinc = r .* sin(inc);
x = rcosinc .* cos(az);
y = rcosinc .* sin(az);