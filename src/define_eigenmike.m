function[mic] = define_eigenmike()
%DEFINE_EIGENMIKE creates a structure representing an EigenMike microphone
%
%[mic] = define_eigenmike()
%
%mic has the following fields:
%              a: radius of sphere in metres
%  sensor_angles: [32 x 2] matrix specifying angle of each transducer in
%                 Daniel's [azimuth, inclination] coordinates system
%        sphType: 'rigid' or 'open'
%
%requires: sphDist from SPHERICAL_TOOLBOX
%

mic.a = 0.042;
[mic.sensor_angle(:,1),mic.sensor_angle(:,2),mic.quad] = sphDist('em32');
mic.sphType = 'rigid';