% O. T. Schmidt (oschmidt@ucsd.edu)
% Last revision: 07-April-2019

function [dS] = trapzWeightsSpherical_2D(lon,lat,R)
    
n_lon       = length(lon);
n_lat       = length(lat);

lon         = linspace(0,360,n_lon+1);  
lon         = lon(1:end-1);
lat         = linspace(-90,90,n_lat);

lon_rad     = lon/360*2*pi;
lat_rad     = lat/360*2*pi;

d_lon       = [lon_rad(1)/2 diff(lon_rad)];

diff_lat    = diff(lat_rad);
diff_lat    = diff_lat(1:end-1);
d_lat       = [diff_lat(1)/2 diff_lat diff_lat(end)/2];

d_lon       = repmat(d_lon',[1 n_lat]);
d_lat       = repmat(d_lat, [n_lon 1]);

% just for meshgrid format
d_lon       = d_lon';
d_lat       = d_lat';

% cos(latitude) since lat \in [-90 90] deg 
dS          = abs(R^2*cos(lat_rad').*d_lon.*d_lat);
