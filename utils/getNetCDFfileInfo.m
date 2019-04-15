% G. Mengaldo (gianmarco.mengaldo@gmail.com)
% O. T. Schmidt (oschmidt@ucsd.edu)
% Last revision: 07-April-2019

function [lon,lat,time,level,fields]=getNetCDFfileInfo(file)

if iscell(file)
    file = char(file(1));
% else
%     file = files;
end

disp(' ')
disp('NetCDF file')
disp('------------------------------------')
disp(['Name                      : ' file])
ncid    = netcdf.open(file, 'NOWRITE');
finfo   = dir(file);
disp(['Size                      : ' num2str(finfo.bytes/1e9,'%.2f GB')])

% Get information about the contents of the file.
[ndims, nvars, ngattributes, unlimdimID] = netcdf.inq(ncid);

% Load longitude, latitude and time
disp(' ')
disp('List of variables')
disp('------------------------------------')
vars{nvars}     = [];
fields{nvars-3} = [];
fieldi          = 1;
level           = 0;
for varid = 0:nvars-1
    [vars{varid+1}, xtype, dimids, numatts] ...
        = netcdf.inqVar(ncid, varid);    
    disp(vars{varid+1})
    if strcmp(vars{varid+1},'longitude')
        longitude   = ncread(file,vars{varid+1});
    elseif strcmp(vars{varid+1},'latitude')
        latitude    = ncread(file,vars{varid+1}); 
    elseif strcmp(vars{varid+1},'time')
        time    = ncread(file,vars{varid+1});
    elseif strcmp(vars{varid+1},'level')
        level   = ncread(file,vars{varid+1});        
    else
        fields{fieldi}  = vars{varid+1};
        fieldi          = fieldi+1;
    end
end
netcdf.close(ncid);

n_lon       = length(longitude);
n_lat       = length(latitude);
n_time      = length(time);
n_level     = length(level);

% cast to float
lon         = cast(longitude,'double');
lat         = cast(latitude,'double');
time        = cast(time,'double');
level       = cast(level,'double');
dt          = time(2)-time(1);
disp(' ')
disp('Dimensions')
disp('------------------------------------')
disp(['Longitude                 : ' num2str(n_lon) ' points in [' num2str(lon(1)) ','  num2str(lon(end)) ']']);
disp(['Latitude                  : ' num2str(n_lat) ' points in [' num2str(lat(1)) ','  num2str(lat(end)) ']']);
disp(['Number of time steps      : ' num2str(n_time)])
disp(['Number of levels          : ' num2str(n_level)])
disp(['Time step                 : ' num2str(dt) ' hours'])
disp(['Start date                : ' datestr(double(time(1))/24 + datenum(1900,1,1))]);
disp(['End date                  : ' datestr(double(time(end))/24 + datenum(1900,1,1))]);

