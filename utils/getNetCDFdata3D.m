% O. T. Schmidt (oschmidt@ucsd.edu)
% Last revision: 07-April-2019

function [dat] = getNetCDFdata3D(files,variable,ti,nt_skip,lvl_start,lvl_skip,replaceNaN)

if iscell(variable)
    for vari = 1:length(variable)
        thisvar         = char(variable(vari));
        if iscell(files)
            file = char(files(vari));
        else
            file = files;
        end
        % [lon, lat, level, time]
        dat(:,:,:,vari)   = ncread(file, thisvar, [1 1 lvl_start nt_skip*(ti-1)+1], [Inf Inf Inf 1], [1 1 lvl_skip 1]);
    end
else
    dat = ncread(files, variable, [1 1 lvl_start nt_skip*(ti-1)+1], [Inf Inf Inf 1], [1 1 lvl_skip 1]);
end
dat = cast(dat, 'single');

if exist('replaceNaN','var')
    dat(isnan(dat))  =  replaceNaN;
end
