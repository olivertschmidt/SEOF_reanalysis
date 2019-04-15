% G. Mengaldo (gianmarco.mengaldo@gmail.com)
% O. T. Schmidt (oschmidt@ucsd.edu)
% Last revision: 07-April-2019

function [dat] = getNetCDFdata(files,variable,ti,nt_skip,replaceNaN)

if iscell(variable)
    for vari = 1:length(variable)
        thisvar         = char(variable(vari));
        if iscell(files)
            file = char(files(vari));
        else
            file = files;
        end
        % [lon, lat, level, time]
        dat(:,:,vari)   = ncread(file, thisvar, [1 1 ,nt_skip*(ti-1)+1], [Inf Inf 1], [1 1 1]);
    end
else
    dat = ncread(files, variable, [1 1 ,nt_skip*(ti-1)+1], [Inf Inf 1], [1 1 1]);
end
dat = cast(dat, 'single');

if exist('replaceNaN','var')
    dat(isnan(dat))  =  replaceNaN;
end
