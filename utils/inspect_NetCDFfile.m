% G. Mengaldo (gianmarco.mengaldo@gmail.com)
% O. T. Schmidt (oschmidt@ucsd.edu), 
% Last revision: 07-April-2019

function [] = inspect_NetCDFfile(file,nt)

addpath utils
load blue2red.mat
load coast.mat

if nargin==1, nt = 25; end

[lon,lat,time,fields]=getNetCDFfileInfo(file);
nLon    = length(lon);
nLat    = length(lat);
nFields = length(fields);

% load data into array
q           = zeros(nLat,nLon,nt,nFields);
for ti=1:nt
    for i=1:nFields
        q(:,:,ti,i) = getNetCDFdata(file,char(fields(i)),ti,1,0).';
    end
end
q   = fftshift(real(q),2);

% animate fields
figure('Color','w')
colormap(blue2red)
load('utils/coast.mat')
nFields     = length(fields);
for ti=1:nt
    for i=1:nFields
        subplot(1,nFields,i)
        pcolor(lon-180,lat,q(:,:,ti,i)); shading interp, axis equal tight
        if ti==1
            ca = caxis; caxis(ca);
            cb = colorbar('southoutside'); cb.Label.String = char(fields(i));
        end
        title(['t=' num2str(time(ti))])
        hold on
        plot(coastlon,coastlat,'k-','LineWidth',0.5)
        xlim([-180 180]); ylim([-90 90])
        set(findall(gcf,'type','Axes'),'layer','top');
        grid on
    end
    drawnow
end


% plot mean and variance
q_mean      = mean(q,3);
q_var       = sum((q-q_mean).^2,3)/(nt-1);
figure('name','mean','Color','w')
colormap(blue2red)
for i=1:nFields
    subplot(2,nFields,i)
    pcolor(lon-180,lat,squeeze(q_mean(:,:,i))); shading interp, axis equal tight
    ca = caxis; caxis(ca);
    cb = colorbar('southoutside'); cb.Label.String = ['mean(' char(fields(i)) ')'];
    hold on
    plot(coastlon,coastlat,'k-','LineWidth',0.5)
    isnan_idx   = isnan(q_mean(:,:,i));
    xlim([-180 180]); ylim([-90 90])
    set(findall(gcf,'type','Axes'),'layer','top');
    grid on
    
    subplot(2,nFields,i+nFields)
    pcolor(lon-180,lat,squeeze(q_var(:,:,i))); shading interp, axis equal tight
    ca = caxis; caxis(ca);
    cb = colorbar('southoutside'); cb.Label.String = ['variance(' char(fields(i)) ')'];
    hold on
    plot(coastlon,coastlat,'k-','LineWidth',0.5)
    isnan_idx   = isnan(q_mean(:,:,i));
    xlim([-180 180]); ylim([-90 90])
    set(findall(gcf,'type','Axes'),'layer','top');
    grid on    
end

%%

% animate fields
figure('Color','w')
colormap(blue2red)
load('utils/coast.mat')
nFields     = length(fields);
for ti=1:nt
    for i=1:nFields
        subplot(1,nFields,i)
        pcolor(lon-180,lat,(q(:,:,ti,i)-q_mean(:,:,i))./q_var(:,:,i)); shading interp, axis equal tight
        if ti==1
            ca = caxis; caxis(min(abs(ca))*[-1 1]);
            cb = colorbar('southoutside'); cb.Label.String = ['(' char(fields(i)) '-mean(' char(fields(i)) '))./variance(' char(fields(i)) ')'];
        end
        title(['t=' num2str(time(ti))])
        hold on
        plot(coastlon,coastlat,'k-','LineWidth',0.5)
        xlim([-180 180]); ylim([-90 90])
        set(findall(gcf,'type','Axes'),'layer','top');
        grid on
    end
    drawnow
end

