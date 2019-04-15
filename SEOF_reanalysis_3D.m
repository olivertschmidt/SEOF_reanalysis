% SPECTRAL EMPIRICAL ORTHOGONAL FUNCTION ANALYSIS OF REANALYSIS DATA
% This script computes the three-dimensional climate pattern from the
% reference paper [1]. 
%
% REFERENCE
% [1] O. T. Schmidt, G. Mengaldo, G. Balsamo and N. P. Wedi
% "Spectral Empirical Orthogonal Function analysis of weather and climate
% data", Monthly Weather Review, 2019
%
%
% O. T. Schmidt (oschmidt@ucsd.edu)
% Last revision: 14-April-2019

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');
set(groot, 'defaultFigureColor','w');
clear   variables
clc
addpath utils
load    blue2red.mat
load    coast.mat
load    stdatmos.mat
render_video        = false;

dt_in_hours         = 744;
period_in_hours     = dt_in_hours*12*12;
overlap_in_percent  = 0;
number_of_dt        = 1332;
file                = 'data/E20C/E20C_MONTHLYMEAN00_1900_2010_U131128_3D.nc';
variable            = 'u'; 
T_approx            = 744;          % approximate period (in days)
mode_idx            = 1;

absIso              = 0.75;         % plot isosurfaces at plus/minus 75% of the max(abs())
lvl_start           = 1;            % data thinning if not enough memory
lvl_skip            = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inspect and load data 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[lon,lat,time,levels,vars]  = getNetCDFfileInfo(file);
lon                         = lon-180;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SEOF setup & computation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% time step
dt_data         = time(2)-time(1);
nt_skip         = dt_in_hours/dt_data;          % set skip to dt_in_hours
time            = time(1:nt_skip:end);
dt              = dt_in_hours;
nFFT            = ceil(period_in_hours/dt_in_hours);
nFreq           = nFFT/2+1;
overlap         = ceil(nFFT*overlap_in_percent/100);

% for large data
opts.savefft    = false;                        % save FFT blocks insteasd of keeping them in memory
opts.deletefft  = true;                         % delete FFT blocks after use
opts.loadfft    = false;                        % check if FFT blocks are already saved
if iscell(file), file_info = dir(char(file(1))); else file_info = dir(file);  end
[~,file_name,~] = fileparts(file_info.name);
opts.savedir    = ['results/' file_name];       % save results to 'results' folder in the current directory
opts.savefreqs  = 1:nFreq;                      % save modes frequencies of indices
opts.nsave      = 3;                            % save the most energetic modes only
opts.nt         = min(number_of_dt,length(time));
opts.mean       = 'blockwise';
opts.normvar    = false;

% trapezoidal quadrature weights for sperical coordinates
dS  = zeros(length(lat),length(lon),length(levels));
disp('WARNING! Pressure level/altitude is not weighted!')
for lvli = 1:length(levels)
    dS(:,:,lvli)    = trapzWeightsSpherical_2D(lon, lat, 1);
end
if iscell(variable)
    dS          = repmat(dS,[1 1 1 length(variable)]);
end

% analysis interval
date_start  = datestr(double(time(1))/24 + datenum(1900,1,1));
date_end    = datestr(double(time(1)+opts.nt*dt_in_hours)/24 + datenum(1900,1,1));

% calculate SEOF
[L,P,f]         = spod(@(ti)getNetCDFdata3D(file,variable,ti,nt_skip,lvl_start,lvl_skip,0), nFFT, dS, overlap, dt, opts);
[nFreq,nModes]  = size(L);

% save info file
if opts.savefft
    saveDir = fullfile(opts.savedir,['nfft' num2str(nFFT) '_novlp' num2str(overlap) '_nblks' num2str(nBlks)]);
    save([saveDir '/info.mat'],'-regexp','^(?!(L|P)$).');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spectrum (different representations and confidence intervals)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f_days          = f*24;
if iscell(variable)
    var_name    = strcat(variable{:});
else
    var_name    = variable;
end

fig_name    = [upper(var_name) '_dt' num2str(dt_in_hours,'%g') 'days_Tblock' num2str(nFFT*dt_in_hours/24,'%g') '_spectrum'];
figure('name',[upper(var_name) ', ' date_start ' to ' date_end ', d$t=' num2str(dt_in_hours) '$ hrs, $T_{\rm blk}=' num2str(nFFT*dt_in_hours/24) '$ days'], ...
       'position',[1 700 570 275],'fileName',fig_name);

subplot(1,3,1)
set(gca, 'ColorOrder', gray(ceil(1.2*nModes))); hold on
plot(1./f_days,L)
xlabel('period [days]'), ylabel('PSD')
set(gca,'XScale','log','YScale','log'), axis tight, xlim(xlim.*[0.75 1.2]), box on
legendCell = [cellstr(num2str((1:nModes)', 'Mode $%g$'))];
legend(legendCell); legend boxoff
set(gca,'XTick',[1 7 30 365 365*5])

subplot(1,3,2)
set(gca, 'ColorOrder', gray(ceil(1.2*nModes))); hold on
L(L==0)     = eps;
plot(1./f_days,L.*f_days')
xlabel('period [days]'), ylabel('premultiplied PSD')
set(gca,'XScale','log','YScale','log','YColor','k'), axis tight, xlim(xlim.*[0.75 1.2]), box on
title([upper(var_name) ', ' date_start ' to ' date_end ', d$t$=' num2str(dt_in_hours) ' hrs, $T_{\rm blk}=$' num2str(nFFT*dt_in_hours/24) ' days']);
set(gca,'XTick',[1 7 30 365 365*5])

subplot(1,3,3)
set(gca, 'ColorOrder', gray(ceil(1.2*nModes))); hold on
L(L==0)     = eps;
plot(1./f_days,L./sum(L,2))
xlabel('period [days]'), ylabel('normalized PSD')
set(gca,'XScale','log','YScale','linear','YColor','k'), axis tight, xlim(xlim.*[0.75 1.2]), box on
ylim([1e-4 1])
set(gca,'XTick',[1 7 30 365 365*5])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Contour plot, Hovmoeller diagrams and d-lambda spectrum 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f_idx           = findnearest(f_days,1/T_approx,1);
period          = 1/f_days(f_idx);
alt             = interp1(stdatmos(:,4)*100,stdatmos(:,1),levels)/1000; % atmospalt(levels*100)/1000;

% 3D mode
if numel(P)==1
    mode        = cast(squeeze(P(f_idx,mode_idx)),'double');
else
    mode        = cast(squeeze(P(f_idx,:,:,:,mode_idx)),'double');
end
[~,~,lvl_idx]   = ind2sub([numel(lat) numel(lon) numel(levels)],find(abs(mode(:))==max(abs(mode(:)))));

fig_name    = [upper(var_name) '_period' num2str(period,'%g') 'days_mode' num2str(mode_idx)];
figure('position',[1 700 570 275],'name',fig_name);
set(gcf,'fileName',get(gcf,'name'))
colormap(blue2red)
load('utils/coast.mat')
% 2D mode at lvl_idx
mode_2D        = permute(squeeze(mode(:,:,lvl_idx)),[2 1]);

% map
subplot(3,3,[1,2,4,5])
pcolor(lon,lat,real(mode_2D)); shading interp, axis equal tight
caxis(max(abs(mode_2D(:)))*[-1 1]);
hold on
plot(coastlon,coastlat,'k-','LineWidth',0.5)
xlim([min(lon) max(lon)]); ylim([min(lat) max(lat)])
set(findall(gcf,'type','Axes'),'layer','top');
ylabel('Latitude');
title([upper(var_name) ', period: ' num2str(period,'%g') ' days, mode ' num2str(mode_idx) ', ' num2str(levels(lvl_idx)) ' hPa'],'Interpreter','latex')
grid on

[max_i,max_j] = ind2sub(size(mode_2D),find(abs(mode_2D)==max(abs(mode_2D(:)))));

plot([1 1]*lon(max_j),lat([1 end]),'m-')
plot(lon([1 end]),[1 1]*lat(max_i),'m-')

% Hovmoeller (lon)
subplot(3,3,[7,8])

t_1period   = linspace(0,period,25);
mode_max    = mode_2D(max_i,:);
phase       = exp(1i*linspace(0,2*pi,25));

pcolor(lon,t_1period,real(phase'*mode_max)); shading interp
xlabel('Longitude'); ylabel('Time [days]');
grid on
caxis(max(abs(mode_2D(:)))*[-1 1]);
daspect([1 5*period/(lon(end)-lon(1)) 1])

% Hovmoeller (lat)
subplot(3,3,[3,6])
mode_max    = mode_2D(:,max_j);

pcolor(t_1period,lat,real(mode_max*phase)); shading interp
xlabel('Time [days]');
grid on
%set(gca,'YTickLabel',[])
cb = colorbar('eastoutside'); caxis(max(abs(mode_2D(:)))*[-1 1]);
xlabel(cb, upper(var_name),'Interpreter','latex')
daspect([5*period/(lon(end)-lon(1)) 1 1])
set(gca,'XAxisLocation','top')

% spectrum
subplot(3,3,9)
set(gca, 'ColorOrder', gray(ceil(1.2*nModes))); hold on
L(L==0)     = eps;
plot(1./f_days,(L(:,1)-L(:,2))./sum(L,2),'b-')
xlabel('period [days]'), ylabel('$\Delta\lambda_{1,2}$')
set(gca,'XScale','log','YScale','linear','YColor','k'), axis tight, xlim(xlim.*[0.75 1.2]), box on
ylim([-0.01 1])
plot([period period],[0 1],'m-')
set(gca,'XTick',[1 7 30 365 365*5])

% Hovmoeller (alt)
fig_name    = [upper(var_name) '3DHoevm_period' num2str(period,'%g') 'days_mode' num2str(mode_idx)];
figure('position',[1 700 570 275],'name',fig_name);
set(gcf,'fileName',get(gcf,'name'))
colormap(blue2red)

mode_max    = permute(squeeze(mode(max_j,max_i,:)),[1 2]);
t_1period   = linspace(0,period,25);
phase       = exp(1i*linspace(0,2*pi,25));
pcolor(t_1period,alt,real(phase'*mode_max.')'); shading interp
ylabel('Altitude [km]'); xlabel('Time [days]');
grid on
daspect([5*period/(alt(1)-alt(end)) 1 1])
cb = colorbar('eastoutside'); caxis(max(abs(mode_2D(:)))*[-1 1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3D Isosurface visualization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if numel(P)==1
    dat        = cast(squeeze(P(f_idx,mode_idx)),'double')';
else
    dat        = squeeze(P(f_idx,:,:,:,mode_idx));
end

if render_video, nt_anim = 60; else nt_anim = 1; end
phase       = exp(1i*linspace(0,2*pi,nt_anim+1));
phase       = phase(1:end-1);

figure
vid_name    = [upper(var_name) '3D_period' num2str(period,'%g') 'days_mode' num2str(mode_idx)];
set(gcf,'fileName',vid_name)

for framei = 1:nt_anim
    cla
    mode_2D = real(phase(framei)*permute(dat,[2 1 3]));
    
    if framei==1, isoVal = absIso*abs(max(mode_2D(:))); end
    [Lon,Lat,Alt]   = meshgrid(lon,lat,alt);
    p = patch(isosurface(Lon,Lat,Alt,mode_2D,isoVal));
    hold off
    isonormals(Lon,Lat,Alt,mode_2D,p)
    set(p,'FaceColor','red','EdgeColor','none','FaceAlpha',0.75,'FaceLighting','phong');
    hold on
    % negative contour
    p = patch(isosurface(Lon,Lat,Alt,mode_2D,-isoVal));
    isonormals(Lon,Lat,Alt,mode_2D,p)
    set(p,'FaceColor','blue','EdgeColor','none','FaceAlpha',0.75,'FaceLighting','phong');
    
    plot(coastlon,coastlat,'k-','LineWidth',1)
    xlim([min(lon) max(lon)]); ylim([min(lat) max(lat)]); zlim([0 max(Alt(:))])
    
    plot3([1 1]*lon(max_j),lat([1 end]),[1 1]*alt(lvl_idx),'m-')
    plot3(lon([1 end]),[1 1]*lat(max_i),[1 1]*alt(lvl_idx),'m-')
    plot3([1 1]*lon(max_j),[1 1]*lat(max_i),alt([1 end]),'m-')
    
    view(3)
    daspect([1 1 0.5])
    xlabel('Longitude'); ylabel('Latitude'); zlabel('Altitude [km]');
    title([upper(var_name) ', period: ' num2str(period,'%g') ' days, mode ' num2str(mode_idx)],'Interpreter','latex')
    box on
    legend(num2str(isoVal,'%.2g'),num2str(-isoVal,'%.2g'),'Location','northwest'); legend boxoff

    drawnow

    if render_video
        print(['./animation/' vid_name '_frame' num2str(framei,'%5.5i.png')],'-dpng','-r200')
    end
end