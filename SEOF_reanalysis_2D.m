% SPECTRAL EMPIRICAL ORTHOGONAL FUNCTION ANALYSIS OF REANALYSIS DATA
% This script computes the two-dimensional climate patterns from the
% reference paper [1]. Before using this script, the reanalysis data has to be
% downloaded. The corresponding download scripts are located in the
% 'data/EI' (ERA-Interim) and 'data/E20C' (ERA 20C) folders. Details on how
% to gain access to the various ECMWF reanalysis data sets can be found under 
% https://confluence.ecmwf.int/display/WEBAPI/Access+ECMWF+Public+Datasets.
% The core routine spod() that performs the SEOF (called SPOD in fluid 
% mechanics) is located in the 'utils' folder. The most recent version can
% be found under
% https://www.mathworks.com/matlabcentral/fileexchange/65683-spectral-proper-orthogonal-decomposition-spod.
% 
% SEOF PARAMETERS
% dt_in_hours         : temoral resolution ('dt' in [1])   
% period_in_hours     : segment length ('T' in [1]);
% overlap_in_percent  : overlap between FFT blocks
% number_of_dt        : number of spanshots ('Inf' for all available data)
% T_approx            : period of the mode to plot in days (closest)
% mode_idx            : mode number ('1' for leading mode)
% lat_min, lat_max, lon_min, lon_max : longitudinal and lateral restriction of the SEOF domain (-Inf, Inf for entire data domain)
%
% EXAMPLE: ENSO (reproduces figures 1 & 2 in [1])
% Choosing pattern = 'ENSO' [default] performs an SEOF analysis for the sea
% surface temperature data in the file 'data/E20C/E20C_AN_1900_2010_SST34128.nc'.
% Once your ECMWF account is set up (follow the instructions in the link above),
% the data file can be downloaded using the corresponding Python script
% 'data/E20C/E20C_AN_1900_2010_SST34128.py'. To perform a 12-year monthly 
% analysis, we let dt_in_hours = 1 month = 24*30 hours and 
% period_in_hours = 12 years = 12*365*24 hours. A clear ENSO signature is 
% found for a period of 6 years. For visulaization purposes, we hence let
% T_approx = 6 years = 2190 days and mode_idx = 1, as we are interested in
% the pattern that contains most of the variance.
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
render_video = false;    % animate mode and save PNGs for video rendering


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choose pattern
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pattern  = 'ENSO';
% 'ENSO'    : El Nino“Southern Oscillation and Pacific decadal oscillation
% 'MJO'     : Madden-Julian Oscillation
% 'QBO'     : Quasi-biennial oscillation
% 'MEI'     : Multivariate ENSO Index (MEI)
% 'BEI'     : Bi-variate index (BEI)




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform SEOF analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch pattern
    case 'ENSO'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ENSO - El Nino“Southern Oscillation%
        % PDO  - Pacific decadal oscillation  %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % SST: Sea surface temperature [K]
        % -> ENSO (730-2555 days)
        % -> PDO (2920-4380 days)
        % data/E20C/E20C_AN_1900_2010_SST34128.py
        file        = 'data/E20C/E20C_AN_1900_2010_SST34128.nc';
        variable    = 'sst';
        % 12-year monthly analysis
        dt_in_hours         = 24*30;        
        period_in_hours     = 12*365*24;
        overlap_in_percent  = 0;
        number_of_dt        = Inf;
        % visualization
        T_approx            = 2190;        % approximate period (in days)
        mode_idx            = 1;
        lat_min             = -Inf;
        lat_max             =  Inf;
        lon_min             = -Inf;
        lon_max             =  Inf;
    case 'MJO'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % MJO - Madden-Julian Oscillation %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % TTR: Top net thermal radiation [W m**-2]
        % -> MJO (30-90 days)
        file        = 'data/EI/EI_FC3_1979_2017_TTR179128.nc';
        variable    = 'ttr';
        % file        = 'data/EI/EI_FC3_1979_2017_TP228128.nc';
        % variable    = 'tp';
        % 2-year twice-daily analysis
        dt_in_hours         = 12;
        period_in_hours     = 1*365*24;
        overlap_in_percent  = 0;
        number_of_dt        = Inf;
        % visualization
        T_approx            = 40.56;        % approximate period (in days)
        mode_idx            = 1;
        lat_min             = -Inf;
        lat_max             =  Inf;
        lon_min             = -Inf;
        lon_max             =  Inf;
    case 'QBO'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % QBO  - Quasi-biennial oscillation %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % MSLP: Mean sea level pressure [Pa]
        % -> QBO (850 months)
        file        = 'data/E20C/E20C_AN_1900_2010_MSLP151128.nc';
        variable    = 'msl';
        % 12-year monthly analysis
        dt_in_hours         = 24*30;
        period_in_hours     = 12*365*24;
        overlap_in_percent  = 0;
        number_of_dt        = Inf;
        % visualization
        T_approx            = 876;        % approximate period (in days)
        mode_idx            = 1;
        lat_min             = -Inf;
        lat_max             =  Inf;
        lon_min             = -Inf;
        lon_max             =  Inf;
    case 'MEI'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Multivariate ENSO Index (MEI) %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        file        = 'data/E20C/E20C_MONTHLYMEAN00_1900_2010_MEI.nc';
        variable    = {'sst';'msl';'tcc';'u10';'v10';'t2m'};
        normalize   = true;
        % 12-year monthly analysis
        dt_in_hours         = 720;
        period_in_hours     = dt_in_hours*12*5;
        overlap_in_percent  = 0;
        number_of_dt        = 1332;
        % visualization
        T_approx            = 876;        % approximate period (in days)
        mode_idx            = 1;
        lat_min             = -Inf;
        lat_max             =  Inf;
        lon_min             = -Inf;
        lon_max             =  Inf;
    case 'BEI'
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        % Bi-variate index (BEI) %
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        % SST and U@300hPa
        file        = {'data/E20C/E20C_AN_1900_2010_U300hPa.nc'; 'data/E20C/E20C_AN_1900_2010_SST34128.nc'};
        variable    = {'u'; 'sst'};
        normalize   = true;
        % 12-year monthly analysis
        dt_in_hours         = 24*30;
        period_in_hours     = 12*365*24;
        overlap_in_percent  = 0;
        number_of_dt        = Inf;
        % visualization
        T_approx            = 2190;        % approximate period (in days)
        mode_idx            = 1;
        lat_min             = -Inf;
        lat_max             =  Inf;
        lon_min             = -Inf;
        lon_max             =  Inf;
    otherwise
        error(['Pattern ' pattern ' not implemented.']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inspect and load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[lon,lat,time,level,vars]   = getNetCDFfileInfo(file);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SEOF setup & computation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% time step
dt_data         = time(2)-time(1);
nt_skip         = round(dt_in_hours/dt_data);   % set skip to dt_in_hours
time            = time(1:nt_skip:end);
dt              = dt_in_hours;
nFFT            = ceil(period_in_hours/dt_in_hours);
nFreq           = nFFT/2+1;
overlap         = ceil(nFFT*overlap_in_percent/100);

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
opts.conflvl    = 0.95;

% variable-wise normalization by variance by variance via weight matrix
if ~exist('normalize','var'), normalize = false; end
if ~iscell(variable), variable    = {variable}; end
variance    = ones(length(variable),1);
if normalize
    disp(' ')
    disp('Normalization by variance')
    disp('------------------------------------')
    for vari= 1:length(variable)
        if iscell(file)
            filei = char(file(vari));
        else
            filei = file;
        end
        var_name    = char(variable(vari));
        disp(['Reading ' var_name ' to compute variance for normalization']);
        dat = ncread(filei, var_name, [1 1 1], [Inf Inf opts.nt], [1 1 1]);
        variance(vari)  = var(dat(~isnan(dat(:))));
    end
end
clear dat

% trapezoidal quadrature weights for sperical coordinates including
% variance normalization
dS              = trapzWeightsSpherical_2D(lon, lat, 1);
dS(lat<=lat_min|lat>=lat_max,:)                 = 0;
if lon_min>=lon_max
    dS(:,(lon-180)<=lon_min&(lon-180)>=lon_max)     = 0;
else
    dS(:,(lon-180)<=lon_min|(lon-180)>=lon_max)     = 0;
end
if iscell(variable)
    dS          = repmat(dS,[1 1 length(variable)]);
    for vari= 1:length(variable)
        dS(:,:,vari)          = dS(:,:,vari)/variance(vari);
    end
else
    dS          = dS/variance;
end

% analysis interval
date_start  = datestr(double(time(1))/24 + datenum(1900,1,1));
date_end    = datestr(double(time(1)+opts.nt*dt_in_hours)/24 + datenum(1900,1,1));

% calculate SEOF
[L,P,f,Lc]      = spod(@(ti)getNetCDFdata(file,variable,ti,nt_skip,0), nFFT, dS, overlap, dt, opts);
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
figure('name',[upper(var_name) ', ' date_start ' to ' date_end ', d$t=' num2str(dt_in_hours) '$ hrs, $T_{\rm blk}=' num2str(nFFT*dt_in_hours/24) '$ days'], 'fileName',fig_name);

subplot(1,3,1)
for mi=1:nModes
    patch('XData',1./[f_days(2:end) fliplr(f_days(2:end))],'YData',[squeeze(Lc(2:end,mi,1)); flipud(squeeze(Lc(2:end,mi,2)))],'EdgeColor','none','FaceColor',0.8*[1 1 1]);
    hold on
end
plot(1./f_days,L,'k-')
xlabel('period [days]'), ylabel('$\lambda$')
set(gca,'XScale','log','YScale','log'), axis tight, xlim(xlim.*[0.75 1.2]), box on
set(gca,'XTick',[1 7 30 365 365*5])

subplot(1,3,2)
set(gca, 'ColorOrder', gray(ceil(1.2*nModes))); hold on
L(L==0)     = eps;
plot(1./f_days,L.*f_days')
xlabel('period [days]'), ylabel('$f\lambda$')
set(gca,'XScale','log','YScale','log','YColor','k'), axis tight, xlim(xlim.*[0.75 1.2]), box on
title([upper(var_name) ', ' date_start ' to ' date_end ', d$t$=' num2str(dt_in_hours) ' hrs, $T_{\rm blk}=$' num2str(nFFT*dt_in_hours/24) ' days']);
set(gca,'XTick',[1 7 30 365 365*5])

subplot(1,3,3)
set(gca, 'ColorOrder', gray(ceil(1.2*nModes))); hold on
L(L==0)     = eps;
plot(1./f_days,L./sum(L,2))
xlabel('period [days]'), ylabel('$\lambda/\Sigma_j \lambda_j$')
set(gca,'XScale','log','YScale','linear','YColor','k'), axis tight, xlim(xlim.*[0.75 1.2]), box on
ylim([1e-4 1])
set(gca,'XTick',[1 7 30 365 365*5])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Contour plot, Hovmoeller diagrams and d-lambda spectrum
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f_idx       = findnearest(f_days,1/T_approx);

for f_idx = f_idx
    period      = 1/f_days(f_idx);
    for vari= 1:length(variable)
        var_name    = char(variable(vari));
        fig_name    = [upper(var_name) '_period' num2str(period,'%g') 'days_mode' num2str(mode_idx)];
        figure('name',fig_name);
        set(gcf,'fileName',get(gcf,'name'))
        colormap(blue2red)
        load('utils/coast.mat')
        if numel(P)==1
            if length(variable)>1
                mode        = cast(squeeze(P(f_idx,vari,mode_idx)),'double')';
            else
                mode        = cast(squeeze(P(f_idx,mode_idx)),'double')';
            end
        else
            if length(variable)>1
                mode        = cast(squeeze(P(f_idx,:,:,vari,mode_idx)),'double')';
            else
                mode        = cast(squeeze(P(f_idx,:,:,mode_idx)),'double')';
            end
        end
        mode = fftshift(mode,2);

        % map
        subplot(3,3,[1,2,4,5])

        pcolor(lon-180,lat,real(mode)); shading interp, axis equal tight
        caxis(max(abs(mode(:)))*[-1 1]);
        hold on
        plot(coastlon,coastlat,'k-','LineWidth',0.5)
        xlim([-180 180]); ylim([-90 90])
        if lon_min<=lon_max
            plot([lon_min lon_min],[lat_min lat_max],'g-','LineWidth',1)
            plot([lon_max lon_max],[lat_min lat_max],'g-','LineWidth',1)
            plot([lon_min lon_max],[lat_min lat_min],'g-','LineWidth',1)
            plot([lon_min lon_max],[lat_max lat_max],'g-','LineWidth',1)
        else
            plot([lon_min lon_min],[lat_min lat_max],'g-','LineWidth',1)
            plot([lon_max lon_max],[lat_min lat_max],'g-','LineWidth',1)
            plot([lon_min 180],[lat_min lat_min],'g-','LineWidth',1)
            plot([lon_min 180],[lat_max lat_max],'g-','LineWidth',1)

            plot([lon_min lon_min],[lat_min lat_max],'g-','LineWidth',1)
            plot([lon_max lon_max],[lat_min lat_max],'g-','LineWidth',1)
            plot([lon_max -180],[lat_min lat_min],'g-','LineWidth',1)
            plot([lon_max -180],[lat_max lat_max],'g-','LineWidth',1)
        end
        set(findall(gcf,'type','Axes'),'layer','top');

        ylabel('Latitude');
        title([upper(var_name) ', period: ' num2str(period,'%g') ' days, mode ' num2str(mode_idx)],'Interpreter','latex')
        grid on

        [max_i,max_j] = ind2sub(size(mode),find(abs(mode)==max(abs(mode(:)))));
        max_i = max_i(1); max_j = max_j(1);

        plot([1 1]*lon(max_j)-180,lat([1 end]),'m-','LineWidth',0.5)
        plot(lon([1 end])-180,[1 1]*lat(max_i),'m-','LineWidth',0.5)

        % Hovmoeller (lon)
        subplot(3,3,[7,8])

        t_1period   = linspace(0,period,25);
        mode_max    = mode(max_i,:);
        phase       = exp(1i*linspace(0,2*pi,25));

        pcolor(lon-180,t_1period,real(phase'*mode_max)); shading interp
        xlabel('Longitude'); ylabel('Time [days]');
        grid on
        daspect([1 5*period/(lon(end)-lon(1)) 1])
        set(gca,'YTick',ylim)
        set(findall(gcf,'type','Axes'),'layer','top');

        % Hovmoeller (lat)
        subplot(3,3,[3,6])
        mode_max    = mode(:,max_j);

        pcolor(t_1period,lat,real(mode_max*phase)); shading interp
        xlabel('Longitude'); xlabel('Time [days]');
        grid on
        cb = colorbar('eastoutside'); caxis(max(abs(mode(:)))*[-1 1]);
        xlabel(cb, upper(var_name),'Interpreter','latex')
        daspect([5*period/(lon(end)-lon(1)) 1 1])
        set(gca,'XAxisLocation','top')
        set(gca,'XTick',xlim)
        set(findall(gcf,'type','Axes'),'layer','top');

        % spectrum
        subplot(3,3,9)
        set(gca, 'ColorOrder', gray(ceil(1.2*nModes))); hold on
        L(L==0)     = eps;
        plot(1./f_days,(L(:,1)-L(:,2))./sum(L,2),'k-','LineWidth',0.5)
        xlabel('period [days]'), ylabel('$\Delta\lambda$')
        set(gca,'XScale','log','YScale','linear','YColor','k'), axis tight, xlim(xlim.*[0.75 1.2]), ylim([0 1]), box on
        plot([period period],[0 1],'m-')
        set(gca,'XTick',[1 7 30 365 365*5])

        drawnow
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Animate and save PNGs for video rendering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if render_video
    nt_anim   = 60;
    f_anim    = nt_anim/period;
    t_anim    = linspace(0,period,nt_anim); t_anim = t_anim(1:end-1);

    dat_max   = max(abs(mode(:)));
    fig_vid   = figure('renderer','opengl','Color','w','position',[1 700 570 275]);
    colormap(blue2red)
    load('utils/coast.mat')
    vid_name    = [upper(var_name) '_period' num2str(period,'%g') 'days_mode' num2str(mode_idx)];

    for ti = [1 1:nt_anim-1]
        figure(fig_vid)
        title(['Period: ' num2str(period) ' hrs, mode #' num2str(mode_idx)])
        phi = exp(-1i*f_days(f_idx)*2*pi*t_anim(ti));
        hold off
        pcolor(lon-180,lat,fftshift(real(phi*mode),2)); shading interp, axis equal tight
        caxis(0.8*dat_max*[-1 1]);
        hold on
        plot(coastlon,coastlat,'k-','LineWidth',0.5)
        xlim([-180 180]); ylim([-90 90])
        set(findall(gcf,'type','Axes'),'layer','top');
        grid on
        xlabel('Longitude'); ylabel('Latitude');
        title([upper(var_name) ', period: ' num2str(period,'%g') ' days, mode ' num2str(mode_idx)],'Interpreter','latex')

        drawnow

        print(['./animation/' vid_name '_frame' num2str(ti,'%5.5i.png')],'-dpng','-r200')
    end
    render_cmd  = ['ffmpeg -f image2 -r 25 -pattern_type glob -i ''' vid_name '_frame*.png'' ' vid_name '.mp4'];
    disp(' ');
    disp('Video rendering command');
    disp('------------------------------------');
    disp(render_cmd);
end
