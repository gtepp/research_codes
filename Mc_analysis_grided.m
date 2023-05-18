%% Script for calculating Mc for a grided area with the maximum curvature method

% Currently loads and uses a matfile with a ZMAP catalog (named a) for each year
% formatted as: [lon, lat, year, month, day, mag, depth, hour, minute]

% Converted from R algorithms available in:
% http://www.corssa.org/export/sites/corssa/.galleries/articles-pdf/Mignan-Woessner-2012-CORSSA-Magnitude-of-completeness.pdf

% can also output movie showing temporal changes

% Created: 3/24/2023
% Updated: 4/4/2023
% By: Gabrielle Tepp, Caltech

%-------------------

%% INPUT PARAMETERS

mbin = 0.1; %Magnitude bin size
nbsample = 200; %Bootstrapping

%select function: maxc (default), gft, mbs, emr (only one method allowed at a time)
mc_func = 'maxc';

y_range = [1932:20:1972,1983:10:2023];%[1932:10:1972,1978:5:2023];%1932:2023; %[1932:5:1972,1975:2023];% list of years to loop over (end year not included)
lat_range = [31.5,37.5]; %[32.74,33.64];% latitude bounds of area
lon_range = [-121.2,-114]; %[-116,-114.8];% longitude bounds of area
grid_sz = 0.5; % in degrees
mineq = 100; % minimum number of events required for grid bin Mc calculation

makemov = 0; % make movie (1=yes)
mname = 'Mc_1932-2023_10-20yrbin.mp4';

%% MAKE GRID

% grid bounds
gr_lat = (lat_range(1)+grid_sz/2):grid_sz:(lat_range(2)-grid_sz/2); 
gr_lon = (lon_range(1)+grid_sz/2):grid_sz:(lon_range(2)-grid_sz/2);

% center points
mp_lat = (gr_lat(1:end-1)+gr_lat(2:end))/2;
mp_lon = (gr_lon(1:end-1)+gr_lon(2:end))/2; 


%% COMPUTE Mc

allMc = nan(length(gr_lat)-1,length(gr_lon)-1,length(y_range)-1,2);
n_eqs = nan(length(gr_lat)-1,length(gr_lon)-1,length(y_range)-1);

for y = 1:length(y_range)-1

    mag = [];
    latlon = [];

    % get data for year bin
    for i = 0:(y_range(y+1)-y_range(y)-1) 

        load(strcat('./zmap/data_files/scsn_',num2str(y_range(y)+i),'.mat')) % matfile with catalog in ZMAP format

        if (y_range(y)+i) < 1990 % for older data
            a(a(:,6)==0,:) = []; % remove rows with M=0 (1976 & 1982 have a lot)
        end

        % remove data outside geographic bounds
        a((a(:,1)<lon_range(1)|a(:,1)>lon_range(2)|a(:,2)<lat_range(1)|a(:,2)>lat_range(2)),:) = [];

        mag = [mag;a(:,6)]; % get magnitudes
        latlon = [latlon;a(:,[1,2])]; % keep locations

    end

    % loop over grid
    for la = 1:length(gr_lat)-1
        for lo = 1:length(gr_lon)-1

            % get mags for current grid bin
            gr_mag = mag((latlon(:,1)>gr_lon(lo)&latlon(:,1)<gr_lon(lo+1)&latlon(:,2)>gr_lat(la)&latlon(:,2)<gr_lat(la+1)),:);

            % if there are EQs in grid bin, calculate Mc
            if length(gr_mag) >= mineq

                %For mbass(), see algorithm Amorese [2007]
                if strcmpi(mc_func,'gft') == 1
                    fname = 'Goodness-of-fit Test (GFT)';
                    Mc_bootstrap = bootstrp(nbsample,@(samps) gft(samps,mbin),gr_mag);

                    allMc(la,lo,y,1) = mean(cell2mat(Mc_bootstrap(:,1)),'omitnan');
                    allMc(la,lo,y,2) = std(cell2mat(Mc_bootstrap(:,1)), 'omitnan');

                elseif strcmpi(mc_func,'mbs') == 1
                    fname = 'Mc by b-val Stability (MBS)';
                    Mc_bootstrap = bootstrp(nbsample,@(samps) mbs(samps,mbin),gr_mag);

                    allMc(la,lo,y,1) = mean(Mc_bootstrap,'omitnan');
                    allMc(la,lo,y,2) = std(Mc_bootstrap, 'omitnan');

                elseif strcmpi(mc_func,'emr') == 1
                    fname = 'Entire Magnitude Range method (EMR)';
                    %when using emr(), the loop may break due to failure of nlsfit(),
                    %in this case use:
                    %Mc_bootstrap[i] = as.numeric(try(emr(sample(mag, replace=TRUE),mbin)$Mc))
                    disp('EMR is not currently working')
                    break;
                else % default to maxc
                    fname = 'Maximum Curvature (MAXC)';
                    [Mc_bootstrap,bs] = bootstrp(nbsample,@(samps) maxc(samps,mbin),gr_mag);

                    allMc(la,lo,y,1) = mean(Mc_bootstrap,'omitnan');
                    allMc(la,lo,y,2) = std(Mc_bootstrap, 'omitnan');
                end

                n_eqs(la,lo,y) = length(gr_mag); % keep number of EQs in year group

            end

            %disp(strcat("Using the ",mc_func{m}," method for ",num2str(y_range(y))," to ",num2str(y_range(y+1))))
            %disp(strcat("Mc (mean): ", num2str(allMc(y,1,m))))
            %disp(strcat("Sigma0 (std. dev.): "," ", num2str(allMc(y,2,m))))
        end
    end

end


%% PLOT RESULTS

if makemov == 1
    % for movie, preallocate
    M(length(y_range)-1) = struct('cdata',[],'colormap',[]);
end

for y = 1:length(y_range)-1

    f = figure;
    f.Position(3:4) = [1600 600];

    ax1 = subplot(1,2,1);

    %usamap(lat_range,lon_range)
    %cali = shaperead('usastatehi', 'UseGeoCoords', true,'Selector',{@(name) strcmpi(name,'California'), 'Name'});
    axesm('mercator','MapLatLimit',lat_range,'MapLonLimit',lon_range,'Grid','on','MeridianLabel','on','MLineLocation',2,...
        'ParallelLabel','on','PLineLocation',2,'FLatLimit',lat_range,'FLonLimit',lon_range,'Frame','on') %Define Map axes
   
    surfm(gr_lat(1:end-1),gr_lon(1:end-1),allMc(:,:,y,1));
    geoshow('usastatehi.shp',"EdgeColor","black",'FaceColor','none')

    colormap(ax1,'parula')
    c = colorbar('eastoutside');
    clim([0,4]);
    c.Label.String = 'Mc';
    fontsize(gcf,16,'points')

    if strcmpi(mc_func,'gft') == 1
        title('Goodness-of-fit Test (GFT)');
    elseif strcmpi(mc_func,'mbs') == 1
        title('Mc by b-val Stability (MBS)');
    elseif strcmpi(mc_func,'emr') == 1
        title('Entire Magnitude Range method (EMR)');
    else % default to maxc
        title('Maximum Curvature (MAXC)');
    end


    % plot number of earthquakes in each year group

    ax2 = subplot(1,2,2);

    axesm('mercator','MapLatLimit',lat_range,'MapLonLimit',lon_range,'Grid','on','MeridianLabel','on','MLineLocation',2,...
        'ParallelLabel','on','PLineLocation',2,'FLatLimit',lat_range,'FLonLimit',lon_range,'Frame','on') %Define Map axes
    %geoshow(cali,"FaceColor","black");

    surfm(gr_lat(1:end-1),gr_lon(1:end-1),n_eqs(:,:,y));
    %geoshow(repmat(mp_lat,1,length(mp_lon)),repmat(mp_lon,1,length(mp_lat)),n_eqs(:,:,y))
    geoshow('usastatehi.shp',"EdgeColor","black",'FaceColor','none')

    title([num2str(y_range(y)),' to ',num2str(y_range(y+1))])
    colormap(ax2,'copper')
    c = colorbar('eastoutside');
    c.Label.String = 'Number of Earthquakes';
    clim([50,1500])
    fontsize(gcf,16,'points')
    grid on;
    box on;

    if makemov == 1
        M(y) = getframe(f); % get frame for movie
    end

end

if makemov == 1
    %movie(M,1,0.2);
    v = VideoWriter(mname,'MPEG-4');
    v.FrameRate = 1;
    open(v)
    writeVideo(v,M)
    close(v)
end

%% MC FUNCTIONS

%FMD
function res = fmd(mag,mbin)

mi = transpose(min(round(mag/mbin)*mbin):mbin:max(round(mag/mbin)*mbin));
nbm = length(mi);
cumnbmag = zeros(nbm,1);

for i = 1:nbm
    cumnbmag(i) = length(find(mag > (mi(i)-mbin/2)));
end

cumnbmagtmp = [cumnbmag;0];
nbmag = abs(diff(cumnbmagtmp));
res = [mi, cumnbmag, nbmag]; %m cum noncum

end


%Maximum Curvature (MAXC) [e.g., Wiemer & Wyss, 2000]
function Mc = maxc(mag,mbin)

FMD = fmd(mag,mbin);
[~,ind] = max(FMD(:,3)); % max noncum value
Mc = FMD(ind,1); % M of max 

end

%Goodness-of-fit test (GFT) [Wiemer & Wyss, 2000]
function params = gft(mag,mbin)

FMD = fmd(mag,mbin);
McBound = maxc(mag,mbin);
Mco = McBound-0.4+(0:14)/10;
R = zeros(15,1);

for i = 1:15
    indmag = find(mag > (Mco(i)-mbin/2));
    b = log10(exp(1))/(mean(mag(indmag))-(Mco(i)-mbin/2));
    a = log10(length(indmag))+b*Mco(i);
    FMDcum_model = 10.^(a-b*FMD(:,1));
    indmi = find(FMD(:,1) >= Mco(i));
    R(i) = sum(abs(FMD(indmi,2)-FMDcum_model(indmi)))/sum(FMD(indmi,2))*100;
    %in Wiemer&Wyss [2000]: 100-R
end

indGFT = find(R <= 5); %95% confidence
if ~isempty(indGFT)
    Mc = Mco(indGFT(1));
    best = "95%";
else
    indGFT = find(R <= 10); %90% confidence
    if ~isempty(indGFT)
        Mc = Mco(indGFT(1));
        best = "90%";
    else
        Mc = McBound;
        best = "MAXC";
    end
end

params = {Mc,best,Mco,R};

end

%Mc by b-val Stability (MBS) [Cao & Gao, 2002]
%Modification with Shi & Bolt [1982] uncertainty [Woesner & Wiemer, 2005]
function [Mc,Mco,bi,unc,bave] = mbs(mag,mbin)

McBound = maxc(mag,mbin);
Mco = McBound-0.7+(0:19)/10;
bi = zeros(20,1); 
unc = zeros(20,1);

for i = 1:20
    indmag = find(mag > Mco(i)-mbin/2);
    nbev = length(indmag);
    bi(i) = log10(exp(1))/(mean(mag(indmag))-(Mco(i)-mbin/2));
    unc(i) = 2.3*bi(i).^2*sqrt(sum((mag(indmag)- mean(mag(indmag))).^2)/(nbev*(nbev-1)));
end

bave = zeros(15,1);

for i = 1:15
    bave(i) = mean(bi(i:i+5));
end

%dbi_old = abs(diff(bi));
%indMBS_old = find(dbi_old <= 0.03);
dbi = abs(bave(1:15)-bi(1:15));
indMBS = find(dbi <= unc(1:15));

if ~isempty(indMBS)
    Mc = Mco(indMBS(1));
else % if no residuals <= uncertainty
    Mc = NaN;
end

end

%Entire Magnitude Range method (EMR) [Woesner & Wiemer, 2005]
% function res = emr(mag,mbin)
% FMD = fmd(mag,mbin)
% nbm = length(FMD$m)
% McMAXC = maxc(mag,mbin)$Mc
% mu = abs(McMAXC/2); sig = abs(McMAXC/4)
% if(mu > 1)mu = abs(McMAXC/10); sig = abs(McMAXC/20)
% McBound = McMAXC
% Mco = McBound-0.3+(seq(9)-1)/10
% params = zeros(9*4); dim(params) = c(9,4) %a, b, mu, sigma
% prob = zeros(9)
% savedmodel = zeros(9*nbm); dim(savedmodel) = c(9,nbm)
% for(i in 1:9){
% indmag = which(mag > Mco[i]-mbin/2)
% nbev = length(indmag)
% b = log10(exp(1))/(mean(mag[indmag])-(Mco[i]-mbin/2))
% a = log10(length(indmag))+b*Mco[i]
% cumN = 10^(a-b*FMD$m)
% params[i,1] = a; params[i,2] = b
% cumNtmp = 10^(a-b*(max(FMD$m)+mbin))
% cumNtmp = c(cumN, cumNtmp)
% N = abs(diff(cumNtmp))
% data = data.frame(N=N, m=FMD$m, Nd=FMD$noncum)
% indLow = which(FMD$m < Mco[i]); indHigh = which(FMD$m >= Mco[i])
% dataTest = data.frame(N=data$N[indLow], m=data$m[indLow], Nd=data$Nd[indLow])
% dataTmp = data.frame(N=data$N[indHigh], m=data$m[indHigh], Nd=data$Nd[indHigh])
% checkNo0 = which(dataTest$Nd != 0)
% dataTest = data.frame(N=dataTest$N[checkNo0], m=dataTest$m[checkNo0],
% Nd=dataTest$Nd[checkNo0])
% %Nmax = max(dataTmp$Nd)
% Nmax = max(dataTest$Nd)
% %Nmax = dataTest$Nd[length(dataTest$Nd)]
% Mmintmp = min(dataTest$m)
% dataTest$Nd = dataTest$Nd/Nmax
% dataTest$m = dataTest$m-Mmintmp
% data4fit = data.frame(N=dataTest$Nd, m=dataTest$m)
% %non-linear least squares fit
% nlsfit = nls(N~pnorm(m, mean=mean, sd=sd), data=data4fit,
% start=list(mean=mu, sd=sig), control=list(maxiter=100, warnOnly = TRUE))
% params[i,3] = coef(nlsfit)["mean"]; params[i,4] = coef(nlsfit)["sd"]
% dataTest$N = pnorm(dataTest$m, mean=coef(nlsfit)["mean"],
% sd=coef(nlsfit)["sd"])*Nmax
% dataTest$m = dataTest$m+Mmintmp
% dataTest$Nd = dataTest$Nd*Nmax
% dataPred = data.frame(N=c(dataTest$N, dataTmp$N), m=c(dataTest$m, dataTmp$m),
% Nd=c(dataTest$Nd, dataTmp$Nd))
% dataPred$N = round(dataPred$N)
% savedmodel[i,c(checkNo0,indHigh)] = dataPred$N
% %CHECK EMR METHOD%
% %pdf(strcat(wd,"plot_NonCumModel_",Mco[i],".pdf", sep=""))
% %plot(dataPred$m, dataPred$Nd, pch=18, xlab="Magnitude",
% ylab="Cumulative Number", log="y")
% %points(dataPred$m, dataPred$N, pch=1)
% %abline(v=Mco[i], lty="dashed")
% %legend("topright", c("Data","EMR model"), cex=0.8, lty=c(0,0), pch=c(18,1))
% %dev.off()
% %write.table(dataPred, file=strcat(wd, "file_NonCumModel_",Mco[i],
% ".txt", sep=""))
% %Logarithm to the basis of 10 of Poisson probability density
% probtmp = zeros(nbm)
% CheckNo0 = which(dataPred$N != 0)
% Pmodel = dataPred$N[CheckNo0]; Pdata = dataPred$Nd[CheckNo0]
% Completeness Magnitude 41
% probtmp[CheckNo0] = 1/log(10)*(-Pmodel+Pdata*log(Pmodel)-lgamma(Pdata+1))
% prob[i] = -sum(probtmp)
% }
% indbestfit = which(prob == min(prob, 'omitnan'))
% res = list(Mc=Mco[indbestfit], a=params[indbestfit,1], b=params[indbestfit,2],
% mu=params[indbestfit,3], sigma=params[indbestfit,4],
% model=savedmodel[indbestfit,], Mco=Mco, prob=prob)
% 
% end
