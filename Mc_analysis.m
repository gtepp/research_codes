%% Calculate Mc with different methods to examine temporal changes in catalog by year(s)

% Currently loads and uses a matfile with a ZMAP catalog (named a) for each year
% formatted as: [lon, lat, year, month, day, mag, depth, hour, minute]

% Converted from R algorithms available in:
% http://www.corssa.org/export/sites/corssa/.galleries/articles-pdf/Mignan-Woessner-2012-CORSSA-Magnitude-of-completeness.pdf

% Methods included (see above document for details):
% Max curvature (maxc)
% Goodness-of-fit for b-value (gft)
% Mc by b-val Stability (mbs)
% Entire Magnitude Range method (emr) - not currently functional (still in R - not converted)

% Created: 3/20/2023
% Updated: 4/4/2023
% By: Gabrielle Tepp, Caltech

%-------------------

%% INPUT PARAMETERS

mbin = 0.1; %Magnitude bin size
nbsample = 200; %Bootstrapping

%select function: maxc (default), gft, mbs, emr (can use multiple - comma separate in cell array)
mc_func = {'maxc','gft','mbs'};

y_range = [1932:5:1972,1975:2023];%1932:2023;% list of years to loop over (end year not included)
lat_range = [32.74,33.64];%[34.4,36];%[35.5,36.5];%[-90,90]; % latitude bounds
lon_range = [-116.2,-114.8];%[-120.5,-118];%[-118,-117];%[-180,180]; % longitude bounds
rm_eqs = 'out'; % events to remove: inside provided bounds ('in') or outside bounds ('out', default)

fdir = './zmap/data_files/scsn_'; % directory and initial filename (searches for '{fdir}{year}.mat')


%% COMPUTE Mc

allMc = nan(length(y_range)-1,3,size(mc_func,2));
n_eqs = nan(length(y_range)-1,1);
lrgM = [];

for y = 1:length(y_range)-1

    mag = [];

    for i = 0:(y_range(y+1)-y_range(y)-1)

        load(strcat(fdir,num2str(y_range(y)+i),'.mat')) % matfile with catalog in ZMAP format

        if (y_range(y)+i) < 1990 % for older data
            a(a(:,6)==0,:) = []; % remove rows with M=0 (1976 & 1982 have a lot)
        end

        % filter data by location
        if strcmpi(rm_eqs,'in') == 1 % remove data inside geographic bounds
            a((a(:,1)>lon_range(1)&a(:,1)<lon_range(2)&a(:,2)>lat_range(1)&a(:,2)<lat_range(2)),:) = [];
        else % remove data outside geographic bounds
            a((a(:,1)<lon_range(1)|a(:,1)>lon_range(2)|a(:,2)<lat_range(1)|a(:,2)>lat_range(2)),:) = [];
        end

        mag = [mag;a(:,6)]; % get magnitudes

        lrgM = [lrgM;a(a(:,6)>=6.5,[3,6])]; % keep year & mag of any large EQs 

    end

    n_eqs(y,:) = length(mag); % keep number of EQs in year group

    for m = 1:size(mc_func,2)

        %For mbass(), see algorithm Amorese [2007]
        if strcmpi(mc_func{m},'gft') == 1
            fname = 'Goodness-of-fit Test (GFT)';
            Mc_bootstrap = bootstrp(nbsample,@(samps) gft(samps,mbin),mag);

            allMc(y,1,m) = mean(cell2mat(Mc_bootstrap(:,1)),'omitnan');
            allMc(y,2,m) = std(cell2mat(Mc_bootstrap(:,1)), 'omitnan');
            allMc(y,3,m) = median(cell2mat(Mc_bootstrap(:,1)),'omitnan');

        elseif strcmpi(mc_func{m},'mbs') == 1
            fname = 'Mc by b-val Stability (MBS)';
            Mc_bootstrap = bootstrp(nbsample,@(samps) mbs(samps,mbin),mag);

            allMc(y,1,m) = mean(Mc_bootstrap,'omitnan');
            allMc(y,2,m) = std(Mc_bootstrap, 'omitnan');
            allMc(y,3,m) = median(Mc_bootstrap,'omitnan');

        elseif strcmpi(mc_func{m},'emr') == 1
            fname = 'Entire Magnitude Range method (EMR)';
            %when using emr(), the loop may break due to failure of nlsfit(),
            %in this case use:
            %Mc_bootstrap[i] = as.numeric(try(emr(sample(mag, replace=TRUE),mbin)$Mc))
            disp('EMR is not currently working')
            break;
        else % default to maxc
            fname = 'Maximum Curvature (MAXC)';
            [Mc_bootstrap,bs] = bootstrp(nbsample,@(samps) maxc(samps,mbin),mag);

            allMc(y,1,m) = mean(Mc_bootstrap,'omitnan');
            allMc(y,2,m) = std(Mc_bootstrap, 'omitnan');
            allMc(y,3,m) = median(Mc_bootstrap,'omitnan');
        end

        disp(strcat("Using the ",mc_func{m}," method for ",num2str(y_range(y))," to ",num2str(y_range(y+1))))
        disp(strcat("Mc (mean): ", num2str(allMc(y,1,m))))
        disp(strcat("Sigma0 (std. dev.): "," ", num2str(allMc(y,2,m))))

    end

end


%% PLOT RESULTS

figure;

for m = 1:size(mc_func,2)

    subplot(size(mc_func,2),1,m);
    hold on;
    
    % plot lines for major EQs
    for e = 1:size(lrgM,1)
        plot([lrgM(e,1),lrgM(e,1)],[0,4.5],'--','color',[0.5,0.5,0.5])
    end

    scatter(y_range(1:end-1),allMc(:,1,m),49,allMc(:,2,m),'filled','markeredgecolor','k')

    ylabel('Magnitude of Completeness (Mc)')
    ylim([0,4.5])
    c = colorbar('eastoutside');
    c.Label.String = 'Standard Deviation';
    fontsize(gcf,16,'points')
    grid on;
    box on;

    if strcmpi(mc_func{m},'gft') == 1
        title('Goodness-of-fit Test (GFT)');
    elseif strcmpi(mc_func{m},'mbs') == 1
        title('Mc by b-val Stability (MBS)');
    elseif strcmpi(mc_func{m},'emr') == 1
        title('Entire Magnitude Range method (EMR)');
    else % default to maxc
        title('Maximum Curvature (MAXC)');
    end

end

% plot number of earthquakes in each year group

figure;
hold on;

% plot lines for major EQs
for e = 1:size(lrgM,1)
    plot([lrgM(e,1),lrgM(e,1)],[0,4.5],'--','color',[0.5,0.5,0.5])
end

scatter(y_range(1:end-1),n_eqs,49,'r','filled','markeredgecolor','k');

ylabel('Number of Earthquakes')
fontsize(gcf,16,'points')
grid on;
box on;


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
