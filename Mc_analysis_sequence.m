%% Script for calculating Mc with different methods for specific earthquake sequences

% Currently loads and uses a matfile with a ZMAP catalog (named a) for each year
% formatted as: [lon, lat, year, month, day, mag, depth, hour, minute]

% Converted from R algorithms available in:
% http://www.corssa.org/export/sites/corssa/.galleries/articles-pdf/Mignan-Woessner-2012-CORSSA-Magnitude-of-completeness.pdf

% Methods included (see above document for details):
% Max curvature (maxc)
% Goodness-of-fit for b-value (gft)
% Mc by b-val Stability (mbs)
% Entire Magnitude Range method (emr) - not currently functional (still in R - not converted)

% Created: 4/4/2023
% Updated: 5/18/2023
% By: Gabrielle Tepp, Caltech

%-------------------

%% INPUT PARAMETERS

mbin = 0.1; %Magnitude bin size
nbsample = 200; %Bootstrapping

%select function: maxc (default), gft, mbs, emr (can use multiple - comma separate in cell array)
mc_func = {'mbs'};%{'maxc','gft','mbs'};

dt_min = '2019-07-04';%'2022-01-01';%'2010-04-04';%'2019-07-04';
dt_max = '2023-02-04';%'2023-01-01';%'2015-04-04';%'2023-01-01';
lat_range = [35.45,36.25];%[-90,90];%[31.75,32.75];%[-90,90]; %[32.74,33.64];%[34.4,36];%[35.5,36.5];% latitude bounds
lon_range = [-118.15,-117.25];%[-180,180];%[-115.63,180];%[-115.8,-114.9];%[-180,180]; %[-116.2,-114.8];%[-120.5,-118];%[-118,-117];% longitude bounds

fdir = './zmap/data_files/scsn_'; % directory and initial filename (searches for '{fdir}{year}.mat')


%% COMPUTE Mc

allMc = nan(3,size(mc_func,2));

alleqs = [];

% get time info
[st_yr,~,~] = datevec(dt_min);
[end_yr,~,~] = datevec(dt_max);

% get events for the desired years
for i = st_yr:end_yr

    load(strcat(fdir,num2str(i),'.mat')) % matfile with catalog in ZMAP format

    if i < 1990 % for older data
        a(a(:,6)==0,:) = []; % remove rows with M=0 (1976 & 1982 have a lot)
    end

    % remove data outside geographic bounds
    a((a(:,1)<lon_range(1)|a(:,1)>=lon_range(2)|a(:,2)<lat_range(1)|a(:,2)>lat_range(2)),:) = [];

    alleqs = [alleqs;a]; % add events to master list

end

% remove events outside time bounds
alleqs(:,10) = datenum([alleqs(:,[3:5,8:9]),zeros(length(alleqs),1)]);
alleqs((alleqs(:,10)<datenum(dt_min)|alleqs(:,10)>datenum(dt_max)),:) = [];

alleqs = sortrows(alleqs,6,'ascend');

% calculate Mc
for m = 1:size(mc_func,2)

    %For mbass(), see algorithm Amorese [2007]
    if strcmpi(mc_func{m},'gft') == 1
        fname = 'Goodness-of-fit Test (GFT)';
        Mc_bootstrap = bootstrp(nbsample,@(samps) gft(samps,mbin),alleqs(:,6));

        allMc(1,m) = mean(cell2mat(Mc_bootstrap(:,1)),'omitnan');
        allMc(2,m) = std(cell2mat(Mc_bootstrap(:,1)), 'omitnan');
        allMc(3,m) = median(cell2mat(Mc_bootstrap(:,1)),'omitnan');

    elseif strcmpi(mc_func{m},'mbs') == 1
        fname = 'Mc by b-val Stability (MBS)';
        Mc_bootstrap = bootstrp(nbsample,@(samps) mbs(samps,mbin),alleqs(:,6)); % outputs: {Mc,Mco,bi,unc,bave};

        allMc(1,m) = mean(cell2mat(Mc_bootstrap(:,1)),'omitnan');
        allMc(2,m) = std(cell2mat(Mc_bootstrap(:,1)), 'omitnan');
        allMc(3,m) = median(cell2mat(Mc_bootstrap(:,1)),'omitnan');

    elseif strcmpi(mc_func{m},'emr') == 1
        fname = 'Entire Magnitude Range method (EMR)';
        %when using emr(), the loop may break due to failure of nlsfit(),
        %in this case use:
        %Mc_bootstrap[i] = as.numeric(try(emr(sample(mag, replace=TRUE),mbin)$Mc))
        disp('EMR is not currently working')
        break;
    else % default to maxc
        fname = 'Maximum Curvature (MAXC)';
        [Mc_bootstrap,bs] = bootstrp(nbsample,@(samps) maxc(samps,mbin),alleqs(:,6));

        allMc(1,m) = mean(Mc_bootstrap,'omitnan');
        allMc(2,m) = std(Mc_bootstrap, 'omitnan');
        allMc(3,m) = median(Mc_bootstrap,'omitnan');
    end

    disp(strcat("Using the ",mc_func{m}," method for ",dt_min," to ",dt_max))
    disp(strcat("Mc (mean): ", num2str(allMc(1,m))))
    disp(strcat("Sigma0 (std. dev.): "," ", num2str(allMc(2,m))))

end

% get FMD (frequency magnitude distribution)
res = fmd(alleqs(:,6),mbin);


%% PLOT RESULTS

figure;
hold on;

yyaxis left
% plot Mc's
for m = 1:size(mc_func,2)
    plot([allMc(1,m),allMc(1,m)],[0,round(max(log(res(:,2))),1,'significant')],'--','color',[0.5,0.5,0.5])
end
scatter(res(:,1),log(res(:,2)),64,'b','filled','markeredgecolor','k')
ylabel('Cumulative Magnitude (log)')

yyaxis right
scatter(res(:,1),res(:,3),64,'r','filled','markeredgecolor','k')
ylabel('Non-cumulative Magnitude')

fontsize(gcf,16,'points')
grid on;
box on;

% plot earthquakes

figure;

geoscatter(alleqs(:,2),alleqs(:,1),49,alleqs(:,6),'filled','markeredgecolor','k');

colormap parula
c = colorbar('eastoutside');
c.Label.String = 'Magnitude';

set(gca,'fontsize',16)
box on;
grid on;



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
function params = mbs(mag,mbin)

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

% put into one variable for output
% params.Mc = Mc;
% params.Mco = Mco;
% params.bi = bi;
% params.unc = unc;
% params.bave = bave;
params = {Mc,Mco,bi,unc,bave};

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
