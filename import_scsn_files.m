%% Import data from text file for zmap format
%
% Script for importing data from the SCSN/SCEDC yearly catalog text files with the following name format:
%
% filename format: path/{year}.catalog
%
% And putting into zmap format:
% 
% a = [lon, lat, year, month, day, mag, depth, hour, minute]
%
% Auto-generated by MATLAB on 16-Mar-2023 15:42:50
% edited by G. Tepp on 4/4/2023

%% Set up the Import Options and import the data

% specify directory
dir = "path to {year}.catalog files";

opts = delimitedTextImportOptions("NumVariables", 14);

% Specify range and delimiter
opts.DataLines = [11, Inf];
opts.Delimiter = " ";

% Specify column names and types
opts.VariableNames = ["YYYMMDD", "HHmmSSss", "ET", "GT", "MAG", "M", "LAT", "LON", "DEPTH", "Q", "EVID", "NPH", "NGRM", "Var14"];
opts.SelectedVariableNames = ["YYYMMDD", "HHmmSSss", "ET", "GT", "MAG", "M", "LAT", "LON", "DEPTH", "Q", "EVID", "NPH", "NGRM"];
opts.VariableTypes = ["string", "string", "string", "string", "double", "categorical", "double", "double", "double", "categorical", "double", "double", "double", "string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.ConsecutiveDelimitersRule = "join";
opts.LeadingDelimitersRule = "ignore";

% Specify variable properties
opts = setvaropts(opts, ["YYYMMDD", "HHmmSSss", "Var14"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["YYYMMDD", "HHmmSSss", "ET", "GT", "M", "Q", "Var14"], "EmptyFieldRule", "auto");

% Import the data

for y = 1932:1999%2000:2023

    tbl = readtable(strcat(dir,num2str(y),".catalog"), opts);

    %% Convert to output type
    DATE = tbl.YYYMMDD;
    TIME = tbl.HHmmSSss;
    ET = tbl.ET;
    GT = tbl.GT;
    MAG = tbl.MAG;
    M = tbl.M;
    LAT = tbl.LAT;
    LON = tbl.LON;
    DEPTH = tbl.DEPTH;
    Q = tbl.Q;
    EVID = tbl.EVID;
    NPH = tbl.NPH;
    NGRM = tbl.NGRM;

    %% Put into ZMAP catalog format

    a = nan(size(DATE,1)-2,9);

    for r = 1:size(DATE,1)-2

        if strcmpi(ET(r),'eq') == 1 && strcmpi(GT(r),'l') == 1 % only keep local earthquakes

        tempD = char(DATE(r));
        tempT = char(TIME(r));

        a(r,:) = [LON(r),LAT(r),str2num(tempD(1:4)),str2num(tempD(6:7)),str2num(tempD(9:10)),MAG(r),DEPTH(r),...
            str2num(tempT(1:2)),str2num(tempT(4:5))];

        end

    end

    a(any(isnan(a), 2), :) = []; % clear NaN rows

    save(strcat('scsn_',num2str(y),'.mat'),'a');

end

%clear all