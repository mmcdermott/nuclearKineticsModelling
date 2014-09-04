%Parameters: 
%  Data Parameters: 
dataFileSuffix = '.csv';
if (exist('dataDir','var') == 0)
    dataDir = '../data/';
end
if (exist('dataFile','var') == 0) 
    dataFileFull = [dataDir '2-MT' dataFileSuffix];
else 
    dataFileFull = [dataDir dataFile dataFileSuffix];
end

% Grabbing Parameters from Output:

modelParams = csvread(dataFileFull, 0,0, [0,0, 0,12]);
F_MT = modelParams(1);
contact_length = modelParams(2);
Eta = modelParams(3);
Mu = modelParams(4);
Vg = modelParams(5);
Vs_c = modelParams(6);
Vs = modelParams(7);
kc = modelParams(8);
kr = modelParams(9);
translation = modelParams(10) == 1;
startX = modelParams(11);
startY = modelParams(12);
startPsi = modelParams(13);

PlottingParams = csvread(dataFileFull, 1, 0, [1,0, 1,6]);

mt_numb_m = PlottingParams(1);
mt_numb_d = PlottingParams(2);
R1_max    = PlottingParams(3);
R2_max    = PlottingParams(4);
Prad      = PlottingParams(5);
Duration  = PlottingParams(6);
Tau       = PlottingParams(7);

numRegions = csvread(dataFileFull, 2, 0, [2,0, 2,0]);
regionAngles = csvread(dataFileFull, 3, 0, [3,0, 3,numRegions]);
regionProbabilities = csvread(dataFileFull, 4, 0, [4, 0, 4, numRegions-1]);
regionForceMultipliers = csvread(dataFileFull, 5, 0, [5, 0, 5, numRegions-1]);

%  Initialization of Plotting Parameters: 
x    = startX;
y    = startY;
psi  = startPsi;

cosinePrt = Prad*cos(psi);
sinePrt   = Prad*sin(psi);

xP_m = x + cosinePrt;
xP_d = x - cosinePrt;
yP_m = y + sinePrt;
yP_d = y - sinePrt;

%  Data Fromat Params:
rowStart = 6;
colStart = 0;
