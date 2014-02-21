meanMx1 = [];
meanMy1 = [];
meanMx2 = [];
meanMy2 = [];
stdvMx1 = [];
stdvMy1 = [];
stdvMx2 = [];
stdvMy2 = [];

dataDir = '../data/histData/centrosomeData/';
dataFile = 'offCenter';
dataFile2 = 'offCenterBanded';
run parameters.m
dataFileFull2 = [dataDir dataFile2 dataFileSuffix];

numRuns = 100;
csvrange = [rowStart colStart rowStart+20*numRuns-1 colStart+6];
DATA1 = csvread(dataFileFull,rowStart,colStart,csvrange);
DATA2 = csvread(dataFileFull2,rowStart,colStart,csvrange);

psi1Temp = abs((DATA1(:,4)-pi/2)*180/pi);
x1Temp = DATA1(:,2);
MTbaseMX1Temp = DATA1(:,5);
MTbaseMY1Temp = DATA1(:,6);
MTbaseDX1Temp = DATA1(:,7);
MTbaseDY1Temp = DATA1(:,8);

psi2Temp = abs((DATA2(:,4)-pi/2)*180/pi);
x2Temp = DATA2(:,2);
MTbaseMX2Temp = DATA2(:,5);
MTbaseMY2Temp = DATA2(:,6);
MTbaseDX2Temp = DATA2(:,7);
MTbaseDY2Temp = DATA2(:,8);

for slice = 1 : 20
    meanMx1(slice) = mean(x1temp(slice:20:size(x1temp)));
    stdvMx1(slice) = std(x1temp(slice:20:size(x1temp)));
    meanMy1(slice) = mean(y1temp(slice:20:size(y1temp)));
    stdvMy1(slice) = std(y1temp(slice:20:size(y1temp)));
    meanMx2(slice) = mean(x2temp(slice:20:size(x2temp)));
    stdvMx2(slice) = std(x2temp(slice:20:size(x2temp)));
    meanMy2(slice) = mean(y2temp(slice:20:size(y2temp)));
    stdvMy2(slice) = std(y2temp(slice:20:size(y2temp)));
end

tData = DATA1(:,1);
t = tData(1:20);

figure;


figure;
