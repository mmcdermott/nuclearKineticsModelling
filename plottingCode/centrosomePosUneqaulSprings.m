meanMx1 = [];
meanMy1 = [];
meanMx2 = [];
meanMy2 = [];
stdvMx1 = [];
stdvMy1 = [];
stdvMx2 = [];
stdvMy2 = [];

dataDir = '../data/histData/centrosomeData/';
dataFile = '20t1OffCenter';
dataFile2 = '20t1OffCenterBanded';
run parameters.m
dataFileFull2 = [dataDir dataFile2 dataFileSuffix];

numRuns = 100;
csvrange = [rowStart colStart rowStart+20*numRuns-1 colStart+7];
DATA1 = csvread(dataFileFull,rowStart,colStart,csvrange);
DATA2 = csvread(dataFileFull2,rowStart,colStart,csvrange);

psi1Temp = abs((DATA1(:,4)-pi/2)*180/pi);
x1Temp = (1/(2*R1_max))*DATA1(:,2);
MTbaseMX1Temp = (1/(2*R1_max))*DATA1(:,5);
MTbaseMY1Temp = (1/(2*R2_max))*DATA1(:,6);
MTbaseDX1Temp = (1/(2*R1_max))*DATA1(:,7);
MTbaseDY1Temp = (1/(2*R2_max))*DATA1(:,8);

psi2Temp = abs((DATA2(:,4)-pi/2)*180/pi);
x2Temp = (1/(2*R1_max))*DATA2(:,2);
MTbaseMX2Temp = (1/(2*R1_max))*DATA2(:,5);
MTbaseMY2Temp = (1/(2*R2_max))*DATA2(:,6);
MTbaseDX2Temp = (1/(2*R1_max))*DATA2(:,7);
MTbaseDY2Temp = (1/(2*R2_max))*DATA2(:,8);

for slice = 1 : 20
    meanMx1(slice) = mean(MTbaseMX1Temp(slice:20:size(MTbaseMX1Temp)));
    stdvMx1(slice) = std(MTbaseMX1Temp(slice:20:size(MTbaseMX1Temp)));
    meanMy1(slice) = mean(MTbaseMY1Temp(slice:20:size(MTbaseMY1Temp)));
    stdvMy1(slice) = std(MTbaseMY1Temp(slice:20:size(MTbaseMY1Temp)));
    meanMx2(slice) = mean(MTbaseMX2Temp(slice:20:size(MTbaseMX2Temp)));
    stdvMx2(slice) = std(MTbaseMX2Temp(slice:20:size(MTbaseMX2Temp)));
    meanMy2(slice) = mean(MTbaseMY2Temp(slice:20:size(MTbaseMY2Temp)));
    stdvMy2(slice) = std(MTbaseMY2Temp(slice:20:size(MTbaseMY2Temp)));

    meanDx1(slice) = mean(MTbaseDX1Temp(slice:20:size(MTbaseDX1Temp)));
    stdvDx1(slice) = std(MTbaseDX1Temp(slice:20:size(MTbaseDX1Temp)));
    meanDy1(slice) = mean(MTbaseDY1Temp(slice:20:size(MTbaseDY1Temp)));
    stdvDy1(slice) = std(MTbaseDY1Temp(slice:20:size(MTbaseDY1Temp)));
    meanDx2(slice) = mean(MTbaseDX2Temp(slice:20:size(MTbaseDX2Temp)));
    stdvDx2(slice) = std(MTbaseDX2Temp(slice:20:size(MTbaseDX2Temp)));
    meanDy2(slice) = mean(MTbaseDY2Temp(slice:20:size(MTbaseDY2Temp)));
    stdvDy2(slice) = std(MTbaseDY2Temp(slice:20:size(MTbaseDY2Temp)));
end

tData = DATA1(:,1);
t = tData(1:20);

figure;
errorbarxy(meanMx1,meanMy1,stdvMx1,stdvMy1,{'ro','r','r'})
hold on;
errorbarxy(meanDx1,meanDy1,stdvDx1,stdvDy1,{'bo','b','b'})
xlabel('Fraction Egg Length in X');
ylabel('Fraction Egg Length in Y');
legend('Centrosome M','Centrosome D (20-fold stronger spring)');
title('Average Centrosome Trajectories (Discrete)--No Band, Unequal Springs','FontSize', 16)

figure;
plot(meanMx1,meanMy1,'ro-')
hold on;
plot(meanDx1,meanDy1,'bo-')
xlabel('Fraction Egg Length in X');
ylabel('Fraction Egg Length in Y');
legend('Centrosome M','Centrosome D (20-fold stronger spring)');
title('Average Centrosome Trajectories (Lined)--No Band, Unequal Springs','FontSize', 16)

figure;
errorbar(meanMx1,stdvMx1,'ro')
hold on;
errorbar(meanDx1,stdvDx1,'bo')
xlabel('Time (minutes)');
ylabel('Fraction Egg Length in X');
legend('Centrosome M','Centrosome D (20-fold stronger spring)');
title('Change in X position--No Band')

figure;
errorbarxy(meanMx2,meanMy2,stdvMx2,stdvMy2,{'ro','r','r'})
hold on;
errorbarxy(meanDx2,meanDy2,stdvDx2,stdvDy2,{'bo','b','b'})
xlabel('Fraction Egg Length in X');
ylabel('Fraction Egg Length in Y');
legend('Centrosome M','Centrosome D (20-fold stronger spring)');
title('Average Centrosome Trajectories (Discrete)--Banded, Unequal Springs','FontSize', 16)

figure;
plot(meanMx2,meanMy2,'ro-')
hold on;
plot(meanDx2,meanDy2,'bo-')
xlabel('Fraction Egg Length in X');
ylabel('Fraction Egg Length in Y');
legend('Centrosome M','Centrosome D (20-fold stronger spring)');
title('Average Centrosome Trajectories (Lined)--Banded, Unequal Springs','FontSize', 16)

figure;
errorbar(meanMx2,stdvMx2,'ro')
hold on;
errorbar(meanDx2,stdvDx2,'bo')
xlabel('Time (minutes)');
ylabel('Fraction Egg Length in X');
legend('Centrosome M','Centrosome D (20-fold stronger spring)');
title('Change in X position--Banded')