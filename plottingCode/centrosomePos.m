run parameters

meanMx1 = [];
meanMy1 = [];
stdvMx1 = [];
stdvMy1 = [];

meanDx1 = [];
meanDy1 = [];
stdvDx1 = [];
stdvDy1 = [];

if (exist('numRuns','var') == 0)
    numRuns = 100;
end

numDataPoints = Duration*2;

csvrange = [rowStart colStart rowStart+numDataPoints*numRuns-1 colStart+7];
DATA = csvread(dataFileFull,rowStart,colStart,csvrange);

psi1Temp = abs((DATA(:,4)-pi/2)*180/pi);
x1Temp = (1/(2*R1_max))*DATA(:,2);
MTbaseMX1Temp = (1/(2*R1_max))*DATA(:,5);
MTbaseMY1Temp = (1/(2*R2_max))*DATA(:,6);
MTbaseDX1Temp = (1/(2*R1_max))*DATA(:,7);
MTbaseDY1Temp = (1/(2*R2_max))*DATA(:,8);

for slice = 1:numDataPoints
    meanMx1(slice) = mean(MTbaseMX1Temp(slice:numDataPoints:size(MTbaseMX1Temp)));
    stdvMx1(slice) = std(MTbaseMX1Temp(slice:numDataPoints:size(MTbaseMX1Temp)));
    meanMy1(slice) = mean(MTbaseMY1Temp(slice:numDataPoints:size(MTbaseMY1Temp)));
    stdvMy1(slice) = std(MTbaseMY1Temp(slice:numDataPoints:size(MTbaseMY1Temp)));

    meanDx1(slice) = mean(MTbaseDX1Temp(slice:numDataPoints:size(MTbaseDX1Temp)));
    stdvDx1(slice) = std(MTbaseDX1Temp(slice:numDataPoints:size(MTbaseDX1Temp)));
    meanDy1(slice) = mean(MTbaseDY1Temp(slice:numDataPoints:size(MTbaseDY1Temp)));
    stdvDy1(slice) = std(MTbaseDY1Temp(slice:numDataPoints:size(MTbaseDY1Temp)));
end

tData = DATA(:,1);
t = tData(1:20);

figure;
errorbarxy(meanMx1,meanMy1,stdvMx1,stdvMy1,{'ro','r','r'})
hold on;
errorbarxy(meanDx1,meanDy1,stdvDx1,stdvDy1,{'bo','b','b'})
xlabel('Fraction Egg Length in X');
ylabel('Fraction Egg Length in Y');
legend('Centrosome M','Centrosome D');
title('Average Centrosome Trajectories (Discrete)','FontSize', 16)

figure;
plot(meanMx1,meanMy1,'ro-')
hold on;
plot(meanDx1,meanDy1,'bo-')
xlabel('Fraction Egg Length in X');
ylabel('Fraction Egg Length in Y');
legend('Centrosome M','Centrosome D');
title('Average Centrosome Trajectories (Lined)','FontSize', 16)

figure;
errorbar(meanMx1,stdvMx1,'ro')
hold on;
errorbar(meanDx1,stdvDx1,'bo')
xlabel('Time (minutes)');
ylabel('Fraction Egg Length in X');
legend('Centrosome M','Centrosome D');
title('Change in X position')