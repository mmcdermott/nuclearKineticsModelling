dataDir = '../data/histData/unequalSprings/';
dataFile = 'offCenter2t1';
run parameters.m

dataFile2 = 'offCenter5t1';
dataFile3 = 'offCenter10t1';
dataFile4 = 'offCenter20t1';
dataFile5 = 'offCenter100t1';
dataFileFull2 = [dataDir dataFile2 dataFileSuffix];
dataFileFull3 = [dataDir dataFile3 dataFileSuffix];
dataFileFull4 = [dataDir dataFile4 dataFileSuffix];
dataFileFull5 = [dataDir dataFile5 dataFileSuffix];

numRuns = 1000;
csvrange = [rowStart colStart rowStart+numRuns-1 colStart+2];
DATA1 = csvread(dataFileFull,rowStart,colStart,csvrange);
DATA2 = csvread(dataFileFull2,rowStart,colStart,csvrange);
DATA3 = csvread(dataFileFull3,rowStart,colStart,csvrange);
DATA4 = csvread(dataFileFull4,rowStart,colStart,csvrange);
DATA5 = csvread(dataFileFull5,rowStart,colStart,csvrange);

psi1 = abs((DATA1(:,3)-pi/2)*180/pi);
x1 = DATA1(:,1);
psi2 = abs((DATA2(:,3)-pi/2)*180/pi);
x2 = DATA2(:,1);
psi3 = abs((DATA3(:,3)-pi/2)*180/pi);
x3 = DATA3(:,1);
psi4 = abs((DATA4(:,3)-pi/2)*180/pi);
x4 = DATA4(:,1);
psi5 = abs((DATA5(:,3)-pi/2)*180/pi);
x5 = DATA5(:,1);

figure;
axes('Linewidth', 3.5);
nhist({psi1,psi2,psi3,psi4,psi5},'legend',{'Difference is 2-Fold', 'Difference is 5-fold', 'Difference is 10-fold', 'Difference is 20-fold', 'Difference is 100-fold'},'box','location','EastOutside','fsize',16,'linewidth',3);
hold on;
StartingAng1 = plot(abs((startPsi-pi/2)*180/pi), 0, '*');
title('Final Pronucleus Angle For Varying Spring Constants')
xlabel('Final Pronucleus Angle');

figure;
axes('Linewidth', 3.5);
nhist({x1,x2,x3,x4,x5},'legend',{'Difference is 2-Fold', 'Difference is 5-fold', 'Difference is 10-fold', 'Difference is 20-fold', 'Difference is 100-fold'},'box','location','EastOutside','fsize',16,'linewidth',3);
hold on;
StartingPos1 = plot(startX, 0, '*');
title('Final Pronucleus Position for Varying Spring Constants')
xlabel('Final Pronucleus x-Position Displacement');