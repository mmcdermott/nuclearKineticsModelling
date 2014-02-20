dataDir = '../data/histData/uncenteredOffAngleTendenciesBanded/';
dataFile = 'pPo8';
run parameters.m
dataFile2 = 'p2Po8';
dataFile3 = 'p3Po8';
dataFile4 = 'p4Po8';
dataFileFull2 = [dataDir dataFile2 dataFileSuffix];
dataFileFull3 = [dataDir dataFile3 dataFileSuffix];
dataFileFull4 = [dataDir dataFile4 dataFileSuffix];

numRuns = 1000;
csvrange = [rowStart colStart rowStart+numRuns-1 colStart+2];
DATA1 = csvread(dataFileFull,rowStart,colStart,csvrange);
DATA2 = csvread(dataFileFull2,rowStart,colStart,csvrange);
DATA3 = csvread(dataFileFull3,rowStart,colStart,csvrange);
DATA4 = csvread(dataFileFull4,rowStart,colStart,csvrange);
psi1 = abs((DATA1(:,3)-pi/2)*180/pi);
psi2 = abs((DATA2(:,3)-pi/2)*180/pi);
psi3 = abs((DATA3(:,3)-pi/2)*180/pi);
psi4 = abs((DATA4(:,3)-pi/2)*180/pi);

psi1Diff = psi1-22.5;
psi2Diff = psi2-45;
psi3Diff = psi3-67.5;
psi4Diff = psi4-90;

figure;
axes('Linewidth', 3.5);
nhist({psi1, psi2, psi3, psi4},'legend',{'Starting pi/8 Off Center', 'Starting 2pi/8 Off Center', 'Starting 3pi/8 Off Center', 'Starting 4pi/8 Off Center'},'box','location','EastOutside','fsize',16,'linewidth',3);
hold on;
StartingAng1 = plot(22.5, 0, '*');
StartingAng2 = plot(45, 0, '*');
StartingAng3 = plot(3*22.5, 0, '*');
StartingAng4 = plot(90, 0, '*');
title('Rotation Stabilization with No Translation')
xlabel('Final Pronucleus Angle');

figure;
axes('Linewidth', 3.5);
nhist({psi1Diff, psi2Diff, psi3Diff, psi4Diff},'legend',{'Starting pi/8 Off Center', 'Starting 2pi/8 Off Center', 'Starting 3pi/8 Off Center', 'Starting 4pi/8 Off Center'},'box','location','EastOutside','fsize',16,'linewidth',3);
hold on;
title('Rotation Stabilization with No Translation')
xlabel('Final Pronucleus Anglular Displacement');