dataDir = '../data/baseTests/';
dataFile = '2p5Offset';
run parameters.m
dataFile2 = '5Offset';
dataFile3 = '7p5Offset';
dataFile4 = '10Offset';
dataFileFull2 = [dataDir dataFile2 dataFileSuffix];
dataFileFull3 = [dataDir dataFile3 dataFileSuffix];
dataFileFull4 = [dataDir dataFile4 dataFileSuffix];

numRuns = 1000;
csvrange = [rowStart colStart rowStart+numRuns-1 colStart+2]; 
DATA1 = csvread(dataFileFull,rowStart,colStart,csvrange);
DATA2 = csvread(dataFileFull2,rowStart,colStart,csvrange);
DATA3 = csvread(dataFileFull3,rowStart,colStart,csvrange);
DATA4 = csvread(dataFileFull4,rowStart,colStart,csvrange);

x1 = DATA1(:,1);
x2 = DATA2(:,1);
x3 = DATA3(:,1);
x4 = DATA4(:,1);

diffX1 = 2.5-x1;
diffX2 = 5-x2;
diffX3 = 7.5-x3;
diffX4 = 10-x4;

figure;
axes('Linewidth', 3.5);
nhist({x1, x2, x3, x4},'legend',{'Starting at 2.5', 'Starting at 5',...
                                 'Starting at 7.5', 'Starting at 10'},...
                               'box','location','EastOutside','fsize',16,...
                               'linewidth',3);
hold on;
xlabel('Final Position')
StartingPosition1 = plot(2.5, 0, '*');
StartingPosition2 = plot(5, 0, '*');
StartingPosition3 = plot(7.5, 0, '*');
StartingPosition4 = plot(10, 0, '*');
title('Central Stable Equilibrium Tendencies','FontSize',18)

figure;
axes('Linewidth', 3.5);
nhist({diffX1, diffX2, diffX3, diffX4},...
      'legend',{'Starting at 2.5', 'Starting at 5','Starting at 7.5',...
                'Starting at 10'},'box','location','EastOutside','fsize',16,...
                'linewidth',3);
hold on;
xlabel('Displacement from the Initial Positions')
title('Central Stable Equilibrium Tendencies - Expanding Returns','FontSize',18)
