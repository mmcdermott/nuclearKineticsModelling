dataDir = '../data/histData/unequalSpringsBanded/';
dataFile = '100t1';
run parameters.m

numRuns = 1000;
csvrange = [rowStart colStart rowStart+numRuns-1 colStart+2];
DATA1 = csvread(dataFileFull,rowStart,colStart,csvrange);

psi1 = abs((DATA1(:,3)-pi/2)*180/pi);
x1 = DATA1(:,1);

figure;
axes('Linewidth', 3.5);
nhist({psi1},'legend',{'Springs constant differ by a factor of 100'},'box','location','EastOutside','fsize',16,'linewidth',3);
hold on;
StartingAng1 = plot(0, 0, '*');
title('Final Pronucleus Angle For Drastically Different Spring Constants')
xlabel('Final Pronucleus Angle');

figure;
axes('Linewidth', 3.5);
nhist({x1},'legend',{'Springs constant differ by a factor of 100'},'box','location','EastOutside','fsize',16,'linewidth',3);
hold on;
StartingPos1 = plot(0, 0, '*');
title('Rotation Stabilization with No Translation')
xlabel('Final Pronucleus x-Position Displacement');