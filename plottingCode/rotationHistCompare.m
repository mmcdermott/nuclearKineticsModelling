% nplot   = 1;
% figpath = ['graphs'];
% 
% run clean_fig_data_folder
if (exist('numRuns') == 0)
    numRuns = 100;
end
dataDir = '../data/hist-data/';
run parameters
dataFileFull2 = [dataDir dataFile2 dataFileSuffix];
csvrange = [rowStart colStart rowStart+numRuns-1 colStart+2]; 
DATA1 = csvread(dataFileFull,rowStart,colStart,csvrange);
DATA2 = csvread(dataFileFull2,rowStart,colStart,csvrange);
psi1 = abs((DATA1(:,3)-pi/2)*180/pi);
psi2 = abs((DATA2(:,3)-pi/2)*180/pi);

figure
dataPsi = nhist({psi1, psi2},'legend',{dataFile, dataFile2},'pdf','box','location','NorthOutside');
hold on;
StartingAng = plot(startPsi, 0, '*');

xlabel('|Pronucleus Rotation| (deg)');
ylabel('PDF of Runs');