% nplot   = 1;
% figpath = ['graphs'];
% 
% run clean_fig_data_folder
if (exist('numRuns','var') == 0)
    numRuns = 100;
end
if (exist('dataDir','var') == 0)
    dataDir = '../data/histData/';
end
run parameters
csvrange = [rowStart colStart rowStart+numRuns-1 colStart+2]; 
DATA = csvread(dataFileFull,rowStart,colStart,csvrange);
x = DATA(:,1);
y = DATA(:,2);
psi = abs((DATA(:,3)-pi/2)*180/pi);

figure
subplot(2,2,1);
dataX = nhist(x,'legend',{dataFile},'pdf','box','location','NorthOutside');
hold on;
StartingX = plot(startX, 0,'*');
axis 'auto x'
xlabel('x-translation (mum)');
ylabel('PDF of Runs');
%legend([StartingX],{'Starting Position'});

subplot(2,2,2);
dataX = nhist(y,'legend',{dataFile},'pdf','box','location','NorthOutside');
hold on;
StartingY = plot(startY, 0,'*');
axis 'auto x'
xlabel('y-translation (mum)');
ylabel('PDF of Runs');
%legend([StartingX],{'Starting Position'});


subplot(2,2,3);
dataPsi = nhist(psi,'legend',{dataFile},'pdf','box','location','NorthOutside');
hold on;
StartingAng = plot(startPsi, 0, '*');
xlabel('|Pronucleus Rotation| (deg)');
ylabel('PDF of Runs');

