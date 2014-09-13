% nplot   = 1;
% figpath = ['graphs'];
% 
% run clean_fig_data_folder
if (exist('numRuns') == 0)
    numRuns = 100;
end
dataDir = '../data/hist-data/';
run parameters
csvrange = [rowStart colStart rowStart+numRuns-1 colStart+2]; 
DATA1 = csvread(dataFileFull,rowStart,colStart,csvrange);
dataFileFull2 = [dataDir dataFile2 dataFileSuffix];
DATA2 = csvread(dataFileFull2,rowStart,colStart,csvrange);
x1 = DATA1(:,1);
y1 = DATA1(:,2);
psi1 = abs((DATA1(:,3)-pi/2)*180/pi);
x2 = DATA2(:,1);
y2 = DATA2(:,2);
psi2 = abs((DATA2(:,3)-pi/2)*180/pi);

figure
subplot(1,2,1);
hold on;
nhist({x1, x2},'legend',{dataFile, dataFile2},'pdf','box','location','NorthOutside');
StartingX = plot(startX, 0,'*');
%legend('Starting Position');
%text(0.1,0.1,dataX);
axis 'auto x'
xlabel('x-translation (mum)');
ylabel('PDF of Runs');


subplot(1,2,2);
dataPsi = nhist({psi1, psi2},'legend',{dataFile, dataFile2},'pdf','box','location','NorthOutside');
hold on;
StartingAng = plot(startPsi, 0, '*');

%axis([0 360 0 1])
%axis 'auto y'

xlabel('|Pronucleus Rotation| (deg)');
ylabel('PDF of Runs');
% %Bands:
% for region = 1:numRegions
%     regionStart = regionAngles(region)*180/pi;
%     regionEnd = regionAngles(region+1)*180/pi;
%     if (regionProbabilities(region) ~= 1)
%         Band=plot([regionStart,regionEnd],[0,0],'Color','b','LineWidth',4);
%     end
% end
% if (exist('Band') == 1)
%     legend([StartingAng Band],{'Starting Angle', '50% Probability Band'});
% else 
%     legend([StartingAng],{'Starting Angle'});
% end
