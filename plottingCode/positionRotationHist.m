% nplot   = 1;
% figpath = ['graphs'];
% 
% run clean_fig_data_folder
if (exist('numRuns') == 0)
    numRuns = 100;
end
dataDir = '../data/histData/';
run parameters
csvrange = [rowStart colStart rowStart+numRuns-1 colStart+2]; 
DATA = csvread(dataFileFull,rowStart,colStart,csvrange);
x = DATA(:,1);
y = DATA(:,2);
psi = abs((DATA(:,3)-pi/2)*180/pi);

figure
subplot(1,2,1);
dataX = nhist(x,'legend',{dataFile},'pdf','box','location','NorthOutside');
hold on;
StartingX = plot(startX, 0,'*');
%text(0.1,0.1,dataX);
axis 'auto x'
xlabel('x-translation (mum)');
ylabel('PDF of Runs');
%legend([StartingX],{'Starting Position'});


subplot(1,2,2);
dataPsi = nhist(psi,'legend',{dataFile},'pdf','box','location','NorthOutside');
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
