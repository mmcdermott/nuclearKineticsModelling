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
psi = DATA(:,3);

figure
subplot(2,1,1);
position = [x, y];
posHist = hist3(position,[18,18]);
posHist2 = posHist';
posHist2(size(posHist,1) + 1, size(posHist,2) + 1) = 0;
colormap(flipud(hot));
xb = linspace(min(position(:,1)),max(position(:,1)),size(posHist,1)+1);
yb = linspace(min(position(:,2)),max(position(:,2)),size(posHist,1)+1);
StartingPos = plot(startX, startY,'*','LineWidth',3);
hold on;
regionPlotter(0,0,R1_max,R2_max,regionAngles,regionProbabilities,regionForceMultipliers);
posHeatMap = pcolor(xb,yb,posHist2);

shading flat;
shading interp;

axis 'equal'
axis([-30,30,-15,15]);
xlabel('Final x-Position mum');
ylabel('Final y-Position mum');
colorbar;
legend([StartingPos],{'Starting Position'});


subplot(2,1,2);
psiRose = rose(psi,350);%,'legend',{dataFile},'pdf','box','location','NorthOutside');
colormap(flipud(hot));
x = get(psiRose, 'XData') ;
y = get(psiRose, 'YData') ;
%p = patch(x, y,'blue') ;
hold on;

startingAng = polar([startPsi,startPsi],[10,-10],'k-');
%startingAng = polar([startPsi,startPsi],[2500,0],'k-');
%startingAng = compass(2500*cos(startPsi), 2500*sin(startPsi), 'k-');

set(startingAng,'LineWidth',2);

xlabel('Final Rotation (deg)');
ylabel('CDF of Runs');
legend([startingAng],{'Starting Orientation'});
