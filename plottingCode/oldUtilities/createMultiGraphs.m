%TODO: Add annotations

nplot   = 1;
figpath = ['graphs'];

run clean_fig_data_folder
dataFileSuffix = '.csv';
dataFile1 = [dataDir '200-MT-1' dataFileSuffix];
dataFile2 = [dataDir '200-MT-2' dataFileSuffix];
dataFile3 = [dataDir '200-MT-3' dataFileSuffix];
dataFile4 = [dataDir '200-MT-4' dataFileSuffix];
dataFile5 = [dataDir '200-MT-5' dataFileSuffix];
dataFile = dataFile1;
run parameters

color1 = 'blue';
color2 = 'green';
color3 = 'red';
color4 = 'black';
color5 = 'yellow';

N = Duration/Tau;
csvrange = [rowStart colStart rowStart+N-1 colStart+3]; 

figure
subplot(2,1,1);
StartingAng = plot(0,startPsi,'*');
hold on;
axis([0 Duration 0 360]);
xlabel('Time (minutes)');
ylabel('|Pronucleus Rotation| (deg)');
%Bands:
useBand = false;
for region = 1:numRegions
    regionStart = regionAngles(region)*180/pi;
    regionEnd = regionAngles(region+1)*180/pi;
    if (regionProbabilities(region) ~= 1)
        useBand = true;
        Band=plot([0,0],[regionStart,regionEnd],'Color','b','LineWidth',4);
        hold on;
    end
end


subplot(2,1,2);
StartingX = plot(0,startX, '*');
hold on;
axis([0 Duration -50 50]);
xlabel('Time (minutes)');
ylabel('x-translation (mum)');

%%%%===========File 1=====================================
DATA = csvread(dataFile1,rowStart,colStart,csvrange);
disp('Read Data 1!');

time = DATA(:,1);
x = DATA(:,2);
y = DATA(:,3);
psi = abs((DATA(:,4)-pi/2)*180/pi);

subplot(2,1,1);
psi1 = plot(time, psi, 'Color', color1);
hold on;

subplot(2,1,2);
x1 = plot(time,x, 'Color', color1);
hold on;

%%%%===========File 2=====================================
DATA = csvread(dataFile2,rowStart,colStart,csvrange);
disp('Read Data 2!');

time = DATA(:,1);
x = DATA(:,2);
y = DATA(:,3);
psi = abs((DATA(:,4)-pi/2)*180/pi);

subplot(2,1,1);
psi2 = plot(time, psi, 'Color', color2);
hold on;

subplot(2,1,2);
x2 = plot(time,x, 'Color', color2);
hold on;

%%%%===========File 3=====================================
DATA = csvread(dataFile3,rowStart,colStart,csvrange);
disp('Read Data 3!');

time = DATA(:,1);
x = DATA(:,2);
y = DATA(:,3);
psi = abs((DATA(:,4)-pi/2)*180/pi);

subplot(2,1,1);
psi3 = plot(time, psi, 'Color', color3);
hold on;

subplot(2,1,2);
x3 = plot(time,x, 'Color', color3);
hold on;

%%%%===========File 4=====================================
DATA = csvread(dataFile4,rowStart,colStart,csvrange);
disp('Read Data 4!');

time = DATA(:,1);
x = DATA(:,2);
y = DATA(:,3);
psi = abs((DATA(:,4)-pi/2)*180/pi);

subplot(2,1,1);
psi4 = plot(time, psi, 'Color', color4);
hold on;

subplot(2,1,2);
x4 = plot(time,x, 'Color', color4);
hold on;

%%%%===========File 5=====================================
DATA = csvread(dataFile1,rowStart,colStart,csvrange);
disp('Read Data 5!');

time = DATA(:,1);
x = DATA(:,2);
y = DATA(:,3);
psi = abs((DATA(:,4)-pi/2)*180/pi);

subplot(2,1,1);
psi5 = plot(time, psi, 'Color', color5);
hold on;

subplot(2,1,2);
x5 = plot(time,x, 'Color', color5);
hold on;
%%%%============Legend===================================
subplot(2,1,1);
title('x-Position of Center and Rotation of Pronucleus');
if (useBand)
  legend([StartingAng psi1 psi2 psi3 psi4 psi5 Band], {'Starting Angle','Run 1','Run 2','Run 3','Run 4','Run 5', '50% Probability Band'});
else
  legend([StartingAng psi1 psi2 psi3 psi4 psi5], {'Starting Angle','Run 1','Run 2','Run 3','Run 4','Run 5'});
end

subplot(2,1,2);
if (useBand)
  legend([StartingX x1 x2 x3 x4 x5 Band], {'Starting Position','Run 1','Run 2','Run 3','Run 4','Run 5', '50% Probability Band'});
else
  legend([StartingX x1 x2 x3 x4 x5], {'Starting Position','Run 1','Run 2','Run 3','Run 4','Run 5'});
end
