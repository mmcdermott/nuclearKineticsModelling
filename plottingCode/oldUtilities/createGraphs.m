%TODO: Add annotations

do2000 = false;
nplot   = 1;
figpath = ['graphs'];

run clean_fig_data_folder
run parameters

color2 = 'blue';
color20 = 'green';
color200 = 'red';
color2000 = 'black';
%color20000 = 'purple';

N = Duration/Tau;
csvrange = [rowStart colStart rowStart+N-1 colStart+3]; 
dataFileSuffix = '-MT.csv';
dataFile2MT = [dataDir '2' dataFileSuffix];
dataFile20MT = [dataDir '20' dataFileSuffix];
dataFile200MT = [dataDir '200' dataFileSuffix];
dataFile2000MT = [dataDir '2000' dataFileSuffix];

figure
subplot(2,1,1);
StartingAng = plot(0,startPsi,'*');
hold on;
axis([0 Duration 0 360]);
xlabel('Time (minutes)');
ylabel('|Pronucleus Rotation| (deg)');
%Bands:
for region = 1:numRegions
    regionStart = regionAngles(region)*180/pi;
    regionEnd = regionAngles(region+1)*180/pi;
    if (regionProbabilities(region) ~= 1)
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

%Old: 
%%Cell: 
%regionPlotter(0,0,R1_max,R2_max,regionAngles - pi/2);
%hold on;
%ellipse(0,0,R1_max,R2_max,'k');

%axis equal;
%xlabel('X Position (mum)');
%ylabel('Y Position (mum)');
%hold on;
%plot(startX,startY,'*');
%hold on;

%%%%===========2 MT=====================================
DATA = csvread(dataFile2MT,rowStart,colStart,csvrange);
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

%%%%===========20 MT=====================================
DATA = csvread(dataFile20MT,rowStart,colStart,csvrange);
disp('Read Data 20!');

x = DATA(:,2);
y = DATA(:,3);
psi = abs((DATA(:,4)-pi/2)*180/pi);


subplot(2,1,1);
psi20 = plot(time, psi, 'Color', color20);
hold on;

subplot(2,1,2);
x20 = plot(time, x, 'Color', color20);
hold on;

%%%%===========200 MT=====================================
DATA = csvread(dataFile200MT,rowStart,colStart,csvrange);
disp('Read Data 200!');

x = DATA(:,2);
y = DATA(:,3);
psi = abs((DATA(:,4)-pi/2)*180/pi);


subplot(2,1,1);
psi200 = plot(time, psi, 'Color', color200);
hold on;

subplot(2,1,2);
x200 = plot(time, x, 'Color', color200);
hold on;

if (do2000)
    %%%===========2000 MT=====================================
    DATA = csvread(dataFile2000MT,rowStart,colStart,csvrange);
    disp('Read Data 2000!');

    x = DATA(:,2);
    y = DATA(:,3);
    psi = abs((DATA(:,4)-pi/2)*180/pi);


    subplot(2,1,1);
    psi2000 = plot(time, psi, 'Color', color2000);
    hold on;

    subplot(2,1,2);
    x2000 = plot(time, x, 'Color', color2000);
    hold on;
end

subplot(2,1,1);
title('x-Position of Center and Rotation of Pronucleus');
legend([StartingAng psi2 psi20 psi200 Band], {'Starting Angle','2 Total MTs', '20 Total MTs', '200 Total MTs', '50% Probability Band'});%, '2000 Total MTs');

subplot(2,1,2);
legend([StartingX x2 x20 x200], {'Starting Position','2 Total MTs', '20 Total MTs', '200 Total MTs'});%, '2000 Total MTs');
