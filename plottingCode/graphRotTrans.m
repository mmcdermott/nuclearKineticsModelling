nplot   = 1;
figpath = ['graphs'];

run clean_fig_data_folder
run parameters

color = 'blue';
N = Duration/Tau;
csvrange = [rowStart colStart rowStart+N-1 colStart+3];


figure
title('Position (x,y) of Center and Rotation of Pronucleus');
subplot(2,1,1);
plot(0,startPsi,'*');
hold on;
axis([0 Duration 0 2*pi]);
xlabel('Time (minutes)');
ylabel('Angle of Pronucleus Rotation (radians)');

subplot(2,1,2);
ellipse(0,0,R1_max,R2_max,'k');
regionPlotter(0,0,R1_max,R2_max,regionAngles);
axis equal;
xlabel('X Position (mum)');
ylabel('Y Position (mum)');
hold on;
plot(startX,startY,'*');
hold on;

DATA = csvread(dataFile,rowStart,colStart,csvrange);
disp('Read Data!');

time = DATA(:,1);
x = DATA(:,2);
y = DATA(:,3);
psi = DATA(:,4);

subplot(2,1,1);
plot(time, psi, 'Color', color);
hold on;

subplot(2,1,2);
plot(x,y, 'Color', color);
hold on;