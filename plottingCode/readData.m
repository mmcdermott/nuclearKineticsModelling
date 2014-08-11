run parameters

nplot   = 1;
figpath = ['figs'];

run clean_fig_data_folder

figure('visible','on');
subplot(2, 2, 1);
hold on;

axis equal
xlim([-R1_max R1_max])
ylim([-R2_max R2_max]) 
xlabel('x (nm)','fontsize',16)
ylabel('y (nm)','fontsize',16) 
% plot boundary and interior points

regionPlotter(0,0,R1_max,R2_max,regionAngles,regionProbabilities,regionForceMultipliers)
ellipse(0,0,R1_max,R2_max,'k')
circle_1(x,y,Prad,'g')
circle_1(xP_m,yP_m,Prad/8,'r')
circle_1(xP_d,yP_d,Prad/8,'b')

title('Initial Configuration','fontsize',16)
file_print=[figpath '/init.jpg'];
print(gcf,'-djpeg',file_print);
pause(0.1)

hold off;

DATA = csvread(dataFileFull,rowStart,colStart);
disp('Read Data!');
t         = 0;
tOld      = 0;
mag       = 0;
magOld    = 0;
theta     = 0;
thetaOld  = 0;
%torque    = 0;
%torqueOld = 0;

time = DATA(:,1);

start_M = 5;
end_M = start_M + 2*mt_numb - 1;
start_D = end_M + 1;
end_D = start_D + 2*mt_numb - 1;

forceStart = end_D + 5;
forceEnd   = end_D + 6;
forceFull  = DATA(:,forceStart:forceEnd);
forceMag   = sqrt(sum(forceFull.^2,2));
forceAng   = atan2(forceFull(:,2),forceFull(:,1));

torqueData = DATA(:,end_D + 9);

startBasePosM = end_D + 10;
endBasePosM   = startBasePosM + 1;
startBasePosD = endBasePosM + 1;
endBasePosD   = startBasePosD + 1;

%plot force:
subplot(2, 2, 2);
hold on;

axisLimForce = max(max(forceMag),-min(forceMag)) + 2;
axis([0 Duration -axisLimForce axisLimForce]);
xlabel('Time (min)');
forceMagP = plot(time, forceMag, 'Color', 'g', 'LineWidth', 1);
forceAngP = plot(time, forceAng, 'Color', 'k', 'LineWidth', 1);
legend([forceMagP forceAngP], {'Force Magnitude (pN)','Force Angle (radians)'});
timeMarkerForce = plot([0 0], [-axisLimForce axisLimForce],'Color','g','LineWidth',2);


title('Force Magnitude & Direction vs time','fontsize',16);

subplot(2,2,3);
polar(0, max(forceMag));
hold on;
forceAngZeros = zeros(1,2*length(forceAng));
forceMagZeros = zeros(1,2*length(forceMag));
forceAngZeros(2:2:end) = forceAng;
forceMagZeros(2:2:end) = forceMag;
tp = polar(forceAngZeros, forceMagZeros,':');
set(tp, 'linewidth', 0.2);
timeMarkerPolar = polar([0 0],[0 0],'g*');
set(timeMarkerPolar, 'linewidth', 3);
legend([tp timeMarkerPolar], {'Force History','Current Net Force'});

subplot(2,2,4);
hold on;
axisLimTorque = max(max(torqueData),-min(torqueData)) + 2;
axis([0 Duration -axisLimTorque axisLimTorque]);
xlabel('Time (min)');
ylabel('Torque (pN mum)');
torqueP = plot(time, torqueData, 'Color', 'r', 'LineWidth',1);
legend([torqueP], {'Torque (pN mum)'});
timeMarkerTorque = plot([0 0], [-axisLimTorque axisLimTorque],'Color','g','LineWidth',2);


step = 1;
for i = 1 : step : length(DATA)
    data = DATA(i,:);
    t = data(1);
    x = data(2);
    y = data(3);
    psi = data(4);
    
    cosinePrt = Prad*cos(psi);
    sinePrt   = Prad*sin(psi);
    xP_m = (x + cosinePrt);
    yP_m = (y + sinePrt);
    xP_d = (x - cosinePrt);
    yP_d = (y - sinePrt);
    
    MT_pos_M = data(start_M:end_M);
    MT_xm = (MT_pos_M(1:2:length(MT_pos_M)))';
    MT_ym = (MT_pos_M(2:2:length(MT_pos_M)))';
    MT_pos_D = data(start_D:end_D);
    MT_xd = (MT_pos_D(1:2:length(MT_pos_M)))';
    MT_yd = (MT_pos_D(2:2:length(MT_pos_M)))';
    
    MT_posXm = [xP_m.*ones(mt_numb,1) MT_xm]';
    MT_posYm = [yP_m.*ones(mt_numb,1) MT_ym]';
    MT_posXd = [xP_d.*ones(mt_numb,1) MT_xd]';
    MT_posYd = [yP_d.*ones(mt_numb,1) MT_yd]';
    
    force_M  = data((end_D+1):(end_D+2));
    force_D  = data((end_D+3):(end_D+4));
    force    = data((end_D+5):(end_D+6));
    
    torque_M = data(end_D + 7);
    torque_D = data(end_D + 8);
    torque   = data(end_D + 9);
    
    startBasePosM = end_D + 10;
    endBasePosM   = startBasePosM + 1;
    startBasePosD = endBasePosM + 1;
    endBasePosD   = startBasePosD + 1;
    basePosM      = data(startBasePosM:endBasePosM);
    basePosD      = data(startBasePosD:endBasePosD);
    
    MT_posXm = [basePosM(1).*ones(mt_numb,1) MT_xm]';
    MT_posYm = [basePosM(2).*ones(mt_numb,1) MT_ym]';
    MT_posXd = [basePosD(1).*ones(mt_numb,1) MT_xd]';
    MT_posYd = [basePosD(2).*ones(mt_numb,1) MT_yd]';

    
    subplot(2, 2, 1);
    cla;
    hold on;
    axis equal
    xlim([-R1_max R1_max])
    ylim([-R2_max R2_max]) 
    xlabel('x (nm)','fontsize',16)
    ylabel('y (nm)','fontsize',16) 
    title(['Time = ' num2str(t) ' of ' num2str(length(DATA)*Tau)],'fontsize',16)

    regionPlotter(0,0,R1_max,R2_max,regionAngles,regionProbabilities,regionForceMultipliers);
    ellipse(0,0,R1_max,R2_max,'k')
    circle_1(x,y,Prad,'g') 
    line([basePosM(1) xP_m], [basePosM(2) yP_m]);
    line([basePosD(1) xP_d], [basePosD(2) yP_d]);
    circle_1(basePosM(1),basePosM(2),Prad/8,'r')
    circle_1(basePosD(1),basePosD(2),Prad/8,'b')

    %Plot the MTs.
    line([xP_m xP_d], [yP_m yP_d]);
    line(MT_posXm,MT_posYm,'Color','r','LineWidth',1)
    line(MT_posXd,MT_posYd,'Color','b','LineWidth',1)
    
    hold off;
    % print
    % file_print=[figpath '/' num2str(nplot/1000) '.jpg'];
    % print(gcf,'-djpeg',file_print);
   
    subplot(2, 2, 2);
    delete(timeMarkerForce);
    timeMarkerForce = plot([t t], [-axisLimForce axisLimForce],'Color','g','LineWidth',2);

    mag      = norm(force);
    theta    = atan2(force(2),force(1));
    
    subplot(2,2,3);
    thetaPolar = [0 theta];
    rho        = [0 mag];
    delete(timeMarkerPolar);
    timeMarkerPolar = polar(thetaPolar, rho, '-g*');
    set(timeMarkerPolar, 'linewidth', 3);
    
    subplot(2,2,4);
    delete(timeMarkerTorque);
    timeMarkerTorque = plot([t t], [-Prad*10*F_MT Prad*10*F_MT],'Color','g','LineWidth',2);

    % print
    file_print=[figpath '/' num2str(i/length(DATA)) '.jpg'];
    print(gcf,'-djpeg',file_print);
end
