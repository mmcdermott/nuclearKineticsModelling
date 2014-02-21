run parameters

nplot   = 1;
figpath = ['figs'];

run clean_fig_data_folder

figure('visible','off');
hold on;

axis equal
xlim([-R1_max R1_max])
ylim([-R2_max R2_max]) 
xlabel('x (nm)','fontsize',16)
ylabel('y (nm)','fontsize',16) 
% plot boundary and interior points

regionPlotter(0,0,R1_max,R2_max,regionAngles,regionProbabilities)
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

startBasePosM = end_D + 10;
endBasePosM   = startBasePosM + 1;
startBasePosD = endBasePosM + 1;
endBasePosD   = startBasePosD + 1;

step = 100;
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

    cla;
    hold on;
    axis equal
    xlim([-R1_max R1_max])
    ylim([-R2_max R2_max]) 
    xlabel('x (nm)','fontsize',16)
    ylabel('y (nm)','fontsize',16) 
    title(['Time = ' num2str(t) ' of ' num2str(length(DATA)*Tau)],'fontsize',16)

    regionPlotter(0,0,R1_max,R2_max,regionAngles,regionProbabilities);
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
   
    % print
    file_print=[figpath '/' num2str(i/length(DATA)) '.jpg'];
    print(gcf,'-djpeg',file_print);
end
