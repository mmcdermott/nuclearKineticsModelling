nplot   = 1;
figpath = ['figs'];

run clean_fig_data_folder

figure
subplot(2, 2, 1);
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
polar(0, 10*F_MT);

subplot(2,2,4);
hold on;
axisLimTorque = max(max(torqueData),-min(torqueData)) + 2;
axis([0 Duration -axisLimTorque axisLimTorque]);
xlabel('Time (min)');
ylabel('Torque (pN mum)');
torqueP = plot(time, torqueData, 'Color', 'r', 'LineWidth',1);
legend([torqueP], {'Torque (pN mum)'});
timeMarkerTorque = plot([0 0], [-axisLimTorque axisLimTorque],'Color','g','LineWidth',2);