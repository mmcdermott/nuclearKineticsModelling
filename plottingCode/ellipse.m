function ellipse(x,y,r1,r2,option)
%x and y are the coordinates of the center of the circle
%r is the radius of the circle
%0.01 is the angle step, bigger values will draw the circle faster but
%you might notice imperfections (not very smooth)
ang=0:0.01:2*pi; 
xp=r1*cos(ang);
yp=r2*sin(ang);
plot(x+xp,y+yp,option);
hold on;
fill(x+xp,y+yp,[0.8,0.8,0.8])
%hold on;
end           
