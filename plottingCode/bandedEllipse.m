function bandedEllipse(x,y,r1,r2,angStart,angEnd,color, width)
    %x and y are the coordinates of the center of the circle
    %r is the radius of the circle
    %0.01 is the angle step, bigger values will draw the circle faster but
    %you might notice imperfections (not very smooth)
    ang=angStart:0.01:angEnd; 
    xp=r1*cos(ang);
    yp=r2*sin(ang);
    plot(x+xp,y+yp,'Color', color, 'LineWidth', width);
    hold on;
end           
