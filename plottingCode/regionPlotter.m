function regionPlotter(x,y,r1,r2,regionAngles,regionProbabilities,regionForceMultipliers)
    %x and y are the coordinates of the center of the circle
    %r is the radius of the circle
    %0.01 is the angle step, bigger values will draw the circle faster but
    %you might notice imperfections (not very smooth)
    numRegions = length(regionAngles)-1;
    colorBand = 'b';
    colorCortex = 'k';
    widthBand = 4;
    widthCortex = 1;
    for (region = 1:numRegions)
        % First we need to convert the angle. We know that a point on the
        % ellipse is describable via 
        % $p = (r1 \cos(\theta'), r2 \sin(\theta'))$. But, we also have that
        % $p = (r \cos(regionAngles(region)), r \sin(regionAngles(region))$.
        % We don't know either $\theta'$ or $r$. But, we see that we can
        % solve for both, here, actually, as this describes the intersection
        % between two lines that intersect exactly once. In particular, we
        % see that $r = r1 \cos(\theta')/\cos(\theta)$ and 
        % $r = r2 \sin(\theta')/\sin(\theta)$. So, 
        % \tan(\theta') = r1/r2 \cotan(\theta). We'll use this to find
        % theta. 
        thetaStart = convertAngle(r1, r2, regionAngles(region));
        thetaEnd = convertAngle(r1, r2, regionAngles(region+1));
        %display(strcat('thetaStart = ', num2str(thetaStart)));
        %display(strcat('thetaEnd = ', num2str(thetaEnd)));
        if (regionProbabilities(region) == 1 && regionForceMultipliers(region) == 1)
            bandedEllipse(x,y,r1,r2,thetaStart,thetaEnd,colorCortex,widthCortex);
        else
            bandedEllipse(x,y,r1,r2,thetaStart,thetaEnd,colorBand,widthBand);
        end
    end
end           
