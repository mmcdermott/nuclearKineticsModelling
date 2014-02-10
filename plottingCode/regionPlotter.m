function regionPlotter(x,y,r1,r2,regionAngles,regionProbabilities)
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
        if (regionProbabilities(region) == 1)
            bandedEllipse(x,y,r1,r2,regionAngles(region),regionAngles(region+1),colorCortex,widthCortex);
        else
            bandedEllipse(x,y,r1,r2,regionAngles(region),regionAngles(region+1),colorBand,widthBand);
        end
    end
end           
