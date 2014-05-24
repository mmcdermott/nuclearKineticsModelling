if (exist('translation', 'var') == 0)
    translation = true;
end
if (exist('motherEnvelope', 'var') == 0)
    motherEnvelope = pi;
end
if (exist('numRuns','var') == 0)
    numRuns = 100;
end
if (exist('fixedWidth','var') == 0)
    fixedWidth = pi/6;
end
dataBaseDir = '../data/rotVsBand/';
if (motherEnvelope == pi)
    dataEnvDir = 'mEnvP/';
elseif (motherEnvelope == 2*pi)
    dataEnvDir = 'mEnv2P/';
end
if (translation)
    dataTransDir = 'translation/';
else
    dataTransDir = 'noTranslation/';
end
dataDir = [dataBaseDir dataTransDir dataEnvDir];
dataFiles = getAllFiles(dataDir);
plotXval = []; plotXratio = []; plotYpsi = []; plotYx = []; plotYy = [];
errorY = []; errorX = []; errorPsi = [];
for i = 1:length(dataFiles)
    dataFile = char(dataFiles(i));
    fileParts = strsplit(dataFile,'/');
    dataFileFull = char(fileParts(6));
    fileNameParts = strsplit(dataFileFull,'.csv');
    %resetting dataFile
    dataFile = char(fileNameParts(1));
    fileNameStripped = strsplit(dataFile,'bw');
    fileNameStripped = strsplit(char(fileNameStripped(2)),'ew');
    bandWidth = char(fileNameStripped(1));
    envWidth = char(fileNameStripped(2));
    
    if (strcmp(envWidth, 'P'))
        envWidth = pi;
    else 
        envWidthParts = strsplit(envWidth,'o');
        envNum = char(envWidthParts(1));
        if (strcmp(envNum, 'P'))
            envNum = pi;
        else
            envNum = 2*pi;
        end
        envDen = str2double(char(envWidthParts(2)));
        envWidth = envNum/envDen;
    end
    if (strcmp(bandWidth, '0') || strcmp(bandWidth, 'P'))
        if (strcmp(bandWidth, '0'))
            bandWidth = 0;
        else 
            bandWidth = pi;
        end
    else 
        bandWidthParts = strsplit(bandWidth,'o');
        bandNum = char(bandWidthParts(1));
        if (strcmp(bandNum, 'P'))
            bandNum = pi;
        else
            bandNum = 2*pi;
        end
        bandDen = str2double(char(bandWidthParts(2)));
        bandWidth = bandNum/bandDen;
    end
    if (envWidth == fixedWidth)
        display(['grabbing Data for file ' dataFileFull]);
        run parameters
        csvrange = [rowStart colStart rowStart+numRuns-1 colStart+2]; 
        DATA = csvread(dataFileFull,rowStart,colStart,csvrange);
        x = DATA(:,1);
        y = DATA(:,2);
        psi = abs((DATA(:,3)-pi/2)*180/pi);
        plotXval = [plotXval bandWidth];
        plotXratio = [plotXratio bandWidth/envWidth];
        plotYpsi = [plotYpsi mean(psi)];
        plotYx = [plotYx mean(x)];
        plotYy = [plotYy mean(y)];
        
        errorPsi = [errorPsi std(psi)];
        errorX = [errorX std(x)];
        errorY = [errorY std(y)];
    end
end
[sortedXval,indices] = sort(plotXval*180/pi);
sortedXratio = plotXratio(indices);
sortedYpsi = plotYpsi(indices);
sortedYx = plotYx(indices);
sortedYy = plotYy(indices);

sortedErrorPsi = errorPsi(indices);
sortedErrorX = errorX(indices);
sortedErrorY = errorY(indices);

figure;
if (translation)
    subplot(3,2,1);
    %text(0,0,'All angles are in degrees');
    clear off;
    errorbar(sortedXval,sortedYpsi,sortedErrorPsi);
    xlabel('Band Width');
    ylabel('Abs(Rotation)');
    subplot(3,2,3);
    errorbar(sortedXval,sortedYx,sortedErrorX);
    xlabel('Band Width (degrees)');
    ylabel('Final x positions (mum)');
    subplot(3,2,5);
    errorbar(sortedXval,sortedYy,sortedErrorY);
    xlabel('Band Width (degrees)');
    ylabel('Final x positions (mum)');

    subplot(3,2,2);
    errorbar(sortedXratio,sortedYpsi,sortedErrorPsi);
    xlabel('(Band Width)/(Daughter Envelope Width)');
    ylabel('Abs(Rotation)');
    subplot(3,2,4);
    errorbar(sortedXratio,sortedYx,sortedErrorX);
    xlabel('(Band Width)/(Daughter Envelope Width)');
    ylabel('Final x positions (mum)');
    subplot(3,2,6);
    errorbar(sortedXratio,sortedYy,sortedErrorY);
    xlabel('(Band Width)/(Daughter Envelope Width)');
    ylabel('Final x positions (mum)');
else
    subplot(1,2,1);
    errorbar(sortedXval,sortedYpsi,sortedErrorPsi);
    xlabel('Band Width');
    ylabel('Abs(Rotation)');
    
    subplot(1,2,2);
    errorbar(sortedXratio,sortedYpsi,sortedErrorPsi);
    xlabel('(Band Width)/(Daughter Envelope Width)');
    ylabel('Abs(Rotation)');
end
