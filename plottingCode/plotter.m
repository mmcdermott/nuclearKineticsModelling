if (exist('translation', 'var') == 0)
    translation = true;
end
if (exist('motherEnvelope', 'var') == 0)
    motherEnvelope = pi;
end
if (exist('numRuns','var') == 0)
    numRuns = 100;
end
if (exist('ratio','var') == 0)
    ratio = true;
end
dataBaseDir = '../data/rotVsBand/springed/';
figsBaseDir = '../../figures/rotVsBand/springed/';
if (motherEnvelope == pi)
    envDir = 'mEnvP/';
elseif (motherEnvelope == 2*pi)
    envDir = 'mEnv2P/';
end
if (translation)
    transDir = 'translation/';
else
    transDir = 'noTranslation/';
end
dataDir = [dataBaseDir transDir envDir];
dataFiles = getAllFiles(dataDir);

fixedWidths = [pi/2 pi/2.5 pi/3 pi/3.5 pi/4 pi/5 pi/6 pi/7 pi/8 pi/10 ...
               pi/12 pi/15 pi/18];
fixedWidthFiles = {'Po2', 'Po2p5', 'Po3', 'Po3p5', 'Po4', 'Po5', ...
                   'Po6', 'Po7', 'Po8', 'Po10', 'Po12', 'Po15', 'Po18'};
for j = 1:length(fixedWidths)
    %%%TODO: This is done stupidly.
    fixedWidth = fixedWidths(j);
    fixedWidthFile = ['ew' char(fixedWidthFiles(j))];
    plotXval = []; plotXratio = []; plotYpsi = []; plotYx = []; 
    plotYy = []; errorY = []; errorX = []; errorPsi = [];
    for i = 1:length(dataFiles)
        dataFile = char(dataFiles(i));
        fileParts = strsplit(dataFile,'/');
        dataFileFull = char(fileParts(7));
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
            if (bandWidth == pi/2)
                continue;
            end
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

    fig = figure('visible', 'off');
    if (translation)
        if (ratio)
            subplot(3,1,1);
            hold on;
            errorbar(sortedXratio,sortedYpsi,sortedErrorPsi);
            plot([sortedXratio(1) sortedXratio(length(sortedXratio))], ...
                 [90 90], 'g--');
            plot([sortedXratio(1) sortedXratio(length(sortedXratio))], ...
                 [180/8 180/8], 'k--');
            xlabel('(Band Width)/(Daughter Envelope Width)');
            ylabel('Abs(Rotation)');
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            subplot(3,1,2);
            hold on;
            errorbar(sortedXratio,sortedYx,sortedErrorX);
            xlabel('(Band Width)/(Daughter Envelope Width)');
            ylabel('Final x positions (mum)');
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            subplot(3,1,3);
            hold on;
            errorbar(sortedXratio,sortedYy,sortedErrorY);
            xlabel('(Band Width)/(Daughter Envelope Width)');
            ylabel('Final x positions (mum)');
        else
            subplot(3,1,1);
            hold on;
            errorbar(sortedXval,sortedYpsi,sortedErrorPsi);
            xlabel('Band Width');
            ylabel('Abs(Rotation)');
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            subplot(3,1,2);
            hold on;
            errorbar(sortedXval,sortedYx,sortedErrorX);
            xlabel('Band Width (degrees)');
            ylabel('Final x positions (mum)');
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            subplot(3,1,3);
            hold on;
            errorbar(sortedXval,sortedYy,sortedErrorY);
            xlabel('Band Width (degrees)');
            ylabel('Final x positions (mum)');
        end
    else
        if (ratio)
            hold on;
            errorbar(sortedXratio,sortedYpsi,sortedErrorPsi);
            xlabel('(Band Width)/(Daughter Envelope Width)');
            ylabel('Abs(Rotation)');
        else
            hold on;
            errorbar(sortedXval,sortedYpsi,sortedErrorPsi);
            xlabel('Band Width');
            ylabel('Abs(Rotation)');
        end
    end
    figureFile = [figsBaseDir transDir envDir fixedWidthFile];
    print(fig,'-dpdf', figureFile);
    %savefig(figureFile);
    saveas(fig,[figureFile '.fig']);
end
