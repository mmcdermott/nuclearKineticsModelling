figflg = exist(figpath);
if figflg == 7
    'deleting old figures'
    cd(figpath)
    delete *.jpg
    cd ..
else
   mkdir(figpath);
end