clear
fileBase="rotVsBand/springed/"
commandBase="./mtKineticModel 100"

mEnvPrefix="const float_T envWidthM = "
translationPrefix="const bool translation = "
envPrefix="const float_T envWidthD = "
bandPrefix="const float_T width = "

mEnvs="2*pi"
translations="true"
envWidths="pi\/2 pi\/2.5 pi\/3 pi\/3.5 pi\/4 pi\/5 pi\/6 pi\/7 pi\/8 pi\/10 pi\/12 pi\/15 pi\/18"
bandWidths="pi\/2 pi\/2.5 pi\/3 pi\/3.5 pi\/4 pi\/5 pi\/6 pi\/7 pi\/8 pi\/10 pi\/12 pi\/15 pi\/18"

mEnvDirs="mEnv2P/"
mEnvDirs=($mEnvDirs)
translationDirs="translation/"
translationDirs=($translationDirs)
envWidthFiles="ewPo2 ewPo2.5 ewPo3 ewPo3.5 ewPo4 ewPo5 ewPo6 ewPo7 ewPo8 ewPo10 ewPo12 ewPo15 ewPo18"
envWidthFiles=($envWidthFiles)
bandWidthFiles="bwPo2 bwPo2.5 bwPo3 bwPo3.5 bwPo4 bwPo5 bwPo6 bwPo7 bwPo8 bwPo10 bwPo12 bwPo15 bwPo18"
bandWidthFiles=($bandWidthFiles)

prevMEnv="2*pi"
prevTranslation="true"
prevEnv="pi\/6"
prevBand="pi\/4"

mEnvCount=0
declare -i mEnvCount
transCount=0
declare -i transCount
bandCount=0
declare -i bandCount
envCount=0
declare -i envCount

for mEnv in $mEnvs; do
  mEnvDir=${mEnvDirs[$mEnvCount]}
  sed -i "s/$mEnvPrefix$prevMEnv;/$mEnvPrefix$mEnv;/g" parameters.cpp
  transCount=0
  declare -i transCount
  for translation in $translations; do
    translationDir=${translationDirs[$transCount]}
    sed -i "s/$translationPrefix$prevTranslation;/$translationPrefix$translation;/g" parameters.cpp
    envCount=0
    declare -i envCount
    for envWidth in $envWidths; do
      envWidthFile=${envWidthFiles[$envCount]}
      sed -i "s/$envPrefix$prevEnv;/$envPrefix$envWidth;/g" parameters.cpp
      bandCount=0
      declare -i bandCount
      for bandWidth in $bandWidths; do
        bandWidthFile=${bandWidthFiles[$bandCount]}
        sed -i "s/$bandPrefix$prevBand;/$bandPrefix$bandWidth;/g" parameters.cpp
        make && $commandBase $fileBase$translationDir$mEnvDir$bandWidthFile$envWidthFile
        prevBand=$bandWidth
        bandCount=`expr $bandCount + 1`
      done
      prevEnv=$envWidth
      envCount=`expr $envCount + 1`
    done
    prevTranslation=$translation
    transCount=`expr $transCount + 1`
  done
  prevMEnv=$mEnv
  mEnvCount=`expr $mEnvCount + 1`
done
