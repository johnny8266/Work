#baseDirCoreSoft=/media/nschmidt/local/coresoftware
#baseDirDetector=/media/nschmidt/local/fun4all_eicdetectors
#baseDirCoreSoft=/home/fbock/eic/developEIC/coresoftware
#baseDirDetector=/home/fbock/eic/developEIC/eicdetectors

baseDirCoreSoft=/vol0/pwang-l/Singularity/repository/coresoftware
baseDirDetector=/vol0/pwang-l/Singularity/repository/fun4all_eicdetectors
export MYINSTALL=/vol0/pwang-l/Singularity/Myinstall/

# full build
source /cvmfs/eic.opensciencegrid.org/gcc-8.3/opt/fun4all/core/bin/eic_setup.sh -n
export MYINSTALL=/vol0/pwang-l/Singularity/Myinstall/
#export MYINSTALL=$HOME/install

cd $MYINSTALL
rm -rf *

# for FLDR2 in g4main g4detectors
# do
#     echo -e "\n\n\n"
#     echo "entering " $FLDR2 " building"
#     echo -e "\n\n\n"
# 
#     cd $baseDirCoreSoft/simulation/g4simulation/$FLDR2
#     mkdir build
#     cd build
#     ../autogen.sh --prefix=$MYINSTALL
#     make -j4 install
# done
# 
# for FLDR in CaloBase trackbase trackbase_historic mvtx intt tpc micromegas
# do
#     echo -e "\n\n\n"
#     echo "entering " $FLDR " building"
#     echo -e "\n\n\n"
# 
#     cd $baseDirCoreSoft/offline/packages/$FLDR
#     mkdir build
#     cd build
#     ../autogen.sh --prefix=$MYINSTALL
#     make -j4 install
# done
# 
# for FLDR2 in g4tpc g4mvtx g4intt g4micromegas g4epd g4bbc g4calo g4vertex g4jets 
# do
#     echo -e "\n\n\n"
#     echo "entering " $FLDR2 " building"
#     echo -e "\n\n\n"
# 
#     cd $baseDirCoreSoft/simulation/g4simulation/$FLDR2
#     mkdir build
#     cd build
#     ../autogen.sh --prefix=$MYINSTALL
#     make -j4 install
# done



for PKW in CaloBase CaloReco ClusterIso
do
    echo -e "\n\n\n"
    echo "entering " $PKW " building"
    echo -e "\n\n\n"

    cd $baseDirCoreSoft/offline/packages/$PKW
    rm -rf build
    mkdir build
    cd build
    ../autogen.sh --prefix=$MYINSTALL
    make -j4 install
done


#for PKW1 in g4calo g4eval g4detectors
for PKW1 in g4calo g4detectors
do
    echo -e "\n\n\n"
    echo "entering " $PKW1 " building"
    echo -e "\n\n\n"

    cd $baseDirCoreSoft/simulation/g4simulation/$PKW1
    rm -rf build
    mkdir build
    cd build
    ../autogen.sh --prefix=$MYINSTALL
    make -j4 install
done


#for PKW2 in g4eiccalos g4drcalo g4beastmagnet g4lblvtx g4rich g4mrich
for PKW2 in g4eiccalos g4drcalo
do
    echo -e "\n\n\n"
    echo "entering " $PKW2 " building"
    echo -e "\n\n\n"
    
    cd $baseDirDetector/simulation/g4simulation/$PKW2
    rm -rf build
    mkdir build
    cd build
    ../autogen.sh --prefix=$MYINSTALL
    make -j4 install
done


cd ~/Singularity/
source /cvmfs/eic.opensciencegrid.org/gcc-8.3/opt/fun4all/core/bin/setup_local.sh $MYINSTALL


# 
# source /cvmfs/eic.opensciencegrid.org/gcc-8.3/opt/fun4all/core/bin/setup_local.sh $MYINSTALL
# 
# 
# for FLDR in KFParticle_sPHENIX trackreco PHTpcTracker tpccalib CaloReco ClusterIso particleflow NodeDump
# do
#     echo -e "\n\n\n"
#     echo "entering " $FLDR " building"
#     echo -e "\n\n\n"
#     cd $baseDirCoreSoft/offline/packages/$FLDR
#     mkdir build
#     cd build
#     ../autogen.sh --prefix=$MYINSTALL
#     make -j4 install
# done


#cd $baseDirCoreSoft

#'
#build order:
#simulation/g4simulation/g4main
#simulation/g4simulation/g4detectors
#offline/packages/CaloBase
#offline/packages/trackbase
#offline/packages/trackbase_historic
#offline/packages/mvtx
#offline/packages/intt
#offline/packages/tpc
#offline/packages/micromegas
#simulation/g4simulation/g4tpc
#simulation/g4simulation/g4mvtx
#simulation/g4simulation/g4intt
#simulation/g4simulation/g4micromegas
#simulation/g4simulation/g4epd
#simulation/g4simulation/g4bbc
#simulation/g4simulation/g4calo
#offline/packages/trigger
#offline/packages/PHGenFitPkg/GenFitExp
#offline/packages/PHGenFitPkg/PHGenFit
#simulation/g4simulation/g4vertex
#simulation/g4simulation/g4jets
#offline/packages/jetbackground
#simulation/g4simulation/g4eval
#simulation/g4simulation/g4trackfastsim
#simulation/g4simulation/g4histos
#offline/packages/KFParticle_sPHENIX
#offline/packages/trackreco
#offline/packages/PHTpcTracker
#offline/packages/tpccalib
#offline/packages/CaloReco
#offline/packages/ClusterIso
#offline/packages/particleflow
#offline/packages/NodeDump
#'
