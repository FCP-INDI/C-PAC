# CPAC Image Resources installation script
# script usage: ./install_resources.sh $FSLDIR
#    with $FSLDIR being your environment variable pointing
#    to your FSL directory

cp -R MNI_3mm $1/data/standard;
cp -R symmetric $1/data/standard;
cp -R tissuepriors/2mm $1/data/standard/tissuepriors;
cp -R tissuepriors/3mm $1/data/standard/tissuepriors;
cp HarvardOxford-lateral-ventricles-thr25-2mm.nii.gz $1/data/atlases/HarvardOxford;
