#!/bin/bash

BCL=/home/iwe20/bcl/bin/bcl-apps.exe
ROSETTA_SOURCE=~/Rosetta/main/source
ROSETTA_APP=rosetta_scripts.default.linuxgccrelease

runs=25
align=true
overwrite='false'

#-----------------------------------------------------------------#

while getopts l:r:a:o: flag
do
    case "${flag}" in
        l) ligname=${OPTARG};;
        r) runs=${OPTARG};;
        a) align=${OPTARG};;
        o) overwrite=${OPTARG};;
    esac
done

#-----------------------------------------------------------------#

echo Create file structure

if [ -d "$ligname" ]; then
    echo "Direcotry $ligname does exist."
    if [ $overwrite = true ] ; then
        echo "Overwrite is allowed. Old files will be deleted..."
    else
        echo "Aborting program. Use option '-o true' to allow overwrite or rename your ligand."
        exit 1
    fi
fi

rm -rf $ligname
mkdir $ligname
if test -f $ligname.sdf; then 
	cp $ligname.sdf $ligname/$ligname.mol
fi
if test -f $ligname.mol; then 
	cp $ligname.mol $ligname/
fi

cd $ligname/


#-----------------------------------------------------------------#

if [ $align = true ] ; then
	echo Align ligand to start fragment
	
	$BCL molecule:AlignToScaffold ../input/lig_model.sdf $ligname.mol $ligname.final.mol &> error.log
else
	cp $ligname.mol $ligname.final.mol
fi

if [ ! -f "$ligname.final.mol" ] ; then
    echo "[ERROR] Alignment failed. Check error.log for details."
    exit 1
fi

#-----------------------------------------------------------------#

echo Start conformer generation...

$BCL molecule:ConformerGenerator -top_models 50 -max_iterations 500 -add_h -ensemble_filenames "$ligname".final.mol -conformers_single_file $ligname.conf.sdf -conformation_comparer SymmetryRMSD 0.25 -cluster &> error.log

if [ ! -s "$ligname.conf.sdf" ] ; then
    echo "[ERROR] Conformer generation failed. Check error.log for details."
    exit 1
fi

#-----------------------------------------------------------------#

echo Create ligand params file

$ROSETTA_SOURCE/scripts/python/public/molfile_to_params.py -n LIG -p LIG --conformers-in-one-file $ligname.conf.sdf &> error.log
# this line is suited for my selfmade rosetta standalone
#$ROSETTA_SOURCE/molfile_to_params.py -n LIG -p LIG --conformers-in-one-file $ligname.conf.sdf &> error.log

if [ ! -f "LIG.pdb" ] ; then
    echo "[ERROR] Param file generation failed. Check error.log for details."
    exit 1
fi

#-----------------------------------------------------------------#

echo Create docking pose

cat ../input/protein_target.pdb LIG.pdb > protein_ligand_complex.pdb

#-----------------------------------------------------------------#

echo Start docking. This might take some time...

$ROSETTA_SOURCE/bin/$ROSETTA_APP @../input/options -nstruct $runs &> error.log

if [ ! -f "score.sc" ] ; then
    echo "[ERROR] Docking failed. Check error.log for details."
    exit 1
fi

echo Done.

#-----------------------------------------------------------------#

echo Prepare result files

sort -n -k 2 score.sc | head -n 10 | awk '{print $2 " " $50 " " $57}' > best_complexes.sc
sort -n -k 50 score.sc | head -n 10 | awk '{print $2 " " $50 " " $57}' > best_bindings.sc

echo
echo ------------------------------------------------
echo Best complex is:
sort -n -k 2 score.sc | head -n 1 | awk '{print $57}'
sort -n -k 2 score.sc | head -n 1 | awk '{print $2}'
echo
echo Best binding is:
sort -n -k 50 score.sc | head -n 1 | awk '{print $57}'
sort -n -k 50 score.sc | head -n 1 | awk '{print $50}'
echo ------------------------------------------------
echo

echo Scores are expected to be negative. The more negative, the better.
echo Find more results in best_complexes.sc and best_bindings.sc in $ligname/
echo 
echo All done.

#-----------------------------------------------------------------#

rm -f error.log
rm -f $ligname.conf.sdf
rm -f $ligname.final.mol
rm -f $ligname.mol
rm -f LIG.params
rm -f LIG.pdb
rm -f LIG_conformers.pdb
rm -f protein_ligand_complex.pdb
rm -f score.sc
