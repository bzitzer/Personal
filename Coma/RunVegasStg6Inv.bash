#!/bin/bash

ConfigPath=/homes/gretel/bzitzer/vegas/configfiles/$1

Opt=" -config=$ConfigPath/Stg6CBG.config -cuts=$ConfigPath/Stg6cuts.config -S6A_Batch=1 "
Opt=" $Opt -S6A_UserDefinedExclusionList=/homes/gretel/bzitzer/veritas/ZitzerExclusionList_wHAWC.txt"
Opt=" $Opt -S6A_ExcludeSource=1"
Opt=" $Opt -RBM_UseModifiedLiMa=1"
echo $Opt
Stg6List=$2

Vegas=~/vegas/vegas_liMaDev/bin

rm -rf Stg6CmdList.txt

cat $Stg6List |
    while read File
    do
#	echo $File
	if [ -f "$File" ]; then
	    echo $File
	    run=${File%'.stg5.root'}
	    run=${File%'.root'}
	    last=${#run}
	    let first=last-5
	    run=${run:first:last}
	    echo $run
	    stg5File=${File%'.root'}'.stg5.root'
	    echo $stg5File > Stage6List$run.txt
	    Cmd="$Vegas/vaStage6 $Opt Stage6List$run.txt"
	    echo $Cmd >> Stg6CmdList.txt
	    echo "mv config/results_s6.root config/results_"$run"_modLiMa_s6.root " >> Stg6CmdList.txt
	fi

    done
