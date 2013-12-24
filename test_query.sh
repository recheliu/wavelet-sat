#!/bin/bash
ROOT_PATH=build/vs2010/x64
NC_PATH=D:/projects/WaveletSAT/build/vs2010/x64/output
if [ "$1" == "F" ] ;
then
    ROOT_PATH="${ROOT_PATH}/SATFile"
    MS="8192 0"
    ZS="1 0"
    EXT="nc4"
    WS="1"
else
    MS="8192"
    ZS="1 0"
    EXT="wnc4"
    if [ "$1" == "W" ] ;
    then
	WS="1 4 16"
    else
	WS="1"
    fi
fi
BINS="32 64 96 128"
DECODER="${ROOT_PATH}/SimpleNDQuery/Release/SimpleNDQuery.exe"
for w in ${WS};
do
    for z in ${ZS};
    do
	for m in ${MS};
	do
	    for d in hydrogenAtom Foot Ocean.20010201 MJO.2008-01-01_00-00-00.QVAPOR;
	    do
		out_of_core="false"
		if [ "$1" == "F" -a "$m" == "0" ] ;
		then
		    out_of_core=true;
		fi

		if [ $out_of_core == true ]
		then
		    nq="128"
		else
		    nq="8192"
		fi
	    
		ARGLIST="--size-of-full-arrays $m --is-testing-query true --n-testing-values ${nq} --query-win-length $w"

		echo  ====================================================
		echo $d: ${ARGLIST}
		for b in ${BINS};
		do
		    TEST_ARGLIST="--nc-filepath ${NC_PATH}/${d}.b_${b}.z_${z}.${EXT} ${ARGLIST}"
		    echo  ${TEST_ARGLIST}
		    do_test=true
		    if [ "$d" == "Foot" -a "$1" == "F" -a "$m" == 8192 ] ;
		    then
			if [ "$b" == "96" -o "$b" == "128" ] ; 
			then
			    do_test=false
			fi
		    fi
		    if [ "$d" == "Ocean.20010201" -a "$1" == "F" -a "$m" == 8192 ] ;
		    then
			if [ "$b" == "128" ] ; 
			then
			    do_test=false
			fi
		    fi
		    if [ $do_test == true ] ;
		    then
			${DECODER} ${TEST_ARGLIST}
		    else
			echo "main:-1,-1,-1,-1,0"
		    fi
		done
	    done
	done
    done
done

############################################################
# Copyright (c) 2013 Teng-Yok Lee
#
# See the file LICENSE.txt for copying permission.
############################################################
