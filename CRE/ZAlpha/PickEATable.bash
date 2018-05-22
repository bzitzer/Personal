#!/bin/bash

# Returns the correct effective area table for vaStage6 
# $1 is the run in question
# $2 is the configuration directory, which has a file called EA.config which has the names of the correct EA files
# ATM21 and ATM22 have been switching on 3/15 and 11/15 every year
# old array files being used before summer 2009
# ANL L2 installed right before run number 58653
# upgraded array started Sept. 2012
# 20130731 - BJZ modified to read EA.config in the configuration directory


run=$1
CONFIGPATH=$2
EAPATH=/raid/reedbuck/bzitzer/tables/stage6

EANameArray=( $( cat "$CONFIGPATH/EA.config" ) ) 
DoCheck=0

if [ "$DoCheck" == "1" ]; then
    if ! test -e $CONFIGPATH/EA.config ; then
	echo EA.config does not exist! 
    fi	

    for i in "${EANameArray[@]}"
      do
      if ! test -e $EAPATH/$i ; then
	  echo "$i does not exist!"
      fi
    done
fi

OAATM21TABLE=$EAPATH/${EANameArray[0]}
OAATM22TABLE=$EAPATH/${EANameArray[1]}
NAATM21TABLE=$EAPATH/${EANameArray[2]}
NAATM22TABLE=$EAPATH/${EANameArray[3]}
UAATM21TABLE=$EAPATH/${EANameArray[4]}
UAATM22TABLE=$EAPATH/${EANameArray[5]}

if test $run -lt 37990 ; then
	table=$OAATM22TABLE
elif test $run -ge 37990 -a $run -lt 39958 ; then
	table=$OAATM21TABLE
elif test $run -ge 39958 -a $run -lt 42971 ; then
	table=$OAATM22TABLE
elif test $run -ge 42971 -a $run -lt 44881 ; then
	table=$OAATM21TABLE
elif test $run -ge 44881 -a $run -lt 46542 ; then
	table=$OAATM22TABLE
elif test $run -ge 46542 -a $run -lt 48283 ; then
	table=$NAATM22TABLE
elif test $run -ge 48283 -a $run -lt 50447 ; then
	table=$NAATM21TABLE
elif test $run -ge 50447 -a $run -lt 53313 ; then
	table=$NAATM22TABLE
elif test $run -ge 53313 -a $run -lt 55608 ; then
	table=$NAATM21TABLE
elif test $run -ge 55608 -a $run -lt 58653 ; then
	table=$NAATM22TABLE
elif test $run -ge 58653 -a $run -lt 61178 ; then
	table=$NAATM21TABLE
elif test $run -ge 61178 -a $run -lt 63408 ; then
	table=$NAATM22TABLE
elif test $run -ge 63408 -a $run -lt 64799 ; then
	table=$UAATM22TABLE
elif test $run -ge 64799 -a $run -lt 67411 ; then
	table=$UAATM21TABLE
elif test $run -ge 67410; then
	table=$UAATM22TABLE	
fi

echo $table 
