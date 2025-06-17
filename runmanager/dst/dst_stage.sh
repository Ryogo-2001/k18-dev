#!/bin/bash

#_____________________________
#    stager for DST
#  AUTHOR: Yoshiyuki NAKADA
#      April 2, 2019
#_____________________________


APP_PATH=$(readlink -f $0)
APP_DIR=$(dirname $APP_PATH)

HSTAGE_DIR=/ghi/fs02/hstage/requests

#__________________________________________________

QUIET=0
while getopts ":q" opt; do
    case $opt in
	q ) QUIET=1 ;;
	* ) echo "ERROR: Invalid option"; exit 1 ;;
    esac
done

shift $(($OPTIND - 1))

#__________________________________________________

if [ $# -ne 1 ]; then
    echo "USAGE: $0 [runlist]"
    exit 0
fi

if ! [ -e "$1" ]; then
    echo "ERROR: No such file > $1"
    exit 1
fi

runlist=$1

tmp=$(basename $runlist)
label=${tmp%.*}

#__________________________________________________

work=""
bin=""
rootin=""
runid=()

#__________________________________________________

while read line; do

    if [ -z "$line" -o "${line:0:1}" = '#' ]; then
        continue
    elif [ "${line:0:5}" = 'work:' ]; then
        tmp=( $line )
        work=${tmp[1]}
        continue
    elif [ "${line:0:4}" = 'bin:' ]; then
        tmp=( $line )
        bin=${tmp[1]}
        continue
    elif [ "${line:0:7}" = 'rootin:' ]; then
        tmp=( $line )
        rootin=${tmp[1]}
        continue
    fi

    words=( $line )
    for item in ${words[@]}; do

	if [ "${item:0:1}" == '#' ]; then
	    break
	fi

	expr $item + 1 >/dev/null 2>&1
	if [ $? -lt 2 ]; then
	    tmp=$(( 10#$item ))
	    runid+=( $(printf %05d $tmp) )
	fi
    done
    
done < $runlist

#__________________________________________________

if [ -z "$work" ] || ! [ -d "$work" ]; then
    echo "ERROR: Invalid file declaration [work: $work]"
    exit 1
fi

#__________________________________________________

cd $work
if [ $? -eq 1 ]; then
    exit 1
fi

if ! [ -x "$bin" ]; then
    if [ -z "$bin" ]; then bin="null"; fi
    echo "ERROR: Invalid file declaration [bin: $bin]"
    exit 1
fi

if ! [ -d "$rootin" ]; then
    if [ -z "$rootin" ]; then rootin="null"; fi
    echo "ERROR: Invalid file declaration [rootin: $rootin]"
    exit 1
fi

#__________________________________________________

element=()
basebin=$(basename $bin)
if [ "$basebin" = "DstPiKAna" ]; then
    element=( "KuramaTracking" "K18Tracking" "Hodoscope" "Easiroc" )
elif [ "$basebin" = "DstPiKCatch" ]; then
    element=( "CatchTracking" "KuramaTracking" "K18Tracking" "Hodoscope" "Easiroc" )
elif [ "$basebin" = "DstKuramaHodoscope" ]; then
    element=( "KuramaTracking" "Hodoscope" )
else
    echo "ERROR: Unsupported DST [$basebin]"
    exit 1
fi

#__________________________________________________

cd $work/$rootin
if [ $? -eq 1 ]; then
    exit 1
fi

fl_exit=0
for id in ${runid[@]}; do
    for item in ${element[@]}; do
	target="run${id}_${item}.root"
	if ! [ -r "$target" ]; then
	    echo "ERROR: Unreachable ROOT input file [$work/$rootin/$target]"
	    fl_exit=1
	fi
    done
done

if [ $fl_exit -eq 1 ]; then
    exit 1
fi

#__________________________________________________

index=0
inpath=()
stagepath=()
for id in ${runid[@]}; do
    for item in ${element[@]}; do
        tmppath=( $(readlink -f $work/$rootin/run${id}_${item}.root) )
        while read status path; do
            case "$status" in
            'G' );;
            'B' );;
            'BP' )
                echo "ERROR: Cannot be staged [$item]"
                ;;
            'H' )
                # od $tmppath >/dev/null 2>&1 &
                stagepath+=( $tmppath )
                ;;
            '*' )  
                echo "ERROR: Unknow status was detected [$item]"
                ;;
            esac
        done < <(ghils $tmppath)
        inpath+=( $tmppath )
    done
done

HSTAGE_PATH=$HSTAGE_DIR/$label.lst.$(date "+%Y%m%d%H%M%S")
for item in ${stagepath[@]}; do
    echo $item >> $HSTAGE_PATH
done

#__________________________________________________

cd $work
if [ $? -eq 1 ]; then
    exit 1
fi

fl_done=0

stat=()

while : ; do
    index=0
    for item in ${inpath[@]}; do
        while read status path; do
            case "$status" in
            'G' )
                stat[$index]=1
                fl_done=$(($fl_done + 1))
                ;;
            'B' )
                stat[$index]=1
                fl_done=$(($fl_done + 1))
                ;;
            'BP' )
                stat[$index]=2
                fl_done=$(($fl_done + 1))
                ;;
            'H' )  stat[$index]=0 ;;
            '*' )  
                echo "ERROR: Unknow status was detected [$item]"
                stat[$index]=-1 ;;
            esac
        done < <(ghils $item)
        index=$(($index + 1))
    done
    if [ $QUIET -eq 0 ]; then
        prog=""
        for item in ${stat[@]}; do
            case $item in
            0 ) prog="${prog}." ;;
            1 ) prog="${prog}!" ;;
            2 ) prog="${prog}x" ;;
            * ) prog="${prog}?" ;;
            esac
        done
        echo -en "executing... [$prog] $fl_done/${#inpath[@]}\r"
    fi
    if [ ${#inpath[@]} -eq $fl_done ]; then break; fi
done

if [ $QUIET -eq 0 ]; then echo -en "\n"; fi
