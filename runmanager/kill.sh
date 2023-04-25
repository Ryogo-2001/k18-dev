#!/bin/sh

bjobs | while read jid buff; do

    if [ $jid == 'JOBID' ]; then
	continue
    fi

    bkill $jid

done

killall hadd
killall run.py
