#!/bin/bash
function run_all_file() {
filename=$1
thread=$2

lines=`cat ${filename} | awk '{print $0}'`
IFS=$'\n' read -rd '' -a cmds <<< "${lines}"

for window in $(seq 0 $(( ${#cmds[@]} - 1 ))); do
	wait=`echo "${window}%${thread}" | bc`
	if [[ ${wait} == 0 ]]; then
		wait
	fi
	${cmds[$window]} &
done
wait
}

if [ -z "$2" ]; then
	thread=1
else
	thread=$2
fi
run_all_file $1 ${thread}