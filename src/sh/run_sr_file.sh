#!/bin/bash
function run_all_file() {
filename=$1

lines=`cat ${filename} | awk '{print $0}'`
IFS=$'\n' read -rd '' -a cmds <<< "${lines}"

for window in $(seq 0 $(( ${#cmds[@]} - 1 ))); do
	${cmds[$window]}
done
}

run_all_file $1