#!/bin/bash
function run_all_cmd() {
cmds=$1

for window in $(seq 0 $(( ${#cmds[@]} - 1 ))); do
	${cmds[$window]} &
done
wait
}

IFS=':' read -r -a cmds <<< "$@"

run_all_cmd "${cmds[@]}"