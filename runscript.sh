#!/bin/bash

test_run()
{
  paramfolder="$1"
  cd "$paramfolder"
  echo python ../main.py "$paramfolder.json"
  cd ..
}
run()
{
  paramfolder="$1"
  cd "$paramfolder"
  python ../main.py "$paramfolder.json"
  cd ..
}

export -f test_run 
export -f run

test_run "$1"
	echo "Enter to continue, Ctrl+c to end"
	read option
	if [[ -z $option ]]; then
	 	run "$1"
	else
		echo "run program denied!"
	fi 