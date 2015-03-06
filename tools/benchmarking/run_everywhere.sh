#!/bin/bash

for node in `cat nodes`
do
	ssh ${node} "cd `pwd`; nohup ./test_cuda.sh >${node}.out 2>${node}.err < /dev/null &"
done
