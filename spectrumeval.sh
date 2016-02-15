#!/bin/bash

exec python spectrumeval.py

echo $*

while read line; do
	echo $line >> spectrumeval.log
	echo '>>>test.json'
done



