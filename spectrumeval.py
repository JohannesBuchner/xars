#!/usr/bin/env python

import sys

logfile = open('spectrumeval.log', 'a')
while True:
	line = sys.stdin.readline()
	if line == '':
		break
	logfile.write(line)
	sys.stdout.write('>>>test.json\n')
	sys.stdout.flush()


