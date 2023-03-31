#!/usr/bin/env python

import os,sys
import time

outputfile = 'output.txt'
if len(sys.argv)>2:
       outputfile = sys.argv[2]

print(outputfile)
with open(outputfile,'w') as f:
       if len(sys.argv)>1:
              f.write('Hello from python tool ' + sys.argv[1])
       else:
              f.write('Hello with no arguments')

time.sleep(60)

