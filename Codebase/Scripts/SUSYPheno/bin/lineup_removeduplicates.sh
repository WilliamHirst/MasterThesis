#!/bin/bash
lineup.py | awk '!x[$0]++' | lineup.py 
#lineup.py | awk_removeduplicates | lineup.py 

