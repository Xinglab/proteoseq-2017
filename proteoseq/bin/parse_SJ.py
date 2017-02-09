#!/usr/bin/python
from optparse import OptionParser
from collections import defaultdict
import os

def main():
	parser = OptionParser('') # parse_SJ.py RNA/GM18486.rna/SJ.out.tab
	(options, args) = parser.parse_args()
	args[0]=os.path.abspath(args[0])
	d = readSJ(args[0])
	with open(args[1], 'w') as out:
		for i in d[1]:
			#if i.find('chr14') == -1:
				# continue
			key3 = d[0][i+ '_3'].keys()
			key3 = map(int, key3)
			key3.sort()
			key5 = d[0][i+ '_5'].keys()
			key5 = map(int, key5)
			key5.sort()
			lines = d[1][i]		
			for j in lines:
				if len(key5) == 1:
					out.write(j+"\t-1\t-1\t-1\t-1\n")
					continue
				ele = j.split()			
				x = leftSearch(int(ele[1]),key5,0,len(key5)-1)
				y = rightSearch(int(ele[2]),key3,0,len(key3)-1)
				out.write(j+"\t"+str(x)+"\t"+str(y)+"\t"+str(int(ele[1])-x)+"\t"+str(y-int(ele[2]))+"\n")


def leftSearch(key, arr, low, high):
	mid = int((high+low)/2)
	if low + 1 == high:
		return arr[low]
	if arr[mid] >= key:
		high = mid
	else:
		low = mid
	#print low,"\t",high
	return leftSearch(key, arr, low, high)

def rightSearch(key, arr, low, high):
	mid = int((high+low)/2)
	if low + 1 == high:
		return arr[high]
	if arr[mid] > key:
		high = mid
	else:
		low = mid
	#print low,"\t",high
	return rightSearch(key, arr, low, high)

def readSJ(file):
	rev_dict = defaultdict(dict)
	data = defaultdict(list)
	with open(file, 'r') as f:
		for line in f:
        		line = line.rstrip()
			ele = line.split()
			rev_dict[ele[0]+"_"+ele[3]+"_3"].update({ele[1]:''})
			rev_dict[ele[0]+"_"+ele[3]+"_5"].update({ele[2]:''})
			data[ele[0]+"_"+ele[3]].append(line)
	return [rev_dict, data]

if __name__=='__main__':
	main()
