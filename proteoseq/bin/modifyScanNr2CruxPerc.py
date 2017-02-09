#!/usr/bin/python

import glob,sys,re
import os.path

def main():
	if len(sys.argv) < 3:
		sys.exit("[ERROR]: "+sys.argv[0] + " cometpindir")

	indir = re.sub(re.compile("\/$"),'',sys.argv[1])
	outdir = sys.argv[2]
	if not os.path.exists(outdir):
		os.makedirs(outdir)
	outtab = outdir + '/' + os.path.basename(indir)+'.2cruxprec'
	i = 1
	tag = 0
	with open(outtab,'w') as out:
		pinfiles = glob.glob(indir + '/*.pin')
		for f in pinfiles:
			#if f != 'GM18486.rna/GM18486_MSB15024_01B.pin':
			#	continue
			prefix = str(i) + '00'
			with open(f,'r') as fh:
				for line in fh:
					line = line.rstrip()
					if line.find('id') == 0 and tag == 0:
						out.write(line+"\n")
						tag = 1
					elif line.find('id') != 0:
						ele = line.split("\t")
						arr = ele[0].split("_")
						arr[len(arr)-3] = prefix + arr[len(arr)-3]
						ele[2] = arr[len(arr)-3]
						ele[0] = '_'.join(arr)
						out.write("\t".join(ele)+"\n")
			i = i + 1

if __name__=='__main__':
	main()
