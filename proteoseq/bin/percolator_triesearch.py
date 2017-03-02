#!/usr/bin/python

from optparse import OptionParser
from collections import defaultdict
import os,sys,warnings,re,glob
import logging
import subprocess
import TrieSearch

def main():
	usage = 'usage: %prog <options> -p junctionPep -c percolatorfile -t tmpdir -b bedtooldir -o outfile'
	parser = OptionParser(usage)
	parser.add_option('-p', dest='pepfile', help='junction pep file [Default %default]')
	parser.add_option('-c', dest='cruxfile', help='percolator.target.peptides.txt from crux [Default %default]')
	parser.add_option('-t', dest='tmpdir',default='tmp',help='tmp dir [Default %default]')
	parser.add_option('-b', dest='bedtooldir',help='bedtool bin directory path')
	parser.add_option('-o', dest='outfile',default='result.txt',help='final output file')
	
	(options, args) = parser.parse_args()

	if options.pepfile is None or options.cruxfile is None or options.bedtooldir is None:
		sys.exit("[ERROR] "+parser.get_usage())
	
	warnings.formatwarning = custom_formatwarning
	
	## start
	# step 1: '''read percolator peptides under FDR 0.05'''
	sample = os.path.basename(options.pepfile)
	warnings.warn('======'+sample)
	warnings.warn('# Obtain percolator result:\t' + options.cruxfile)
	
	pephash,firstAA = {},{}
	with open(options.cruxfile,'r') as f:
		for line in f:
			ele = line.rstrip().split("\t")
			if ele[0] != 'PSMId' and float(ele[2]) < 0.05:
				pep = re.sub(re.compile("\.|\*|-"),'',ele[4])
				pephash[pep] = ele[3]
				firstAA[str.upper(pep)[0]] = ''

	# ste 2: ''' store junction pep into dict'''
	warnings.warn("# Read Junction Pep fasta:\t" + options.pepfile)
	seq = defaultdict(str)
	head = ''
	with open(options.pepfile,'r') as f:
		for line in f:
			if line[0] == '>':
				head = line.rstrip()[1:]
			else: seq[head] += line.rstrip()

	# step 3: ''' Re-map percolator peptides to Junction Pep'''
	warnings.warn('# Re-map percolator peptides to Junction Pep')
	info = search(pephash,seq,firstAA)

	# step 4: ''' FDR filter'''
	warnings.warn('# Get peptide under fdr')
	allPepHash = defaultdict(dict)
	novelChrPepHash = {}
	for peptide in info:
		ks = info[peptide]
		tag = 0
		for k in ks:
			head = k.split("\t")[1]
			if head.find('startTrim:') == -1: # if peptide could be mapped to uniprot, means annotated peptide
				tag = 1
		if tag == 0:
			novelChrPepHash[peptide] = pephash[peptide]

	warnings.warn("# Num of pep from novel junction:\t %s" % len(novelChrPepHash))
	warnings.warn("# Num of pep from annot junction:\t %s" % (len(pephash) - len(novelChrPepHash)))

	result = ['# Peptide\tProtein\tPEP\tKnown/Novel']
	for peptide in pephash.keys():
		if peptide in novelChrPepHash: continue
		ks = info[peptide]
		for k in ks:
			s = k+"\t"+pephash[peptide] + "\t1"
			result.append(s)
	
	localResult = localFDR(novelChrPepHash)
	arr = re.split(re.compile(";"),localResult)
	for a in arr:
		ele = a.split("\t")
		ks = info[ele[0]]
		for k in ks:
			s = k + "\t"+pephash[ele[0]] + "\t0"
			result.append(s)
	
	# step 5: ''' output to files'''
	with open(options.outfile,'w') as fout:
		fout.write('======'+sample+'\n')
		fout.write('# Obtain percolator result from ' + options.cruxfile + '\n')
		fout.write('# Read Junction Pep fasta:\t' + options.pepfile + '\n')
		for r in result:
			fout.write(r+'\n')

	## end
def localFDR(chrPepHash,FDR=0.05):
	if len(chrPepHash) == 0: return ''
	uniq = sorted(list(set(chrPepHash.values())))
	warnings.warn("# Before FDR %f filter:\t%d" %(FDR,len(chrPepHash)))
	result,n = '',0
	for s in uniq:
		sums,content,num,total = 0.0,'',0,len(chrPepHash.keys())
		for l in chrPepHash.keys():
			pep = chrPepHash[l]
			if pep <= s:
				sums += float(pep)
				num += 1
				content += str(l)+"\t"+str(pep)+";"
		fdr = sums / num
		if fdr > FDR: break
		n += 1
		result = content
	warnings.warn( "# After FDR %f filter:\t%d" %(FDR,n))
	return result

def search(pephash,seqs,firstAA):
	info = defaultdict(dict)
	trieTree = TrieSearch.TrieSearch(pephash.keys())
	trieTree.make_trie()
	for aa in firstAA:
		warnings.warn("## AA:\t"+aa)
		result = trieTree.doSearch(seqs,[aa])
		for r in result:
			ele = r.split('\t')
			pep,ids = ele[1],ele[2].split(';')
			for i in ids:
				info[pep]['\t'.join([pep,i])] = ''
	return info
	
def custom_formatwarning(msg, *a):
	# ignore everything except the message
	return str(msg) + '\n'

if __name__ == '__main__':
	main()
