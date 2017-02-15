#!/usr/bin/python

from optparse import OptionParser
from collections import defaultdict
import os,sys,warnings,re,glob
import logging
import subprocess
import threading

# control for multi-thread
tindex = 0

mutex = threading.Lock()

# result (global variables used in multi-thread)
chrPepHash = {}
infoArr = []

def main():
	usage = 'usage: %prog <options> -p junctionPep -c percolatorfile -e Alu.unique.bed -j SJdir -t tmpdir -n threadNum -d bedtooldir'
	parser = OptionParser(usage)
	parser.add_option('-p', dest='pepfile', help='junction pep file [Default %default]')
	parser.add_option('-c', dest='cruxfile', help='percolator.target.peptides.txt from crux [Default %default]')
	parser.add_option('-e', dest='exonfile',default='None', help='Ensembl_Alu_25bp_0.5.unique.sorted.bed [Default %default]')
	parser.add_option('-j', dest='sjfile', help='SJ.out.tab file from STAR [Default %default]')
	parser.add_option('-t', dest='tmpdir',default='tmp',help='tmp dir [Default %default]')
	parser.add_option('-n', dest='threadNum',type='int',default=1,help='thread num [Default %default]')
	parser.add_option('-d', dest='bedtooldir',help='bedtool bin directory path')
	
	(options, args) = parser.parse_args()

	if options.sjfile is None or options.pepfile is None or options.cruxfile is None or options.exonfile is None or options.bedtooldir is None:
		sys.exit("[ERROR] "+parser.get_usage())
	
	warnings.formatwarning = custom_formatwarning
	
	## start
	# 1
	sample = os.path.basename(options.pepfile)
	print '======'+sample
	print '# obtain percolator result from '+sample
	print '# bedfile:',options.exonfile
	
	peparr = []
	pephash = {}
	seqhash = {}
	n = 0
	with open(options.cruxfile,'r') as f:
		for line in f:
			ele = line.rstrip().split("\t")
			if ele[0] != 'PSMId' and float(ele[2]) < 0.05:
				peparr.append(ele[3]+"\t"+ele[4])
				pep = re.sub(re.compile("\.|\*|-"),'',ele[4])
				pephash[pep] = ele[3]
				seqhash["seq_"+str(n)] = pep
				n += 1
	# 2
	junctionHash = {}
	print '# building index for junction SJ.out.tab'
	SJfiles = [options.sjfile]
	for sjf in SJfiles:
		with open(sjf,'r') as f:
			for line in f:
				ele = line.rstrip().split("\t")
				strand = '+' if ele[3] == '1'else '-'
				junction = ele[0]+'_'+str(int(ele[1])-1)+'_'+str(int(ele[2])+1) + '_' + strand
				junctionHash[junction] = ele[5]
	# 3
	print "# building index for chr.. seq"
	lenSeqHash = defaultdict(dict)
	hid = ''
	with open(options.pepfile,'r') as f:
		for line in f:
			if line[0] == '>': hid = line.rstrip().replace(">",'')
			else: lenSeqHash[len(line.rstrip())][hid] = ''
	
	arrSeq = []
	arrHash = {}
	for k in sorted(lenSeqHash.keys()):
		arrHash[k] = len(arrSeq) + 1 - 1
		arrSeq.append(lenSeqHash[k].keys())

	# 4

	# 5
	print "# read junction pep and write to bed"
	seq = {}
	head = ''
	if not os.path.exists(options.tmpdir): os.makedirs(options.tmpdir)
	with open(options.tmpdir+'/'+sample+'.bed','w') as fw:
		with open(options.pepfile,'r') as f:
			for line in f:
				if line[0] == '>':
					head = line.rstrip()[1:]
					arr = head.split("_")
					chrom,strand = arr[0],arr[4]
					left,jl = arr[5].split(',')[1:3]
					jr,right = arr[6].split(',')[0:2]
					fw.write("\t".join([chrom,str(int(left)-1),jl,head+'_left','0',strand,"\n"]))
					fw.write("\t".join([chrom,str(int(jr)-1),right,head+'_right','0',strand,"\n"]))
				else: seq[head] = line.rstrip()
	# 6

	# 7
	print '# build index2 for chr.. seq'
	tmplen = defaultdict(lambda: defaultdict(lambda: defaultdict(str)))
	tmpfirst = dict([(k[0], '') for k in pephash.keys()])
	for k in seq.keys():
		s = seq[k].upper()
		for f in tmpfirst.keys():
			i = s.find(f)
			if i == -1: continue
			l = len(s) - i
			tmplen[f][l][k] = ''

	arrSeq2 = defaultdict(list)
	arrStart = defaultdict(dict)
	for f in tmplen.keys():
		for l in sorted(tmplen[f].keys()):
			if f in arrSeq2:
				arrStart[f][l] = len(arrSeq2[f]) + 1- 1
				arrSeq2[f] += tmplen[f][l].keys()
			else:
				arrStart[f][l] = 0
				arrSeq2[f] = tmplen[f][l].keys()

	# 8
	print '# multiple thread to get junction peptides'
	thread_list = []
	for i in xrange(options.threadNum):
		sthread = threading.Thread(target = run, args = (str(i),pephash,arrHash,arrStart,arrSeq2,seq))
		sthread.setDaemon(True)
		sthread.start()
		thread_list.append(sthread)
	for i in xrange(options.threadNum):
		thread_list[i].join()
	
	info = defaultdict(list)
	for l in infoArr:
		ele = l.split("\t")
		info[ele[0]].append(l)

	# 9
	print '# get peptide under fdr'
	print '# peptide\tstart\tend\tAlu/HSE_exon\ttag\tid\tseq'
	novelChrPepHash = {}
	for peptide in info.keys():
		ks = info[peptide]
		for k in ks:
			head = k.split("\t")[1]
			arr = head.split("_")
			junc = '_'.join(arr[0:3]+arr[4:5])
			if junc not in junctionHash:
				sys.exit("[ERROR]:\t junction '%s' not exists in junction files\n" % options.sjfile)
			if junctionHash[junc] == 0:
				novelChrPepHash[peptide] = chrPepHash[peptide]

	warnings.warn("# num of pep from novel junction:\t %s" % len(novelChrPepHash))
	warnings.warn("# num of pep from annot junction:\t %s" % (len(chrPepHash) - len(novelChrPepHash)))

	
	for peptide in chrPepHash.keys():
		#if peptide in novelChrPepHash: continue
		ks = info[peptide]
		for k in ks:
			print k+"\t"+chrPepHash[peptide] + "\t1"
	
	localResult = localFDR(novelChrPepHash)
	arr = re.split(re.compile(";"),localResult)
	for a in arr:
		ele = a.split("\t")
		ks = info[ele[0]]
		for k in ks:
			print k + "\t"+chrPepHash[ele[0]] + "\t0"
	
	## end


def aaindex(start,end):
	length = abs(end - start)+1
	index = length / 3 if length % 3 == 0 else (length - length % 3) / 3 + 1
	return index

def localFDR(chrPepHash,FDR=0.05):
	if len(chrPepHash) == 0: return ''
	uniq = sorted(list(set(chrPepHash.values)))
	warnings.warn("# Before FDR %f filter:\t%d" %(FDR,len(chrPepHash)))
	result,n = '',0
	for s in uniq:
		sums,content,num,total = 0,'',0,len(chrPepHash.keys())
		for l in chrPepHash.keys():
			pep = chrPepHash[l]
			if pep <= s:
				sums += pep
				num += 1
				content += str(l)+"\t"+str(pep)+";"
		fdr = sums / num
		if fdr > FDR: break
		n += 1
		result = content
	warnings.warn( "# After FDR %f filter:\t%d" %(FDR,n))
	return result

def run(tNum,pephash,arrHash,arrStart,arrSeq2,seq):
	global tindex,chrPepHash,infoArr
	arrPepHash = sorted(pephash.keys())
	j = 0
	while 1:
		if mutex.acquire(1):
			if tindex == len(arrPepHash):
				mutex.release()
				break
			j = tindex
			#warnings.warn("Thread-%s\t%d" %(tNum,j))
			tindex += 1
			mutex.release()
		pepj = arrPepHash[j]
		lenj = len(pepj)
		index = arrHash[lenj] if lenj in arrHash else 0
		firstAA = pepj[0]
		iStart = arrStart[firstAA][lenj] if firstAA in arrStart and lenj in arrStart[firstAA] else 0
		iEnd = len(arrSeq2[firstAA]) - 1
		arr = arrSeq2[firstAA]
		#warnings.warn("Thread-%s\t%d\t%d\t%s\t%d\t%d\t%s\t%d" %(tNum,j,lenj,firstAA,iStart,iEnd,pepj,lenj))
		for i in arr[iStart:iEnd+1]:
			index = seq[i].upper().find(pepj.upper())
			if index != -1:
				chrPepHash[pepj] = pephash[pepj]
				s = '\t'.join([pepj,i])
				#s = "\t".join([pepj,str(start),str(end),str(k),tag,str(i),seq[i],exons,str(startPos),str(endPos)])
				infoArr.append(s)
		j += 1
	return (chrPepHash,infoArr)
			
def getindex():
	global indexThread
	indexThread += 1
	return indexThread


def custom_formatwarning(msg, *a):
	# ignore everything except the message
	return str(msg) + '\n'

if __name__ == '__main__':
	main()
