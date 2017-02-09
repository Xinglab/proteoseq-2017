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
	usage = 'usage: %prog <options> -p junctionPep -c cruxfile -e Alu.unique.bed -j SJdir -t tmpdir -n threadNum -d bedtooldir'
	parser = OptionParser(usage)
	parser.add_option('-p', dest='pepfile', help='junction pep file [Default %default]')
	parser.add_option('-c', dest='cruxfile', help='percolator.target.peptides.txt from crux [Default %default]')
	parser.add_option('-e', dest='exonfile', help='Ensembl_Alu_25bp_0.5.unique.sorted.bed [Default %default]')
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
	print '# read bed file and generate position Hash'
	exonPosHash = defaultdict(dict)
	with open(options.exonfile,'r') as f:
		for line in f:
			ele = line.rstrip().split("\t")
			strline = line.rstrip().replace("\t","_")
			s1,s2 = ele[0]+'_'+ele[5]+'exonLeft',int(ele[1])+1
			s3,s4 = ele[0]+'_'+ele[5]+'exonRight',int(ele[2])
			if s2 in exonPosHash[s1]:
				exonPosHash[s1][s2].append(strline)
			else:
				exonPosHash[s1][s2] = [strline]

			if s4 in exonPosHash[s3]:
				exonPosHash[s3][s4].append(strline)
			else:
				exonPosHash[s3][s4] = [strline]

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
	print '# bedtools to find Exon location'
	exonloc = defaultdict(dict)
	p = subprocess.Popen(options.bedtooldir + 'bedtools coverage -a '+options.tmpdir+'/'+sample+'.bed -b '+ options.exonfile + ' -s', shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
	for line in p.stdout.readlines():
		ele = line.rstrip().split("\t")
		if ele[7] == '0': continue
		strand = ele[5]
		loc = ele[3].split("_")[-1]
		myid = re.sub(re.compile('_left|_right'),'',ele[3])
		sequence = seq[myid]
		ja,jb,trimaa,orf = myid.split("_")[5:9]
		trimaa = int(trimaa.replace('startTrim:',''))
		orf = int(orf.replace('ORF:',''))
		a1,a2,a3 = [int(x) for x in ja.split(',')]
		b1,b2,b3 = [int(x) for x in jb.split(',')]
		startPos = a1 + trimaa*3 + orf if strand == '+' else b3 - trimaa*3 - orf
		if strand == '+' and trimaa*3 >= (a3 - (a1 + orf) + 1): startPos = b1 + trimaa*3 - (a3 - (a1 + orf) + 1)
		if strand == '-' and trimaa*3 >= (b3 - orf - b1 + 1): startPos = a3 - (trimaa*3 - (b3 - orf - b1 + 1))
		if strand == '+':
			if loc == 'left':
				if a3 >= startPos:
					s = aaindex(startPos, a3) - 1
					exonloc[myid]['0-'+str(s)] = strand+"\t"+loc
				else:
					exonloc[myid]['NA-NA'] = strand+"\t"+loc
			else:
				if a3 >= startPos or startPos == b1:
					s = 0 if startPos == b1 else aaindex(startPos, a3+1) - 1
					exonloc[myid][str(s)+'-'+str(len(sequence)-1)] = strand +"\t"+loc
				else:
					exonloc[myid]['NA-NA'] = strand+"\t"+loc
		else:
			if loc == 'left':
				if b1 <= startPos or startPos == a3:
					s = 0 if startPos == a3 else aaindex(startPos, b1-1) - 1
					exonloc[myid][str(s)+'-'+str(len(sequence)-1)] = strand +"\t"+loc
				else:
					exonloc[myid]['NA-NA'] = strand+"\t"+loc
			else:
				if b1 <= startPos:
					s = aaindex(startPos, b1) - 1
					exonloc[myid]['0-'+str(s)] = strand+"\t"+loc
				else:
					exonloc[myid]['NA-NA'] = strand+"\t"+loc
	retval = p.wait()

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
		sthread = threading.Thread(target = run, args = (str(i),pephash,arrHash,arrStart,arrSeq2,seq,exonloc,exonPosHash))
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
			head = k.split("\t")[5]
			arr = head.split("_")
			junc = '_'.join(arr[0:3]+arr[4:5])
			if junc not in junctionHash:
				sys.exit("[ERROR]:\t junction '%s' not exists in junction files\n" % options.sjfile)
			if junctionHash[junc] == 0:
				novelChrPepHash[peptide] = chrPepHash[peptide]

	warnings.warn("# num of pep from novel junction:\t %s" % len(novelChrPepHash))
	warnings.warn("# num of pep from annot junction:\t %s" % (len(chrPepHash) - len(novelChrPepHash)))

	
	for peptide in chrPepHash.keys():
		if peptide in novelChrPepHash: continue
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

def run(tNum,pephash,arrHash,arrStart,arrSeq2,seq,exonloc,exonPosHash):
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
				start,end = index, len(pepj) + index - 1
				k = ';'.join(exonloc[i].keys())
				tag = 'N'
				for k2 in exonloc[i].keys():
					if k2.find('NA') != -1: continue
					s,e = [int(x) for x in k2.split("-")]
					mini = min([s,e,start,end])
					maxi = max([s,e,start,end])
					if abs(end-start) + 1 + abs(e-s) + 1 > maxi-mini+1: tag = 'Y'
				ele = i.split('_')
				chrom,jl,jr,strand,annol,annor,startTrim,ORF = ele[0],int(ele[1]),int(ele[2]),ele[4],ele[5],ele[6],ele[7],ele[8]
				startTrim = startTrim.replace("startTrim:",'')
				ORF = int(ORF.replace('ORF:',''))
				startTrim = int(startTrim) + start
				l = int(annol.split(',')[0])
				r = int(annol.split(',')[2])
				startPos = l+startTrim*3+ORF if strand == '+' else r-startTrim*3-ORF
				if strand == '+' and startPos > jl:
					startPos = jr + startPos -jl - 1
				if strand == '-' and startPos < jr:
					startPos = jl - (jr - startPos) + 1
				endPos = startPos-1 + len(pepj)*3 if strand == '+' else startPos+1-len(pepj)*3
				if strand == '+' and startPos <= jl and endPos > jl:
					endPos = jr + endPos - jl - 1
				if strand == '-' and startPos >= jr and endPos < jr:
					endPos = jl - (jr - endPos) + 1
				exons = ''
				if jr in exonPosHash[chrom+"_"+strand+"exonLeft"]:
					arr = exonPosHash[chrom+"_"+strand+"exonLeft"][jr]
					for a in arr:
						minLeft,maxRight = a.split("_")[1:3]
						minLeft = int(minLeft)+1
						maxRight = int(maxRight)
						p = sorted([maxRight, minLeft, startPos, endPos])
						dis = abs(maxRight-minLeft)+1+abs(endPos-startPos)+1
						overlapStatus = 1 if abs(p[3]-p[0]) + 1 < dis else 0
						if maxRight >= startPos and maxRight >= endPos and overlapStatus == 1: exons = a + ';' + exons
				if jl in exonPosHash[chrom+"_"+strand+"exonLeft"]:
					arr = exonPosHash[chrom+"_"+strand+"exonLeft"][jl]
					for a in arr:
						minLeft,maxRight = a.split("_")[1:3]
						minLeft = int(minLeft)+1
						maxRight = int(maxRight)
						p = sorted([maxRight, minLeft, startPos, endPos])
						dis = abs(maxRight-minLeft)+1+abs(endPos-startPos)+1
						overlapStatus = 1 if abs(p[3]-p[0]) + 1 < dis else 0
						if minLeft <= startPos and maxLeft <= endPos and overlapStatus == 1: exons = a + ';' + exons

				chrPepHash[pepj] = pephash[pepj]
				s = "\t".join([pepj,str(start),str(end),str(k),tag,str(i),seq[i],exons,str(startPos),str(endPos)])
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
