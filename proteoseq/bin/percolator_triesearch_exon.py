#!/usr/bin/python

from optparse import OptionParser
from collections import defaultdict
import os,sys,warnings,re,glob
import logging
import subprocess
import TrieSearch

def main():
	usage = 'usage: %prog <options> -p junctionPep -c cruxfile -e Alu.unique.bed -j SJdir -t tmpdir -b bedtooldir -o outfile'
	parser = OptionParser(usage)
	parser.add_option('-p', dest='pepfile', help='junction pep file [Default %default]')
	parser.add_option('-c', dest='cruxfile', help='percolator.target.peptides.txt from crux [Default %default]')
	parser.add_option('-e', dest='exonfile', help='Ensembl_Alu_25bp_0.5.unique.sorted.bed [Default %default]')
	parser.add_option('-j', dest='sjfile', help='SJ.out.tab file from STAR [Default %default]')
	parser.add_option('-t', dest='tmpdir',default='tmp',help='tmp dir [Default %default]')
	parser.add_option('-b', dest='bedtooldir',help='bedtool bin directory path')
	parser.add_option('-o', dest='outfile',default='result.txt',help='final output file')
	
	(options, args) = parser.parse_args()

	if options.sjfile is None or options.pepfile is None or options.cruxfile is None or options.exonfile is None or options.bedtooldir is None:
		sys.exit("[ERROR] "+parser.get_usage())
	
	warnings.formatwarning = custom_formatwarning
	
	## start
	# step 1: '''read percolator peptides under FDR 0.05'''
	sample = os.path.basename(options.pepfile)
	warnings.warn('======'+sample)
	warnings.warn('# Percolator result file:\t' + options.cruxfile)
	warnings.warn('# Exon bedfile:\t' + options.exonfile)
	
	pephash,firstAA = {},{}
	with open(options.cruxfile,'r') as f:
		for line in f:
			ele = line.rstrip().split("\t")
			if ele[0] != 'PSMId' and float(ele[2]) < 0.05:
				pep = re.sub(re.compile("\.|\*|-"),'',ele[4])
				pephash[pep] = ele[3]
				firstAA[str.upper(pep)[0]] = ''
	
	
	# step 2: '''stor junction coordinate into dict for novel or known junction judge'''
	junctionHash = {}
	warnings.warn('# Building index for junction from:\t' + options.sjfile)
	SJfiles = [options.sjfile]
	for sjf in SJfiles:
		with open(sjf,'r') as f:
			for line in f:
				ele = line.rstrip().split("\t")
				strand = '+' if ele[3] == '1'else '-'
				junction = ele[0]+'_'+str(int(ele[1])-1)+'_'+str(int(ele[2])+1) + '_' + strand
				junctionHash[junction] = ele[5]

	# step 3: '''store exon coordinate into dict for further mapping'''
	warnings.warn('# Read bed file and generate position Hash')
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

	# ste 4: ''' store junction pep into dict'''
	warnings.warn("# Read junction pep and write to bed")
	seq,head = {},''
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

	# step 5: ''' bedtool to find the location of exon'''
	warnings.warn('# Bedtools to find Exon location')
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

	# step 6: ''' Re-map percolator peptides to Junction Pep'''
	warnings.warn('# Re-map percolator peptides to Junction Pep')
	chrPepHash,info = search(pephash, seq, exonloc, exonPosHash, firstAA)
	
	# step 7: ''' FDR filter'''
	warnings.warn('# Get peptide under fdr')
	novelChrPepHash = {}
	for peptide in info:
		ks = info[peptide]
		for k in ks:
			head = k.split("\t")[5]
			arr = head.split("_")
			junc = '_'.join(arr[0:3]+arr[4:5])
			if junc not in junctionHash:
				sys.exit("[ERROR]:\t junction '%s' not exists in junction files\n" % options.sjfile)
			if junctionHash[junc] == '0':
				#print '--',junc,junctionHash[junc]
				novelChrPepHash[peptide] = chrPepHash[peptide]

	warnings.warn("# Num of pep from novel junction:\t %s" % len(novelChrPepHash))
	warnings.warn("# Num of pep from annot junction:\t %s" % (len(chrPepHash) - len(novelChrPepHash)))

	result = ['# Peptide\tStart\tEnd\tregion of exon\ttag\tJunction ID\tJunction Pep\tExon\tExon S\tExon E\tPEP\tKnown/Novel']
	
	for peptide in chrPepHash.keys():
		if peptide in novelChrPepHash: continue
		ks = info[peptide]
		for k in ks:
			s = k +"\t" + chrPepHash[peptide] + "\t1"
			result.append(s)

	localResult = localFDR(novelChrPepHash)
	arr = re.split(re.compile(";"),localResult)
	for a in arr:
		ele = a.split("\t")
		ks = info[ele[0]]
		for k in ks:
			s = k + "\t" + chrPepHash[ele[0]] + "\t0"
			result.append(s)
	
	# step 8: ''' output to files'''
	with open(options.outfile,'w') as fout:
		fout.write('======'+sample+'\n')
		fout.write('# Percolator result file:\t' + options.cruxfile + '\n')
		fout.write('# Exon bedfile:\t' + options.exonfile + '\n')
		fout.write('# Junction file:\t' + options.sjfile + '\n')
		for r in result:
			fout.write(r+'\n')
	## end

def localFDR(chrPepHash,FDR=0.01):
	''' do local FDR filter'''
	if len(chrPepHash) == 0: return ''
	pepvalues = [float(x) for x in chrPepHash.values()]
	uniq = sorted(list(set(pepvalues)))
	warnings.warn("# Before FDR %f filter:\t%d" %(FDR,len(chrPepHash)))
	result,n = '',0
	for s in uniq:
		sums,content,num,total = 0.0,'',0,len(chrPepHash.keys())
		for l in chrPepHash.keys():
			pep = float(chrPepHash[l])
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

def aaindex(start,end):
	length = abs(end - start)+1
	index = length / 3 if length % 3 == 0 else (length - length % 3) / 3 + 1
	return index

def search(pephash, seqs, exonloc, exonPosHash, firstAA):
        info = defaultdict(dict)
	chrPepHash = dict()
        trieTree = TrieSearch.TrieSearch(pephash.keys())
        trieTree.make_trie()
        for aa in firstAA:
                warnings.warn("## AA:\t"+aa)
                result = trieTree.doSearch(seqs,[aa])
                for r in result:
                        ele = r.split('\t')
                        pep,ids = ele[1],ele[2].split(';')
                        for i in ids:
				index = seqs[i].upper().find(pep.upper())
				start,end = index, len(pep) + index - 1
				if i not in exonloc:
					warnings.warn('[ERROR]: \'' + i + '\' \'' +aa +'\' not in exonloc. [Check the exon file input]')
					continue
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
				endPos = startPos-1 + len(pep)*3 if strand == '+' else startPos+1-len(pep)*3
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
				chrPepHash[pep] = pephash[pep]
				s = "\t".join([pep,str(start),str(end),str(k),tag,str(i),seqs[i],exons,str(startPos),str(endPos)])
				info[pep][s] = ''
	return (chrPepHash,info)	
	
def custom_formatwarning(msg, *a):
	# ignore everything except the message
	return str(msg) + '\n'

if __name__ == '__main__':
	main()
