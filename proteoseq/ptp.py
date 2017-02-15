#!/usr/bin/python

from optparse import OptionParser
from collections import defaultdict
import os,sys,random,datetime,warnings,re
import logging
import subprocess
import ConfigParser

OUTDIR,BINDIR,CHROMS,WINE,COMETEXE,COMETPAR,CRUX,BEDTOOLDIR,PERCOLATOR = 'outdir','','','','','','','',''

def main():
	global OUTDIR,BINDIR,CHROMS,WINE,COMETEXE,COMETPAR,CRUX,BEDTOOLDIR,PERCOLATOR
	## read global parameters
	config = ConfigParser.ConfigParser()
	config.readfp(open('config.ini',"rb"))
	BINDIR = config.get('global','BINDIR')
	CHROMS = config.get('global','CHROMS')
	WINE = config.get('global','WINE')
	COMETPAR = config.get('global','COMETPAR')
	COMETEXE = config.get('global','COMETEXE')
	CRUX = config.get('global','CRUX')
	BEDTOOLDIR = config.get('global','BEDTOOLDIR')
	PERCOLATOR = config.get('global','PERCOLATOR')

	## read pipeline parameters
	usage = 'usage: %prog -b Aligned.out.sorted.bam -j SJ.tab.out -p proteomicsdir -e HSExonfile/None -o outdir --l 66 --g genome_file --min-junc-reads 2 --trim-RK False'
	parser = OptionParser(usage)
	# necessary parameters
	parser.add_option('-b','--bamfile', dest='bamfile', help='bam file from STAR [Default %default]')
	parser.add_option('-j','--sjfile', dest='sjfile', help='SJ.tab.out file from STAR [Default %default]')
	parser.add_option('-p','--rawdir', dest='rawdir', help='proteomics dir (format: raw/mzXML) [Default %default]')
	parser.add_option('-e','--exonfile', dest='exonfile',default='None', help='exons file (bed format) or use None to search all junctions peptide [Default %default]')
	parser.add_option('-o', dest='outdir',default=OUTDIR, help='Output dir filename [Default %default]')	
	# parameters for translation
        parser.add_option('--l', dest='flank', type='int', default=66, help='Extend flanking junction ends by this number of bp [Default %default]')
        parser.add_option('--g', dest='genome_file', default=CHROMS, help='genomic fasta directory (by chromosomes) [Default %default]')
        parser.add_option('--min-junc-reads', dest='min_junc_reads', default=2, type='int', help='Minimum number of reads required spanning the junction [Default %default]')
	(options, args) = parser.parse_args()
	# check parameters
	if options.sjfile is None or options.bamfile is None or options.rawdir is None:
		sys.exit("[ERROR] "+parser.get_usage())
	if options.exonfile != 'None' and not os.path.exists(options.exonfile):
		sys.exit("[ERROR] Please input exon file or use '-e None'\n"+parser.get_usage())
	if options.outdir is not None:
		OUTDIR = re.sub(re.compile("/$"),"",options.outdir)

	# mkdir 'outdir' if not exists
	if not os.path.exists(OUTDIR):
		os.makedirs(OUTDIR)

	# warning and logging
	warnings.formatwarning = custom_formatwarning
	logging.basicConfig(filename=OUTDIR+'/pipeline.log', level=logging.INFO)

	## start
	# 1. parse SJ.tab.out file
	warnings.warn("## starting parse junction files")
	sjfilename = parseSJ(options.sjfile)
	warnings.warn('## output file name: %s' % sjfilename.replace('.SJ',''))
	# 2. translate junctions to peptide
	warnings.warn("## starting translating into junction peptides")
	fastaname = translate(options.bamfile,OUTDIR+"/"+sjfilename,options.exonfile,options.flank,options.genome_file,options.min_junc_reads,OUTDIR+"/"+sjfilename.replace('SJ','fa'))
	# 3. merge peptides to uniprot database
	customdb = mergePeps2Database(fastaname)
	# 4. database search using comet
	cometoutdir = databaseSearch(options.rawdir, OUTDIR+"/"+customdb)
	# 5. percolator (or crux percolator)
	percolatorfile = percolator(cometoutdir)
	# 6. post-filter(extract peptide mapped Alu/HSE exons, and FDR filter)
	postPercolatorFilter(OUTDIR+'/'+fastaname,percolatorfile,options.exonfile,OUTDIR+"/"+sjfilename)

def test(outdir,outfile):
	global OUTDIR,BINDIR,CHROMS,WINE,COMETEXE,COMETPAR,CRUX,BEDTOOLDIR,PERCOLATOR
	## read global parameters
	config = ConfigParser.ConfigParser()
	config.readfp(open('config.ini',"rb"))
	BINDIR = config.get('global','BINDIR')
	CHROMS = config.get('global','CHROMS')
	WINE = config.get('global','WINE')
	COMETPAR = config.get('global','COMETPAR')
	COMETEXE = config.get('global','COMETEXE')
	CRUX = config.get('global','CRUX')
	BEDTOOLDIR = config.get('global','BEDTOOLDIR')
	PERCOLATOR = config.get('global','PERCOLATOR')
	OUTDIR = outdir

	fastaname = outdir+'/'+outfile+'.fa'
	#sjfile = 'data/SJ_out/LCLs/GM18486.rna.SJ'
	sjfile = outdir+'/'+outfile+'.SJ'
	customdb = outdir+'/merge_'+outfile+'.fa'
	cometoutdir = 'comet_'+outfile
	#percolatorCrux(cometoutdir)
	#percolator(cometoutdir)
	percolatorfile = outdir + '/percolator_' + outfile + '/percolator.target.peptides.txt'
	exonfile = 'data/Ensembl_Alu_25bp_0.5.unique.sorted.bed'
	postPercolatorFilter(fastaname,percolatorfile,exonfile,sjfile)
	
def parseSJ(SJfile):
	if os.path.exists(SJfile) == False:
		sys.exit("[ERROR]: Junction file not exists!\n")
	outfile = datetime.datetime.now().strftime('%Y%m%d%H%M%S') + "-" + str(random.randint(1000,9999)) + ".SJ"
	os.system("python "+BINDIR+"/"+"parse_SJ.py " + SJfile + " " + OUTDIR+ "/" + outfile)
	return os.path.basename(outfile)

def translate(bamfile,sjfile,exonfiles,flank,genome_file,min_junc_reads,outfile):
	if os.path.exists(bamfile) == False:
		sys.exit("[ERROR]: Bam file not exists!\n")
	scriptname = 'translateJunc-star.py' if exonfiles != 'None' else 'translateJunc-star-allSJ.py'
	cmd = "python "+BINDIR + "/" + scriptname + " -o " + outfile + " -l " + str(flank) + " --min-junc-reads=" + str(min_junc_reads) + " -g " + genome_file + " " + bamfile + " " + sjfile + " " + exonfiles
	os.system(cmd)
	return os.path.basename(outfile)

def mergePeps2Database(fastafile):
	os.system("cat data/UP000005640_9606_additional_cdhit1.fasta " + OUTDIR+"/"+fastafile + ">" + OUTDIR + "/merge_" + fastafile)
	return os.path.basename(OUTDIR + "/merge_" + fastafile)

def databaseSearch(rawdir, database):
	warnings.warn('## database search start')
	COMETOUTDIR = OUTDIR + "/comet_" + re.sub(re.compile("merge_|\.fa$"),"",os.path.basename(database))
	RAWDIR = re.sub(re.compile("/$"),"",rawdir)
	mylogger = logging.getLogger("comet")

	allfiles = os.listdir(RAWDIR)
	rawfiles = [x for x in allfiles if re.search('\.raw|\.mzXML',x) is not None]
	if not os.path.exists(COMETOUTDIR): os.makedirs(COMETOUTDIR)
	for i in rawfiles[0:]:
		warnings.warn("\t# raw file:\t%s" % i)
		inf = RAWDIR + '/' + i
		outf = COMETOUTDIR + '/' + re.sub(re.compile("\..*$"),"",i)
		cmd = "WINEDEBUG=fixme-all,err-all " + WINE +" "+ COMETEXE +" -P"+ COMETPAR+" -D"+database+" -N"+outf+" "+inf
		warnings.warn("\t"+cmd)
		p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
		for line in p.stdout.readlines():
			mylogger.info(line.rstrip())
		retval = p.wait()
	return os.path.basename(COMETOUTDIR)

def percolatorCrux(cometdir):
	cometfiles = os.listdir(OUTDIR+'/'+cometdir)
	cruxPrecTmp = OUTDIR+'/'+cometdir.replace('comet','cruxPrecTmp')	
	cruxoutdir = OUTDIR+'/'+cometdir.replace('comet','cruxoutput')
	if not os.path.exists(cruxoutdir):
                os.makedirs(cruxoutdir)
	os.system('python ' + BINDIR+'/modifyScanNr2CruxPerc.py'+' '+OUTDIR+'/'+cometdir +' '+cruxPrecTmp)
	os.system(CRUX +' percolator --train-fdr 0.05 --test-fdr 0.05 --overwrite T --output-dir ' + cruxoutdir + ' ' + cruxPrecTmp + '/' + cometdir + '.2cruxprec')
	return cruxoutdir + '/percolator.target.peptides.txt'

def percolator(cometdir):
	cometfiles = os.listdir(OUTDIR+'/'+cometdir)
	PercolatorOutDir = OUTDIR+'/'+cometdir.replace('comet','percolator')
	if not os.path.exists(PercolatorOutDir):
		os.makedirs(PercolatorOutDir)
	os.system('cat ' + OUTDIR+'/'+cometdir + '/*.pin > ' + PercolatorOutDir + '/percolatorTmp.pin')
	# merge multiple pin files to single one to input to percolator
	data = []
	n = 0
	with open(PercolatorOutDir+'/percolatorTmp.pin', 'r') as f:
		for line in f:
			line = line.rstrip('\r\n')
			if n == 0:
				data.append(line)
				n += 1
			if line.find('id') != 0:
				data.append(line)
	with open(PercolatorOutDir+'/percolatorTmp.pin', 'w') as fout:
		for d in data:
			fout.write(d+'\n')
	os.system(PERCOLATOR + ' ' + PercolatorOutDir+'/percolatorTmp.pin -t 0.05 -F 0.05 -m ' + PercolatorOutDir + '/percolator.target.psm.txt -r '+ PercolatorOutDir +'/percolator.target.peptides.txt')
	return PercolatorOutDir + '/percolator.target.peptides.txt'

def postPercolatorFilter(fastaname,percolatorfile,exonfille,sjfile, threadNum=1):
	tmpdir = OUTDIR+'/tmp'
	scriptname = 'percolator_test2parellel.py' if exonfille != 'None' else 'percolator_directoutput.py'
	#os.system(BINDIR+'/cruxpep_percolator_test2parellel.py -p %s -c %s -e %s -j %s -t %s -n %d -d %s' % (fastaname,percolatorfile,exonfille,sjfile,tmpdir,threadNum,BEDTOOLDIR))
	os.system('python ' + BINDIR+'/' + scriptname + ' -p %s -c %s -e %s -j %s -t %s -n %d -d %s' % (fastaname,percolatorfile,exonfille,sjfile,tmpdir,threadNum,BEDTOOLDIR))

def custom_formatwarning(msg, *a):
	# ignore everything except the message
	return str(msg) + '\n'

if __name__ == '__main__':
	main()
	#test(sys.argv[1],sys.argv[2])
