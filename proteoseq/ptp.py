#!/usr/bin/python

from optparse import OptionParser
from collections import defaultdict
import os,sys,random,datetime,warnings,re,datetime,glob
import logging
import subprocess
import ConfigParser

OUTDIR,BINDIR,CHROMS,WINE,COMETEXE,COMETPAR,CRUX,BEDTOOLDIR,PERCOLATOR = 'outdir','','','','','','','',''

def main(options):
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
	
	warnings.formatwarning = custom_formatwarning
        logging.basicConfig(filename=OUTDIR + '/pipeline.log', level=logging.INFO)	

	# print parameters
	printParameters(options, OUTDIR)

	## start
	outfilename = os.path.basename( glob.glob(OUTDIR+'/*.fa')[0].replace('.fa','')) if len(glob.glob(OUTDIR+'/*.fa')) != 0 else ''
	# 1. parse SJ.tab.out file
	sjfilename = parseSJ(options.sjfile) if options.step == 1  else outfilename + '.SJ'
	# 2. translate junctions to peptide
	fastaname = translate(options.bamfile, OUTDIR+"/"+sjfilename, options.exonfile, options.flank, options.genome_file, options.min_junc_reads, OUTDIR+"/"+sjfilename.replace('SJ','fa'), OUTDIR) if options.step <= 2 else outfilename + '.fa'
	if not os.path.exists(OUTDIR + '/' + fastaname):
		sys.exit('[ERROR]:\tTranslation failed, please alloc more memory')
	# 3. merge peptides to uniprot database
	customdb = mergePeps2Database(fastaname,options.database) if options.step <= 3 else 'merge_' + outfilename + '.fa'
	# 4. database search using comet
	cometoutdir = databaseSearch(options.rawdir, OUTDIR+"/"+customdb) if options.step <= 4 else 'comet_'+outfilename
	# 5. percolator (or crux percolator)
	percolatorfile = percolator(cometoutdir) if options.step <= 5 else 'percolator_'+outfilename
	# 6. post-filter(extract peptide mapped Alu/HSE exons, and FDR filter)
	resultfile = OUTDIR+"/"+sjfilename.replace('.SJ','.result.txt')
	if options.step <= 6:
		postPercolatorFilter(OUTDIR+'/'+fastaname,percolatorfile,options.exonfile,OUTDIR+"/"+sjfilename, resultfile)

def printParameters(options, outdir):
	time = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
	warnings.warn('# Parameters\t'+str(time))
	warnings.warn('@ Bam file:\t-b '+ options.bamfile)
	warnings.warn('@ Junction file:\t-j '+ options.sjfile)
	warnings.warn('@ Proteom dir:\t-p '+ options.rawdir)
	warnings.warn('@ Exon bed file:\t-e '+ options.exonfile)
	warnings.warn('@ Protein Database:\t-d '+ options.database)
	warnings.warn('@ Genome dir:\t-g '+ options.genome_file)
	warnings.warn('@ Output dir:\t-o '+ outdir)
	warnings.warn('@ Flanking region:\t--l ' + str(options.flank))
	warnings.warn('@ Min reads:\t--min-junc-reads ' + str(options.min_junc_reads))
	warnings.warn('@ Start step:\t--step ' + str(options.step))
	warnings.warn('# ##########\n')
	
def parseSJ(SJfile, outfilename = ''):
	warnings.warn('## Step 1:\tparse junction file start')
	if os.path.exists(SJfile) == False:
		sys.exit("[ERROR]: Junction file not exists!\n")
	outfile = datetime.datetime.now().strftime('%Y%m%d%H%M%S') + "-" + str(random.randint(1000,9999)) + ".SJ"
	if outfilename != '': outfile = outfilename
	os.system("python "+BINDIR+"/"+"parse_SJ.py " + SJfile + " " + OUTDIR+ "/" + outfile)
	return os.path.basename(outfile)

def translate(bamfile,sjfile,exonfiles,flank,genome_file,min_junc_reads,outfile,outdir):
	warnings.warn('## Step 2:\ttranslation start')
	if os.path.exists(bamfile) == False:
		sys.exit("[ERROR]: Bam file not exists!\n")
	scriptname = 'translateJunc-star.py' if exonfiles != 'None' else 'translateJunc-star-allSJ.py'
	cmd = "python "+BINDIR + "/" + scriptname + " -o " + outfile + " -l " + str(flank) + " --min-junc-reads=" + str(min_junc_reads) + " -g " + genome_file + " -G " + outdir + " " + bamfile + " " + sjfile + " " + exonfiles
	os.system(cmd)
	return os.path.basename(outfile)

def mergePeps2Database(fastafile,database):
	warnings.warn('## Step 3:\tmerge database')
	os.system("cat "+database+' ' + OUTDIR+"/"+fastafile + ">" + OUTDIR + "/merge_" + fastafile)
	return os.path.basename(OUTDIR + "/merge_" + fastafile)

def databaseSearch(rawdir, database):
	warnings.warn('## Step 4:\tcomet search start')
	COMETOUTDIR = OUTDIR + "/comet_" + re.sub(re.compile("merge_|\.fa$"),"",os.path.basename(database))
	RAWDIR = re.sub(re.compile("/$"),"",rawdir)
	mylogger = logging.getLogger("comet")

	allfiles = os.listdir(RAWDIR)
	rawfiles = [x for x in allfiles if re.search('\.raw|\.mzXML|\.mzML|\.mzML\.gz|\.mzXML\.gz',x,re.I) is not None]
	if not os.path.exists(COMETOUTDIR): os.makedirs(COMETOUTDIR)
	for i in rawfiles[0:]:
		warnings.warn("\t# raw file:\t%s" % i)
		inf = RAWDIR + '/' + i
		outf = COMETOUTDIR + '/' + re.sub(re.compile("\..*$"),"",i)
		cmd = COMETEXE.replace('win64','linux') +" -P"+ COMETPAR+" -D"+database+" -N"+outf+" "+inf
		if re.search(r'raw$',i,re.I):
			cmd = "WINEDEBUG=fixme-all,err-all " + WINE +" "+ COMETEXE +" -P"+ COMETPAR+" -D"+database+" -N"+outf+" "+inf
		if re.search(r'.gz',i):
			tmpdir = OUTDIR + '/tmp'
			if not os.path.exists(tmpdir): os.makedirs(tmpdir)
			inf_no_gz = tmpdir+'/'+re.sub(r'\.gz$','',i)
			os.system('gunzip -c ' + inf + ' > ' + inf_no_gz)
			cmd = COMETEXE.replace('win64','linux') +" -P"+ COMETPAR+" -D"+database+" -N"+outf+" "+inf_no_gz
		warnings.warn("\t"+cmd)
		p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
		for line in p.stdout.readlines():
			mylogger.info(line.rstrip())
		retval = p.wait()
	if os.path.exists(OUTDIR + '/tmp'): os.system('rm -r ' + OUTDIR + '/tmp')
	return os.path.basename(COMETOUTDIR)

def percolatorCrux(cometdir):
	warnings.warn('## Step 5:\tpercolator search start')
	cometfiles = os.listdir(OUTDIR+'/'+cometdir)
	cruxPrecTmp = OUTDIR+'/'+cometdir.replace('comet','cruxPrecTmp')	
	cruxoutdir = OUTDIR+'/'+cometdir.replace('comet','cruxoutput')
	if not os.path.exists(cruxoutdir):
                os.makedirs(cruxoutdir)
	os.system('python ' + BINDIR+'/modifyScanNr2CruxPerc.py'+' '+OUTDIR+'/'+cometdir +' '+cruxPrecTmp)
	os.system(CRUX +' percolator --train-fdr 0.05 --test-fdr 0.05 --overwrite T --output-dir ' + cruxoutdir + ' ' + cruxPrecTmp + '/' + cometdir + '.2cruxprec')
	return cruxoutdir + '/percolator.target.peptides.txt'

def percolator(cometdir):
	warnings.warn('## Step 5:\tpercolator search start')
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

def postPercolatorFilter(fastaname,percolatorfile,exonfille,sjfile, outfile, threadNum=1):
	warnings.warn('## Step 6:\tpercolator filter search start')
	tmpdir = OUTDIR+'/tmp'
	if exonfille != 'None':
		os.system('python ' + BINDIR+'/percolator_triesearch_exon.py  -p %s -c %s -e %s -j %s -t %s -b %s -o %s' % (fastaname, percolatorfile, exonfille, sjfile, tmpdir, BEDTOOLDIR, outfile))
	else:
		basename = os.path.basename(fastaname)
		mergefasta = fastaname.replace(basename,'merge_'+basename)
		os.system('python ' + BINDIR+'/percolator_triesearch.py  -p %s -c %s -t %s -b %s -o %s' % (mergefasta, percolatorfile, tmpdir, BEDTOOLDIR, outfile))
		

def custom_formatwarning(msg, *a):
	# ignore everything except the message
	return str(msg) + '\n'

def getoptions():
	global OUTDIR
	config = ConfigParser.ConfigParser()
        config.readfp(open('config.ini',"rb"))
        CHROMS = config.get('global','CHROMS')
	## read pipeline parameters
	usage = 'usage: %prog -b Aligned.out.sorted.bam -j SJ.tab.out -p proteomicsdir -e HSExonfile/None -d dabatase -g genome_file -o outdir --l 66 --min-junc-reads 2 --trim-RK False --step 1/2/3/4/5/6'
	parser = OptionParser(usage)
	# necessary parameters
	parser.add_option('-b','--bamfile', dest='bamfile', help='bam file from STAR [Default %default]')
	parser.add_option('-j','--sjfile', dest='sjfile', help='SJ.tab.out file from STAR [Default %default]')
	parser.add_option('-p','--rawdir', dest='rawdir', help='proteomics dir (format: raw/mzXML) [Default %default]')
	parser.add_option('-e','--exonfile', dest='exonfile',default='None', help='exons file (bed format) or use None to search all junctions peptide [Default %default]')
	parser.add_option('-d','--database', dest='database',default='data/UP000005640_9606_additional_cdhit1.fasta', help='fasta file to perform the MS search[Default %default]')
	parser.add_option('-o', dest='outdir',default=OUTDIR, help='Output dir filename [Default %default]')	
	# parameters for translation
        parser.add_option('--l', dest='flank', type='int', default=66, help='Extend flanking junction ends by this number of bp [Default %default]')
        parser.add_option('-g', dest='genome_file', default=CHROMS, help='genomic fasta directory (by chromosomes) [Default %default]')
        parser.add_option('--min-junc-reads', dest='min_junc_reads', default=2, type='int', help='Minimum number of reads required spanning the junction [Default %default]')
        parser.add_option('--step', dest='step', default = 1, type='int', help='Start from certain step. 1: start from parse SJ.out.tab. 2: start from translation. 3: start from database merge. 4: start from comet. 5: start from percolator. 6: start from percolatorfilter [Default %default]')
	(options, args) = parser.parse_args()
	# check parameters
	if options.sjfile is None or options.bamfile is None or options.rawdir is None:
		sys.exit("[ERROR] "+parser.get_usage())
	if options.exonfile != 'None' and not os.path.exists(options.exonfile):
			sys.exit("[ERROR] Please input exon file or use '-e None'\n"+parser.get_usage())
	if not os.path.exists(options.database):
			sys.exit("[ERROR] Please input database file '-d database'\n"+parser.get_usage())
	if options.outdir is not None:
		OUTDIR = re.sub(re.compile("/$"),"",options.outdir)

	# mkdir 'outdir' if not exists
	if not os.path.exists(OUTDIR):
		os.makedirs(OUTDIR)
	return options, args

if __name__ == '__main__':
	options,args = getoptions()
	main(options)
