#! /u/home/y/ybwang/python

from optparse import OptionParser
import logging
import os,sys,random,datetime,warnings,re
import subprocess

def main():
	usage = 'usage: %prog --homedir /u/home/y/ybwang --install installdir'
	parser = OptionParser(usage)
	parser.add_option('--homedir', dest='homedir', help='HOME directory [Default %default]')
	parser.add_option('--install', dest='install', help='INSTALL directory [Default %default]')
	(options, args) = parser.parse_args()
	# check parameters
	if options.homedir is None or options.install is None:
		sys.exit("[ERROR] "+parser.get_usage())
	# warning information	
	warnings.formatwarning = custom_formatwarning
	logging.basicConfig(filename='pipeline.log', level=logging.INFO)

	## star install
	# 1. copy wine directory to home directory
	warnings.warn('[INSTALL] copy wine dir to home directory: ' + options.homedir)
	HOMEDIR = re.sub(r'\/$','',options.homedir)
	INSTALLDIR = re.sub(r'\/$','',options.install)
	os.system('tar -zxvf wine_dir.tgz -C ' + HOMEDIR)
	homedir1 = HOMEDIR.replace('/','\\')
	cmd1 = 'sed -i s\'/homedirXXXXXA/' + homedir1 +'/g\' ' + HOMEDIR + '/.wine/*.reg'
	cmd2 = 'sed -i s\'/homedirXXXXXB/' + options.homedir +'/g\' ' + HOMEDIR + '/.wine/*.reg'
	os.system(cmd1)
	os.system(cmd2)
	
	# 2. copy files to install directory
	if not os.path.exists(INSTALLDIR):
		os.makedirs(INSTALLDIR)
	warnings.warn('[INSTALL] copy files to install dir: ' + options.install)
	os.system('cp -rf proteoseq/* ' + INSTALLDIR)
	
	# 3. ln wine
	warnings.warn('[INSTALL] ln wine to install dir')
	os.system('ln -s ' + HOMEDIR + '/wine-1.6.2/bin/wine64 ' + INSTALLDIR + '/bin/wine')
	
	# 4. change config.ini file
	warnings.warn('[INSTALL] generate configure file')
	with open(INSTALLDIR + '/config.ini','w') as fout:
		fout.write('[global]\n')
		fout.write('BINDIR = ' + INSTALLDIR + '/bin\n')
		fout.write('CHROMS = ' + '/u/home/f/frankwoe/nobackup/hg19/hg19_by_chrom/\n')
		fout.write('WINE = '+ INSTALLDIR + '/bin/wine\n')
		fout.write('COMETEXE = ' + INSTALLDIR + '/bin/comet/comet.2015025.win64.exe\n')
		fout.write('COMETPAR = ' + INSTALLDIR + '/bin/comet/comet.params.high-low\n')
		fout.write('CRUX = ' + 	INSTALLDIR + '/bin/crux\n')
		fout.write('BEDTOOLDIR = ' + INSTALLDIR + '/bin/bedtools/\n')
		fout.write('PERCOLATOR = ' + INSTALLDIR + '/bin/percolator-2.08/bin/percolator')

def custom_formatwarning(msg, *a):
	return str(msg) + '\n'

if __name__ == '__main__':
	main()
	
	


