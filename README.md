# proteoseq
##Proteotranscriptome pipeline
RNA-seq and proteomics integrated pipeline for identification of peptide evidence of Alu exon or Human specific exons
<br />


## Test (2017-02-08)
		# download the install files
		git clone https://github.com/Xinglab/proteoseq.git

		# install by running the 'install file'
		./install --homedir /u/home/y/ybwang --install directorToInstall
		Example:
		./install --homedir /u/home/y/ybwang --install ~/scratch/ptp

		# run the pipeline to test
		Usage: ptp.py -b Aligned.out.sorted.bam -j SJ.tab.out -p proteomicsdir -e HSExonfile -o outdir
		Example:
		python ./ptp.py -b GM18486.rna/Aligned.out.sorted.bam -j GM18486.rna/GM18486.rna.SJ -p Yoav_Gilad_proteom/GM18486/ -e data/Ensembl_Alu_25bp_0.5.unique.sorted.bed

## Test (2016-12-08)
		# copy wine，comet to user's home directory
		cd ~
		cp -r /u/home/y/ybwang/wine-1.6.2 ./
		cp -r /u/home/y/ybwang/.wine ./
		cp -r /u/home/y/ybwang/.local ./
		cd ~/.wine
		sed -i 's/\\y\\\\ybwang/\\p\\\\panyang/g' *.reg
		sed -i 's/ybwang/panyang/g' *.reg
		
		cd ~
		cp -r /u/home/y/ybwang/comet ./
		ln -s ~/wine-1.6.2/bin/win64 ~/comet/wine
		
		# test whether wine and comet could run normally, if the output is the usage of comet tool, then success
		WINEDEBUG=fixme-all,err-all ~/comet/wine ~/comet/bin/comet.2015025.win64.exe
		
		# if the command could run normally
		# and if crux tool installed successed，then run pipeline to test
		ln -s /u/home/y/ybwang/nobackup-yxing-PROJECT/ProteoTranscriptomePipeline/ ybwangdir
		python ./ptp.py -b ybwangdir/data/RNA/GM18486.rna/Aligned.out.sorted.bam -j ybwangdir/data/SJ_out/LCLs/GM18486.rna.SJ -p /u/home/y/ybwang/nobackup-yxing/data/Yoav_Gilad_proteom/GM18486 -e ybwangdir/data/Ensembl_Alu_25bp_0.5.unique.sorted.bed -o output

##Dependencies

##Build

##Usage

##Output
