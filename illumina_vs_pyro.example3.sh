
#(c) matteo.ramazzotti@unifi.it - example script for riboFrame
# this script draws the full riboFrame pipeline for the HMP sample SRS017139 form reads download to comparison of 
# abundance profiles between Illumina WGS riboFrame processed data and pyrosequencing data.
# WARNING: serveral Gb of data are taken from the network!!!
# riboTrap.pl, riboMap.pl, hmmsearch must be in path. Please point classifier.jar in appropriate folder according to user's setup

#sample runs at SRA can be found with
#perl SRS2SRR.pl SRS017139
#Sample SRS017139: 4 experiments (25166 25165 21118 21107 )
#ID: 25166 -> SRX: SRX023295 -> SRS: SRS017139 -> SRR: SRR059903 -> Deep WGS Genome Sequencing of HMP clinical samples
#DOWNLOAD: ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR059/SRR059903/SRR059903.sra
#ID: 25165 -> SRX: SRX023294 -> SRS: SRS017139 -> SRR: SRR059902 -> Deep WGS Genome Sequencing of HMP clinical samples
#DOWNLOAD: ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR059/SRR059902/SRR059902.sra
#ID: 21118 -> SRX: SRX019695 -> SRS: SRS017139 -> SRR: SRR042005 -> HMP 16S sequencing on 454
#DOWNLOAD: ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR042/SRR042005/SRR042005.sra
#ID: 21107 -> SRX: SRX019684 -> SRS: SRS017139 -> SRR: SRR041115 -> HMP 16S sequencing on 454
#DOWNLOAD: ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR041/SRR041115/SRR041115.sra


#SRA reports that for WGS "The SRA run(s) below have been pre-filtered by NCBI to remove contaminating human sequence."
#hmpdacc actually have indexed the WGS data, so SRA is not needed for WGS...

mkdir SRS017139

wget -O SRS017139/SRS017139.tar.bz2 http://downloads.hmpdacc.org/data/Illumina/supragingival_plaque/SRS017139.tar.bz2
tar -xjvf SRS017139/SRS017139.tar.bz2

egrep -A1 -e "^@.+/[12]\$" SRS017139/SRS017139.denovo_duplicates_marked.trimmed.1.fastq | tr "@" ">" | grep -v "^--" > SRS017139/ilmn.1.fasta
egrep -A1 -e "^@.+/[12]\$" SRS017139/SRS017139.denovo_duplicates_marked.trimmed.2.fastq | tr "@" ">" | grep -v "^--" > SRS017139/ilmn.2.fasta
egrep -A1 -e "^@.+/[1]\$" SRS017139/SRS017139.denovo_duplicates_marked.trimmed.singleton.fastq | tr "@" ">" | grep -v "^--" >> SRS017139/ilmn.1.fasta
egrep -A1 -e "^@.+/[2]\$" SRS017139/SRS017139.denovo_duplicates_marked.trimmed.singleton.fastq | tr "@" ">" | grep -v "^--" >> SRS017139/ilmn.2.fasta

rm SRS017139/*.fastq
 
#sample SRS017139
hmmsearch -E 0.00001 --domtblout SRS017139/ilmn.1.fwd.bact.ribosomal.table --noali --cpu 2 -o /dev/null hmms/16s_bact_for3.hmm SRS017139/ilmn.1.fasta &
hmmsearch -E 0.00001 --domtblout SRS017139/ilmn.1.rev.bact.ribosomal.table --noali --cpu 2 -o /dev/null hmms/16s_bact_rev3.hmm SRS017139/ilmn.1.fasta &
hmmsearch -E 0.00001 --domtblout SRS017139/ilmn.1.fwd.arch.ribosomal.table --noali --cpu 2 -o /dev/null hmms/16s_arch_for3.hmm SRS017139/ilmn.1.fasta 
hmmsearch -E 0.00001 --domtblout SRS017139/ilmn.1.rev.arch.ribosomal.table --noali --cpu 2 -o /dev/null hmms/16s_arch_rev3.hmm SRS017139/ilmn.1.fasta 
hmmsearch -E 0.00001 --domtblout SRS017139/ilmn.2.fwd.bact.ribosomal.table --noali --cpu 2 -o /dev/null hmms/16s_bact_for3.hmm SRS017139/ilmn.2.fasta &
hmmsearch -E 0.00001 --domtblout SRS017139/ilmn.2.rev.bact.ribosomal.table --noali --cpu 2 -o /dev/null hmms/16s_bact_rev3.hmm SRS017139/ilmn.2.fasta &
hmmsearch -E 0.00001 --domtblout SRS017139/ilmn.2.fwd.arch.ribosomal.table --noali --cpu 2 -o /dev/null hmms/16s_arch_for3.hmm SRS017139/ilmn.2.fasta 
hmmsearch -E 0.00001 --domtblout SRS017139/ilmn.2.rev.arch.ribosomal.table --noali --cpu 2 -o /dev/null hmms/16s_arch_rev3.hmm SRS017139/ilmn.2.fasta 

riboTrap.pl SRS017139/ilmn

java -Xmx1g -jar rdp_classifier-2.3.jar -q SRS017139/ilmn.16S.fasta -o SRS017139/ilmn.16S.rdp

riboMap.pl file=SRS017139/ilmn.16S.rdp var=full conf=0.8 cross=over percmin=1 covplot=0 abuplot=0 out=SRS017139/ilmn_full_count
riboMap.pl file=SRS017139/ilmn.16S.rdp var=all conf=0.8 cross=over percmin=1 covplot=0 abuplot=0 out=SRS017139/ilmn_all_count
riboMap.pl file=SRS017139/ilmn.16S.rdp var=V1-V3 conf=0.8 cross=over percmin=1 covplot=0 abuplot=0 out=SRS017139/ilmn_v1v3_count
riboMap.pl file=SRS017139/ilmn.16S.rdp var=V3-V5 conf=0.8 cross=over percmin=1 covplot=0 abuplot=0 out=SRS017139/ilmn_v3v5_count

#16S data are from SRA searching for SRS017139 (4 entries), then searching in the page the 454 multiplex associated to that sample and then using the filter to identify the correct run 
wget -O SRS017139/v3v5_1.sra ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR042/SRR042005/SRR042005.sra
wget -O SRS017139/v3v5_2.sra ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR041/SRR041115/SRR041115.sra

fastq-dump --fasta -O SRS017139 SRS017139/v3v5_1.sra
fastq-dump --fasta -O SRS017139 SRS017139/v3v5_2.sra

rm SRS017139/*.sra

cat SRS017139/v3v5_1.fasta SRS017139/v3v5_2.fasta > SRS017139/pyro_v3v5.fasta
rm SRS017139/v3v5*

java -Xmx1g -jar rdp_classifier-2.3.jar -q SRS017139/pyro_v3v5.fasta -o SRS017139/pyro_v3v5.rdp

riboMap.pl file=SRS017139/pyro_v3v5.rdp ori=454 var=full conf=0.8 cross=over percmin=1 covplot=0 abuplot=0 out=SRS017139/pyro_v3v5_count

sed 's/xxxxx/SRS017139/g' autoplot.general.R > SRS017139/autoplot.SRS017139.R
nohup R --slave --vanilla < SRS017139/autoplot.SRS017139.R
mv nohup.out SRS017139/results.txt
