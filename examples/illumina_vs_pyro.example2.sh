
#(c) matteo.ramazzotti@unifi.it - example script for riboFrame
# this script draws the full riboFrame pipeline for the HMP sample SRS011084 form reads download to comparison of 
# abundance profiles between Illumina WGS riboFrame processed data and pyrosequencing data.
# WARNING: serveral Gb of data are taken from the network!!!
# riboTrap.pl, riboMap.pl, hmmsearch must be in path. Please point classifier.jar in appropriate folder according to user's setup

#sample runs at SRA can be found with
#perl SRS2SRR.pl SRS011084
#Sample SRS011084: 8 experiments (26989 26988 26790 26789 22138 22137 22119 22118 )
#ID: 26989 -> SRX: SRX024715 -> SRS: SRS011084 -> SRR: SRR062103 -> Deep WGS Genome Sequencing of HMP clinical samples
#DOWNLOAD: ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR062/SRR062103/SRR062103.sra
#ID: 26988 -> SRX: SRX024714 -> SRS: SRS011084 -> SRR: SRR062102 -> Deep WGS Genome Sequencing of HMP clinical samples
#DOWNLOAD: ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR062/SRR062102/SRR062102.sra
#ID: 26790 -> SRX: SRX024516 -> SRS: SRS011084 -> SRR: SRR061904 -> Deep WGS Genome Sequencing of HMP clinical samples
#DOWNLOAD: ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR061/SRR061904/SRR061904.sra
#ID: 26789 -> SRX: SRX024515 -> SRS: SRS011084 -> SRR: SRR061903 -> Deep WGS Genome Sequencing of HMP clinical samples
#DOWNLOAD: ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR061/SRR061903/SRR061903.sra
#ID: 22138 -> SRX: SRX020622 -> SRS: SRS011084 -> SRR: SRR057013 -> Survey of multiple body sites
#DOWNLOAD: ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR057/SRR057013/SRR057013.sra
#ID: 22137 -> SRX: SRX020621 -> SRS: SRS011084 -> SRR: SRR056932 -> Survey of multiple body sites
#DOWNLOAD: ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR056/SRR056932/SRR056932.sra
#ID: 22119 -> SRX: SRX020603 -> SRS: SRS011084 -> SRR: SRR055794 -> Survey of multiple body sites
#DOWNLOAD: ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR055/SRR055794/SRR055794.sra
#ID: 22118 -> SRX: SRX020602 -> SRS: SRS011084 -> SRR: SRR055711 -> Survey of multiple body sites
#DOWNLOAD: ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR055/SRR055711/SRR055711.sra

#SRA reports that for WGS "The SRA run(s) below have been pre-filtered by NCBI to remove contaminating human sequence."
#hmpdacc actually have indexed the WGS data, so SRA is not needed for WGS...

mkdir SRS011084

wget -O SRS011084/SRS011084.tar.bz2 http://downloads.hmpdacc.org/data/Illumina/stool/SRS011084.tar.bz2
tar -xjvf SRS011084/SRS011084.tar.bz2

egrep -A1 -e "^@.+/[12]\$" SRS011084/SRS011084.denovo_duplicates_marked.trimmed.1.fastq | tr "@" ">" | grep -v "^--" > SRS011084/ilmn.1.fasta
egrep -A1 -e "^@.+/[12]\$" SRS011084/SRS011084.denovo_duplicates_marked.trimmed.2.fastq | tr "@" ">" | grep -v "^--" > SRS011084/ilmn.2.fasta
egrep -A1 -e "^@.+/[1]\$" SRS011084/SRS011084.denovo_duplicates_marked.trimmed.singleton.fastq | tr "@" ">" | grep -v "^--" >> SRS011084/ilmn.1.fasta
egrep -A1 -e "^@.+/[2]\$" SRS011084/SRS011084.denovo_duplicates_marked.trimmed.singleton.fastq | tr "@" ">" | grep -v "^--" >> SRS011084/ilmn.2.fasta

rm SRS011084/*.fastq
 
#sample SRS011084
hmmsearch -E 0.00001 --domtblout SRS011084/ilmn.1.fwd.bact.ribosomal.table --noali --cpu 2 -o /dev/null hmms/16s_bact_for3.hmm SRS011084/ilmn.1.fasta &
hmmsearch -E 0.00001 --domtblout SRS011084/ilmn.1.rev.bact.ribosomal.table --noali --cpu 2 -o /dev/null hmms/16s_bact_rev3.hmm SRS011084/ilmn.1.fasta &
hmmsearch -E 0.00001 --domtblout SRS011084/ilmn.1.fwd.arch.ribosomal.table --noali --cpu 2 -o /dev/null hmms/16s_arch_for3.hmm SRS011084/ilmn.1.fasta &
hmmsearch -E 0.00001 --domtblout SRS011084/ilmn.1.rev.arch.ribosomal.table --noali --cpu 2 -o /dev/null hmms/16s_arch_rev3.hmm SRS011084/ilmn.1.fasta &
hmmsearch -E 0.00001 --domtblout SRS011084/ilmn.2.fwd.bact.ribosomal.table --noali --cpu 2 -o /dev/null hmms/16s_bact_for3.hmm SRS011084/ilmn.2.fasta &
hmmsearch -E 0.00001 --domtblout SRS011084/ilmn.2.rev.bact.ribosomal.table --noali --cpu 2 -o /dev/null hmms/16s_bact_rev3.hmm SRS011084/ilmn.2.fasta &
hmmsearch -E 0.00001 --domtblout SRS011084/ilmn.2.fwd.arch.ribosomal.table --noali --cpu 2 -o /dev/null hmms/16s_arch_for3.hmm SRS011084/ilmn.2.fasta &
hmmsearch -E 0.00001 --domtblout SRS011084/ilmn.2.rev.arch.ribosomal.table --noali --cpu 2 -o /dev/null hmms/16s_arch_rev3.hmm SRS011084/ilmn.2.fasta &

riboTrap.pl SRS011084/ilmn

java -Xmx1g -jar classifier.jar -q SRS011084/ilmn.16S.fasta -o SRS011084/ilmn.16S.rdp

riboMap.pl file=SRS011084/ilmn.16S.rdp var=full conf=0.8 cross=over percmin=1 covplot=0 abuplot=0 out=SRS011084/ilmn_full_count
riboMap.pl file=SRS011084/ilmn.16S.rdp var=all conf=0.8 cross=over percmin=1 covplot=0 abuplot=0 out=SRS011084/ilmn_all_count
riboMap.pl file=SRS011084/ilmn.16S.rdp var=V1-V3 conf=0.8 cross=over percmin=1 covplot=0 abuplot=0 out=SRS011084/ilmn_v1v3_count
riboMap.pl file=SRS011084/ilmn.16S.rdp var=V3-V5 conf=0.8 cross=over percmin=1 covplot=0a abuplot=0 out=SRS011084/ilmn_v3v5_count

#16S data are from SRA searching for SRS011084 (4 entries), then searching in the page the 454 multiplex associated to that sample and then using the filter to identify the corect run 
wget -O SRS011084/v1v3_1.sra ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR056/SRR056932/SRR056932.sra 	
wget -O SRS011084/v1v3_2.sra ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR055/SRR055794/SRR055794.sra

wget -O SRS011084/v3v5_1.sra ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR055/SRR055711/SRR055711.sra
wget -O SRS011084/v3v5_2.sra ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR057/SRR057013/SRR057013.sra

fastq-dump --fasta -O SRS011084 SRS011084/v1v3_1.sra
fastq-dump --fasta -O SRS011084 SRS011084/v1v3_2.sra

fastq-dump --fasta -O SRS011084 SRS011084/v3v5_1.sra
fastq-dump --fasta -O SRS011084 SRS011084/v3v5_2.sra

rm SRS011084/*.sra

cat SRS011084/v3v5_1.fasta SRS011084/v3v5_2.fasta > SRS011084/pyro_v3v5.fasta
cat SRS011084/v1v3_1.fasta SRS011084/v1v3_2.fasta > SRS011084/pyro_v1v3.fasta
rm SRS011061/v3v5* SRS011061/v1v3*

java -Xmx1g -jar classifier.jar -q SRS011084/pyro_v1v3.fasta -o SRS011084/pyro_v1v3.rdp
java -Xmx1g -jar classifier.jar -q SRS011084/pyro_v3v5.fasta -o SRS011084/pyro_v3v5.rdp

riboMap.pl file=SRS011084/pyro_v1v3.rdp ori=454 var=full conf=0.8 cross=over percmin=1 covplot=0 abuplot=0 out=SRS011084/pyro_v1v3_count
riboMap.pl file=SRS011084/pyro_v3v5.rdp ori=454 var=full conf=0.8 cross=over percmin=1 covplot=0 abuplot=0 out=SRS011084/pyro_v3v5_count

sed 's/xxxxx/SRS011084/g' autoplot.general.R > SRS011084/autoplot.SRS011084.R
nohup R --slave --vanilla < SRS011084/autoplot.SRS011084.R
mv nohup.out SRS011084/results.txt
