#(c) matteo.ramazzotti@unifi.it - example script for riboFrame
# this script draws the full riboFrame pipeline for the HMP sample SRS011061 form reads download to comparison of 
# abundance profiles between Illumina WGS riboFrame processed data and pyrosequencing data.
# WARNING: serveral Gb of data are taken from the network!!!
# riboTrap.pl, riboMap.pl, hmmsearch must be in path. Please point classifier.jar in appropriate folder according to user's setup

mkdir SRS011061
#sample runs at SRA can be found with
#perl SRS2SRR.pl SRS011061
#Sample SRS011061: 6 experiments (25437 25436 22138 22137 22119 22118 )
#ID: 25437 -> SRX: SRX023566 -> SRS: SRS011061 -> SRR: SRR060371 -> Deep WGS Genome Sequencing of HMP clinical samples
#DOWNLOAD: ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR060/SRR060371/SRR060371.sra
#ID: 25436 -> SRX: SRX023565 -> SRS: SRS011061 -> SRR: SRR060370 -> Deep WGS Genome Sequencing of HMP clinical samples
#DOWNLOAD: ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR060/SRR060370/SRR060370.sra
#ID: 22138 -> SRX: SRX020622 -> SRS: SRS011061 -> SRR: SRR056969 -> Survey of multiple body sites
#DOWNLOAD: ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR056/SRR056969/SRR056969.sra
#ID: 22137 -> SRX: SRX020621 -> SRS: SRS011061 -> SRR: SRR056886 -> Survey of multiple body sites
#DOWNLOAD: ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR056/SRR056886/SRR056886.sra
#ID: 22119 -> SRX: SRX020603 -> SRS: SRS011061 -> SRR: SRR055748 -> Survey of multiple body sites
#DOWNLOAD: ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR055/SRR055748/SRR055748.sra
#ID: 22118 -> SRX: SRX020602 -> SRS: SRS011061 -> SRR: SRR055668 -> Survey of multiple body sites
#DOWNLOAD: ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR055/SRR055668/SRR055668.sra


#SRA reports that for WGS "The SRA run(s) below have been pre-filtered by NCBI to remove contaminating human sequence."
#hmpdacc actually have indexed the WGS data, so SRA is not needed for WGS...
mkdir SRS011061
wget -O SRS011061/SRS011061.tar.bz2 http://downloads.hmpdacc.org/data/Illumina/stool/SRS011061.tar.bz2
tar -xjvf SRS011061/SRS011061.tar.bz2

egrep -A1 -e "^@.+/[12]\$" SRS011061/SRS011061.denovo_duplicates_marked.trimmed.1.fastq | tr "@" ">" | grep -v "^--" > SRS011061/ilmn.1.fasta
egrep -A1 -e "^@.+/[12]\$" SRS011061/SRS011061.denovo_duplicates_marked.trimmed.2.fastq | tr "@" ">" | grep -v "^--" > SRS011061/ilmn.2.fasta
egrep -A1 -e "^@.+/[1]\$" SRS011061/SRS011061.denovo_duplicates_marked.trimmed.singleton.fastq | tr "@" ">" | grep -v "^--" >> SRS011061/ilmn.1.fasta
egrep -A1 -e "^@.+/[2]\$" SRS011061/SRS011061.denovo_duplicates_marked.trimmed.singleton.fastq | tr "@" ">" | grep -v "^--" >> SRS011061/ilmn.2.fasta

rm SRS011061/*.fastq
 
hmmsearch -E 0.00001 --domtblout SRS011061/ilmn.1.fwd.bact.ribosomal.table --noali --cpu 2 -o /dev/null hmms/16s_bact_for3.hmm SRS011061/ilmn.1.fasta
hmmsearch -E 0.00001 --domtblout SRS011061/ilmn.1.rev.bact.ribosomal.table --noali --cpu 2 -o /dev/null hmms/16s_bact_rev3.hmm SRS011061/ilmn.1.fasta
hmmsearch -E 0.00001 --domtblout SRS011061/ilmn.1.fwd.arch.ribosomal.table --noali --cpu 2 -o /dev/null hmms/16s_arch_for3.hmm SRS011061/ilmn.1.fasta 
hmmsearch -E 0.00001 --domtblout SRS011061/ilmn.1.rev.arch.ribosomal.table --noali --cpu 2 -o /dev/null hmms/16s_arch_rev3.hmm SRS011061/ilmn.1.fasta 
hmmsearch -E 0.00001 --domtblout SRS011061/ilmn.2.fwd.bact.ribosomal.table --noali --cpu 2 -o /dev/null hmms/16s_bact_for3.hmm SRS011061/ilmn.2.fasta
hmmsearch -E 0.00001 --domtblout SRS011061/ilmn.2.rev.bact.ribosomal.table --noali --cpu 2 -o /dev/null hmms/16s_bact_rev3.hmm SRS011061/ilmn.2.fasta
hmmsearch -E 0.00001 --domtblout SRS011061/ilmn.2.fwd.arch.ribosomal.table --noali --cpu 2 -o /dev/null hmms/16s_arch_for3.hmm SRS011061/ilmn.2.fasta 
hmmsearch -E 0.00001 --domtblout SRS011061/ilmn.2.rev.arch.ribosomal.table --noali --cpu 2 -o /dev/null hmms/16s_arch_rev3.hmm SRS011061/ilmn.2.fasta 

riboTrap.pl SRS011061/ilmn

java -Xmx1g -jar classifier.jar -q SRS011061/ilmn.16S.fasta -o SRS011061/ilmn.16S.rdp

riboMap.pl file=SRS011061/ilmn.16S.rdp var=full conf=0.8 cross=over percmin=1 covplot=0 abuplot=0 out=SRS011061/ilmn_full_count
riboMap.pl file=SRS011061/ilmn.16S.rdp var=all conf=0.8 cross=over percmin=1 covplot=0 abuplot=0 out=SRS011061/ilmn_all_count
riboMap.pl file=SRS011061/ilmn.16S.rdp var=V1-V3 conf=0.8 cross=over percmin=1 covplot=0 abuplot=0 out=SRS011061/ilmn_v1v3_count
riboMap.pl file=SRS011061/ilmn.16S.rdp var=V3-V5 conf=0.8 cross=over percmin=1 covplot=0 abuplot=0 out=SRS011061/ilmn_v3v5_count

#16S data are from SRA searching for SRS011061 (4 entries), then searching in the page the 454 multiplex associated to that sample and then using the filter to identify the corect run 
wget -O SRS011061/v1v3_1.sra ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR056/SRR056886/SRR056886.sra
wget -O SRS011061/v1v3_2.sra ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR045/SRR055748/SRR055748.sra

wget -O SRS011061/v3v5_1.sra ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR056/SRR056969/SRR056969.sra
wget -O SRS011061/v3v5_2.sra ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR055/SRR055668/SRR055668.sra

fastq-dump --fasta -O SRS011061 SRS011061/v1v3_1.sra
fastq-dump --fasta -O SRS011061 SRS011061/v1v3_2.sra

fastq-dump --fasta -O SRS011061 SRS011061/v3v5_1.sra
fastq-dump --fasta -O SRS011061 SRS011061/v3v5_2.sra

rm SRS011061/*.sra

cat SRS011061/v1v3_1.fasta SRS011061/v1v3_2.fasta > SRS011061/pyro_v1v3.fasta
cat SRS011061/v3v5_1.fasta SRS011061/v3v5_2.fasta > SRS011061/pyro_v3v5.fasta
rm SRS011061/v3v5* SRS011061/v1v3*

java -Xmx1g -jar classifier.jar -q SRS011061/pyro_v1v3.fasta -o SRS011061/pyro_v1v3.rdp
java -Xmx1g -jar classifier.jar -q SRS011061/pyro_v3v5.fasta -o SRS011061/pyro_v3v5.rdp
rm SRS011061/*.fasta

riboMap.pl file=SRS011061/pyro_v1v3.rdp ori=454 var=full conf=0.8 cross=over percmin=0.5 covplot=0 abuplot=0 out=SRS011061/pyro_v1v3_count
riboMap.pl file=SRS011061/pyro_v3v5.rdp ori=454 var=full conf=0.8 cross=over percmin=0.5 covplot=0 abuplot=0 out=SRS011061/pyro_v3v5_count

sed 's/xxxxx/SRS011061/g' compare.general.R > SRS011061/compare.SRS011061.R
nohup R --slave --vanilla < SRS011061/compare.SRS011061.R
mv nohup.out SRS011061/results.txt

