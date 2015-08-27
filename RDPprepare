#!/usr/bin/perl
#(c) matteo.ramazzotti@unifi.it - part of riboFrame
#usage: zcat current_Bacteria_unaligned.gb.gz current_Archaea_unaligned.gb.gz | perl RDPprepare.pl
#the gz files must have been present on the disk (downloaded from RDP prokect in unaligned genbank format)
#if not, please download from the two following
#https://rdp.cme.msu.edu/download/current_Bacteria_unaligned.gb.gz
#https://rdp.cme.msu.edu/download/current_Archaea_unaligned.gb.gz

print STDERR "\n - Reading RDP database -\n\n";
open (DB,">RDPgb-fasta.txt");
$start = 0;
$stop = 1;
$miss = 0;
$avail = 0;
while(<STDIN>) {
	$start = 1 if ($_ =~ /^LOCUS/);
	$start = 0 if ($_ =~ /^\/\// );
	$stop = 0 if ($_ =~ /^LOCUS/);
	$stop = 1 if ($_ =~ /^\/\//);
	$rep = '' if ($_ =~ /^LOCUS/);
	&analyze if ($stop); 
	if ($start) {
		$rep .= $_;
	}
}
print STDERR "\n\n - Creating RDP index -\n";
open (DESC,">RDPgb-index.txt");
foreach (keys %name) {
	print DESC "GEN: $_\n";
	print DESC "LIN: ",$lin{$_},"\n";
	print DESC "SPE: ",$name{$_},"\n";
	print DESC "IDS: ",$id{$_},"\n";
}
close DESC;
close DB;


sub analyze {
	### USED BY INDEX TO PARSE RDPgb ###
	$tot++;
	undef $id;
	undef $name;
	undef $lineage;
	undef $taxid;
	undef $seq;
	$rep =~ /SOURCE      (.+?)\n/;
	$name = $1;
	$perc = sprintf("%02d",int($avail/$tot*100));
	$rep =~ /LOCUS\s+(.+?)\s/;
	$id = $1;
	$rep =~ /ORGANISM.+?\n(.+?)REFERENCE/s;
	$lineage = $1;
	$lineage =~ s/ {2,}//g;
	$lineage =~ s/\n//g;
	$lineage =~ s/\t//g;
	$lineage =~ s/\"//g;
	$lineage =~ s/\.$//;
	return if ($lineage && $lineage =~ /unclassified_bacteria/i);
	return if ($name && $name =~ /uncultured/i);
#	return if ($lineage && $lineage =~ /archaea/i);
	$avail++;
	$rep =~ /taxon:(.+?)\"/;
	$taxid = $1;
	$rep =~ /ORIGIN\n(.+)\s/s;
	$seq = $1;
	$seq =~ s/[\d\s]//g;
	$seq = uc($seq);
	if (!$id || !$name || !$lineage || !$taxid || !$seq || $seq =~ /n/i) {
		$miss++;
		return;
	}
	@tmp = split(/;/,$lineage);
	$tmp[$#tmp] =~ s/^ //;
	$name{$tmp[$#tmp]} .= $name."|";
	$id{$tmp[$#tmp]} .= $id."|";
	$lin{$tmp[$#tmp]} = $lineage;
	print DB ">$id\n$seq\n"; 
	print STDERR "\r   Collecting sequence $avail of $tot ($perc%), $miss errors "; 
}
