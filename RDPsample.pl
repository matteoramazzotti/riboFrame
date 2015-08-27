#!/usr/bin/perl
#(c) matteo.ramazzotti@unifi.it - part of riboFrame

if ($ARGV[0] eq 'h' || $ARGV[0] eq 'help' || ! $ARGV[0]) {
	print "USAGE: RDPsampler.pl length=100 size=all [sel=] [excl=] [16Smin=]\n\n";
	print "       lenght: length of the read\n";
	print "       size  : number of species (reads) per genus\n";
	print "       sel   : list of genera (or a file)\n";
	print "       excl  : a string that remove genera if matched\n";
	print "       incl  : a string that remove all non-matching genera (case insensitive) if matched\n";
	print "       16Smin: minimal length of the original 16S gene\n\n";
	exit;
}
$len16S = 0;
$length = 100;
$size = 'all';

foreach(@ARGV) {
	($what,$is) = split (/=/,$_);
	$length = $is if ($what eq 'length');
	$size = $is if ($what eq 'size');
	$sel = $is if ($what eq 'sel');
	$incl = $is if ($what eq 'incl');
	$excl = $is if ($what eq 'excl');
	$len16S = $is if ($what eq '16Smin');
}

if ($sel) {
	if (-e $sel) {
		open(IN, $sel);
		while(<IN>) {
			chomp;
			push(@sel,$_);
		}
		close IN;
	} else {
		@sel = split (/ /,$sel);
	}
	%sel = map {$_ => 1} @sel;
}
print STDERR "Loading index...";
open(IN,"/home/matteo/Desktop/riboFrame/16Ssampling/rdp/RDPgb-index.txt") or die "No index available";
while($line = <IN>) {
	chomp $line;
	$line =~ s/[\n\r]//g;
	if ($line =~ /GEN: (.+)/) {
		$name = $1;
		$name =~ s/[ \.:]/_/g;
	}
	if ($line =~ /SPE: (.+)/) {
		$spe = $1;
		@{$spen{$name}} = split (/\|/,$spe);
	}
	if ($line =~ /IDS: (.+)/) {
		$ids = $1;
		@{$spei{$name}} = split (/\|/,$ids);
#		print "TOT: ", $#{$spei{$name}},"\n";
		foreach $i (0..$#{$spei{$name}}) {
			$spen{$name}[$i] =~ s/[ \.:]/_/g;
			$name{$spei{$name}[$i]} = $spen{$name}[$i];
			$genus{$spei{$name}[$i]} = $name;
		}
	}
}
close IN;
print STDERR "done.\nLoading sequences...";
#print scalar keys %spen," - ",scalar keys %spei," - ",scalar keys %name,"\n";

open(IN,"/home/matteo/Desktop/riboFrame/16Ssampling/rdp/RDPgb-fasta.txt") or die "No sequence";
while($line = <IN>) {
	chomp $line;
	$line =~ s/[\n\r]//g;
	if ($line =~ /^>/) {
		$line =~ s/>//;
		$name = $line;
	} else {
		$seq{$name} .= $line;
	}
}
close IN;

$| = 1;
@genera =sort keys %spei if !@sel;
@genera = @sel if @sel;
$cnt = 0;
print STDERR "done.\nExtracting reads...";
foreach $gen (@genera) {
	$genus_all{$gen} = 1;
	next if ($excl && $gen =~ /$excl/i);
	next if ($incl && $gen !~ /$incl/i);
	$genus_ok1{$gen} = 1;
	#filters out short 16S genes if requested
	if ($len16S) {
#		print STDERR $gen, scalar @{$spei{$gen}}," -> ";
		@ok = ();
		foreach $s (@{$spei{$gen}}) {
			$species_all{$s} = 1;
			next if (length($seq{$s}) < $len16S);
			push(@ok,$s);
			$species_ok{$s} = 1;
		}
#		print STDERR scalar @ok,"\n";
		next if (scalar @ok == 0);
		$genus_ok2{$gen} = 1;
		@{$spei{$gen}} = @ok;
	}
	#select species among current genus
	@list = ();
	if ($size eq 'all') {
		@list = @{$spei{$gen}};
	} 
	else {
		$i = 0;
		while ($i < $size) {
			$r = int(rand(scalar @{$spei{$gen}}));
#			print STDERR "$i) $gen -> $spei{$gen}[$r] -> ",length($seq{$spei{$gen}[$r]})," > $len16S\n"; 
			push(@list,$spei{$gen}[$r]);
			$i++;
		}
	}
	foreach $id (@list) {
#		next if ($len16S && length($seq{$id}) < $len16S);
		$st = int(rand(length($seq{$id})-$length));
		print ">",$genus{$id}."@".$name{$id},".",length($seq{$id}),".$length.$st:",($st+$length),"\n";
		print substr($seq{$id},$st,$length),"\n";
		$cnt++;
#		print STDERR "\r$gen: creating read $cnt";
	}
}
print STDERR "done.\n\n";
print STDERR "GENERA:\n  TOTAL: ", scalar keys %genus_all,"\n  NAME OK:", scalar keys %genus_ok1,"\n  LENGTH OK:", scalar keys %genus_ok2,"\n\n";
print STDERR "SPECIES:\n  TOTAL: ", scalar keys %species_all,"\n  LENGTH OK: ",scalar keys %species_ok,"\n\n";
print STDERR "READS:\n  TOTAL: $cnt\n  AVG x GENUS: ",int($cnt / scalar keys %genus_ok2),"\n\n";
 
print STDERR "\n";
