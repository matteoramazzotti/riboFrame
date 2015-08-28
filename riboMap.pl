
#!/usr/bin/perl
# (c) matteo.ramazzotti@unifi.it
# part of the riboFrame v1.0 package
$ver = '20150820';
if (!$ARGV[0]) {
print "\n----------- riboMap ----------\n------- riboFrame v1.0 -------\n------- rev.  $ver --------\n(c) matteo.ramazzotti\@unifi.it\n\n";
print "USAGE: riboMap.pl file=infile [conf=0] [var=all] [tol=0] [cross=cont] [abuplot=0] [covplot=0] [percmode=present|all] [percmin=0] out=name\n\n";
print "file     : a .RDPClassified.txt file (mandatory).\n";
print "conf     : a confidence threshold below which the rank is not counted.\n";
print "var      : a region to be targeted (accepted values are e.g. \"V1-V2\" or \"V8,V9\", full or all). Default is all\n";
print "tol      : if unspecified, the boundaries of a variable region are exactly defined\n";
print "         : if an integer, the boundary of a region are expanded by the specified bps.\n";
print "         : if an integer followed by the % sign, the boundary of a region are expanded of the specfied percent bps of the read length.\n";
print "cross    : if set to \"in\", the read must entirely be in a variable region (considering tol)\n";
print "         : if set to \"cont\", the read must contain a variable region (considering tol)\n";
print "         : if set to \"any\", the read can overlap or bein overlapped by a variable region (considering tol)\n";
print "percmode : if set to all the percentage is computed using all reads as total\n";
print "         : if set to present only the actually selected reads are used\n";
print "percmin  : if specified, output will be restricted to ranks higher than the specified % value\n"; 
print "pairmark : if paired reads are present, this will specify the mark used for detecting them (e.g. if read/1 and read/2 >= pairmark=\"/\"\n"; 
print "out      : a basename for the output files. If absent, no output/plotting is written and only dumb statistics is reported on screen.\n";
print "abuplot  : if set to 0, only text files are written as output. No plots.\n";
print "covplot  : if set to 1, a covergae plot is written.\n";
exit;
}

&initialize;

if ($ARGV[0] eq 'sel') {
print "\n----------- riboMap ----------\n------- riboFrame v1.0 -------\n(c) matteo.ramazzotti\@unifi.it\n------------------------------\n\n";
print <<EOX; 
example usage of the var command
                  | ..|V1|....|V2|....|V3|....|V4|....|V5|....|V6|....|V7|....|V8|....|V9|.. |
var=all	          |   ----    ----    ----    ----    ----    ----    ----    ----    ----   | 
var=full          | ------------------------------------------------------------------------ |
var=V1,V5,V7      |   ----                            ----            ----                   |
var=V3-V5         |                   --------------------                                   |
var=100-900	      |    (~) ---------------------------------- (~)                            |
 
var=V3 cross=in   |                   ----                                                   |
var=V3 cross=cont |                 --------                                                 |
var=V3 cross=any  |                ----  ----                                                |

EOX
;

exit;
}
if ($ARGV[0] eq 'vars') {
	foreach $v (0..$#V) {
		print "$V[$v]: $Vs[$v] - $Ve[$v]\n";
	}
	exit;
}


###### EVALUATION OF COMMAND LINE ARGUMENTS ######

# see default value in the initiaize sub at the bottom

foreach(@ARGV) {
	chomp $_;
	$_ =~ s/ //g;
	($what,$is) = split (/=/,$_);
	$file = $is if ($what eq 'file');
	$out = $is if ($what eq 'out');
	$thr = $is if ($what eq 'thr');
	$var = $is if ($what eq 'var');
	$tol = $is if ($what eq 'tol');
	$cross = $is if ($what eq 'cross');
	$abuplot = $is if ($what eq 'abuplot');
	$skip = $is if ($what eq 'skip');
	$percmode = $is if ($what eq 'percmode');
	$percmin = $is if ($what eq 'percmin');
	$covplot = $is if ($what eq 'covplot');
	$pairmark = $is if ($what eq 'pairmark');
	$ori = $is if ($what eq 'ori');
}
$pairmark =~ s/([ \\\|\(\)\[\]\{\}\^\$\*\+\?])/\\$1/g;
$pairmark .= "[12]";

&decode_regions;

##### A FIRST SUMMARY WITH CONFIGURATION ######### 

print STDERR "\nConfiguration\n--------------------------------------------------\n";
print STDERR "INFILE    : $file\n";
print STDERR "REGION    : ", join " ",(sort keys %valid),"\n" if ($var);
print STDERR "CROSS     : $cross\n" if ($cross);
print STDERR "TOLERANCE : $tol\n" if ($tol);
print STDERR "MIN CONF  : $thr\n" if ($thr);
print STDERR "MIN PERC  : $percmin\n" if ($percmin);

$R = `which R`;
print STDERR "R is not in you path, no plot will be produced.\n" if (!$R);
$abuplot = 0 if (!$R);
$covplot = 0 if (!$R);
print STDERR "\n";

open (IN,"$file") or die "No file $file";
print STDERR "Loading ranks from $file...";
$tot = 0;
$invalid = 0;
$valid = 0;
$pair = 0;
%cov = ();
while($line = <IN>) {
	chomp $line;
	next if ($line !~ /^\w/);
	$line =~ s/\"//g;
	$line =~ s/$pairmark// if ($pairmark);
	@tmp = split (/\t+/,$line);
	$st = 0;
	$st = 1 if ($tmp[1] eq '-');
	$name = shift @tmp;
	$name .=":2000:70-98" if ($ori eq '454'); # a fake V1 in case of 16S pyrosequencing reads, that is specified using skip=1 at command line
	@tmp1 = split (/\:/,$name); #  the last two fields of the read name are the length and position of the read
	$pos = pop @tmp1;
	$len = pop @tmp1;
	$region = '';
	$region = 'full';
	$region = localize($pos,$len,$tol); #tries to determine if a read is in variable region (among those under selection, see decode_regions)
	if ($covplot) {
		($st,$en) = split (/\./,$pos);
		for ($i=$st;$i<=$en;$i++) {
			$cov1{$i}++;
			$cov2{$i}++ if ($valid{$region});
		}
	}
	$read = join ":",@tmp1;
	$n = $read.".1" if (!$cnt{$read}); #the pairing status is based on the read name count, not on the name itself ('cause diffent ways exist to label pairs)
	$n = $read.".2" if ($cnt{$read});
	if ($valid{$region}) {
		push(@reads,$read) if (!$cnt{$read}); # reads are collected with name, not with pair information so that pairing is evaluated later
		$exist{$n} = 1;
	}
	$cnt{$read}++;# a read is considered to be paired if its pair is also mapped, no matter the region, because it is mapped on the 16S gene
	for ($l=$st;$l<=$#tmp;$l+=3) { #Archaea	domain	1.0
		next if (!$ranks{$tmp[$l+1]}); # only the 6 main ranks are used
		$all{$tmp[$l+1].":".$tmp[$l]} = 1; # this is just for counting
		if ($valid{$region}) { # this is to save some memory
			$r{$tmp[$l+1]."@".$n} = $tmp[$l];  # rank name
			$c{$tmp[$l+1]."@".$n} = $tmp[$l+2]; #confidence score
			$selected{$tmp[$l+1].":".$tmp[$l]} = 1; # this is just for counting
		}
	}
	$tot++;
	$invalid++ if (!$valid{$region});
	$valid++ if ($valid{$region});
	$pair++ if ($cnt{$read} > 1 && $valid{$region});
}
close IN;
print STDERR "done.\n\n";

print STDERR "Reads selection\n------------------------------------r--------------\n";
print STDERR "Total      : ", $tot,"\n";
print STDERR "Discarded  : ", $invalid,"\n";
print STDERR "Available  : ", $valid,"\n";
print STDERR "Singletons : ", ($valid-$pair),"\n";
print STDERR "Pairs      : ", $pair,"\n\n";

print STDERR "No valid reads found, cannot continue!\n" if ($tot == 0);
exit if ($tot == 0);

########## ESTIMATION OF CONFIDENCES AND ABUNDANCES ###################
%above_c = ();
%tot = ();
%abund = ();
%present = ();
%missing = ();
foreach $rank (@ranks) {
	foreach $name (@reads) {
		READ1:
		if ($exist{$name.".1"}) {
			goto READ2 if ($r{$rank."@".$name.".1"} !~ /\w/); # some genus does not have a full domain -> genus classification, so some rank misses 
			$delta = 0.5 if ($exist{$name.".2"});
			$delta = 1 if (!$exist{$name.".2"});
			$tot{"$rank:".$r{$rank."@".$name.".1"}} += $delta;
			if ($c{$rank."@".$name.".1"} >= $thr) {
				$above_c{$rank.":".$r{$rank."@".$name.".1"}} = 1; # this is just for counting
				$abund{"$rank:".$r{$rank."@".$name.".1"}} += $delta;
				$present{$rank} += $delta;
			}
			else {
				$missing{"$rank:".$r{$rank."@".$name.".1"}} += $delta;
			}
		}
		READ2:	
		if ($exist{$name.".2"}) {
			next if ($r{$rank."@".$name.".2"} !~ /\w/); # some genus does not have a full domain -> genus classification, so some rank misses 
			$delta = 0.5 if ($exist{$name.".1"});
			$delta = 1 if (!$exist{$name.".1"});
			$tot{"$rank:".$r{$rank."@".$name.".2"}} += $delta;
			if ($c{$rank."@".$name.".2"} >= $thr) {
				$above_c{$rank.":".$r{$rank."@".$name.".2"}} = 1; # this is just for counting
				$abund{"$rank:".$r{$rank."@".$name.".2"}} += $delta;
				$present{$rank} += $delta;
			}
			else {
				$missing{"$rank:".$r{$rank."@".$name.".2"}} += $delta;
			}
		} 
	}
}
print STDERR "\n";

########## COMPILATION OF THE ABUNDANCES ###################
%o = ();
%oc = ();
%above_p = ();
foreach $rank (@ranks) {
	$o{$rank} = "Name\tTot\tMis\tCount\tPerc\n"; #this initializes output content 
	foreach $name (sort keys %abund) {
		next if ($name !~ /$rank/);
		$label = $name;
		$label =~ s/$rank://;
		$perc = $abund{$name}/$tot*100 if ($percmode eq 'all');
		$perc = $abund{$name}/$present{$rank}*100 if ($percmode eq 'present');
		next if ($percmin && $perc < $percmin); # if this threshold is not met, the output will not be written 
		$above_p{$rank.":".$label} = 1;
		$missing{$name} = 0 if (!$missing{$name});
		$o{$rank} .= $label."\t".$tot{$name}."\t".$missing{$name}."\t".$abund{$name}."\t".$perc."\n"; # the true output is here
		$oc{$rank}++;
	}
}

&report;
&covplot if($covplot);
print STDERR "\nNo output file specified, no output written.\n" if(!$out);
&write_out if($out);


sub write_out {
	print STDERR "\n";
	foreach $rank (@ranks) {
		print STDERR "$rank: $oc{$rank} entries written to $out.$rank.cnt." if ($oc{$rank});
		print STDERR "$rank: no entries available." if (!$oc{$rank});
		next if (!$oc{$rank});
		open(OUT,">$out.$rank.cnt");
		print OUT $o{$rank};
		close OUT;
		if ($abuplot) {
			print STDERR " Plotting...";
			open(OUT,">Rscript");
			print OUT "data<-read.table(file=\"$out.$rank.cnt\", sep=\"\\t\", header=1)\n";
			print OUT "size<-(dim(data)[1]-1)/3\nif(size < 2) {size<-2}\n";
			print OUT "max<-max(nchar(as.vector(data\$Name)))/1.8\n";
   			print OUT "pdf(file=\"$out.$rank.cnt.pdf\", width=size, height=10, bg=\"white\")\n";
			print OUT "par(mar=c(max,4.1,4.1,2.1))\n";
			print OUT "barplot(data\$Perc, names=data\$Name, las=2)\n";
			print OUT "dev.off()\n";
			close OUT;
			`R --slave --vanilla < Rscript`;
			unlink "Rscript";
			print STDERR "done.";
		}
		print STDERR "\n";
	}
}

sub covplot {
	print STDERR "\nWriting covplot to $out.coverage.pdf...";
	open (OUT,">$out.coverage");
	print OUT "pos\tall\tvar\n";
	foreach $p (1..1515) {
		$cov1{$p} = 0 if (!$cov1{$p}); # full set of reads
		$cov2{$p} = 0 if (!$cov2{$p}); # selected set of reads
		print OUT $p,"\t",$cov1{$p},,"\t",$cov2{$p},"\n";
	}
	close OUT;
	open(OUT,">Rscript");
	print OUT 'varst<-c(69,137,433,576,822,986,1117,1243,1435)',"\n";
	print OUT 'varen<-c(99,242,497,682,879,1043,1173,1294,1465)',"\n";
	print OUT 'labels=c("V1","V2","V3","V4","V5","V6","V7","V8","V9")',"\n";
	print OUT "coverage<-read.table(file=\"$out.coverage\",sep=\"\t\", header=T)\n";
	print OUT "pdf(file=\"$out.coverage.pdf\", width=10, height=5)\n";
	print OUT 'plot(coverage$pos,coverage$all,type="l",xlab="Position", ylab="Coverage", xaxt="n")',"\n";
	print OUT 'lines(coverage$pos,coverage$var,col="red")',"\n";
	print OUT 'abline(v=c(varen,varst),lty=5)',"\n";
	print OUT 'axis(3,at=(varen-varst)/2+varst,labels=labels)',"\n";
	print OUT 'axis(1,at=seq(1,1550,by=100)-1, labels=seq(1,1550,by=100)-1)',"\n";
	print OUT 'rect(varst,min(coverage$all),varen,max(coverage$all), col=rgb(1,0,0,0.3), border=NA)',"\n";
	print OUT 'dev.off()',"\n";
	close OUT;
	`R --slave --vanilla < Rscript`;
	unlink "$out.coverage";
	unlink "Rscript";
	print STDERR "done.\n";
}

print STDERR "\n";

sub localize {
	my $pos = shift;
	my $len = shift;
	my $tol = shift;
	my ($st,$en) = split (/\./,$pos);
	$tol = int($len/100*$tol) if ($tol =~ /%/);
	foreach $ind (0..$#V) {
		$ok = 0;
		$ok = 'in' if ($st >= $Vs[$ind] && $en <= $Ve[$ind]); # the read is exactly contained in one of the variable region
		$ok = 'cont' if ($st <= $Vs[$ind] && $en >= $Ve[$ind]); # the read contains one of the variable region
		$ok = '5over' if ($st >= $Vs[$ind]-$tol && $en >= $Vs[$ind] && $en <= $Ve[$ind]); # the read 5' overhang a variable region (tol defines the degree of overhang)
		$ok = '3over' if ($st <= $Ve[$ind] && $st >= $Vs[$ind] && $en <= $Ve[$ind]+$tol); # the read 5' overhang a variable region (tol defines the degree of overhang)
		$ok = 0 if ($cross eq 'in' && $ok ne 'in');
		$ok = 0 if ($cross eq 'cont' && $ok =~ /over/);
		return $V[$ind] if ($ok); # returns the variable region the read is placed in if criteria are met
	}
}

##### DECODING OF REGION SELECTION ####

sub decode_regions {	
	if ($var eq 'all') { # all variable regions used
		%valid = map {$_ => 1} @V;	
	}
	if ($var eq 'full') { # the full gene is used
		$valid{$var}= 1;
		@V = ($var);
		@Vs =(1);
		@Ve = (1515);
		$Vs{$var} = 1;
		$Ve{$var} = 1515;
	}
	if ($var =~ /V/i && $var !~ /,/ && $var !~ /-/) { # only the specified variable region used
		$valid{$var} = 1;
		@V = ($var);
		@Vs = $Vs{$var};
		@Ve = $Ve{$var};
	}
	if ($var =~ /,/) { # some (comma separated) specified variable regions used
		@tmp = split (/,/,$var);
		%valid = map {$_ => 1} @tmp;	
	}
	if ($var =~ /-/) { # the variable regions spanning from X1 to X2 used 
		$min = 1000000;
		$max = 0;
		@tmp = split (/-/,$var); # e.g. V5-V6
		if ($var =~ /V/i) {
			foreach $r (@tmp) {
				$min = $Vs{$r} if ($Vs{$r} < $min);
				$min = $Ve{$r} if ($Ve{$r} < $min);
				$max = $Vs{$r} if ($Vs{$r} > $max);
				$max = $Ve{$r} if ($Ve{$r} > $max);
			}
			$valid{$var} = 1; # the key here is "V5-V6"
			@V = ($var);
			@Vs = ($min);
			@Ve = ($max);
			$Vs{$var} = $min;
			$Ve{$var} = $max;
		} else {
			$var = 'user';
			$valid{$var} = 1;
			@V = ($var);
			@Vs = ($tmp[0]);
			@Ve = ($tmp[1]);
			$Vs{$var} = $tmp[0];
			$Ve{$var} = $$tmp[1];
		}
	}
}

sub report {
	print STDERR "Ranks selection\n--------------------------------------------------\n\n";
	print STDERR "Rank\tTotal\tSel\tConfOK\tPercOK\n"; 
	foreach $rank (@ranks) {
		$cnta = 0;
		foreach $n (sort keys %all) {
			$cnta++ if ($n =~ /^$rank:/);
		}
		$cntv = 0;
		foreach $n (sort keys %selected) {
			$cntv++ if ($n =~ /^$rank:/);
#			$c1{$n} = 1 if ($n =~ /^class:/);
		}
		$cntc = 0;
		foreach $n (sort keys %above_c) {
			$cntc++ if ($n =~ /^$rank:/);
#			$c2{$n} = 1 if ($n =~ /^class:/);
		}
		$cntp = 0;
		foreach $n (sort keys %above_p) {
			$cntp++ if ($n =~ /^$rank:/);
		}
		print ucfirst($rank),"\t$cnta\t$cntv\t$cntc\t$cntp\n";
	}
#	print STDERR "----\n",(join "\n", sort keys %c1),"\n\n",(join "\n", sort keys %c2),"\n\n";
}

##### INITIALIZE SOME VARIABLE ####

sub initialize { 
	@ranks = qw/domain phylum class order family genus/;
	%ranks = map {$_ => 1} @ranks;	
	@V = qw/V1 V2 V3 V4 V5 V6 V7 V8 V9/;
	@Vs = qw/69 137 433 576 822 986 1117 1243 1435/;
	@Ve = qw/99 242 497 682 879 1043 1173 1294 1465/;
	foreach $v (0..$#V) {
		$Vs{$V[$v]} = $Vs[$v];
		$Ve{$V[$v]} = $Ve[$v];
	}
	#defaults
	$thr = 0;
	$var = 'all';
	$tol = 0;
	$cross = 'over';
	$percmode = 'present';
	$percmin = 0;
	$pairmark = "";
	$abuplot = 1;
}
