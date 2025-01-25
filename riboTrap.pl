#!/usr/bin/perl
# (c) matteo.ramazzotti@unifi.it
# part of the riboFrame v1.0 package

$ver = '20191116';

# WARNINNG: THIS SCRIPT ASSUMES THAT SEQUENCES ARE WRITTEN IN JUST 1 LINE AFTER THE HEADER
if (!$ARGV[0]) {
	print "\n----------- riboTrap ---------\n------- riboFrame v1.0 -------\n------- rev.  $ver --------\n(c) matteo.ramazzotti\@unifi.it\n\n";
	print "USAGE: riboTrap.pl sample1 [sample2] [sample3] ... [sampleN] [nopair]\n\n";
	exit;
}

@samples = @ARGV;
$nopair = 1 if ($samples[$#samples] eq 'nopair');
pop @samples if ($samples[$#samples] eq 'nopair');

foreach $sam (@samples) {
	%r1bf = ();     %r1br = ();     %r2bf = ();     %r2br = (); %rb = (); %flagb = (); %seq = (); %cnta = ();
	%r1af = ();     %r1ar = ();     %r2af = ();     %r2ar = (); %ra = (); %flagb = (); %seq = (); %cntb = ();
	$base1 = "$sam" if ($nopair);
	$base1 = "$sam.1" if (!$nopair);
	foreach $beast (qw/bact arch/) {
		print STDERR "\nParsing predictions of $beast ribosome on $sam\n";
		print STDERR "  Loading pair 1 fwd $beast..." if (!$nopair);
		print STDERR "  Loading fwd $beast..." if ($nopair);
#		print STDERR "$sam/$sam.1.fwd.$beast.ribosomal.table\n";	
		open(IN,"$base1.fwd.$beast.ribosomal.table") or die "no 1.fwd.$beast file";
		while($line = <IN>) {
			next if ($line =~ /#/);
			@tmp = split (/\s+/,$line);
			next if ($tmp[18]-$tmp[17] < 60);
			$tmp[0] =~ s/\/\d$//; # to cope with pair annotation 
			$r1bf{$tmp[0]} = 1 if ($beast eq 'bact');
			$rb{$tmp[0]} = 1  if ($beast eq 'bact');
			$r1af{$tmp[0]} = 1 if ($beast eq 'arch');
			$ra{$tmp[0]} = 1  if ($beast eq 'arch');
			$r1fposb{$tmp[0]} = $tmp[15].".".$tmp[16] if ($beast eq 'bact');
			$r1fposa{$tmp[0]} = $tmp[15].".".$tmp[16] if ($beast eq 'arch');
#			print STDERR $tmp[15].".".$tmp[16]," is $tmp[0] = ",$r1fposb{$tmp[0]}," = ",$r1fposa{$tmp[0]}  if ($tmp[0] eq 'HWUSI-EAS574_103088163:6:26:6229:19857');
		}
		close IN;

		print STDERR "done.\n  Loading pair 1 rev $beast..." if (!$nopair);
		print STDERR "done.\n  Loading rev $beast..." if ($nopair);
		open(IN,"$base1.rev.$beast.ribosomal.table") or die "no 1.rev.$beast file";
		while ($line = <IN>) {
			next if ($line =~ /#/);
			@tmp = split (/\s+/,$line);
			next if ($tmp[18]-$tmp[17] < 60);
			$tmp[0] =~ s/\/\d$//; # to cope with pair annotation 
			$r1br{$tmp[0]} = 1  if ($beast eq 'bact');
			$rb{$tmp[0]} = 1 if ($beast eq 'bact');
			$r1ar{$tmp[0]} = 1 if ($beast eq 'arch');
			$ra{$tmp[0]} = 1  if ($beast eq 'arch');
			$r1rposb{$tmp[0]} = $tmp[15].".".$tmp[16] if ($beast eq 'bact');
			$r1rposa{$tmp[0]} = $tmp[15].".".$tmp[16] if ($beast eq 'arch');
#			print STDERR $tmp[15].".".$tmp[16]," is $tmp[0] = ",$r1rposb{$tmp[0]}," = ",$r1rposa{$tmp[0]}  if ($tmp[0] eq 'HWUSI-EAS574_103088163:6:26:6229:19857');
		}
		close IN;
		print STDERR "done.\n";
		if (!$nopair) {
				print STDERR "  Loading pair 2 fwd $beast...";
				$base2 = "$sam" if ($nopair);
				$base2 = "$sam.2" if (!$nopair);
				open(IN,"$base2.fwd.$beast.ribosomal.table") or die "no 2.fwd.$beast file";
				while($line = <IN>) {
					next if ($line =~ /#/);
					@tmp = split (/\s+/,$line); # tmp[0] is the read name
					next if ($tmp[18]-$tmp[17] < 60);
					$tmp[0] =~ s/\/\d$//; # to cope with pair annotation 
					$r2bf{$tmp[0]} = 1 if ($beast eq 'bact');
					$rb{$tmp[0]} = 1 if ($beast eq 'bact');
					$r2af{$tmp[0]} = 1 if ($beast eq 'arch');
					$ra{$tmp[0]} = 1  if ($beast eq 'arch');
					$r2fposb{$tmp[0]} = $tmp[15].".".$tmp[16] if ($beast eq 'bact');
					$r2fposa{$tmp[0]} = $tmp[15].".".$tmp[16] if ($beast eq 'arch');
					print STDERR $tmp[15].".".$tmp[16]," is $tmp[0] = ",$r2fposb{$tmp[0]}," = ",$r2fposa{$tmp[0]}  if ($tmp[0] eq 'HWUSI-EAS574_103088163:6:26:6229:19857');
				}
				close IN;

				print STDERR "done.\n  Loading pair 2 rev $beast...";
				open(IN,"$base2.rev.$beast.ribosomal.table") or die "no 2.rev.$beast file";
				while($line = <IN>) {
					next if ($line =~ /#/);
					@tmp = split (/\s+/,$line); # tmp[0] is the read name
					next if ($tmp[18]-$tmp[17] < 60);
					$tmp[0] =~ s/\/\d$//; # to cope with pair annotation 
					$r2br{$tmp[0]} = 1 if ($beast eq 'bact');
					$rb{$tmp[0]} = 1 if ($beast eq 'bact');
					$r2ar{$tmp[0]} = 1 if ($beast eq 'arch');
					$ra{$tmp[0]} = 1  if ($beast eq 'arch');
					$r2rposb{$tmp[0]} = $tmp[15].".".$tmp[16] if ($beast eq 'bact');
					$r2rposa{$tmp[0]} = $tmp[15].".".$tmp[16] if ($beast eq 'arch');
					print STDERR $tmp[15].".".$tmp[16]," is $tmp[0] = ",$r2rposb{$tmp[0]}," = ",$r2rposa{$tmp[0]} if ($tmp[0] eq 'HWUSI-EAS574_103088163:6:26:6229:19857');
				}
				close IN;
				print STDERR "done.\n";
			}
		}
		#  pair1--> <--pair2   ==    out
		#    F			F			trash, it is impossible
		#    R			R			trash, it is impossible
		#	 F			R			1f2r
		#    R			F			1r2f
		#    F			o			1f
		#    R			o			1r
		#    o			F			2f
		#    o			R			2r
	print STDERR "\nFlagging bacterial sequences...";
	$err = 0;
	foreach $r (keys %rb) { # see the scheme above
		$err++ if (($r1bf{$r} && $r2bf{$r}) || ($r1br{$r} && $r2br{$r}));
		next if (($r1bf{$r} && $r2bf{$r}) || ($r1br{$r} && $r2br{$r}));
		$flagb{$r} = "1f2r" if ($r1bf{$r} && $r2br{$r});
		$flagb{$r} = "1r2f" if ($r1br{$r} && $r2bf{$r});
		$flagb{$r} = "1f" if ($r1bf{$r} && !$r2br{$r});
		$flagb{$r} = "1r" if ($r1br{$r} && !$r2bf{$r});
		$flagb{$r} = "2f" if (!$r1br{$r} && $r2bf{$r});
		$flagb{$r} = "2r" if (!$r1bf{$r} && $r2br{$r});
#		print STDERR "\n$r $flagb{$r} $r2rposb{$r}\n" if ($r eq 'HWUSI-EAS574_103088163:6:26:6229:19857');
		$cntb{$flagb{$r}}++;
	}
	print STDERR "done.\n";
	print STDERR "   reads: ",scalar keys %rb, "\n";
	print STDERR "   dir/rev errors $err\n";
	print STDERR "   flags: ",scalar keys %flagb,"\n";
	print STDERR "   1f: $cntb{'1f'}, 1r: $cntb{'1r'}\n";
	print STDERR "   2f: $cntb{'2f'}, 2r: $cntb{'2r'}\n";
	print STDERR "   1r2f: $cntb{'1r2f'}, 1f2r: $cntb{'1f2r'}\n";

	print STDERR "\nFlagging archaeal sequences...";
	$err = 0;
	foreach $r (keys %ra) { # see the scheme above
		$err++ if (($r1af{$r} && $r2af{$r}) || ($r1ar{$r} && $r2ar{$r}));
		next if (($r1af{$r} && $r1ar{$r}) || ($r2af{$r} && $r2ar{$r}));
		$flaga{$r} = "1f2r" if ($r1af{$r} && $r2ar{$r});
		$flaga{$r} = "1r2f" if ($r1ar{$r} && $r2af{$r});
		$flaga{$r} = "1f" if ($r1af{$r} && !$r2ar{$r});
		$flaga{$r} = "1r" if ($r1ar{$r} && !$r2af{$r});
		$flaga{$r} = "2f" if (!$r1ar{$r} && $r2af{$r});
		$flaga{$r} = "2r" if (!$r1af{$r} && $r2ar{$r});
		$cnta{$flaga{$r}}++;
	}
	print STDERR "done.\n";
	print STDERR "   reads: ",scalar keys %ra, "\n";
	print STDERR "   dir/rev errors $err\n";
	print STDERR "   flags: ",scalar keys %flaga,"\n";
	print STDERR "   1f: $cnta{'1f'}, 1r: $cnta{'1r'}\n";
	print STDERR "   2f: $cnta{'2f'}, 2r: $cnta{'2r'}\n";
	print STDERR "   1r2f: $cnta{'1r2f'}, 1f2r: $cnta{'1f2r'}\n";

	print STDERR "\nLoading pair1 sequences from $base1.fasta\n";
	$cnt1a = 0; $cnt1b = 0;
	open(IN,"$base1.fasta");
	$t = 0;
	while(<IN>) {
		chomp $_;
		if ($_ =~ />/) {
			$t++;
			$name = $_;
			$name =~ s/>//;
			$name =~ s/ .+//; # this is because hmmsearch cuts fasta headers at the first space as blast...
			$name =~ s/\/1//;
			$name =~ s/ +$//;
			$ok = 0;
			next if (!$flagb{$name} && !$flaga{$name});
#			next if (($flagb{$name} && $flaga{$name})); # the same read is attributed to both archea and bacteria, this is not good...
			$ok = 1;
			next;
		} 
		if ($ok) { 
			if ($flagb{$name} =~ /1f/) {
#				($st,$en) = split (/\./,$r1fposb{$name});
#				$st--;
#				$post = 1515-$st-length($_);
				$seq1b{$name} = $_;
				$len1b{$name} = length($_);
#				$seq1bfull{$name} = "-"x$st.$_."-"x$post;# if ($seqb{$name}); 
				$posb{"1f@".$name} = $r1fposb{$name};
				$position1b{$name} = $r1fposb{$name};
#				print "1f @ $name: ",$posb{"1f@".$name},"\n";
#				<STDIN>;
				$cnt1b++;
				print STDERR "\n",$name.": b1f:".$r1fposb{$name}.":",$position1b{$name},"\n" if ($name eq 'HWUSI-EAS574_103088163:6:26:6229:19857');
			} 
			if ($flagb{$name} =~ /1r/) {
#			   ($st,$en) = split (/\./,$r1rposb{$name});
#				$st--;
#				$post = 1515-$st-length($_);
				$_ = reverse($_);
				$_ =~ tr/ATGC/TACG/;
				$seq1b{$name} .= $_;
				$len1b{$name} = length($_);
#			   	$seq1bfull{$name} = "-"x$st.$_."-"x$post;# if ($seqb{$name}); 
				$posb{"1r@".$name} = $r1rposb{$name};  
				$position1b{$name} = $r1rposb{$name};
				$cnt1b++;
				print STDERR "\n",$name.": b1r:".$r1rposb{$name}.":",$position2b{$name},"\n" if ($name eq 'HWUSI-EAS574_103088163:6:26:6229:19857');
			}
			if ($flaga{$name} =~ /1f/) {
#			   ($st,$en) = split (/\./,$r1fposa{$name});
#				$st--;
#				$post = 1465-$st-length($_);
				$seq1a{$name} = $_;
				$len1a{$name} = length($_);
#				$seq1afull{$name} = "-"x$st.$_."-"x$post;# if ($seqb{$name}); 
				$posa{"1f@".$name} = $r1fposa{$name};  
				$position1a{$name} = $r1fposa{$name};
				$cnt1a++;
				print STDERR "\n",$name.": a1f".$r1fposa{$name}.":",$position1a{$name},"\n" if ($name eq 'HWUSI-EAS574_103088163:6:26:6229:19857');
			} 
			if ($flaga{$name} =~ /1r/) {
#			   ($st,$en) = split (/\./,$r1rposa{$name});
#				$st--;
#				$post = 1465-$st-length($_);
				$_ = reverse($_);
				$_ =~ tr/ATGC/TACG/;
				$seq1a{$name} = $_;
				$len1a{$name} = length($_);
#				$seq1afull{$name} = "-"x$st.$_."-"x$post;# if ($seqb{$name}); 
				$posa{"1r@".$name} = $r1rposa{$name};  
				$position1a{$name} = $r1rposa{$name};
				$cnt1a++;
				print STDERR "\n",$name.": a1r".$r1rposa{$name}.":",$position1a{$name},"\n" if ($name eq 'HWUSI-EAS574_103088163:6:26:6229:19857');
			}
			#print STDERR "\r  $cnt1b + $cnt1a / $t good sequences\n";
		}
	}
	close IN;
	print STDERR "\n  done. $cnt1b bacterial and $cnt1a archeal valid reads.\n";

	if (!$nopair) {
		print STDERR "\nLoading pair2 sequences from $base2.fasta\n";
	#	open(IN,"/media/SATA1/metagenomics/reads/$sam"."2.fasta");
		open(IN,"$base2.fasta");
		$cnt2a = 0; $cnt2b = 0;
		$t = 0;
		while(<IN>) {
			chomp $_;
			if ($_ =~ />/) {
				$t++;
				$name = $_;
				$name =~ s/>//;
				$name =~ s/\/2//;
				$name =~ s/ +$//;
				$ok = 0;
				next if (!$flagb{$name} && !$flaga{$name});
	#			next if (($flagb{$name} && $flaga{$name})); # the same read is attributed to both archea and bacteria, this is not good...
				$ok = 1;
				next;
			}
			if ($ok) {
				if ($flagb{$name} =~ /2f/) {
	#				($st,$en) = split (/\./,$r2fposb{$name});
	#				$st--;
	#				$post = 1515-$st-length($_);
					$seq2b{$name} = $_;# if (!$seqb{$name});
					$len2b{$name} = length($_);
	#				$seq2bfull{$name} = "-"x$st.$_."-"x$post;# if ($seqb{$name}); 
					$posb{"2f@".$name} = $r2fposb{$name};  
					$position2b{$name} = $r2fposb{$name};
					print STDERR "\n",$name.": b2f:".$r2fposb{$name}.":",$position2b{$name},"\n" if ($name eq 'HWUSI-EAS574_103088163:6:26:6229:19857');
					$cnt2b++;
				}
				if ($flagb{$name} =~ /2r/) {
	#				($st,$en) = split (/\./,$r2rposb{$name});
	#				$st--;
	#				$post = 1515-$st-length($_);
					$_ = reverse($_);
					$_ =~ tr/ATGC/TACG/;
					$seq2b{$name} = $_;# if (!$seqb{$name});
					$len2b{$name} = length($_);
	#				$seq2bfull{$name} = "-"x$st.$_."-"x$post;# if ($seqb{$name}); 
	#				$seq2b{$name} .= "NNNNNNNNNN".$_ if ($seqb{$name}); 
					$posb{"2r@".$name} = $r2rposb{$name};  
					$position2b{$name} = $r2rposb{$name};
					print STDERR "\n",$name.": b2r:".$r2rposb{$name}.":",$position2b{$name},"\n" if ($name eq 'HWUSI-EAS574_103088163:6:26:6229:19857');
					$cnt2b++;
				}
				if ($flaga{$name} =~ /2f/) {
	#				($st,$en) = split (/\./,$r2fposa{$name});
	#				$st--;
	#				$post = 1465-$st-length($_);
					$seq2a{$name} = $_;# if (!$seqa{$name});
					$len2a{$name} = length($_);
	#				$seq2afull{$name} = "-"x$st.$_."-"x$post;# if ($seqb{$name}); 
	#				$seq2a{$name} .= "NNNNNNNNNN".$_ if ($seqa{$name}); 
					$posa{"2f@".$name} = $r2fposa{$name};  
					$position2a{$name} = $r2fposa{$name};
					print STDERR "\n",$name.": a2f:".$r2fposa{$name}.":",$position2b{$name},"\n" if ($name eq 'HWUSI-EAS574_103088163:6:26:6229:19857');
					$cnt2a++;
				}
				if ($flaga{$name} =~ /2r/) {
	#				($st,$en) = split (/\./,$r2rposa{$name});
	#				$st--;
	#				$post = 1465-$st-length($_);
					$_ = reverse($_);
					$_ =~ tr/ATGC/TACG/;
					$seq2a{$name} = $_;# if (!$seqa{$name});
					$len2a{$name} = length($_);
	#				$seq2afull{$name} = "-"x$st.$_."-"x$post;# if ($seqb{$name}); 
	#				$seq2a{$name} .= "NNNNNNNNNN".$_ if ($seqa{$name}); 
					$posa{"2r@".$name} = $r2rposa{$name};  
					$position2a{$name} = $r2rposa{$name};
					print STDERR "\n",$name.": a2r:".$r2rposa{$name}.":",$position2b{$name},"\n" if ($name eq 'HWUSI-EAS574_103088163:6:26:6229:19857');
					$cnt2a++;
				}
				#print STDERR "\r  $cnt2b + $cnt2a / $t good sequences";
			}
		}
		close IN;
		print STDERR "\n  done. $cnt2b bacterial and $cnt2a archeal valid reads.\n";
	}
	$goodb = $cnt1b + $cnt2b;
	$gooda = $cnt1a + $cnt2a;

#	$sam = 'tmp';
	open(OUTB, ">$sam.bact.16S.fasta");
	open(OUTA, ">$sam.arch.16S.fasta");
 	open(OUTC, ">$sam.16S.fasta");
#	open(OUTD, ">$sam.16S.full.fasta");

	open(COVB, ">$sam.bact.16S.coverage");
	open(COVA, ">$sam.arch.16S.coverage");
	open(COVC, ">$sam.16S.coverage");

	print STDERR "\nWriting $goodb bacterial sequences to $sam.bact.16S.fasta...";
	%used1 = ();
	%used2 = ();
	foreach $name (keys %seq1b) {
		$namea = $name;
		$namea .= ":".$len1b{$name}.":".$position1b{$name};
		print "\n$name:$namea\n" if (!$position1b{$name});
#		<STDIN> if (!$position1b{$name});
		print OUTB ">",$namea,"\n",$seq1b{$name},"\n";
		print OUTC ">",$namea,"\n",$seq1b{$name},"\n";
#		print OUTD ">",$name,"\n",$seq1bfull{$name},"\n";
		push (@pos,$position1b{$name});
		$used1{$name} = 1;
	}
	foreach $name (keys %seq2b) {
		$namea = $name;
		$namea .= ":".$len2b{$name}.":".$position2b{$name};
		print "\n$name:$namea\n" if (!$position2b{$name});
		<STDIN> if (!$position2b{$name});
		print OUTB ">",$namea,"\n",$seq2b{$name},"\n";
		print OUTC ">",$namea,"\n",$seq2b{$name},"\n";
#		print OUTD ">",$name,"\n",$seq2bfull{$name},"\n";
		push (@pos,$position2b{$name});
		$used2{$name} = 1;
	}
	close OUTB;
	print STDERR "done.\n";

	print STDERR "Writing $gooda archeal sequences to $sam.arch.16S.fasta...";
	foreach $name (keys %seq1a) {
		$namea = $name;
		$namea .= ":".$len1a{$name}.":".$position1a{$name};
		print "\n$name:$namea\n" if (!$position1a{$name});
		<STDIN> if (!$position1a{$name});
		print OUTA ">",$namea,"\n",$seq1a{$name},"\n";
		print OUTC ">",$namea,"\n",$seq1a{$name},"\n" if (!$used1{$name});
#		print OUTD ">",$name,"\n",$seq1afull{$name},"\n"  if (!$used1{$name});
		push (@pos,$position1a{$name}) if (!$used1{$name});
		$used1{$name} = 1;
 
 	}
	foreach $name (keys %seq2a) {
		$namea = $name;
		$namea .= ":".$len2a{$name}.":".$position2a{$name};
		print "\n$name:$namea\n" if (!$position2a{$name});
		<STDIN> if (!$position2a{$name});
		print OUTA ">",$namea,"\n",$seq2a{$name},"\n";
		print OUTC ">",$namea,"\n",$seq2a{$name},"\n" if (!$used2{$name});
#		print OUTD ">",$name,"\n",$seq2afull{$name},"\n"  if (!$used2{$name});
		push (@pos,$position2a{$name}) if (!$used2{$name});
		$used2{$name} = 1;
  	}
	close OUTB;
	close OUTC;
#	close OUTD;
	$merge = (scalar keys %used1) + (scalar keys %used2);
	print STDERR "done.\n\n", scalar keys %used1," pair1  + ", scalar keys %used2, " pair2 = $merge reads also written in merged file $sam.16S.fasta\n\n";

	print STDERR " Writing ribosome coverage table to $sam.coverage...";
	%cnta = ();
	foreach $r (@pos) {
		($st,$en) = split (/\./,$r);
		for ($pos=$st;$pos<=$en;$pos++) {
			$cnta{$pos}++;
		}
	}
	print COVC "pos\tcnt\n";
	foreach $pos (1..1515) {
		$cnta{$pos} = 0 if (!$cnta{$pos});
		print COVC $pos,"\t",$cnta{$pos},"\n";
	}
	print STDERR "done.\n";
	close COVC;

	print STDERR " Writing bacterial ribosome coverage table to $sam.bact.coverage...";
#	$okb = 1;
	%cnt1 = ();
	%cnt2 = ();
	%cnt = ();
	foreach $name (keys %posb) {
#		print "$name: ", $posb{"1f@".$name}," - ",$posb{"1r@".$name}," - ", $posb{"2f@".$name}," - ",$posb{"2r@".$name},"\n";
#		$w = <STDIN> if (!$okb);
#		$okb = 1 if ($w =~ /\w/);
		$name =~ s/^...//;
		if ($posb{"1f@".$name}) {
			($st,$en) = split (/\./,$posb{"1f@".$name});
#			print "1f: $st -> $en\n";
			for ($pos=$st;$pos<=$en;$pos++) {
				$cnt1{$pos}++;
				$cnt{$pos}++;
			}
		}
		if ($posb{"1r@".$name}) {
			($st,$en) = split (/\./,$posb{"1r@".$name});
#			print "1r: $st -> $en\n";
			for ($pos=$st;$pos<=$en;$pos++) {
				$cnt1{$pos}++;
			   $cnt{$pos}++;
 
			}
		}
		if ($posb{"2f@".$name}) {
			($st,$en) = split (/\./,$posb{"2f@".$name});
#			print "2f: $st -> $en\n";
			for ($pos=$st;$pos<=$en;$pos++) {
				$cnt2{$pos}++;
			   $cnt{$pos}++;
 			}
		}
		if ($posb{"2r@".$name}) {
			($st,$en) = split (/\./,$posb{"2r@".$name});
#			print "2r: $st -> $en\n";
			for ($pos=$st;$pos<=$en;$pos++) {
				$cnt2{$pos}++;
			   $cnt{$pos}++;
 			}
		}
	}
	print COVB "pos\tread1\tread2\tboth\n";
	foreach $pos (1..1515) {
		$cnt1{$pos} = 0 if (!$cnt1{$pos});
		$cnt2{$pos} = 0 if (!$cnt2{$pos});
		print COVB $pos,"\t",$cnt1{$pos},"\t",$cnt2{$pos},"\t",$cnt{$pos},"\n";
	}
	print STDERR "done.\n";
	close COVB;

	print STDERR " Writing archeal ribosome coverage table to $sam.arch.coverage...";
#	$oka = 1;
	%cnt1 = ();
	%cnt2 = ();
		%cnt = ();
 	foreach $name (keys %posa) {
#		print "$name: ", $posa{"1f@".$name}," - ",$posa{"1r@".$name}," - ", $posa{"2f@".$name}," - ",$posa{"2r@".$name},"\n";
#		$w = <STDIN> if (!$oka);
#		$oka = 1 if ($w =~ /\w/);
		$name =~ s/^...//;
		if ($posa{"1f@".$name}) {
			($st,$en) = split (/\./,$posa{"1f@".$name});
#			print "1f: $st -> $en\n";
			for ($pos=$st;$pos<=$en;$pos++) {
				$cnt1{$pos}++;
			   $cnt{$pos}++;
 			}
		}
		if ($posa{"1r@".$name}) {
			($st,$en) = split (/\./,$posa{"1r@".$name});
#			print "1r: $st -> $en\n";
			for ($pos=$st;$pos<=$en;$pos++) {
				$cnt1{$pos}++;
			   $cnt{$pos}++;
 			}
		}
		if ($posa{"2f@".$name}) {
			($st,$en) = split (/\./,$posa{"2f@".$name});
#			print "2f: $st -> $en\n";
			for ($pos=$st;$pos<=$en;$pos++) {
				$cnt2{$pos}++;
			   $cnt{$pos}++;
 			}
		}
		if ($posa{"2r@".$name}) {
			($st,$en) = split (/\./,$posa{"2r@".$name});
#			print "2r: $st -> $en\n";
			for ($pos=$st;$pos<=$en;$pos++) {
				$cnt2{$pos}++;
			   $cnt{$pos}++;
 			}
		}
	}
	print COVA "pos\tread1\tread2\tboth\n";
	foreach $pos (1..1465) {
		$cnt1{$pos} = 0 if (!$cnt1{$pos});
		$cnt2{$pos} = 0 if (!$cnt2{$pos});
		print COVA $pos,"\t",$cnt1{$pos},"\t",$cnt2{$pos},"\t",$cnt{$pos},"\n";
	}
	print STDERR "done.\n";

	close OUT;
	close COVA;
	close COVB;
}
