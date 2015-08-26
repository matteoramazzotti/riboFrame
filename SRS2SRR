#!/usr/bin/perl
#(c) matteo.ramazzotti@unifi.it - component of the riboFrame project
use LWP::Simple;
if (!$ARGV[1]) { #sra is searched by sample id (SRS) and runs have to be searched in the experiments
    $v = get("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term=$ARGV[0]");
    while ($v =~ m|<Id>(\d+)</Id>|g) {
        push(@ids,$1);
    }
} else { #sra is searched by experiment and the runs associated to the sample id have to be returned
    $v = get("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term=$ARGV[1]");
    $v =~ m|<Id>(\d+)</Id>|;
    @ids = ($1);
}

print "Sample $ARGV[0]: ",scalar @ids, " experiments (",join " ", @ids,")\n";

foreach $id (@ids) {
    $v = get("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id=$id");
    $v =~ s/[\n\r]//g;
    $v =~ m|<EXPERIMENT_PACKAGE_SET>.+?<PRIMARY_ID>(.+?)</PRIMARY_ID>|;
    $SRX = $1;
    $v =~ m|<TITLE>(.+?)</TITLE>|;
    $title = $1;
    print "ID: $id -> SRX: $SRX -> ";
    $SRS = '';
    $SRR = '';
    while($v =~ m|<RUN alias.+?<PRIMARY_ID>(.+?)</PRIMARY_ID>.+?<Pool>.+?<PRIMARY_ID>(.+?)</PRIMARY_ID>.+?</Pool>.+?</RUN>|g) {
        $SRR = $1;
        $SRS = $2;
        print "SRS: $SRS -> SRR: $SRR -> $title\n" if ($SRS eq $ARGV[0]);
        print "DOWNLOAD: ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/",substr($SRR,0,6),"/$SRR/$SRR.sra\n"  if ($SRS eq $ARGV[0]);

    }
    if (!$SRS) {
        while($v =~ m|<RUN .+?>.+?<PRIMARY_ID>(.+?)</PRIMARY_ID>.+?<Pool>.+?<PRIMARY_ID>(.+?)</PRIMARY_ID>.+?</Pool>.+?</RUN>|g) {
            $SRR = $1;
            $SRS = $2;
            print "SRS: $SRS -> SRR: $SRR -> $title\n" if ($SRS eq $ARGV[0]);
            print "DOWNLOAD: ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/",substr($SRR,0,6),"/$SRR/$SRR.sra\n"  if ($SRS eq $ARGV[0]);
        }
    }
}
