#! /usr/bin/perl

use strict;
use Getopt::Long;

my $sga = "";
my $freebayes = "";
my $min_quality = 30;

GetOptions("min-quality=i" => \$min_quality);

my @files = @ARGV;
my @callers = ("sga", "freebayes");

die("Exactly two VCF files must be provided") unless scalar(@files) == 2;

my $vcfHash;

foreach my $f (@files) {
    # What center are these calls from?
    my $caller;
    my $caller_id = -1;
    for(my $i = 0; $i < scalar(@callers); $i++) {

        if(index($f, $callers[$i]) != -1) {
            $caller_id = $i;
            last;
        }
    }

    die("Caller not found") if($caller_id == -1);
    loadVCF($f, $caller_id, $callers[$caller_id] eq "freebayes");
}

# print the vcf header
print "##fileformat=VCFv4.1\n";
print "##fileDate=20140116\n";
print "##source=union calls from passed freebayes-snvs and all SGA calls\n";
print "##reference=/u/jsimpson/simpsonlab/data/references/hs37d5.fa\n";
print "##FILTER=<ID=NormalEvidence,Description=\"There is evidence for the SGA call in the aligned reads.\">\n";
print "##FILTER=<ID=dbSNP,Description=\"The variant exists in dbSNP.\">\n";
print "##FILTER=<ID=LowComplexity,Description=\"The variant falls in a region of low complexity sequence.\">\n";
print "##FILTER=<ID=Homopolymer,Description=\"The variant is in a homopolymer is length at least 7.\">\n";
print "##FILTER=<ID=StrandBias,Description=\"The variant is not well represented by both sequencing strands.\">\n";
print "##FILTER=<ID=NoAltEvidence,Description=\"The alignments do not have evidence of assembled variant (SNVs only).\">\n";
print "##FILTER=<ID=LowAlleleFreq,Description=\"The variant is present in a low proportion of reads.\">\n";
print "##FILTER=<ID=LowQuality,Description=\"The quality scores of the assembly variants are low.\">\n";
print "##FILTER=<ID=LowNormalDepth,Description=\"The control sample has low depth in this region.\">\n";
print "##INFO=<ID=CALLER,Number=.,Type=String,Description=\"Which calling program(s) found this variant.\">\n";
print "##INFO=<ID=VT,Number=1,Type=String,Description=\"Whether this variant is SOMATIC or GERMLINE\">\n";
print "##INFO=<ID=TYPE,Number=A,Type=String,Description=\"The type of allele, either snp, mnp, ins, del, or complex.\">\n";
print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
print join("\t", "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "NORMAL", "TUMOR") . "\n";

foreach my $k (keys %{$vcfHash})
{
    my ($chr, $pos, $ref, $alt) = split(';', $k);
    my $passed_in_any = 0;
    my @called;
    my %status;

    for(my $i = 0; $i < scalar(@callers); $i++) {
        if($vcfHash->{$k}->[$i] ne "") {
            push @called, $callers[$i];
            my @f = split(' ', $vcfHash->{$k}->[$i]);
            
            $status{$callers[$i]} = $f[6];
            $passed_in_any = 1 if($f[6] eq "PASS" || $f[6] eq ".");
        }
    }


    my $out_status = "FAIL";

    # Pass the variant if either set passed it, otherwise use the SGA status
    if($passed_in_any) {
        $out_status = "PASS";
    } elsif(defined($status{"sga"})) {
        $out_status = $status{"sga"};
    }
    my $info = "CALLER=" . join(",", @called) . ";VT=SOMATIC";
    print join("\t", ($chr, $pos, ".", $ref, $alt, ".", $out_status, $info, "GT", "0/0", "0/1")) . "\n";
}

sub loadVCF
{
    my($filename, $idx, $is_freebayes) = @_;

    if ($filename =~ /\.gz$/) {
        open(F, "gunzip -c $filename |") || die "can't open pipe to $filename";
    } else {
        open(F, $filename) || die "can't open $filename";
    }

    my $n_read = 0;
    my $n_kept = 0;
    my $n_lq = 0;
    my $n_graph_fail = 0;
    my $n_not_snv = 0;

    while(<F>)
    {
        chomp;
        next if /^#/;
        next if /^GL/;
        
        $n_read++;
        my @fields = split(' ');
        my $chr = $fields[0];
        my $pos = $fields[1];
        my $ref = $fields[3];
        my $alt = $fields[4];

        # Hard filter freebayes calls
        if($is_freebayes) {
            my $not_snv = length($ref) > 1 || length($alt) > 1;
            my $is_low_quality = $fields[5] < $min_quality;
            my $not_graph_somatic = (index($fields[7], "KmerClassification=SOMATIC") == -1);

            $n_not_snv += $not_snv;
            $n_lq += $is_low_quality;
            $n_graph_fail += $not_graph_somatic;

            next if $not_snv || $is_low_quality || $not_graph_somatic;
        }

        my $tags = $fields[7];
        my $key = join(";", ($chr, $pos, $ref, $alt));
        $vcfHash->{$key}->[$idx] = $_;
        $n_kept++;
    }

    printf STDERR ("Caller: %s variants: %d kept: %d (LQ: %d COMPLEX: %d GRAPH-FAIL: %d)\n", 
        $callers[$idx], $n_read, $n_kept,  $n_lq, $n_not_snv, $n_graph_fail)
}
