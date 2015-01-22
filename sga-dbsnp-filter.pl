#! /usr/bin/perl
# Perform filtering of SGA variant calls

use strict;
use Getopt::Long;
use File::Basename;
use IPC::Cmd qw[can_run];

my $dbsnp_path = ""; # Filter variants against dbSNP at the given directory
my $cosmic_path = ""; # Filter variants against cosmic at the given directory
my $passed_only = 0;
my $min_af = 0;
my $tumor_bam = "";
my $normal_bam = "";
my $outname = "";
my $extra_dir = ""; # where all the dependencies live

GetOptions("dbsnp=s"       => \$dbsnp_path,
           "cosmic=s"      => \$cosmic_path,
           "extra-dir=s"   => \$extra_dir);

# Set paths to dependencies
my $samtools_bin = "$extra_dir/samtools/samtools";
my $vcfsort_bin = "$extra_dir/vcflib/bin/vcfsort";
my $bgzip_bin = "$extra_dir/htslib/bgzip";
my $bcftools_bin = "$extra_dir/bcftools/bcftools";

# check dependencies
check_prerequisites($samtools_bin,
                    $vcfsort_bin,
                    $bgzip_bin,
                    $bcftools_bin);

my $input_file = $ARGV[0];
if($input_file eq "") {
    print STDERR "Error no input given\n";
    exit(1);
}

my %non_dbsnp_sites;
my %cosmic_sites;

my $using_dbsnp = 0;
my $has_index = 0;
if($dbsnp_path ne "") {
    index_file($input_file);
    $has_index = 1;
    load_nondbsnp_sites($input_file, $dbsnp_path);
    $using_dbsnp = 1;
}

# whitelist cosmic sites if dbsnp is used
my $using_cosmic = 0;
if($using_dbsnp and $cosmic_path ne "") {
    load_cosmic_sites($input_file, $cosmic_path);
    $using_cosmic = 1;
}

if ($has_index == 1) {
    cleanup_index_file($input_file);
}

perform_filter($input_file);

sub perform_filter
{
    my($in) = @_;

    open(IN, $in) || die("Cannot open $in");
    my $total_out = 0;
    my %seen_hash;

    while(<IN>) {
        # Print header lines
        if(/^#/) {
            print;
            next;
        }

        chomp;
        
        # If this call is exactly the same as a previous call, ignore it
        my @f = split;
        my $key = sprintf("%s.%d.%s.%s", $f[0], $f[1], $f[3], $f[4]);
        if(defined($seen_hash{$key})) {
            next;
        } else {
            $seen_hash{$key} = 1;
        }

        if(defined($non_dbsnp_sites{$key}) or defined($cosmic_sites{$key})) {
            $total_out++;
        } else {
            my $ft = "dbSNP";
            if($f[6] eq "PASS" || $f[6] eq ".") {
                $f[6] = $ft;
            } else {
                $f[6] .= ";$ft";
            }
        }

        print join("\t", @f) . "\n";
    }
    
    close(IN);

    print STDERR "$total_out calls not in dbsnp\n";
    return $outname;
}

# Build a hash of calls that are NOT present in dbsnp
sub load_nondbsnp_sites
{
    my($in, $path) = @_;
    
    # find the non-dbsnp sites
    open(SITES, "$bcftools_bin isec -C $in.tmp.sorted.vcf.gz $path |") || die("dbsnp: bcftools isec failed");
    while(<SITES>) {
        chomp;
        my @fields = split;
        my $site_key = "$fields[0].$fields[1].$fields[2].$fields[3]";
        $non_dbsnp_sites{$site_key} = 1;
    }
    close(SITES);
}

# Build a hash of calls that ARE present in cosmic
sub load_cosmic_sites
{
    my($in, $path) = @_;
    
    # find the cosmic sites
    opendir(DIR, $path) || die("Cannot open cosmic directory $path");
    my @vcffiles = grep(/.vcf$/ || /.vcf.gz$/, readdir(DIR));
    closedir(DIR);

    foreach my $cosmicfile (@vcffiles) {
        open(SITES, "$bcftools_bin isec -n=2 $in.tmp.sorted.vcf.gz $path/$cosmicfile |") || die("cosmic: bcftools isec failed");
        while(<SITES>) {
            chomp;
            my @fields = split;
            my $site_key = "$fields[0].$fields[1].$fields[2].$fields[3]";
            $cosmic_sites{$site_key} = 1;
        }
        close(SITES);
    }
    my $nkeys = scalar keys %cosmic_sites;
    print STDERR "$nkeys keys in cosmic hash\n";
}

# check that each program can be run
sub check_prerequisites
{
    my (@programs) = @_;

    foreach my $program (@programs) {
        if(!can_run($program)) {
            print STDERR "Error: could not find program $program.\n";
            print STDERR "Please set the --extra option to the sga-extra directory.\n";
            exit(1);
        }
    }
}

sub index_file
{
    my ($in) = @_;

    # sort, bgzip and tabix the file
    system("$vcfsort_bin $in > $in.tmp.sorted.vcf");
    system("$bgzip_bin -f $in.tmp.sorted.vcf");
    system("$bcftools_bin index -f $in.tmp.sorted.vcf.gz");
}

sub cleanup_index_file
{
    my ($in) = @_;

    # Cleanup tmp
    unlink("$in.tmp.sorted.vcf");
    unlink("$in.tmp.sorted.vcf.gz");
    unlink("$in.tmp.sorted.vcf.gz.tbi");
    unlink("$in.tmp.sorted.vcf.gz.csi");
}
