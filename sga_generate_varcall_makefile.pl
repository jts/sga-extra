#! /usr/bin/perl
use strict;
use Getopt::Long;
use POSIX;
use IPC::Cmd qw[can_run];
use File::Basename;

my $reads_fofn = "";
my $base_fofn = "";
my $variant_bam = "";
my $base_bam = "";
my $reference_file;
my $sga_bin = "sga";
my $sga_extra_dir = "";
my $vcflib = "";
my $pp_opt = '-q 20 --discard-quality';
my $min_dbg_opt = 2;
my $min_count_opt = 4;
my $k_opt = 55;
my $dbsnp = "";
my $project_name = "test";
my $variant_bp_estimate = 0;
my $base_bp_estimate = 0;

my @chromosomes = (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 'X', 'Y');

GetOptions("project=s"        => \$project_name,
           "reads-fofn=s"     => \$reads_fofn,
           "base-fofn=s"      => \$base_fofn,
           "variant-bam=s"    => \$variant_bam,
           "base-bam=s"       => \$base_bam,
           "reference-file=s" => \$reference_file,
           "sga=s"            => \$sga_bin,
           "sga-extra-dir=s"  => \$sga_extra_dir,
           "dbsnp=s"          => \$dbsnp,
           "variant-bp=i"     => \$variant_bp_estimate,
           "base-bp=i"        => \$base_bp_estimate);


die("A reference is required (--reference)\n") if($reference_file eq "");

die("--sga-extra-dir must be provided\n") if($sga_extra_dir eq "");

# Set paths to executable
my $samtools_bin = "$sga_extra_dir/samtools/samtools";
my $freebayes_bin = "$sga_extra_dir/freebayes/bin/freebayes";
my $bam2fastq_bin = "$sga_extra_dir/bam2fastq/bam2fastq";
my $sga_dbsnp_script = "$sga_extra_dir/sga-dbsnp-filter.pl";
my $sga_fb_merge_script = "$sga_extra_dir/sga_freebayes_merge.pl";
my $bcftools_bin = "$sga_extra_dir/bcftools/bcftools";
my $vcflib_bin_path = "$sga_extra_dir/vcflib/bin/";

# check prerequisites
check_prerequisites();

# determine whether chromosomes are named "chr1" or "1"
my $chr_prefix = get_chromosome_prefix($reference_file);

# write out the preamble
print_preamble();

# Write the input files as make variables to have shorter command lines later
print_filepaths();

print_separator("sga calling");

my $variant_pp_file = "";
my $base_pp_file = "";

# Set up the input filenames and write the sga preprocess rules
if($reads_fofn ne "") {

    $variant_pp_file = print_input_and_preprocess_fastq("variant", $reads_fofn);
    $base_pp_file = print_input_and_preprocess_fastq("base", $base_fofn) if($base_fofn ne "");

} elsif($variant_bam ne "") {

    $variant_pp_file = print_input_and_preprocess_bam("variant", $variant_bam);
    $base_pp_file = print_input_and_preprocess_bam("base", $base_bam) if($base_bam ne "");

} else {
    die("No input sources provided\n");
}

# Estimate the number of bases if it was not provided
if($variant_bp_estimate == 0) {
    if($reads_fofn ne "") {
        $variant_bp_estimate = estimate_bases_from_fastq_fofn($reads_fofn);
    } else {
        $variant_bp_estimate = estimate_bases_from_bam($variant_bam);
    }   
}

if($base_bp_estimate == 0) {
    if($base_fofn ne "") {
        $base_bp_estimate = estimate_bases_from_fastq_fofn($base_fofn);
    } elsif($base_bam ne "") {
        $base_bp_estimate = estimate_bases_from_bam($base_bam);
    }
}

# write the sga index rule
my $max_bp = ($variant_bp_estimate > $base_bp_estimate ? $variant_bp_estimate : $base_bp_estimate);
print_index($max_bp);

# write the graph-diff rule
my $sga_somatic = print_graph_diff($project_name, $variant_pp_file, $base_pp_file, $variant_bp_estimate, $base_bp_estimate);

# write the preqc rule
#print_preqc();

print_separator("freebayes calling");

# write the freebayes calling rule
my($fb_somatic, $fb_germline) = print_freebayes_calls($project_name, $variant_bam, $base_bam);

print_separator("filtering");

# write the freebayes filtering rule
my $fb_graph_filtered = print_freebayes_filtering($project_name, $fb_somatic, $fb_germline, 
                                                  $variant_pp_file, $base_pp_file,
                                                  $variant_bp_estimate, $base_bp_estimate);


print_variant_postprocess($project_name);

my $variant_suffix = ".leftalign.filters.dbsnp.vcf";
my $sga_filters = basename($sga_somatic, ".vcf") . $variant_suffix;
my $fb_filters = basename($fb_graph_filtered, ".vcf") . $variant_suffix;

print_vcf_merge($project_name, $sga_filters, $fb_filters);

exit(0);

sub print_input_and_preprocess_fastq
{
    my($name, $fofn) = @_;

    open(F, $fofn) || die("Cannot read $fofn");

    # We assume the files are separated by whitespace
    # in the file-of-filenames
    my @files;
    while(<F>) {
        chomp;
        push @files, split(' ', $_);
    }

    # Classify files as pe/se
    my @pe_files;
    my @se_files;
    my %h_files;

    map { $h_files{$_} = 1 } @files;

    foreach my $f (@files) {
        my $p = "";

        if ($f =~ /R1/) {
            # Check for a file with s/R1/R2/
            ($p = $f) =~ s/R1/R2/;
            if(defined($h_files{$p})) {
                push @pe_files, ($f, $p);
            } else {
                die("Pair file not found for $f\n");
            }
        } elsif($f !~ /R2/) {
            push @se_files, $f;
        }
    }

    my $pe_var = uc($name) . "_PE";
    my $se_var = uc($name) . "_SE";
    my $outfile = "$name.fastq.gz";

    if(scalar(@pe_files) > 0) {
        printf("$pe_var=");
        print join(" \\\n\t", @pe_files) . "\n";
    }
    
    if(scalar(@se_files) > 0) {
        printf("\n$se_var=");
        print join(" \\\n\t", @se_files) . "\n";
    }

    die("No input files\n") if(scalar(@pe_files) == 0 && scalar(@se_files) == 0);

    printf("\n");
    printf("# Preprocess the $name input files\n");
    printf("%s:\n", $outfile);
    printf("\trm -f \$@\n\n");

    my $pp_fmt_str = "\t\$(SGA) preprocess %s \$(%s) | gzip >> \$@\n";
    
    if(scalar(@pe_files) > 0) {
        printf($pp_fmt_str, "$pp_opt --pe-mode 1", $pe_var);
    } 

    if(scalar(@se_files) > 0) {
        printf($pp_fmt_str, $pp_opt, $se_var);
    }

    printf "\n";
    return $outfile;
}

sub print_input_and_preprocess_bam
{
    my($name, $bam) = @_;
    die("No bam provided") if $bam eq "";

    my $outfile = "$name.fastq.gz";

    printf("\n");
    printf("# Preprocess the $name input bam files\n");
    printf("%s:\n", $outfile);
    printf("\trm -f \$@\n\n");

    my $pp_str = "\t$bam2fastq_bin --pairs-to-stdout $bam | \$(SGA) preprocess --pe-mode 2 - | gzip >> \$@";
    printf "$pp_str\n";
    return $outfile;
}

sub print_index
{
    my($max_file_size) = @_;

    # sga index in ropebwt mode uses ~0.33 bytes per input base
    # we add an extra bit of overhead to avoid having unusual
    # jobs fail.
    my $max_memory = ($max_file_size * 0.4) / 1000000000;

    my $resource = get_resource_string($max_memory, 4);

    printf("\n");
    printf("# Build an FM-index for reads\n");
    print("%.bwt: %.fastq.gz\n");
    print("\t$resource ");
    print("/usr/bin/time -v \$(SGA) index -a ropebwt -t 4 --no-reverse --no-sai \$<\n");

    $resource = get_resource_string($max_memory, 8);
    print("%.sai: %.fastq.gz %.bwt\n");
    print("\t$resource ");
    print("/usr/bin/time -v \$(SGA) gen-ssa --sai-only -t 8 \$<\n");
}

sub print_graph_diff
{
    my($name, $var_file, $base_file, $var_size, $base_size) = @_;

    my $base_opt = "";
    if($base_file ne "") {
        $base_opt = "-b $base_file";
    }

    my $max_memory = ($var_size + $base_size) * 0.4 / 1000000000;
    my $resource = get_resource_string($max_memory, 8);

    my $gd_fmt_str = "\$(SGA) graph-diff -p $name.sga " .
                     "--min-dbg-count $min_dbg_opt -k $k_opt -x $min_count_opt ".
                     "-t 8 -a debruijn -r $var_file $base_opt --ref \$(REFERENCE)";
    
    my @var_index = get_index_names($var_file);
    my @base_index = get_index_names($base_file);
    
    my $sga_out = "$name.sga.calls.vcf";

    print "\n# Make SGA calls\n";
    print "$sga_out: @var_index @base_index\n";
    print "\t$resource ";
    print "/usr/bin/time -v $gd_fmt_str\n";

    return $sga_out;
}

sub print_freebayes_calls
{
    my($name, $var_bam, $base_bam) = @_;
    
    # The structure of the file names
    my $chr_vcf_fmt = "%s.freebayes.%s.calls.vcf";

    # Extract the sample names from the bam
    my $var_sample_name = get_sample_name($var_bam);
    my $base_sample_name = get_sample_name($base_bam);

    # Raw calls per chromosome
    my $calling_memory = 4;
    my $resource = get_resource_string($calling_memory, 1);
    my $output_pattern = sprintf($chr_vcf_fmt, $name, "%");
    
    my $fb_fmt_str = "\$(FB) -r \$* -f \$(REFERENCE) --pooled-discrete --pooled-continuous " .
                     "--min-alternate-fraction 0.1 --allele-balance-priors-off --genotype-qualities \$(VARIANT_BAM) \$(BASE_BAM)";

    print "\n# Run freebayes on each chromosome independently\n";
    print "$output_pattern:\n";
    print "\t$resource /usr/bin/time -v $fb_fmt_str > \$@\n\n";

    # Make the per-chromosome file names
    my $stride = 4;
    my $chr_file_string = "";
    for(my $i = 0; $i < scalar(@chromosomes); $i += $stride) {
        $chr_file_string .= " \\\n\t" unless $chr_file_string eq "";
        
        my $end = $i + $stride - 1;
        if($end >= scalar(@chromosomes)) {
            $end = scalar(@chromosomes) - 1;
        }
        my @sub = @chromosomes[$i .. $end];
        $chr_file_string .= join(" ", map{ sprintf($chr_vcf_fmt, $name, $chr_prefix . $_) } @sub);
    }

    # Merge per-chromosome calls
    print "FREEBAYES_PER_CHR=" . $chr_file_string . "\n";
    print "\n# Merge freebayes calls and tag with somatic status\n";
    print "$name.freebayes.allchr.vcf: \$(FREEBAYES_PER_CHR)\n";
    print "\t$vcflib_bin_path/vcfcombine \$^ | $vcflib_bin_path/vcfbreakmulti | $vcflib_bin_path/vcfsamplediff -s VT $base_sample_name $var_sample_name - > \$@\n";

    # Split calls
    my $fb_somatic = "$name.freebayes.allchr.somatic.vcf";
    my $fb_germline = "$name.freebayes.allchr.germline.vcf";

    print "\n# Split freebayes calls by somatic status\n";
    print "$fb_somatic: $name.freebayes.allchr.vcf\n";
    print "\tcat \$< | awk '\$\$1 ~ /#/ || \$\$0 ~ /somatic/' > \$@\n";
    print "$fb_germline: $name.freebayes.allchr.vcf\n";
    print "\tcat \$< | awk '\$\$1 ~ /#/ || \$\$0 ~ /germline/ || \$\$0 ~ /reversion/' > \$@\n";
    
    return ($fb_somatic, $fb_germline);
}

sub print_freebayes_filtering
{
    my($name, $fb_somatic, $fb_germline, $var_file, $base_file, $var_size, $base_size) = @_;
    
    my $max_memory = ($var_size + $base_size) * 0.4 / 1000000000;
    my $resource = get_resource_string($max_memory, 8);

    # Filter freebayes against the assembly graph
    my $gc_fmt_str = "\$(SGA) graph-concordance --ref \$(REFERENCE) " .
                     "-r $var_file -b $base_file " .
                     "-g $fb_germline $fb_somatic";

    my $out = "$name.freebayes.allchr.somatic.graph.vcf";
    my @var_index = get_index_names($var_file);
    my @base_index = get_index_names($base_file);
    
    print "\n# Filter freebayes calls against the assembly graph\n";
    print "$out: $fb_somatic $fb_germline @var_index @base_index\n";
    print "\t$resource ";
    print "/usr/bin/time -v $gc_fmt_str 2> graphfilter.err > \$@\n"; 
    return $out;
}

sub print_variant_postprocess
{
    my($name) = @_;

    print "\n# Left align SGA calls and remove calls on unplaced chromosomes\n";
    print "%.leftalign.vcf: %.vcf\n";
    print "\t$bcftools_bin norm -f \$(REFERENCE) \$< | awk '\$\$1 ~ /#/ || \$\$1 !~ /hs37d5/'> \$@\n";

    print "\n# Filter calls\n";
    my $max_memory = 16;
    my $resource = get_resource_string($max_memory, 1);

    my $filter_str = "--min-var-dp 4 --min-af 0.1";

    print "%.filters.vcf: %.vcf \$(VARIANT_BAM) \$(BASE_BAM)\n";
    print "\t$resource ";
    print "\t\$(SGA) somatic-variant-filters $filter_str --tumor \$(VARIANT_BAM) --normal \$(BASE_BAM) --reference \$(REFERENCE) \$< > \$@\n";

    print "\n# Mark dbSNP variants\n";
    print "%.dbsnp.vcf: %.vcf\n";
    print "\t$sga_dbsnp_script --extra $sga_extra_dir --dbsnp $dbsnp \$< > \$@\n";
}

# Merge the results of sga and freebayes
sub print_vcf_merge
{
    my($name, $sga_vcf, $freebayes_vcf) = @_;
    my $mergedname = "$name.merged.vcf";
    printf("\n");
    printf("# Merge the results of SGA and freebayes\n");
    print("$mergedname: $sga_vcf $freebayes_vcf\n");
    printf("\t$sga_fb_merge_script \$^ > \$@\n");

    my $mergedpassedname = "$name.merged.passed.vcf";
    printf("\n");
    printf("# Subset to passed-only sites\n");
    print("$mergedpassedname: $mergedname\n");
    printf("\tcat \$^ | awk '\$\$1 ~ /#/ || \$\$7 ~ /PASS/' > \$@\n");
}

sub print_preqc
{
    printf("\n");
    printf("# Generate a quality report\n");
    print("%.preqc: %.bwt\n");
    printf("\t\$(SGA) preqc -t 8 \$*.fastq.gz > \$@\n");
}

sub print_preamble
{
    print_separator("Program paths");
    printf("SHELL=/bin/bash -o pipefail\n");
    printf("SGA=%s\n", $sga_bin);
    printf("FB=%s\n", $freebayes_bin);
    printf("SAMTOOLS=%s\n", $samtools_bin);
    printf("# do not delete intermediate files\n");
    printf(".SECONDARY:\n");
}

sub print_filepaths
{
    print_separator("File paths");
    printf("REFERENCE=%s\n", $reference_file);
    printf("VARIANT_BAM=%s\n", $variant_bam);
    printf("BASE_BAM=%s\n", $base_bam);
}

sub print_separator
{
    my($text) = @_;
    print "\n" . join("\n", ("############################", "# $text", "############################")) . "\n";
}

sub check_prerequisites
{
    my @programs = ($sga_bin, 
                    $freebayes_bin,
                    $sga_fb_merge_script, 
                    "$vcflib_bin_path/vcfcombine",
                    $bcftools_bin,
                    $samtools_bin,
                 );

    foreach my $program (@programs) {
        if(!can_run($program)) {
            print STDERR "Could not find program $program. Please install it or update your PATH.\n";
            exit(1);
        }
    }
}

# convert filenames reads.fastq.gz -> reads.bwt reads.sai
sub get_index_names
{
    my($f) = @_;
    (my $b = $f) =~ s/fastq.gz$/bwt/;
    (my $s = $f) =~ s/fastq.gz$/sai/;
    return ($b, $s);
}

# Parse the SM tag from a sam header
sub get_sample_name
{
    my($bam_file) = @_;
    my $sample_name = "";
    open(F, "samtools view -H $bam_file |");
    while(<F>) {
        next unless /^\@RG.*SM:(\S+)/;

        $sample_name = $1 if($sample_name eq "");
        die("Error: multiple sample names in bam file $bam_file ($1 != $sample_name)") if($1 ne $sample_name);
    }
    return $sample_name;
}

# 
sub get_chromosome_prefix
{
    my($reference) = @_;
    open(F, $reference) || die("Cannot open reference $reference");

    # Search for chromosome 1 in file 
    while(<F>) {
        next unless />(\S+)/;
        return "chr" if($1 eq "chr1");
        return "" if($1 eq "1");
    }
    
    die("Could not determine chromosome prefix");
}

# Estimate the total number of bases from a collection of fastq files
sub estimate_bases_from_fastq_fofn
{
    my($fofn) = @_;

    my $sum_gzipped_filesize = 0;
    my $sum_text_filesize = 0;

    open(F, $fofn) || die("Cannot open $fofn");
    while(<F>) {
        chomp;

        $sum_text_filesize += -s $_ if(/\.fastq$/);
        $sum_gzipped_filesize += -s $_ if(/\.fastq.gz$/);
    }

    die("Estimating size is currently only supported from gzipped fastq") if($sum_text_filesize > 0);

    # Gzipped fastq compressed to approximately 1 byte per base
    return $sum_gzipped_filesize;
}

# Estimate the total number of bases from a bam file
sub estimate_bases_from_bam
{
    my($bamfile) = @_;

    my $size = -s $bamfile;

    # Bam requires about 1 byte/base
    return $size;
}


# Make a formatted resource string for SGE
sub get_resource_string
{
    my($memory_gb, $threads) = @_;

    # Memory is requested per-thread
    $memory_gb = ceil($memory_gb / $threads);

    my $thread_req = ($threads > 1) ? "-pe smp $threads" : "";
    return sprintf("SGE_RREQ=\"-l h_vmem=%sG -l h_stack=32M $thread_req\"", $memory_gb, $threads);
}
