#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use IO::File;
use Data::Dumper;

my $gtf_file = '';
my $kallisto_transcript_matrix_in = '';
my $kallisto_transcript_matrix_out = '';

GetOptions('gtf_file=s'=>\$gtf_file, 
    	   'kallisto_transcript_matrix_in=s'=>\$kallisto_transcript_matrix_in, 
	       'kallisto_transcript_matrix_out=s'=>\$kallisto_transcript_matrix_out);


unless($gtf_file && $kallisto_transcript_matrix_in && $kallisto_transcript_matrix_out){
  print "\n\nRequired parameters missing\n\n";
  print "Usage:  kallisto_gene_matrix.pl --gtf_file=chr22_with_ERCC92.gtf  --kallisto_transcript_matrix_in=transcript_tpms_all_samples.tsv --kallisto_transcript_matrix_out=gene_tpms_all_samples.tsv\n\n";
  exit();
}

die "\n\nCould not find GTF file: $gtf_file\n\n" unless(-e $gtf_file);
die "\n\nCould not find input kallisto transcript matrix file: $kallisto_transcript_matrix_in\n\n" unless(-e $kallisto_transcript_matrix_in); 

# Build a map of transcript to gene IDs
my %trans;
my $gtf_fh = IO::File->new($gtf_file, 'r');

while (my $gtf_line = $gtf_fh->getline) {
  chomp($gtf_line);
  my @gtf_entry = split("\t", $gtf_line);
  
  # Process only lines that describe transcripts
  next unless $gtf_entry[2] eq 'transcript';
  
  my $g_id = '';
  my $t_id = '';
  
  # Extract the gene_id and transcript_id using regular expressions that include periods
  if ($gtf_entry[8] =~ /gene_id\s+"([\w\.]+)"/) {
    $g_id = $1;
  }
  
  if ($gtf_entry[8] =~ /transcript_id\s+"([\w\.]+)"/) {
    $t_id = $1;
  }
  
  # If either ID is not found, throw an error
  die "\n\nCould not identify gene and transcript id in GTF transcript line:\n$gtf_line\n\n" 
    unless ($g_id && $t_id);
  
  # Map the transcript ID to the gene ID
  $trans{$t_id}{g_id} = $g_id;
}

$gtf_fh->close;


#Parse the kallisto input file, determine gene IDs for each transcript, and calculate sum TPM values
my %gene_exp;
my %gene_length;
my %samples;
my $ki_fh = IO::File->new($kallisto_transcript_matrix_in, 'r');
my $header = '';
my $h = 0;

while (my $ki_line = $ki_fh->getline) {
    $h++;
    chomp($ki_line);
    my @ki_entry = split("\t", $ki_line);
    my $s = 0;

    if ($h == 1) {
        $header = $ki_line;
        my $first_col = shift @ki_entry;
        my $second_col = shift @ki_entry;
        foreach my $sample (@ki_entry) {
            $s++;
            $samples{$s}{name} = $sample;
        }
        next;
    }

    my $trans_id = shift @ki_entry;
    my $length = shift @ki_entry;
    my $gene_id;

    # Parse gene_id from trans_id formatted as ENST00000456328.2|ENSG00000290825.1|...
    if ($trans_id =~ /^[^\|]+\|([^\|]+)/) {
        $gene_id = $1;
    } else {
        die "\n\nCould not identify gene id from trans id: $trans_id\n\n";
    }

    # Accumulate expression values
    $s = 0;
    foreach my $value (@ki_entry) {
        $s++;
        $gene_exp{$gene_id}{$s} += $value;
    }

    # Update gene length
    if ($gene_length{$gene_id}) {
        $gene_length{$gene_id} = $length if ($length > $gene_length{$gene_id});
    } else {
        $gene_length{$gene_id} = $length;
    }
}

$ki_fh->close;

my $ko_fh = IO::File->new($kallisto_transcript_matrix_out, 'w');
unless ($ko_fh) {
    die('Failed to open file: ' . $kallisto_transcript_matrix_out);
}

# Print header
print $ko_fh "$header\n";

# Print the results
foreach my $gene_id (sort { $a cmp $b } keys %gene_exp) {
    print $ko_fh "$gene_id\t$gene_length{$gene_id}\t";
    my @vals;
    foreach my $s (sort { $a <=> $b } keys %samples) {
        push(@vals, $gene_exp{$gene_id}{$s});
    }
    my $val_string = join("\t", @vals);
    print $ko_fh "$val_string\n";
}

$ko_fh->close;

exit;
