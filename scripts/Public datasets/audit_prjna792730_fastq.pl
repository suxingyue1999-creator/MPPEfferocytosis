#!/usr/bin/env perl

use strict;
use warnings;
use Digest::MD5 qw(md5_hex);
use File::Basename qw(basename);
use IO::Uncompress::Gunzip qw($GunzipError);

my ($manifest, $fastq_dir, $output) = @ARGV;
die "Usage: audit_prjna792730_fastq.pl <ENA manifest with MD5> <fastq directory> <output TSV>\n"
  unless defined $output;

open my $manifest_fh, '<', $manifest or die "Cannot read $manifest: $!\n";
my $header = <$manifest_fh>;
die "Manifest has no header: $manifest\n" unless defined $header;

my %sample_for_run = (
    SRR17345233 => 'C1-0h', SRR17345232 => 'C2-0h', SRR17345231 => 'C3-0h',
    SRR17345230 => 'T1-12h', SRR17345229 => 'T2-12h', SRR17345228 => 'T3-12h',
    SRR17345227 => 'T1-24h', SRR17345226 => 'T2-24h', SRR17345225 => 'T3-24h',
);

my @runs;
while (my $line = <$manifest_fh>) {
    chomp $line;
    my @field = split /\t/, $line, -1;
    my ($run, undef, undef, $bytes, $md5, $read_count) = @field;
    die "Unexpected manifest row for $run\n" unless exists $sample_for_run{$run};
    my @bytes = split /;/, $bytes;
    my @md5 = split /;/, $md5;
    die "Incomplete ENA metadata for $run\n" unless @bytes == 2 && @md5 == 2;
    push @runs, {
        run => $run,
        sample => $sample_for_run{$run},
        expected_bytes => \@bytes,
        expected_md5 => \@md5,
        expected_records => $read_count,
    };
}
close $manifest_fh;

die "The manifest must contain nine runs\n" unless @runs == 9;

sub file_md5 {
    my ($path) = @_;
    open my $fh, '<:raw', $path or die "Cannot read $path: $!\n";
    my $ctx = Digest::MD5->new;
    $ctx->addfile($fh);
    close $fh;
    return $ctx->hexdigest;
}

sub next_fastq_record {
    my ($fh) = @_;
    my $header = <$fh>;
    return undef unless defined $header;
    my $sequence = <$fh>;
    my $plus = <$fh>;
    my $quality = <$fh>;
    return { malformed => 1, id => undef } unless defined $sequence && defined $plus && defined $quality;
    chomp($header, $sequence, $plus, $quality);
    my $malformed = ($header !~ /^@/) || ($plus !~ /^\+/) || (length($sequence) != length($quality));
    my ($id) = $header =~ /^@(\S+)/;
    return { malformed => $malformed, id => $id };
}

sub inspect_pair {
    my ($r1, $r2) = @_;
    my $r1_fh = IO::Uncompress::Gunzip->new($r1)
      or die "Cannot decompress $r1: $GunzipError\n";
    my $r2_fh = IO::Uncompress::Gunzip->new($r2)
      or die "Cannot decompress $r2: $GunzipError\n";

    my ($records_1, $records_2, $malformed_1, $malformed_2, $mismatches) = (0, 0, 0, 0, 0);
    while (1) {
        my $record_1 = next_fastq_record($r1_fh);
        my $record_2 = next_fastq_record($r2_fh);
        last unless defined $record_1 || defined $record_2;
        if (defined $record_1) {
            ++$records_1;
            ++$malformed_1 if $record_1->{malformed};
        }
        if (defined $record_2) {
            ++$records_2;
            ++$malformed_2 if $record_2->{malformed};
        }
        if (!defined $record_1 || !defined $record_2 || $record_1->{malformed} || $record_2->{malformed}
            || !defined $record_1->{id} || !defined $record_2->{id} || $record_1->{id} ne $record_2->{id}) {
            ++$mismatches;
        }
    }
    close $r1_fh or die "Decompression error in $r1: $GunzipError\n";
    close $r2_fh or die "Decompression error in $r2: $GunzipError\n";
    return ($records_1, $records_2, $malformed_1, $malformed_2, $mismatches);
}

open my $output_fh, '>', $output or die "Cannot write $output: $!\n";
print {$output_fh} join("\t", qw(accession run sample read file expected_bytes observed_bytes expected_md5 observed_md5 gzip_test fastq_records expected_fastq_records malformed_records paired_id_mismatches status)), "\n";

for my $entry (@runs) {
    my $r1 = "$fastq_dir/$entry->{run}_1.fastq.gz";
    my $r2 = "$fastq_dir/$entry->{run}_2.fastq.gz";
    die "Missing FASTQ file: $r1\n" unless -f $r1;
    die "Missing FASTQ file: $r2\n" unless -f $r2;

    my ($records_1, $records_2, $malformed_1, $malformed_2, $mismatches) = inspect_pair($r1, $r2);
    for my $read_index (0, 1) {
        my $path = $read_index == 0 ? $r1 : $r2;
        my $observed_bytes = -s $path;
        my $observed_md5 = file_md5($path);
        my $records = $read_index == 0 ? $records_1 : $records_2;
        my $malformed = $read_index == 0 ? $malformed_1 : $malformed_2;
        my $status = ($observed_bytes == $entry->{expected_bytes}->[$read_index]
            && $observed_md5 eq $entry->{expected_md5}->[$read_index]
            && $records == $entry->{expected_records}
            && $malformed == 0
            && $mismatches == 0) ? 'PASS' : 'FAIL';
        print {$output_fh} join("\t",
            'PRJNA792730', $entry->{run}, $entry->{sample}, 'R' . ($read_index + 1),
            basename($path), $entry->{expected_bytes}->[$read_index], $observed_bytes,
            $entry->{expected_md5}->[$read_index], $observed_md5, 'PASS', $records,
            $entry->{expected_records}, $malformed, $mismatches, $status
        ), "\n";
    }
}
close $output_fh;
