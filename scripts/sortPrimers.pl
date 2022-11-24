use warnings;
no warnings ('uninitialized', 'substr');
use Cwd;
use Time::Piece;
use Getopt::Long;

my $sample = "sample";
my $work_dir = "results";
my $prefix = "pcr";
my $primer_seq = "primer";

GetOptions(
    'sample=s'    => \$sample,
    'work_dir=s'     => \$work_dir,
    'primer=s' => \$prefix,
    'primer_seq=s'     => \$primer_seq
) or print "Invalid options\n";

my $localtime = localtime;
print "Primer sorting: starting at $localtime\n\n";
print "$sample is searched for $prefix primer using $primer_seq sequences\n\n";


my %forward_primer = ();
my %reverse_primer = ();
my %forward_primer_revcomp = ();
my %reverse_primer_revcomp = ();
my $count_for = 0;
my $count_rev = 0;

open (IN, "$primer_seq") or print "Cannot read/find $primer_seq\n";
while(<IN>){
	chomp $_;
	if ($_ =~ /^>(\S+)/){
		$id = $1;
	}
	else{
		$seq = $_;
		$seq_revcomp = "";
		if ($id =~ /for/i){
			$count_for++;
			$forward_primer{$count_for} = $seq;
			$seq_revcomp = reverse $seq;
			$seq_revcomp =~ tr/ATGCatgc/TACGtacg/;
			$forward_primer_revcomp{$count_for} = $seq_revcomp;
		} 
		elsif ($id =~ /rev/i){
			$count_rev++;
			$reverse_primer{$count_rev} = $seq;
			$seq_revcomp = reverse $seq;
			$seq_revcomp =~ tr/ATGCatgc/TACGtacg/;
			$reverse_primer_revcomp{$count_rev} = $seq_revcomp;
		}
	}
}

print "Forward primer...\n";
foreach my $p (sort {$a <=> $b} keys %forward_primer){
	print "$p $forward_primer{$p} $forward_primer_revcomp{$p}\n";
}

print "Reverse primer....\n";
foreach my $p (sort {$a <=> $b} keys %reverse_primer){
	print "$p $reverse_primer{$p} $reverse_primer_revcomp{$p}\n";
}
my %unique_sequences = ();
my %primer_group = ();
my $total_primer_reads_flash = 0;
my $total_primer_reads_single = 0;
my $total_unique_seq = 0;
my $id = "";
my $sequence = "";
my $trimmed_seq = "";
my $flag = 0;

sub check_primer{
	$sequence = $_[0];
	$flag = 0;
	LOOP1: foreach my $primer1 (keys %forward_primer){
		if ($sequence =~ /^(\w+|)$forward_primer{$primer1}(\w+|)/) {
			$pre_for = $1;
			foreach my $primer2 (keys %reverse_primer_revcomp){
				if ($2 =~ /(\w+|)$reverse_primer_revcomp{$primer2}(\w+|)$/) {
					if (length($1) > 0){
						if (exists $unique_sequences{$1}){
							$unique_sequences{$1} = $unique_sequences{$1} + 1;
						}
						else {
							$unique_sequences{$1} = 1;
						}
						$flag = 1;
						$primer_group{"$sample\t$pre_for\t$forward_primer{$primer1}\t$reverse_primer_revcomp{$primer2}\t$2\t+"}++;
						$trimmed_seq = $1;
						last LOOP1;
					}
				}
			}
		}
	}
	if ($flag eq 0){
		LOOP2: foreach my $primer1 (keys %reverse_primer){
			if ($sequence =~ /^(\w+|)$reverse_primer{$primer1}(\w+|)/) {
				$pre_for = $1;
				foreach my $primer2 (keys %forward_primer_revcomp){
					if ($2 =~ /(\w+|)$forward_primer_revcomp{$primer2}(\w+|)$/) {
						if (length($1) > 0){
							my $seq_revcomp = reverse $1;
							$seq_revcomp =~ tr/ATGCatgc/TACGtacg/;
							if (exists $unique_sequences{$seq_revcomp}){
								$unique_sequences{$seq_revcomp} = $unique_sequences{$seq_revcomp} + 1;
							} else {
								$unique_sequences{$seq_revcomp} = 1;
							}
							$flag = 1;
							$primer_group{"$sample\t$pre_for\t$forward_primer_revcomp{$primer2}\t$reverse_primer{$primer1}\t$2\t-"}++;
							$trimmed_seq = $seq_revcomp;
							last LOOP2;
						}
					}
				}
			}
		}
	}
	if ($flag eq 1){
		return 1;
	}
	else{
		return 0;
	}
}

my $flash_reads = "$work_dir/$sample.extendedFrags.fastq";
my $single_reads = "$work_dir/$sample.trimmed_singles.fastq";
# open (OUT_PRIMERS, ">$work_dir/$sample.$prefix.fasta") or die "Cannot write $work_dir/$sample/$sample.$prefix.fasta\n";
my %all_len = ();
my %flash_len = ();
my %single_len = ();

if ( -e $flash_reads ){
	open(IN, "$flash_reads");
	$line = 0;
	while(<IN>){
		chomp $_;
		if ($line eq 0){
			$id = $_;
			$line = 1;
			next;
		}
		if ($line eq 1){
			$sequence = $_;
			my $flag_primer = check_primer($sequence);
			if ($flag_primer eq 1){
				$total_primer_reads_flash++;
				$len = length($trimmed_seq);
				$all_len{$len}++;
				$flash_len{$len}++;
				# print OUT_PRIMERS "$id/$sample/extended/$len\n$trimmed_seq\n";
			}
			$line = 2;
			next;
		}
		if ($line eq 2){
			$line = 3;
			next;
		}
		if ($line eq 3){
			$line = 0;
			next;
		}
	}
}
else{
	print "Cannot read $flash_reads.\n";
}
		
if ( -e $single_reads ){
	open(IN, "$single_reads");
	$line = 0;

	while(<IN>){
		chomp $_;
		if ($line eq 0){
			$id = $_;
			$line = 1;
			next;
		}
		if ($line eq 1){
			$sequence = $_;
			my $flag_primer = check_primer($sequence);
			if ($flag_primer eq 1){
				$total_primer_reads_single++;
				$len = length($trimmed_seq);
				$all_len{$len}++;
				$single_len{$len}++;
				# print OUT_PRIMERS "$id/$sample/single/$len\n$trimmed_seq\n";
			}
			$line = 2;
			next;
		}
		if ($line eq 2){
			$line = 3;
			next;
		}
		if ($line eq 3){
			$line = 0;
			next;
		}
	}
}
else{
	print "Cannot read $single_reads.\n";
}
# close(OUT_PRIMERS);

		
print "Writting lengh histogram...";
open (OUT_PRIMERS_LEN, ">$work_dir/$sample.seq.len_hist.tsv") or die "Cannot write $work_dir/$sample/$sample.$prefix.seq.len_hist.tsv\n";
print OUT_PRIMERS_LEN "Sequence Length\tOverlapping reads\tSingle reads\n";
foreach my $len (sort{$a <=> $b} keys %all_len){
	print OUT_PRIMERS_LEN "$len\t$flash_len{$len}\t$single_len{$len}\n";
}

print "Writting unique variants log $work_dir/$sample.primer.info.txt...\n";
open (INFO, ">$work_dir/$sample.primer.info.txt") or die "Cannot write $work_dir/$sample.$prefix.primer.info.txt\n";
print INFO "Primer\tSample\tSequence before Forward primer\tForward primer\tReverse primer\tSequence after Reverse primer\tOrientation\tread counts\n";
foreach $primer (sort {$primer_group{$b} <=> $primer_group{$a}} keys %primer_group){
	print INFO "$prefix\t$sample\t$primer\t$primer_group{$primer}\n";
}
close (INFO);


my $total_primer_reads = $total_primer_reads_single + $total_primer_reads_flash;

print "Total sequences with $prefix primer in overlapped reads: $total_primer_reads_flash\n";
print "Total sequences with $prefix primer in Single reads: $total_primer_reads_single\n";
print "Total sequences with $prefix : $total_primer_reads\n";





#chimera detection algorithm
sub possible_chimaera {
	my $sequence1 = $_[0];
	my $sequence2 = $_[1];
	my $candidate = $_[2];
	$candidate = uc($candidate);
	$sequence1 = uc($sequence1);
	$sequence2 = uc($sequence2);
	my %base1 = ();
	my %base2 = ();

	for (my $i = 0; $i < length($sequence1); $i++) {
		$base1{$i} = substr($sequence1, $i, 1 );
	}
	for (my $i = 0; $i < length($sequence2); $i++) {
		$base2{$i} = substr($sequence2, $i, 1 );
	}
	$len_diff2 = length($candidate) - length($sequence2);

	my $stop_base = 0;
	my $stop_base1 = 0;

	LOOP1: for (my $i = 0; $i < length($candidate); $i++) {
		my $base = substr($candidate, $i, 1);
		if ($base1{$i} ne $base) {
			$stop_base = $i;
			last LOOP1;
		}
	}
	LOOP2: for (my $i = (length($candidate) - 1); $i >= 0; $i--) {
		my $base = substr($candidate, $i, 1);
		if ($base2{$i - $len_diff2} ne $base) {
			$stop_base1 = $i;
			last LOOP2;
		}
	}
	if ($stop_base1 > 0 and $stop_base > 0 and $stop_base1 < $stop_base){
		return 1;
	}
	else{
		return 0;
	}
}


my $count_clusters = 0;
my $count_singles = 0;
my $count = 0;
my $count_variants_reads = 0;
my $count_chimera_reads = 0;
my $count_selected_reads = 0;
my $count_low_reads = 0;
my $count_gap_reads = 0;
my $count_variants = 0;
my $count_chimera = 0;
my $count_selected = 0;
my $count_low = 0;
my $count_gap = 0;
my $count_to_finish = 0;

open (FASTA_OUT, ">$work_dir/$sample.clusters.fasta") or die "Cannot write $work_dir/$sample/$sample.$prefix.clusters.fasta\n";
open (FASTA_OUT_SELECTED, ">$work_dir/$sample.clusters.filtered.fasta") or die "Cannot write $work_dir/$sample/$sample.$prefix.clusters.filtered.fasta\n";
open (SINGLES, ">$work_dir/$sample.singletons.fasta") or die "Cannot write $work_dir/$sample/$sample.$prefix.singletons.fasta\n";
open (INFO, ">$work_dir/$sample.clusters.details.tsv") or die "Cannot write $work_dir/$sample.clusters.details.tsv\n";
print INFO "Cluster\tType\tDetail\tSequence\n";

if($total_primer_reads > 0){
	foreach $seq (sort {$unique_sequences{$b} <=> $unique_sequences{$a}} keys %unique_sequences){
		$count++;
		$len = length($seq);
		# print "$unique_sequences{$seq} and $total_primer_reads\n";
		my $per = sprintf("%.2f", (($unique_sequences{$seq}/$total_primer_reads)*100));
		if ($unique_sequences{$seq} eq 1){
			$count_singles++;
			print SINGLES ">$count:$unique_sequences{$seq}:$per\n$seq\n";
		}
		else{
			$count_clusters++;
			my $flag_pcr_error = 0;
			my $flag_chimera;
			my $mismatch_counts = 0;
			my $mismatch = "";
			my $mismatch_with = "";
			my $fc = 0;
			my $flag_gap = 0;
			if ($count_clusters <= 200){
				print "Checking $count_clusters/$unique_sequences{$seq}/$per: ";
				my @nucl = split //, $seq;
				LOOPN: for($i = 0; $i < scalar(@nucl); $i++){
					if ($nucl[$i] !~ /[A|T|C|G]/i){
						$flag_gap = 1;
						last LOOPN;
					}
				}
				if ($flag_gap eq 1){
					$count_gap_reads = $count_gap_reads + $unique_sequences{$seq};
					$count_gap++;
					$details = "ambiguous_basecall\t";
					print "ambiguous_basecall\n";
					print FASTA_OUT ">$count_clusters/$unique_sequences{$seq}/$per/$len/ambiguous_basecall\n$seq\n";
					print FASTA_OUT_SELECTED ">$count_clusters/$unique_sequences{$seq}/$per/$len/ambiguous_basecall\n$seq\n";
					$ids{$seq} = "$count_clusters/$unique_sequences{$seq}/$per/$len/ambiguous_basecall";
				}
				else{
					LOOP7: foreach $seq1 (sort {$unique_sequences{$b} <=> $unique_sequences{$a}} keys %unique_sequences){
						if ($unique_sequences{$seq1} > $unique_sequences{$seq}){
							$len1 = length($seq1);
							if ($len1 eq $len1){
								$fc = sprintf('%.2f', $unique_sequences{$seq1} / $unique_sequences{$seq});
								$mismatch_with = $ids{$seq1};
								# print "PCR of $mismatch_with\n";
								$mismatch_counts = 0;
								$mismatch = undef;
								my @nucl = split //, $seq;
								my @nucl1 = split //, $seq1;
								for (my $i = 0; $i < scalar(@nucl); $i++){
									if ( $nucl[$i] ne $nucl1[$i] ){
										$mismatch_counts++;
										my $j = $i + 1;
										$mismatch = "$j/$nucl[$i]->$nucl1[$i]";
									}
								}
								if ($mismatch_counts eq 1){
									last LOOP7;
								}
							}
						}
						else{
							last LOOP7;
						}
					}
					# print "$mismatch_counts\n";
					if ($mismatch_counts eq 1){
						$flag_pcr_error = 1;
						$details = "1bpVariant\t$mismatch_with".":$fc:$mismatch";
						print "1bpVariant\t$mismatch_with".":$fc:$mismatch\n";
						print FASTA_OUT ">$count_clusters/$unique_sequences{$seq}/$per/$len/1bpVariant\n$seq\n";
						$ids{$seq} = "$count_clusters/$unique_sequences{$seq}/$per/$len/1bpVariant";
						$count_variants++;
						$count_variants_reads = $count_variants_reads + $unique_sequences{$seq};
					}

					if ($flag_pcr_error eq 0){
						$flag_chimera = 0;
						# print "Checking chimera...\n";
						LOOP6: foreach $seq1 (sort {$unique_sequences{$b} <=> $unique_sequences{$a}} keys %unique_sequences){
							if ($seq ne $seq1){
								if ($unique_sequences{$seq1} > $unique_sequences{$seq}){
									LOOP8: foreach $seq2 (sort {$unique_sequences{$b} <=> $unique_sequences{$a}} keys %unique_sequences){
										if ($seq ne $seq2 and $seq1 ne $seq2){
											if ($unique_sequences{$seq2} > $unique_sequences{$seq}){
												$check = possible_chimaera($seq1, $seq2, $seq);
												if ($check == 1){
													$flag_chimera = 1;
													$details{$seq} = "chimera\t$ids{$seq1}-$ids{$seq2}";
													print "chimera: $ids{$seq1}-$ids{$seq2}\n";
													print FASTA_OUT ">$count_clusters/$unique_sequences{$seq}/$per/$len/chimera\n$seq\n";
													$ids{$seq} = "$count_clusters/$unique_sequences{$seq}/$per/$len/chimera";
													$count_chimera++;
													$count_chimera_reads = $count_chimera_reads + $unique_sequences{$seq};
													last LOOP6;
												}
											}
											else{
												last LOOP8;
											}
										}
									}
								}
								else{
									last LOOP6;
								}
							}
						}
						if ($flag_chimera == 0){
							print "Filtered\n";
							print FASTA_OUT ">$count_clusters/$unique_sequences{$seq}/$per/$len/filtered\n$seq\n";
							print FASTA_OUT_SELECTED ">$count_clusters/$unique_sequences{$seq}/$per/$len/filtered\n$seq\n";
							$ids{$seq} = "$count_clusters/$unique_sequences{$seq}/$per/$len/filtered";
							$details = "filtered\t-";
							$count_selected++;
							$count_selected_reads = $count_selected_reads + $unique_sequences{$seq};
						}
					}
				}
			}
			else{
				$count_low++;
				$count_low_reads = $count_low_reads + $unique_sequences{$seq};
				# print "Not checked\n";
				$details = "notChecked\t-";
				print FASTA_OUT ">$count_clusters/$unique_sequences{$seq}/$per/$len/notChecked\n$seq\n";
				$ids{$seq} = "$count_clusters/$unique_sequences{$seq}/$per/$len/notChecked";
			}
			print INFO "$ids{$seq}\t$details\t$seq\n";
		}
	}
	close (FASTA_OUT);
	close (FASTA_OUT_SELECTED);
	close (SINGLES);  	
}

print "Total cluster: $count_clusters\n";
print "Total singletons: $count_singles\n";

open (STATS, ">$work_dir/$sample.clusters.stats.tsv") or die "Cannot write $work_dir/$sample.clusters.stats.tsv\n";
print STATS "$sample\t$total_primer_reads\t$total_primer_reads_flash\t$total_primer_reads_single\t$count_clusters\t$count_singles\t$count_chimera_reads\t$count_variants_reads\t$count_low_reads\t$count_gap_reads\t$count_selected_reads\t$count_chimera\t$count_variants\t$count_low\t$count_gap\t$count_selected\n";
close(STATS);


$localtime = localtime;
print "\n\nPrimer sorting: finished at $localtime\n";







