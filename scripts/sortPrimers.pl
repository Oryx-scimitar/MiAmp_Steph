use warnings;
no warnings ('uninitialized', 'substr');
use Cwd;
# use IPC::Cmd qw[can_run run];
# use Getopt::Long;
# use Bio::SeqIO;
# use Bio::Root::Exception;
# use Bio::PrimarySeq;
# use Bio::Tools::IUPAC;

my $sample = shift;
my $work_dir = shift;
my $prefix = shift;
my $primers = shift;

print "$sample\t$prefix\n\n";

my %forward_primer = ();
my %reverse_primer = ();
my %forward_primer_revcomp = ();
my %reverse_primer_revcomp = ();
my $count_for = 0;
my $count_rev = 0;

open (IN, "$primers") or print "Cannot read/find $primers\n";
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

# my $in = Bio::SeqIO -> new( -file => "$primers", -format => 'fasta');
# while ( $seq = $in->next_seq() ) {
# 	# Create a sequence with degenerate residues
# 	my $ambiseq = Bio::PrimarySeq->new(-seq => $seq->seq(), -alphabet => 'dna');

# 	# Create all possible non-degenerate sequences
# 	my $iupac = Bio::Tools::IUPAC->new(-seq => $ambiseq);
# 	while ($uniqueseq = $iupac->next_seq()) {
# 		# process the unique Bio::Seq object.
# 	 	my $revcomp = $uniqueseq->revcom;
# 		my $id = $seq->id;
# 		if ($id =~ /for/i){
# 			$for++;
# 			$forward_primer{$for} = $uniqueseq->seq();
# 			$forward_primer_revcomp{$for} = $revcomp->seq();
# 		} elsif ($id =~ /rev/i){
# 			$rev++;
# 			$reverse_primer_revcomp{$rev} = $revcomp->seq();
# 			$reverse_primer{$rev} = $uniqueseq->seq();
# 		} else {
# 			die "$id in $primers is not valid primer id Read the user manual to see the correct format of file and Try again\n";
# 		}
# 	}
# }
print "Forward primers...\n";
foreach my $p (sort {$a <=> $b} keys %forward_primer){
	print "$p $forward_primer{$p} $forward_primer_revcomp{$p}\n";
}

print "Reverse primers....\n";
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
			foreach my $primer2 (keys %reverse_primer_revcomp){
				if ($2 =~ /(\w+|)$reverse_primer_revcomp{$primer2}(\w+|)$/) {
					if (length($1) > 50){
						if (exists $unique_sequences{$1}){
							$unique_sequences{$1} = $unique_sequences{$1} + 1;
						}
						else {
							$unique_sequences{$1} = 1;
						}
						$flag = 1;
						$primer_group{"$sample\t$forward_primer{$primer1}\t$reverse_primer_revcomp{$primer2}\t+"}++;
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
				foreach my $primer2 (keys %forward_primer_revcomp){
					if ($2 =~ /(\w+|)$forward_primer_revcomp{$primer2}(\w+|)$/) {
						if (length($1) > 50){
							my $seq_revcomp = reverse $1;
							$seq_revcomp =~ tr/ATGCatgc/TACGtacg/;
							if (exists $unique_sequences{$seq_revcomp}){
								$unique_sequences{$seq_revcomp} = $unique_sequences{$seq_revcomp} + 1;
							} else {
								$unique_sequences{$seq_revcomp} = 1;
							}
							$flag = 1;
							$primer_group{"$sample\t$forward_primer_revcomp{$primer2}\t$reverse_primer{$primer1}\t-"}++;
							$trimmed_seq = $seq_revcomp;
							# print "original: $1\nconverted: $seq_revcomp\n";
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
open (OUT_PRIMERS, ">$work_dir/$sample.$prefix.fasta") or die "Cannot write $work_dir/$sample/$sample.$prefix.fasta\n";
open (OUT_PRIMERS_LEN, ">$work_dir/$sample.$prefix.len.hist") or die "Cannot write $work_dir/$sample/$sample.$prefix.fasta\n";
my %flash_len = ();


if ( -e $flash_reads ){
	# $in = Bio::SeqIO -> new( -file => "$flash_reads", -format => 'fastq');
	# print "Identifiying primers and writing unique variants from $flash_reads...\n\n";
	# LOOP5: while ( $seq = $in->next_seq() ) {
	# 	$sequence = $seq->seq();
	# 	$id = $seq->id;
	# 	my $flag_primer = check_primer($sequence);
	# 	if ($flag_primer eq 1){
	# 		$total_primer_reads_flash++;
	# 		print OUT_PRIMERS ">$id/$sample/extended\n$trimmed_seq\n";
	# 	}
	# }

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
				$total_primer_reads++;
				$len = length($trimmed_seq);
				$flash_len{$len}++;
				print OUT_PRIMERS "$id/$sample/extended/$len\n$trimmed_seq\n";
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
		

foreach my $len (sort{$a <=> $b} keys %flash_len){
	print OUT_PRIMERS_LEN "$len\t$flash_len{$len}\n";
}
open (INFO, ">$work_dir/$sample.$prefix.primers.info.txt") or die "Cannot write $work_dir/$sample.$prefix.primers.info.txt\n";
print "Writting unique variants log $work_dir/$sample.primers.info.txt...\n";
foreach $primer (sort {$primer_group{$b} <=> $primer_group{$a}} keys %primer_group){
	print INFO "$primer\t$primer_group{$primer}\n";
}
close (INFO);

my $count_clusters = 0;
my $count_singles = 0;
my $count = 0;

if($total_primer_reads > 0){
	
	open (OUT_UNIQUE, ">$work_dir/$sample.$prefix.clusters.fasta") or die "Cannot write $work_dir/$sample/$sample.$prefix.clusters.fasta\n";
	open (SINGLES, ">$work_dir/$sample.$prefix.singletons.fasta") or die "Cannot write $work_dir/$sample/$sample.$prefix.singletons.fasta\n";
	
	foreach $seq (sort {$unique_sequences{$b} <=> $unique_sequences{$a}} keys %unique_sequences){
		$count++;
		my $per = sprintf("%.5f", ($unique_sequences{$seq} / $total_primer_reads) * 100);
		if ($unique_sequences{$seq} eq 1){
			$count_singles++;
			print SINGLES ">$count:$unique_sequences{$seq}:$per\n$seq\n";
		}
		else{
			$count_clusters++;
			print OUT_UNIQUE ">$count:$unique_sequences{$seq}:$per\n$seq\n";
		}
	}

	close (OUT_UNIQUE);
	close (SINGLES);  	
}

print "Total sequences with $prefix primers in overlapped reads: $total_primer_reads\n";
print "Total cluster: $count_clusters\n";
print "Total singletons: $count_singles\n";

open (STAT, ">", "$work_dir/$sample.$prefix.sort.stats.txt");
print STAT "$sample\t$prefix\t$total_primer_reads\t$count_clusters\t$count_singles\n";
close (STAT);





