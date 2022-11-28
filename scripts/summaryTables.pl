use warnings;
no warnings ('uninitialized', 'substr');

my $samplesheet = shift;

open (OUT_SUMMARY, ">summary_read_counts.tsv");
open (LIST_SELECTED, ">summary.txt");
open (LIST_DISCARDED, ">summary.discarded.txt");

print OUT_SUMMARY "Library ID\tSample ID\tPCR primers\tTotal sequencing reads\tQuality passed paired-end reads\tQuality passed single-end read1\tQuality passed single-end read2\tQuality failed reads\tOverlapped pair-end reads\tNon-overlapped reads\tPCR primers in overlapped reads\tPCR primers in single reads\tTotal Reads with PCR primers\tReads failed to identify primers\tSingle copy sequences\tNumber of OTU with 100\% identity\tChimeric reads\t1bp mismatch with higher abundant sequence\t<50bp short sequences\tLow read counts\tAmbiguous basecall\tFiltered reads\tExpected filtered\tUnexpected filtered\tExpected removed\tUnexpected removed\tUnknown\n";

print LIST_SELECTED "Library_ID\tSample_ID\tPCR_primers\tTotal_sequencing_reads\tReads with PCR primers\n";
print LIST_DISCARDED "Library_ID\tSample_ID\tPCR_primers\tTotal_sequencing_reads\tReads with PCR primers\n";

open (SAMPLESHEET, "$samplesheet") or die "Cannot open uploaded $samplesheet. Try again.\n";
while(<SAMPLESHEET>){
	chomp $_;
	@words = split("\t", $_);
	my $sample = $words[1];
	my $primer = $words[2];
	my $library_id = $words[0];
	my $folder = $words[0];
	# print "$sample\t$primer\n";
	print OUT_SUMMARY "$library_id\t$sample\t$primer";

	# $file1 = "total_reads.txt";
	# my $total_reads = 0;
	# if (-e $file1){
	# 	open(IN, "$file1");
	# 	LOOP0: while (<IN>){
	# 		chomp $_;
	# 		@words = split("\t", $_);
	# 		if ($library_id eq $words[1]){
	# 			$total_reads = $words[2];
	# 			last LOOP0;
	# 		}
	# 	}
	# 	close(IN);
	# }
	# else{
	# 	print "Cannot find $file1\n";
	# }

	# print OUT_SUMMARY "\t$total_reads";
	$sickle_log = "results/$folder/$folder.trimming_by_sickle.log";
	my $total_reads = 0;
	my $trimmed_paired = 0;
	my $overlapped_paired = 0;

	if (-e $sickle_log){
		open (IN, "$sickle_log");
		LOOP1: while (<IN>){
			chomp $_;
			if (/^FastQ paired records kept: (\d+) \((\d+) pairs/){
				$trimmed_paired = $2;
				
				$total_reads = $total_reads + $2;
			}
			elsif (/^FastQ single records kept: (\d+) \(from PE1: (\d+), from PE2: (\d+)/){
				$total_reads = $total_reads + $1;
				$trimmed_single_r1 = $2;
				$trimmed_single_r2 = $3;
			}
			elsif (/^FastQ paired records discarded: (\d+) \((\d+) pairs/){
				$total_reads = $total_reads + $2;
			}
		}
		close(IN);
		$removed_by_trimming = $total_reads - ($trimmed_paired + $trimmed_single_r1 + $trimmed_single_r2);

		$removed_by_trimming_per = sprintf('%.1f', (($removed_by_trimming/$total_reads)*100));
		$trimmed_paired_per = sprintf('%.1f', (($trimmed_paired/$total_reads)*100));
		$trimmed_single_r1_per = sprintf('%.1f', (($trimmed_single_r1/$total_reads)*100));
		$trimmed_single_r2_per = sprintf('%.1f', (($trimmed_single_r2/$total_reads)*100));
	}
	else{
		print "Cannot find $sickle_log\n";
	}

	print "$trimmed_paired\n";
	print OUT_SUMMARY "\t$total_reads\t$trimmed_paired ($trimmed_paired_per\%)\t$trimmed_single_r1 ($trimmed_single_r1_per\%)\t$trimmed_single_r2 ($trimmed_single_r2_per\%)\t$removed_by_trimming($removed_by_trimming_per\%)";

	$flash_log = "results/$folder/$folder.overlap_by_flash.log";
	if (-e $flash_log){
		open (IN, "$flash_log");
		LOOP2: while (<IN>){
			chomp $_;
			if (/Combined pairs:\s+(\d+)/){
				$overlapped_paired = $1;
				last LOOP2;
			}
		}
		close(IN);
		$removed_by_overlap = $trimmed_paired - $overlapped_paired;

		$overlapped_paired_per = sprintf('%.1f', (($overlapped_paired/$total_reads)*100));
		$removed_by_overlap_per = sprintf('%.1f', (($removed_by_overlap/$total_reads)*100));
		
	}
	else{
		print "Cannot find $flash_log\n";
	}
	
	print OUT_SUMMARY "\t$overlapped_paired ($overlapped_paired_per\%)\t$removed_by_overlap ($removed_by_overlap_per\%)";

	$filter_log = "results/$folder/$folder.clusters.stats.tsv";
	if (-e $filter_log){
		open (IN, "$filter_log");
		while (<IN>){
			chomp $_;
			@words = split("\t", $_);
			$total_primer = $words[1];
			$primers_overlap = $words[2];
			$primers_single = $words[3];
			$clusters = $words[4];
			$singles = $words[5];
			$chimera_reads = $words[6];
			$variants_reads = $words[7];
			$short_reads = $words[8];
			$low_reads = $words[9];
			$gap_reads = $words[10];
			$selected_reads = $words[11];

			$chimera = $words[12];
			$variants = $words[13];
			$short = $words[14];
			$low = $words[15];
			$gap = $words[16];
			$selected = $words[17];
		}
		close(IN);
		$no_primer = ($overlapped_paired + $trimmed_single_r1 + $trimmed_single_r2) - ($total_primer);
		$total_primer_per = sprintf('%.1f', (($total_primer/$total_reads)*100));
		$primers_overlap_per = sprintf('%.1f', (($primers_overlap/$total_reads)*100));
		$primers_single_per = sprintf('%.1f', (($primers_single/$total_reads)*100));
		$no_primer_per = sprintf('%.1f', (($no_primer/$total_reads)*100));
		$singles_per = sprintf('%.1f', (($singles/$total_reads)*100));
		$chimera_reads_per = sprintf('%.1f', (($chimera_reads/$total_reads)*100));
		$variants_reads_per = sprintf('%.1f', (($variants_reads/$total_reads)*100));
		$short_reads_per = sprintf('%.1f', (($short_reads/$total_reads)*100));
		$low_reads_per = sprintf('%.1f', (($low_reads/$total_reads)*100));
		$gap_reads_per = sprintf('%.1f', (($gap_reads/$total_reads)*100));
		$selected_reads_per = sprintf('%.1f', (($selected_reads/$total_reads)*100));
	}
	else{
		print "Cannot find $filter_log\n";
	}

	print OUT_SUMMARY "\t$primers_overlap ($primers_overlap_per\%)\t$primers_single ($primers_single_per\%)\t$total_primer ($total_primer_per\%)\t$no_primer ($no_primer_per\%)\t$singles ($singles_per\%)\t$clusters";
	print OUT_SUMMARY "\t$chimera_reads ($chimera_reads_per\%) [$chimera]";
	print OUT_SUMMARY "\t$variants_reads ($variants_reads_per\%) [$variants]";
	print OUT_SUMMARY "\t$short_reads ($short_reads_per\%) [$short]";
	print OUT_SUMMARY "\t$low_reads ($low_reads_per\%) [$low]";
	print OUT_SUMMARY "\t$gap_reads ($gap_reads_per\%) [$gap]";
	print OUT_SUMMARY "\t$selected_reads ($selected_reads_per\%) [$selected]";

	$in_cluster = "results/$folder/$folder.database.details.tsv";

	open(SUMMARY_SELECTED, ">summary/$sample.$primer.summary.selected.txt");
	open(SUMMARY_DISCARDED, ">summary/$sample.$primer.summary.discarded.txt");
	open(FASTA_SELECTED, ">summary/$sample.$primer.selected.fasta");
	open(FASTA_DISCARDED, ">summary/$sample.$primer.discarded.fasta");

	print SUMMARY_SELECTED "Sample\tLibrary\tPCR\tsequence_id\tReference_id\tPercent_identity\tQuery_cover\tAlignment_length\tErrors\tGaps\tQuery length\tQuery_start-Query_end\tSpecies\tTaxonomy\tsequence\n";
	print SUMMARY_DISCARDED "Sample\tLibrary\tPCR\tsequence_id\tReference_id\tPercent_identity\tQuery_cover\tAlignment_length\tErrors\tGaps\tQuery length\tQuery_start-Query_end\tSpecies\tTaxonomy\tsequence\n";
	print LIST_SELECTED "$sample\t$library_id\t$primer\t$total_reads\t$total_primer";
	print LIST_DISCARDED "$sample\t$library_id\t$primer\t$total_reads\t$total_primer";

	$exp_above_cutoff = 0;
	$exp_below_cutoff = 0;
	$unexp_above_cutoff = 0;
	$unexp_below_cutoff = 0;
	$unknown = 0;

	$exp_above_cutoff_reads = 0;
	$exp_below_cutoff_reads = 0;
	$unexp_above_cutoff_reads = 0;
	$unexp_below_cutoff_reads = 0;
	$unknown_reads = 0;

	if ($primer eq "AnEh"){
		$cutoff = 1000;
		$search1 = "Anaplasma";
		$search2 = "Anaplasma";
		$min_len = 160;
		$max_len = 163;
		$identity = 99;
	}
	elsif ($primer eq "ThBa"){
		$cutoff = 500;
		$search1 = "Theileria";
		$search2 = "Babesia";
		$min_len = 300;
		$max_len = 500;
		$identity = 98;
	}
	elsif ($primer eq "Tryp"){
		$cutoff = 500;
		$search1 = "Trypanosoma";
		$search2 = "Trypanosoma";
		$min_len = 136;
		$max_len = 164;
		$identity = 96;
	}

	open(TAX, "$in_cluster");
	while(<TAX>){
		chomp $_;
		@words = split("\t", $_);
		$mapping = $_;
		# $sample/$count_clusters/$unique_sequences{$seq}/$per/$len/ambiguous_basecall
		@info = split (/\//, $words[0]);
		if ($info[5] =~ /filtered/ or $info[5] =~ /ambiguous/){
			if ($words[1] ne "unknown"){
				if ($info[2] > $cutoff and $words[3] > $identity and $words[4] =~ /100/ and ($info[4] >= $min_len and $info[4] <= $max_len)){
					if ($words[2] =~ /$search1/ or $words[2] =~ /$search2/){
						$exp_above_cutoff++;
						$exp_above_cutoff_reads = $exp_above_cutoff_reads + $info[2];
						print SUMMARY_SELECTED "$sample\t$library_id\t$primer\t$mapping\n";
						print FASTA_SELECTED ">$words[1] $words[2] $info[2] $info[3] $words[3] $words[8] $sample:$library_id:$primer\n$words[10]\n";
						print LIST_SELECTED "\t$words[2]($info[2])";
					}
					else{
						$unexp_above_cutoff++;
						$unexp_above_cutoff_reads = $unexp_above_cutoff_reads + $info[2];
						print SUMMARY_SELECTED "$sample\t$library_id\t$primer\t$mapping\n";
						print FASTA_SELECTED ">$words[1] $words[2] $info[2] $info[3] $words[3] $words[8] $sample:$library_id:$primer\n$words[10]\n";
						print LIST_SELECTED "\t$words[2]($info[2])";
					}
				}
				else{
					if ($words[2] =~ /$search1/ or $words[2] =~ /$search2/){
						$exp_below_cutoff++;
						$exp_below_cutoff_reads = $exp_below_cutoff_reads + $info[2];
						print SUMMARY_DISCARDED "$sample\t$library_id\t$primer\t$mapping\n";
						print FASTA_DISCARDED ">$words[1] $words[2] $info[2] $info[3] $words[3] $words[8] $sample:$library_id:$primer\n$words[10]\n";
						print LIST_DISCARDED "\t$words[2]($info[2])";
					}
					else{
						$unexp_below_cutoff++;
						$unexp_below_cutoff_reads = $unexp_below_cutoff_reads + $info[2];
						print SUMMARY_DISCARDED "$sample\t$library_id\t$primer\t$mapping\n";
						print FASTA_DISCARDED ">$words[1] $words[2] $info[2] $info[3] $words[3] $words[8] $sample:$library_id:$primer\n$words[10]\n";
						print LIST_DISCARDED "\t$words[2]($info[2])";
					}
				}
			}
			else{
				$unknown++;
				$unknown_reads = $unknown_reads + $info[2];
				print SUMMARY_DISCARDED "$sample\t$library_id\t$primer\t$mapping\n";
				print FASTA_DISCARDED ">Unknown - $info[2] $info[3] - - $sample:$library_id:$primer\n$words[10]\n";
				print LIST_DISCARDED "\tunknown($info[2])";
			}
		}
	}

	$exp_above_cutoff_reads_per = sprintf('%.2f', (($exp_above_cutoff_reads/$total_reads)*100));
	$unexp_above_cutoff_reads_per = sprintf('%.2f', (($unexp_above_cutoff_reads/$total_reads)*100));
	$exp_below_cutoff_reads_per = sprintf('%.2f', (($exp_below_cutoff_reads/$total_reads)*100));
	$unexp_below_cutoff_reads_per = sprintf('%.2f', (($unexp_below_cutoff_reads/$total_reads)*100));
	$unknown_per = sprintf('%.2f', (($unknown_reads/$total_reads)*100));

	print OUT_SUMMARY "\t$exp_above_cutoff_reads ($exp_above_cutoff_reads_per) [$exp_above_cutoff]";
	print OUT_SUMMARY "\t$unexp_above_cutoff_reads ($unexp_above_cutoff_reads_per) [$unexp_above_cutoff]";
	print OUT_SUMMARY "\t$exp_below_cutoff_reads ($exp_below_cutoff_reads_per) [$exp_below_cutoff]";
	print OUT_SUMMARY "\t$unexp_below_cutoff_reads ($unexp_below_cutoff_reads_per) [$unexp_below_cutoff]";
	print OUT_SUMMARY "\t$unknown_reads ($unknown_per) [$unknown]";

	print OUT_SUMMARY "\n";
	print SUMMARY_SELECTED "\n";
	print SUMMARY_DISCARDED "\n";
	print LIST_SELECTED "\n";
	print LIST_DISCARDED "\n";
}


