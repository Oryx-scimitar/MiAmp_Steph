use warnings;
no warnings ('uninitialized', 'substr');

my $samplesheet = shift;
my $run = shift;

open (OUT_SUMMARY, ">$run.display_table.txt");
open (OUT_READS_DIST, ">$run.reads_distribution.txt");
open (OUT_COUNTS, ">$run.summary_read_counts.txt");
open (OUT_PERCENT, ">$run.summary_read_percent.txt");
open (SUMMARY, ">$run.summary.txt");
open (SUMMARY1, ">$run.summary.discarded.txt");

print OUT_SUMMARY "Library ID\tSample ID\tPCR primers\tTotal sequencing reads\tQuality passed paired-end reads\tQuality passed single-end read1\tQuality passed single-end read2\tQuality failed reads\tOverlapped pair-end reads\tNon-overlapped reads\tReads with PCR primers\tPCR primers in overlapped reads\tPCR primers in single reads\tReads failed to identify primers\tSingle copy sequences\tNumner of OTU with 100% identity\tChimeric reads\t1bp mismatch with higher abundant sequence\tLow read counts\tAmbiguous basecall\tFiltered reads\tExpected filtered\tUnexpected filtered\tExpected removed\tUnexpected removed\tUnknown\n";
print OUT_READS_DIST "Sample ID\tQuality failed reads\tNon-overlapped reads\tReads failed to identify primers\tSingle copy sequences\tChimeric reads\t1bp mismatch with higher abundant sequence\tLow read counts\tReads with Ambiguous basecalls\tExpected filtered\tUnexpected filtered\tExpected removed\tUnexpected removed\tUnknown\n";

print OUT_COUNTS "Library_ID\tSample_ID\tPCR_primers\tTotal_sequencing_reads\tQuality_passed_paired-end_reads\tQuality_passed_single-end_read1\tQuality_passed_single-end_read2\tQuality_failed_reads\tOverlapped_pair-end_reads\tNon-overlapped_reads\tReads with PCR primers\tPCR_primers_in_overlapped_reads\tPCR_primers_in_single_reads\tReads_failed_to_identify_primers\tSingle_copy_sequences\tChimeric_reads\t1bp_mismatch_with_higher_abundant_sequence\tLow_read_counts\tAmbiguous_basecall\tExpected filtered\tUnexpected filtered\tExpected removed\tUnexpected removed\tUnknown\n";
print OUT_PERCENT "Library_ID\tSample_ID\tPCR_primers\tTotal_sequencing_reads\tQuality_passed_paired-end_reads\tQuality_passed_single-end_read1\tQuality_passed_single-end_read2\tQuality_failed_reads\tOverlapped_pair-end_reads\tNon-overlapped_reads\tReads with PCR primers\tPCR_primers_in_overlapped_reads\tPCR_primers_in_single_reads\tReads_failed_to_identify_primers\tSingle_copy_sequences\tChimeric_reads\t1bp_mismatch_with_higher_abundant_sequence\tLow_read_counts\tAmbiguous_basecall\tExpected filtered\tUnexpected filtered\tExpected removed\tUnexpected removed\tUnknown\n";
print SUMMARY "Library_ID\tSample_ID\tPCR_primers\tTotal_sequencing_reads\tReads with PCR primers\n";
print SUMMARY1 "Library_ID\tSample_ID\tPCR_primers\tTotal_sequencing_reads\tReads with PCR primers\n";
open (SAMPLESHEET, "$samplesheet") or die "Cannot open uploaded $samplesheet. Try again.\n";
while(<SAMPLESHEET>){
	chomp $_;
	@words = split("\t", $_);
	my $sample = $words[2];
	my $primer = $words[3];
	my $library_id = $words[1];
	my $folder = $words[0];
	# print "$sample\t$primer\n";
	print OUT_SUMMARY "$library_id\t$sample\t$primer";
	print OUT_READS_DIST "$library_id";
	print OUT_COUNTS "$library_id\t$sample\t$primer";
	print OUT_PERCENT "$library_id\t$sample\t$primer";

	$file1 = "total_reads.txt";
	my $total_reads = 0;
	if (-e $file1){
		open(IN, "$file1");
		LOOP0: while (<IN>){
			chomp $_;
			@words = split("\t", $_);
			if ($library_id eq $words[0]){
				$total_reads = $words[1];
				last LOOP0;
			}
		}
		close(IN);
	}
	else{
		print "Cannot find $file1\n";
	}

	print OUT_SUMMARY "\t$total_reads";
	print OUT_COUNTS "\t$total_reads";
	print OUT_PERCENT "\t$total_reads";

	$sickle_log = "results/$folder/$folder.trimming_by_sickle.log";
	if (-e $sickle_log){
		open (IN, "$sickle_log");
		LOOP1: while (<IN>){
			chomp $_;
			if (/^FastQ paired records kept: (\d+) \((\d+) pairs/){
				$trimmed_paired = $2;
			}
			elsif (/^FastQ single records kept: (\d+) \(from PE1: (\d+), from PE2: (\d+)/){
				$trimmed_single_r1 = $2;
				$trimmed_single_r2 = $3;
				last LOOP1;
			}
		}
		close(IN);
		$removed_by_trimming = $total_reads - ($trimmed_paired + $trimmed_single_r1 + $trimmed_single_r2);

		$removed_by_trimming_per = sprintf('%.2f', (($removed_by_trimming/$total_reads)*100));
		$trimmed_paired_per = sprintf('%.2f', (($trimmed_paired/$total_reads)*100));
		$trimmed_single_r1_per = sprintf('%.2f', (($trimmed_single_r1/$total_reads)*100));
		$trimmed_single_r2_per = sprintf('%.2f', (($trimmed_single_r2/$total_reads)*100));
	}
	else{
		print "Cannot find $sickle_log\n";
	}

	print OUT_SUMMARY "\t$trimmed_paired ($trimmed_paired_per%)\t$trimmed_single_r1 ($trimmed_single_r1_per%)\t$trimmed_single_r2 ($trimmed_single_r2_per%)\t$removed_by_trimming($removed_by_trimming_per%)";
	print OUT_COUNTS "\t$trimmed_paired\t$trimmed_single_r1\t$trimmed_single_r2\t$removed_by_trimming";
	print OUT_PERCENT "\t$trimmed_paired_per\t$trimmed_single_r1_per\t$trimmed_single_r2_per\t$removed_by_trimming_per";

	print OUT_READS_DIST "\t$removed_by_trimming";

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

		$overlapped_paired_per = sprintf('%.2f', (($overlapped_paired/$total_reads)*100));
		$removed_by_overlap_per = sprintf('%.2f', (($removed_by_overlap/$total_reads)*100));

	}
	else{
		print "Cannot find $flash_log\n";
	}

	print OUT_SUMMARY "\t$overlapped_paired ($overlapped_paired_per%)\t$removed_by_overlap ($removed_by_overlap_per%)";
	print OUT_COUNTS "\t$overlapped_paired\t$removed_by_overlap";
	print OUT_PERCENT "\t$overlapped_paired_per\t$removed_by_overlap_per";

	print OUT_READS_DIST "\t$removed_by_overlap";

	$sort_primer = "results/$folder/$folder.$primer.sort.stats.txt";
	if (-e $sort_primer){
		open (IN, "$sort_primer");
		while (<IN>){
			chomp $_;
			@words = split("\t", $_);
			$total_primer = $words[2];
			$primers_overlap = $words[3];
			$primers_single = $words[4];
			$clusters = $words[5];
			$singles = $words[6];
		}
		close(IN);
		$no_primer = ($overlapped_paired + $trimmed_single_r1 + $trimmed_single_r2) - ($total_primer);

		$total_primer_per = sprintf('%.2f', (($total_primer/$total_reads)*100));
		$primers_overlap_per = sprintf('%.2f', (($primers_overlap/$total_reads)*100));
		$primers_single_per = sprintf('%.2f', (($primers_single/$total_reads)*100));
		$no_primer_per = sprintf('%.2f', (($no_primer/$total_reads)*100));
		$singles_per = sprintf('%.2f', (($singles/$total_reads)*100));
	}
	else{
		print "Cannot find $sort_primer\n";
	}

	print OUT_SUMMARY "\t$total_primer ($total_primer_per%)\t$primers_overlap ($primers_overlap_per%)\t$primers_single ($primers_single_per%)\t$no_primer ($no_primer_per%)\t$singles ($singles_per%)\t$clusters";
	print OUT_COUNTS "\t$total_primer\t$primers_overlap\t$primers_single\t$no_primer\t$singles";
	print OUT_PERCENT "\t$total_primer\t$primers_overlap_per\t$primers_single_per\t$no_primer_per\t$singles_per";


	print OUT_READS_DIST "\t$no_primer\t$singles";

	# print STATS "$sample\t$count_chimera_reads\t$count_variants_reads\t$count_low_reads\t$count_gap_reads\t$count_selected_reads\t$count_chimera\t$count_variants\t$count_low\t$count_gap\t$count_selected\n";

	$filter_log = "results/$folder/$folder.$primer.clusters.stats.txt";
	if (-e $filter_log){
		open (IN, "$filter_log");
		while (<IN>){
			chomp $_;
			@words = split("\t", $_);
			$chimera_reads = $words[1];
			$variants_reads = $words[2];
			$low_reads = $words[3];
			$gap_reads = $words[4];
			$selected_reads = $words[5];

			$chimera = $words[6];
			$variants = $words[7];
			$low = $words[8];
			$gap = $words[9];
			$selected = $words[10];
		}
		close(IN);

		$chimera_reads_per = sprintf('%.2f', (($chimera_reads/$total_reads)*100));
		$variants_reads_per = sprintf('%.2f', (($variants_reads/$total_reads)*100));
		$low_reads_per = sprintf('%.2f', (($low_reads/$total_reads)*100));
		$gap_reads_per = sprintf('%.2f', (($gap_reads/$total_reads)*100));
		$selected_reads_per = sprintf('%.2f', (($selected_reads/$total_reads)*100));
	}
	else{
		print "Cannot find $filter_log\n";
	}

	print OUT_SUMMARY "\t$chimera_reads ($chimera_reads_per%) [$chimera]";
	print OUT_SUMMARY "\t$variants_reads ($variants_reads_per%) [$variants]";
	print OUT_SUMMARY "\t$low_reads ($low_reads_per%) [$low]";
	print OUT_SUMMARY "\t$gap_reads ($gap_reads_per%) [$gap]";
	print OUT_SUMMARY "\t$selected_reads ($selected_reads_per%) [$selected]";

	print OUT_COUNTS "\t$chimera_reads";
	print OUT_COUNTS "\t$variants_reads";
	print OUT_COUNTS "\t$low_reads";
	print OUT_COUNTS "\t$gap_reads";

	print OUT_PERCENT "\t$chimera_reads_per";
	print OUT_PERCENT "\t$variants_reads_per";
	print OUT_PERCENT "\t$low_reads_per";
	print OUT_PERCENT "\t$gap_reads_per";

	print OUT_READS_DIST "\t$chimera_reads\t$variants_reads\t$low_reads\t$gap_reads";




	open(SUMMARY_FULL, ">summary/$sample.$library_id.$primer.summary.full.txt");
	open(SUMMARY_SHORT, ">summary/$sample.$library_id.$primer.summary.short.txt");
	open(SUMMARY_SELECTED, ">summary/$sample.$library_id.$primer.summary.filtered.txt");
	open(SUMMARY_DISCARDED, ">summary/$sample.$library_id.$primer.summary.discarded.txt");
	open(FASTA1, ">summary/$sample.$library_id.$primer.filtered.fasta");
	open(FASTA2, ">summary/$sample.$library_id.$primer.discarded.fasta");

	print SUMMARY_FULL "Sample\tLibrary\tPCR\tCluster ID\tSeq type\tExpected pathogen\tRead counts\tRead percent\tIdentity\tQuery cover\tAmbiguous nucleotide:Mismatch:Gap\tSize\tNCBI\tSpecies\tTaxonomy\tSequence\n";
	print SUMMARY_SHORT "Sample\tLibrary\tPCR\tCluster_ID\tSeq_type\tExpected_pathogen\tRead_counts\tRead_percent\tIdentity\tQuery_cover\tAmbiguous_nucleotide:Mismatch:Gap\tSize\tNCBI\tSpecies\tTaxonomy\tSequence\n";
	print SUMMARY_SELECTED "Sample\tLibrary\tPCR\tCluster_ID\tSeq_type\tExpected_pathogen\tRead_counts\tRead_percent\tIdentity\tQuery_cover\tAmbiguous_nucleotide:Mismatch:Gap\tSize\tNCBI\tSpecies\tTaxonomy\tSequence\n";
	print SUMMARY_DISCARDED "Sample\tLibrary\tPCR\tCluster_ID\tSeq_type\tExpected_pathogen\tRead_counts\tRead_percent\tIdentity\tQuery_cover\tAmbiguous_nucleotide:Mismatch:Gap\tSize\tNCBI\tSpecies\tTaxonomy\tSequence\n";
	print SUMMARY "$sample\t$library_id\t$primer\t$total_reads\t$total_primer";
	print SUMMARY1 "$sample\t$library_id\t$primer\t$total_reads\t$total_primer";

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
		$search = "Anaplasma";
		$min_len = 160;
		$max_len = 163;
		$identity = 99;
	}
	elsif ($primer eq "ThBa"){
		$cutoff = 500;
		$search = "Piroplasmorida";
		$min_len = 329;
		$max_len = 378;
		# $min_len = 337;
		# $max_len = 337;
		$identity = 98;
	}
	elsif ($primer eq "Tryp"){
		$cutoff = 500;
		$search = "Trypanosoma";
		$min_len = 136;
		$max_len = 164;
		$identity = 96;
	}
	$in_cluster = "results/$folder/$folder.$primer.clusters.blast.details.txt";
	open(TAX, "$in_cluster");
	while(<TAX>){
		chomp $_;
		@words = split("\t", $_);
		@info = split (/\//, $words[0]);
		if ($info[4] =~ /filtered/ or $info[4] =~ /ambiguous/){
			@sps = split(",", $words[7]);
			@tax = split(",", $words[8]);
			@ncbi = split(",", $words[1]);
			@ncbi_s = sort { lc($a) cmp lc($b) } @ncbi;
			if (scalar(@sps) > 1){
				$temp = scalar(@sps) - 1;
				$species = "$sps[0]($ncbi[0])($temp)";
				$taxonomy = "$tax[0]($temp)";
				$nc = "$ncbi_s[0]($temp)";
			}
			else{
				$species = "$words[7]($words[1])";
				$taxonomy = $words[8];
				$nc = $words[1];
			}
			if ($info[4] =~ /ambiguous/){
				$ambiguous = "yes";
			}
			else{
				$ambiguous = "no";
			}
			$len = length($words[9]);
			# print "$info[1] $cutoff $words[2] $identity $words[3] $info[3] $min_len $info[3] $max_len\n";
			if ($words[1] ne "unknown"){
				if ($info[1] > $cutoff and $words[2] > $identity and $words[3] =~ /100/ and ($info[3] >= $min_len and $info[3] <= $max_len)){
					if ($_ =~ /$search/){
						$exp_above_cutoff++;
						$exp_above_cutoff_reads = $exp_above_cutoff_reads + $info[1];
						print SUMMARY_FULL "$sample\t$library_id\t$primer\t$words[0]\tfiltered\tyes\t$info[1]\t$info[2]\t$words[2]\t$words[3]\t$ambiguous:$words[4]\t$len\t$words[1]\t$words[7]\t$words[8]\t$words[9]\n";
						print SUMMARY_SHORT "$sample\t$library_id\t$primer\t$words[0]\tfiltered\tyes\t$info[1]\t$info[2]\t$words[2]\t$words[3]\t$ambiguous:$words[4]\t$len\t$nc\t$species\t$taxonomy\t$words[9]\n";
						print SUMMARY_SELECTED "$sample\t$library_id\t$primer\t$words[0]\tfiltered\tyes\t$info[1]\t$info[2]\t$words[2]\t$words[3]\t$ambiguous:$words[4]\t$len\t$nc\t$species\t$taxonomy\t$words[9]\n";
						print FASTA1 ">$nc $species $info[1] $info[2] $words[2] $len $sample:$library_id:$primer\n$words[9]\n";
						print SUMMARY "\t$species($info[1])";
					}
					else{
						$unexp_above_cutoff++;
						$unexp_above_cutoff_reads = $unexp_above_cutoff_reads + $info[1];
						print SUMMARY_FULL "$sample\t$library_id\t$primer\t$words[0]\tfiltered\tno\t$info[1]\t$info[2]\t$words[2]\t$words[3]\t$ambiguous:$words[4]\t$len\t$words[1]\t$words[7]\t$words[8]\t$words[9]\n";
						print SUMMARY_SHORT "$sample\t$library_id\t$primer\t$words[0]\tfiltered\tno\t$info[1]\t$info[2]\t$words[2]\t$words[3]\t$ambiguous:$words[4]\t$len\t$nc\t$species\t$taxonomy\t$words[9]\n";
						print SUMMARY_SELECTED "$sample\t$library_id\t$primer\t$words[0]\tfiltered\tno\t$info[1]\t$info[2]\t$words[2]\t$words[3]\t$ambiguous:$words[4]\t$len\t$nc\t$species\t$taxonomy\t$words[9]\n";
						print FASTA1 ">$species $info[1] $info[2] $words[2] $len $sample:$library_id:$primer\n$words[9]\n";
						print SUMMARY "\t$species($info[1])";
					}
				}
				else{
					if ($_ =~ /$search/){
						$exp_below_cutoff++;
						$exp_below_cutoff_reads = $exp_below_cutoff_reads + $info[1];
						print SUMMARY_FULL "$sample\t$library_id\t$primer\t$words[0]\tremoved\tyes\t$info[1]\t$info[2]\t$words[2]\t$words[3]\t$ambiguous:$words[4]\t$len\t$words[1]\t$words[7]\t$words[8]\t$words[9]\n";
						print SUMMARY_SHORT "$sample\t$library_id\t$primer\t$words[0]\tremoved\tyes\t$info[1]\t$info[2]\t$words[2]\t$words[3]\t$ambiguous:$words[4]\t$len\t$nc\t$species\t$taxonomy\t$words[9]\n";
						print SUMMARY_DISCARDED "$sample\t$library_id\t$primer\t$words[0]\tremoved\tyes\t$info[1]\t$info[2]\t$words[2]\t$words[3]\t$ambiguous:$words[4]\t$len\t$nc\t$species\t$taxonomy\t$words[9]\n";
						print FASTA2 ">$species $info[1] $info[2] $words[2] $len $sample:$library_id:$primer\n$words[9]\n";
						print SUMMARY1 "\t$species($info[1])";
					}
					else{
						$unexp_below_cutoff++;
						$unexp_below_cutoff_reads = $unexp_below_cutoff_reads + $info[1];
						print SUMMARY_FULL "$sample\t$library_id\t$primer\t$words[0]\tremoved\tno\t$info[1]\t$info[2]\t$words[2]\t$words[3]\t$ambiguous:$words[4]\t$len\t$words[1]\t$words[7]\t$words[8]\t$words[9]\n";
						print SUMMARY_SHORT "$sample\t$library_id\t$primer\t$words[0]\tremoved\tno\t$info[1]\t$info[2]\t$words[2]\t$words[3]\t$ambiguous:$words[4]\t$len\t$nc\t$species\t$taxonomy\t$words[9]\n";
						print SUMMARY_DISCARDED "$sample\t$library_id\t$primer\t$words[0]\tremoved\tno\t$info[1]\t$info[2]\t$words[2]\t$words[3]\t$ambiguous:$words[4]\t$len\t$nc\t$species\t$taxonomy\t$words[9]\n";
						print FASTA2 ">$species $info[1] $info[2] $words[2] $len $sample:$library_id:$primer\n$words[9]\n";
						print SUMMARY1 "\t$species($info[1])";
					}
				}
			}
			else{
				$unknown++;
				$unknown_reads = $unknown_reads + $info[1];
				print SUMMARY_FULL "$sample\t$library_id\t$primer\t$words[0]\tremoved\tno\t$info[1]\t$info[2]\t$words[2]\t$words[3]\t$ambiguous:$words[4]\t$len\t$words[1]\tunknown\ttunknown\t$words[7]\n";
				print SUMMARY_SHORT "$sample\t$library_id\t$primer\t$words[0]\tremoved\tno\t$info[1]\t$info[2]\t$words[2]\t$words[3]\t$ambiguous:$words[4]\t$len\t$nc\tunknown\ttunknown\t$words[7]\n";
				print SUMMARY_DISCARDED "$sample\t$library_id\t$primer\t$words[0]\tremoved\tno\t$info[1]\t$info[2]\t$words[2]\t$words[3]\t$ambiguous:$words[4]\t$len\t$words[1]\ttunknown\ttunknown\t$words[7]\n";
				print FASTA2 ">Unknown $info[1] $info[2] $words[2] $len $sample:$library_id:$primer\n$words[7]\n";
				print SUMMARY1 "\tunknown($info[1])";
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

	print OUT_COUNTS "\t$exp_above_cutoff_reads\t$unexp_above_cutoff_reads\t$exp_below_cutoff_reads\t$unexp_below_cutoff_reads\t$unknown_reads";
	print OUT_PERCENT "\t$exp_above_cutoff_reads_per\t$unexp_above_cutoff_reads_per\t$exp_below_cutoff_reads_per\t$unexp_below_cutoff_reads_per\t$unknown_per";
	print OUT_READS_DIST "\t$exp_above_cutoff_reads\t$unexp_above_cutoff_reads\t$exp_below_cutoff_reads\t$unexp_below_cutoff_reads\t$unknown_reads";

	print OUT_SUMMARY "\n";
	print OUT_COUNTS "\n";
	print OUT_PERCENT "\n";
	print OUT_READS_DIST "\n";
	print SUMMARY "\n";
	print SUMMARY1 "\n";
}
