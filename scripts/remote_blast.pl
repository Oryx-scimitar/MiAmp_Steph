use warnings;
no warnings ('uninitialized', 'substr');
use Cwd;
use IPC::Cmd qw[can_run run];
use Getopt::Long;

my $samplesheet = shift;

# system ("cat results/*/*.clusters.filtered.fasta > fasta/all_clusters.fasta");
# system ("blastn -remote -db nr -query fasta/all_clusters.fasta -outfmt '6 qseqid sacc pident length qlen qstart qend slen sstart send mismatch gapopen evalue bitscore sblastname ssciname scomname stitle' -out fasta/all_clusters.blast -max_target_seqs 1")

my $id = "";
my %refs = ();
my %identity = ();
my %error = ();
my %gap = ();
my %match = ();
my %q_len = ();
my %q_start = ();
my %q_end = ();
my %r_len = ();
my %r_start = ();
my %r_end = ();
my %references = ();

$blast_file = "fasta/all_clusters.blast";
if (-e $blast_file) {
	open (BLAST, "fasta/all_clusters.blast");
	while (<BLAST>){
		chomp $_;
		my @words = split("\t", $_);
# RED000938/6/4845/1.81/368/ambiguous_basecall	JN572700	99.728	368	368	1	368	1580	419	786	1	0	0.0	676	N/A	N/A	N/A	Theileria sp. B15a 18S ribosomal RNA gene, partial sequence
		if ($words[0] ne $id){
			$id = $words[0];
			$query_cover{$id} = sprintf('%.2f', (($words[3]/$words[4])*100));
			$refs{$id} = $words[1];
			$identity{$id} = $words[2];
		    $error{$id} = $words[10];
		    $gap{$id} = $words[11];
		    $match{$id} = $words[3];
		    $q_len{$id} = $words[4];
		    $q_start{$id} = $words[5];
		    $q_end{$id} = $words[6];
		    $species{$id} = $words[17];
		}
	}
	close(BLAST);
}
else{
	print "Cannot find $blast_file\n\n";
}

open (SAMPLESHEET, "$samplesheet") or die "Cannot open uploaded $samplesheet. Try again.\n";
while(<SAMPLESHEET>){
	chomp $_;
	@words = split("\t", $_);
	my $sample = $words[0];
	open (OUT, ">summary/$sample.$words[2].nr.blast.tsv");
	print OUT "sample\tcluster_id\treference\tspecies\tpercent_identity\tquery_cover\tmatch_len\tmismatch\tgap\tquery_len\tquery_start_end\n";

	$fasta_file = "results/$sample/$sample.clusters.filtered.fasta";
	if (-e $fasta_file) {
		open(IN, "results/$sample/$sample.clusters.filtered.fasta");
		while(<IN>){
			chomp $_;
			if (/^>(\S+)/){
				$id = $1;
			}
			else{
				if (exists $refs{$id}){
					print OUT "$sample\t$id\t$refs{$id}\t$species{$id}\t$identity{$id}\t$query_cover{$id}\t$match{$id}\t$error{$id}\t$gap{$id}\t$q_len{$id}\t$q_start{$id}-$q_end{$id}\n";
				}
			}
		}
	}
	else{
		print "Cannot find $fasta_file\n";
	}
	close(OUT);
}



