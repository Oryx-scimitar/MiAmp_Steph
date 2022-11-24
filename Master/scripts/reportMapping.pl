use warnings;
no warnings ('uninitialized', 'substr');
use Cwd;
use IPC::Cmd qw[can_run run];
use Getopt::Long;

my $sample = "sample";
my $work_dir = "results";
my $blast_file = "";
my $in_taxonomy_details = "";
my $primer = "test";
my $output = "test";

GetOptions(
    'sample=s'    => \$sample,
    'work_dir=s'     => \$work_dir, 
    'blast_file=s'     => \$blast_file,
    # 'taxonomy=s'	=> \$in_taxonomy_details,
    'output=s'	=> \$output,
    'primer=s'	=> \$fc
) or print "Invalid options\n";

# my %taxonomy = ();
# my %species = ();
# if ( -e $in_taxonomy_details ){
# 	open(TAX, "$in_taxonomy_details");
# 	while(<TAX>){
# 		chomp $_;
# 		@words = split("\t", $_);
# 		$taxonomy{$words[0]} = $words[1];
# 		$species{$words[0]} = $words[2];
# 	} 
# }
# else{
# 	print "Cannot open $in_taxonomy_details\n";
# }

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

if (-e $blast_file) {
	open (BLAST, "$blast_file");
	print "Looking for known Alleles from reference...\n";
	while (<BLAST>){
		chomp $_;
		my @words = split("\t", $_);
		# @ref = split(/\./, $words[1]);
		if ($words[0] ne $id){
			$id = $words[0];
			$query_cover{$id} = sprintf('%.2f', (($words[3]/$words[4])*100));
			print "$id\t$query_cover\t$words[2]\t$words[1]\n";
			$refs{$id} = $words[1];
			$identity{$id} = $words[2];
		    $error{$id} = $words[10];
		    $gap{$id} = $words[11];
		    $match{$id} = $words[3];
		    $q_len{$id} = $words[4];
		    $q_start{$id} = $words[5];
		    $q_end{$id} = $words[6];
		    $r_len{$id} = $words[7];
		    $r_start{$id} = $words[8];
		    $r_end{$id} = $words[9];
		}
		elsif ($identity{$id} eq $words[2] and $match{$id} eq $words[3]){
			$refs{$id} = $refs{$id}.",".$words[1];
		}
	}
	close(BLAST);
}
else{
	print "Cannot find $blast_file\n\n";
}

open (LOG, ">$output") or print "Cannot write $output\n";
print LOG "sequence_id\tReference_id\tPercent_identity\tQuery_cover\tAlignment_length\tErrors\tGaps\tQuery length\tQuery_start-Query_end\tsequence\n";

if (-e $blast_file) {
	print "$work_dir/$sample.clusters.filtered.fasta\n";
	open(IN, "$work_dir/$sample.clusters.filtered.fasta");
	while(<IN>){
		chomp $_;
		if (/^>(\S+)/){
			$id = $1;
		}
		else{
			$sequence = $_;
			my @info = split(":", $id);

			if (exists $refs{$id}){
				@order_ref = split(",", $refs{$id});
				@ordered_ref = sort { $a cmp $b } @order_ref;
				$to_write_ref = "";
				$to_write_taxonomy = "";
				$to_write_species = "";
				foreach $ref (@ordered_ref){
					$to_write_ref = $to_write_ref."$ref,";
					$to_write_taxonomy = $to_write_taxonomy."$taxonomy{$ref}, ";
					$to_write_species = $to_write_species."$species{$ref}, ";
				}
				chop $to_write_ref;
				chop $to_write_ref;
				chop $to_write_species;
				chop $to_write_species;
				chop $to_write_taxonomy;
				chop $to_write_taxonomy;

				print LOG "$id\t$to_write_ref\t$identity{$id}\t$query_cover{$id}\t$match{$id}\t$error{$id}\t$gap{$id}\t$q_len{$id}\t$q_start{$id}-$q_end{$id}\t$sequence\n";
			}
			else{
				print LOG "$id\tunknown\t-\t-\t-\t-\t-\t".length($sequence)."\t-\t$sequence\n";
			}
		}
	}
}

