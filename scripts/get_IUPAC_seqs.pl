use Bio::SeqIO;
use Bio::Root::Exception;
use Bio::PrimarySeq;
use Bio::Tools::IUPAC;

$primers = shift;

my $in = Bio::SeqIO -> new( -file => "$primers", -format => 'fasta');
while ( $seq = $in->next_seq() ) {
	# Create a sequence with degenerate residues
	my $ambiseq = Bio::PrimarySeq->new(-seq => $seq->seq(), -alphabet => 'dna');
	$id = $seq->id;
	# Create all possible non-degenerate sequences
	my $iupac = Bio::Tools::IUPAC->new(-seq => $ambiseq);
	while ($uniqueseq = $iupac->next_seq()) {
		$sequence = $uniqueseq->seq();
		print ">$id.1\n$sequence\n";
	}
}