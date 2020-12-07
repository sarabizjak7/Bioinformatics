import pickle
from os import path
from typing import Tuple, Generator, List

from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def load(organism_id: str) -> SeqRecord:
    """Load the NCBI record, use cached files if possible."""
    if not path.exists(path.join("data", f"{organism_id}.pkl.gz")):
        with Entrez.efetch(db="nucleotide", rettype="gb", id=organism_id) as handle:
            record = SeqIO.read(handle, "gb")
            with open(path.join("data", f"{organism_id}.pkl.gz"), "wb") as f:
                pickle.dump(record, f)
    else:
        with open(path.join("data", f"{organism_id}.pkl.gz"), "rb") as f:
            record = pickle.load(f)

    return record


def codons(seq: str) -> Generator[str, None, None]:
    """Walk along the string, three nucleotides at a time. Cut off excess."""
    for i in range(0, len(seq) - 2, 3):
        yield seq[i:i + 3]


def extract_gt_orfs(record, start_codons, stop_codons, validate_cds=True, verbose=False):
    """Extract the ground truth ORFs as indicated by the NCBI annotator in the
    gene coding regions (CDS regins) of the genome.

    Parameters
    ----------
    record: SeqRecord
    start_codons: List[str]
    stop_codons: List[str]
    validate_cds: bool
        Filter out NCBI provided ORFs that do not fit our ORF criteria.
    verbose: bool

    Returns
    -------
    List[Tuple[int, int, int]]
        tuples of form (strand, start_loc, stop_loc). Strand should be either 1
        for reference strand and -1 for reverse complement.

    """
    cds_regions = [f for f in record.features if f.type == "CDS"]

    orfs = []
    for region in cds_regions:
        loc = region.location
        seq = record.seq[loc.start.position:loc.end.position]
        if region.strand == -1:
            seq = seq.reverse_complement()
            
        if not validate_cds:
            orfs.append((region.strand, loc.start.position, loc.end.position))
            continue

        try:
            assert seq[:3] in start_codons, "Start codon not found!"
            assert seq[-3:] in stop_codons, "Stop codon not found!"
            # Make sure there are no stop codons in the middle of the sequence
            for codon in codons(seq[3:-3]):
                assert (
                    codon not in stop_codons
                ), f"Stop codon {codon} found in the middle of the sequence!"

            # The CDS looks fine, add it to the ORFs
            orfs.append((region.strand, loc.start.position, loc.end.position))

        except AssertionError as ex:
            if verbose:
                print(
                    "Skipped CDS at region [%d - %d] on strand %d"
                    % (loc.start.position, loc.end.position, region.strand)
                )
                print("\t", str(ex))

    return orfs


####################################################################################################################################################
####################################################################################################################################################
####################################################################################################################################################


def find_orfs(sequence, start_codons, stop_codons):
    """Find possible ORF candidates in a single reading frame.

    Parameters
    ----------
    sequence: Seq
    start_codons: List[str]
    stop_codons: List[str]

    Returns
    -------
    List[Tuple[int, int]]
        tuples of form (start_loc, stop_loc)

    """
    
    ORF_candidates = []
    start_pos = []
    stop_pos = []
    for i in range(0, len(sequence), 3):
        codon = sequence[i: i + 3]
        if codon in start_codons:
            start_pos.append(i)
        if codon in stop_codons:
            stop_pos.append(i + 3) # Stop codon is included
            
    #print("START:" + str(start_pos))
    #print("STOP:" + str(stop_pos))

    # Append cantidates for ORFS in a prepared list
    stop_p = 0
    for start in start_pos:
        if start >= stop_p :
            for stop in stop_pos:   
                if stop > start:
                    # Update the stop position
                    stop_p = stop
                    ORF_candidates.append((start, stop))
                    break
            
    return ORF_candidates
        

####################################################################################################################################################

    
def find_all_orfs(sequence, start_codons, stop_codons):
    """Find ALL the possible ORF candidates in the sequence using all six
    reading frames.

    Parameters
    ----------
    sequence: Seq
    start_codons: List[str]
    stop_codons: List[str]

    Returns
    -------
    List[Tuple[int, int, int]]
        tuples of form (strand, start_loc, stop_loc). Strand should be either 1
        for reference strand and -1 for reverse complement.

    """
    all_ORF_candidates = []

    # Reference strand (starting at pos 0 (sequence), starting at pos 1, starting at pos 2)
    # Find ORFs for each one of three sequences
    strand0 = find_orfs(sequence, start_codons, stop_codons)
    strand1 = find_orfs(sequence[1:], start_codons, stop_codons)
    strand2 = find_orfs(sequence[2:], start_codons, stop_codons)
    #print(strand0)
    #print(strand1)
    #print(strand2)
    
    # Add ORFs to a list with all candidates -- tuples with 1 
    [all_ORF_candidates.append((1, start0, stop0)) for (start0, stop0) in strand0]
    [all_ORF_candidates.append((1, start1 + 1, stop1 + 1)) for (start1, stop1) in strand1]
    [all_ORF_candidates.append((1, start2 + 2, stop2 + 2)) for (start2, stop2) in strand2]

    #########################

    rev_seq = sequence.reverse_complement()

    # Reversed stand (starting at pos 0 (sequence), starting at pos 1, starting at pos 2)
    # Find ORFs for each one of three reversed sequences ([::-1])
    rev_strand0 = find_orfs(rev_seq, start_codons, stop_codons)
    rev_strand1 = find_orfs(rev_seq[1:], start_codons, stop_codons)
    rev_strand2 = find_orfs(rev_seq[2:], start_codons, stop_codons)
    #print(rev_strain0)
    #print(rev_strain1)
    #print(rev_strain2)
    
    # Add ORFs to a list with all candidates -- tuples with -1
    # Indexes as in reference strand -> len(seq) - index_of_reversed_strand 
    [all_ORF_candidates.append((-1, len(sequence) - stop0, len(sequence) - start0)) for (start0, stop0) in rev_strand0]
    [all_ORF_candidates.append((-1, len(sequence) - (stop1 + 1), len(sequence) - (start1 + 1))) for (start1, stop1) in rev_strand1]
    [all_ORF_candidates.append((-1, len(sequence) - (stop2 + 2), len(sequence) - (start2 + 2))) for (start2, stop2) in rev_strand2]

    return all_ORF_candidates


####################################################################################################################################################


def translate_to_protein(seq):
    """Translate a nucleotide sequence into a protein sequence.

    Parameters
    ----------
    seq: str

    Returns
    -------
    str
        The translated protein sequence.

    """

    # Dictionary
    
    A = ["GCT", "GCC", "GCA", "GCG"]
    C = ["TGT", "TGC"]
    D = ["GAT", "GAC"]
    E = ["GAA", "GAG"]
    F = ["TTT", "TTC"]
    G = ["GGT", "GGC", "GGA", "GGG"]
    H = ["CAT", "CAC"]
    I = ["ATT", "ATC", "ATA"]
    K = ["AAA", "AAG"]
    L = ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"]
    M = ["ATG"]
    N = ["AAT", "AAC"]
    P = ["CCT", "CCC", "CCA", "CCG"]
    Q = ["CAA", "CAG"]
    R = ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"]
    S = ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"]
    T = ["ACT", "ACC", "ACA", "ACG"]
    V = ["GTT", "GTC", "GTA", "GTG"]
    W = ["TGG"]
    Y = ["TAT", "TAC"]

    # Translation
    
    protein = ""
    for i in range(0, len(seq), 3):
        codon = seq[i: i + 3]
        if codon in A:
            protein = protein + "A"
        if codon in C:
            protein = protein + "C"
        if codon in D:
            protein = protein + "D"
        if codon in E:
            protein = protein + "E"
        if codon in F:
            protein = protein + "F"
        if codon in G:
            protein = protein + "G"
        if codon in H:
            protein = protein + "H"
        if codon in I:
            protein = protein + "I"
        if codon in K:
            protein = protein + "K"
        if codon in L:
            protein = protein + "L"
        if codon in M:
            protein = protein + "M"
        if codon in N:
            protein = protein + "N"
        if codon in P:
            protein = protein + "P"
        if codon in Q:
            protein = protein + "Q"
        if codon in R:
            protein = protein + "R"
        if codon in S:
            protein = protein + "S"
        if codon in T:
            protein = protein + "T"
        if codon in V:
            protein = protein + "V"
        if codon in W:
            protein = protein + "W"
        if codon in Y:
            protein = protein + "Y"
        
            
    return protein












