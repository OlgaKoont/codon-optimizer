"""
FASTA file handler for reading and writing sequences.
"""

from typing import List, Dict, Tuple, Optional
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class FastaHandler:
    """Handler for reading and writing FASTA files."""
    
    @staticmethod
    def read_fasta(file_path: str) -> List[SeqRecord]:
        """
        Read sequences from a FASTA file.
        
        Args:
            file_path: Path to FASTA file
            
        Returns:
            List of SeqRecord objects
        """
        file_path = Path(file_path)
        if not file_path.exists():
            raise FileNotFoundError(f"FASTA file not found: {file_path}")
        
        records = list(SeqIO.parse(file_path, "fasta"))
        if not records:
            raise ValueError(f"No sequences found in {file_path}")
        
        return records
    
    @staticmethod
    def read_single_sequence(file_path: str) -> Tuple[str, str]:
        """
        Read a single sequence from FASTA file.
        
        Args:
            file_path: Path to FASTA file
            
        Returns:
            Tuple of (sequence_id, sequence_string)
        """
        records = FastaHandler.read_fasta(file_path)
        if len(records) > 1:
            raise ValueError(f"Expected single sequence, found {len(records)}")
        
        record = records[0]
        return record.id, str(record.seq)
    
    @staticmethod
    def write_fasta(sequences: List[Tuple[str, str]], output_path: str):
        """
        Write sequences to a FASTA file.
        
        Args:
            sequences: List of (sequence_id, sequence_string) tuples
            output_path: Path to output FASTA file
        """
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        records = []
        for seq_id, seq_string in sequences:
            record = SeqRecord(Seq(seq_string), id=seq_id, description="")
            records.append(record)
        
        SeqIO.write(records, output_path, "fasta")
    
    @staticmethod
    def write_single_sequence(sequence_id: str, sequence: str, output_path: str):
        """
        Write a single sequence to FASTA file.
        
        Args:
            sequence_id: Sequence identifier
            sequence: Sequence string
            output_path: Path to output FASTA file
        """
        FastaHandler.write_fasta([(sequence_id, sequence)], output_path)
    
    @staticmethod
    def validate_sequence(sequence: str, sequence_type: str = "DNA") -> bool:
        """
        Validate sequence contains only valid nucleotides or amino acids.
        
        Args:
            sequence: Sequence string
            sequence_type: "DNA", "RNA", or "PROTEIN"
            
        Returns:
            True if valid, False otherwise
        """
        sequence = sequence.upper().replace(" ", "").replace("\n", "")
        
        if sequence_type == "DNA":
            valid_chars = set("ATCGN-")
        elif sequence_type == "RNA":
            valid_chars = set("AUCGN-")
        elif sequence_type == "PROTEIN":
            valid_chars = set("ACDEFGHIKLMNPQRSTVWYXZ*-")
        else:
            raise ValueError(f"Unknown sequence_type: {sequence_type}")
        
        return all(c in valid_chars for c in sequence)
    
    @staticmethod
    def read_protein_sequence(file_path: str) -> Tuple[str, str]:
        """
        Read a protein sequence from FASTA file.
        
        Args:
            file_path: Path to FASTA file with protein sequence
            
        Returns:
            Tuple of (sequence_id, protein_sequence_string)
        """
        seq_id, sequence = FastaHandler.read_single_sequence(file_path)
        
        # Validate it's a protein sequence
        if not FastaHandler.validate_sequence(sequence, "PROTEIN"):
            raise ValueError(f"Sequence in {file_path} does not appear to be a protein sequence")
        
        return seq_id, sequence
    
    @staticmethod
    def detect_sequence_type(sequence: str) -> str:
        """
        Detect if sequence is DNA, RNA, or protein.
        
        Args:
            sequence: Sequence string
            
        Returns:
            "DNA", "RNA", or "PROTEIN"
        """
        sequence = sequence.upper().replace(" ", "").replace("\n", "")
        
        # Check for protein characters
        protein_chars = set("ACDEFGHIKLMNPQRSTVWYXZ*-")
        if any(c in protein_chars for c in sequence):
            # Check if mostly protein
            protein_count = sum(1 for c in sequence if c in protein_chars)
            if protein_count > len(sequence) * 0.5:
                return "PROTEIN"
        
        # Check for RNA
        if 'U' in sequence and 'T' not in sequence:
            return "RNA"
        
        # Default to DNA
        return "DNA"
    
    @staticmethod
    def translate_dna_to_protein(dna_sequence: str) -> str:
        """
        Translate DNA sequence to protein sequence.
        
        Args:
            dna_sequence: DNA sequence string
            
        Returns:
            Protein sequence string
        """
        seq = Seq(dna_sequence.upper())
        protein = seq.translate()
        return str(protein)
    
    @staticmethod
    def get_orf_sequence(sequence: str, start_codon: str = "ATG") -> Optional[Tuple[int, int, str]]:
        """
        Extract ORF from sequence starting with start codon.
        
        Args:
            sequence: DNA sequence
            start_codon: Start codon (default: ATG)
            
        Returns:
            Tuple of (start_pos, end_pos, orf_sequence) or None if not found
        """
        sequence = sequence.upper()
        start_pos = sequence.find(start_codon)
        
        if start_pos == -1:
            return None
        
        # Find stop codons
        stop_codons = ["TAA", "TAG", "TGA"]
        orf_end = len(sequence)
        
        for i in range(start_pos + 3, len(sequence) - 2, 3):
            codon = sequence[i:i+3]
            if codon in stop_codons:
                orf_end = i + 3
                break
        
        orf_sequence = sequence[start_pos:orf_end]
        return (start_pos, orf_end, orf_sequence)



