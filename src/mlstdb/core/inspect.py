# SPDX-FileCopyrightText: 2025-present Himal Shrestha <stha.himal2007@gmail.com>
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Core functions for mlstdb inspect command.
Handles phylogenetic analysis of MLST schemes including tree construction,
distance matrix calculation, and outlier detection.

Uses MAFFT for fast multiple sequence alignment and FastTree for rapid
phylogenetic tree construction.
"""

import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
from Bio import SeqIO, AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator


class MLSTInspectError(Exception):
    """Custom exception for MLST inspect operations."""
    pass


def check_tools_installed() -> Tuple[bool, bool]:
    """
    Check if MAFFT and FastTree are installed and accessible.
    
    Returns:
        Tuple of (mafft_installed, fasttree_installed)
    """
    mafft_installed = False
    fasttree_installed = False
    
    try:
        result = subprocess.run(
            ["mafft", "--version"],
            capture_output=True,
            text=True,
            check=False
        )
        mafft_installed = result.returncode == 0
    except FileNotFoundError:
        pass
    
    try:
        result = subprocess.run(
            ["fasttree"],
            capture_output=True,
            text=True,
            check=False
        )
        # FastTree returns non-zero when called without args, but stderr shows version
        fasttree_installed = "FastTree" in result.stderr or result.returncode == 0
    except FileNotFoundError:
        pass
    
    return mafft_installed, fasttree_installed


def read_profile(profile_file: Path) -> Tuple[List[str], Dict[str, Dict[str, str]]]:
    """
    Read the MLST profile file and extract allele information.
    
    Args:
        profile_file: Path to the scheme profile file (e.g., klebsiella.txt)
    
    Returns:
        Tuple containing:
        - List of allele names (column headers)
        - Dictionary mapping ST to allele numbers {ST: {allele_name: allele_number}}
    
    Raises:
        MLSTInspectError: If file cannot be read or parsed
    """
    try:
        with open(profile_file, 'r') as f:
            lines = f.readlines()
        
        if not lines:
            raise MLSTInspectError(f"Profile file is empty: {profile_file}")
        
        # Parse header
        header = lines[0].strip().split('\t')
        if header[0] != 'ST':
            raise MLSTInspectError(f"Invalid profile format. Expected 'ST' as first column.")
        
        allele_names = header[1:]  # Skip 'ST' column
        
        # Parse profiles
        profiles = {}
        for line in lines[1:]:
            if not line.strip():
                continue
            
            parts = line.strip().split('\t')
            st = parts[0]
            alleles = parts[1:]
            
            if len(alleles) != len(allele_names):
                continue  # Skip malformed lines
            
            profiles[st] = {
                allele_names[i]: alleles[i] 
                for i in range(len(allele_names))
            }
        
        return allele_names, profiles
    
    except Exception as e:
        raise MLSTInspectError(f"Error reading profile file {profile_file}: {e}")


def read_allele_sequences(allele_file: Path) -> Dict[str, str]:
    """
    Read allele sequences from a FASTA file.
    
    Args:
        allele_file: Path to the allele FASTA file (e.g., gapA.tfa)
    
    Returns:
        Dictionary mapping allele identifiers to sequences {allele_id: sequence}
    
    Raises:
        MLSTInspectError: If file cannot be read
    """
    try:
        sequences = {}
        for record in SeqIO.parse(allele_file, "fasta"):
            sequences[record.id] = str(record.seq)
        return sequences
    except Exception as e:
        raise MLSTInspectError(f"Error reading allele file {allele_file}: {e}")


def run_mafft_alignment(input_fasta: Path, output_fasta: Path, threads: int = 1) -> None:
    """
    Run MAFFT multiple sequence alignment (much faster than Clustal Omega).
    
    Args:
        input_fasta: Path to input FASTA file
        output_fasta: Path to output aligned FASTA file
        threads: Number of threads to use
    
    Raises:
        MLSTInspectError: If alignment fails
    """
    try:
        # Use --auto for automatic algorithm selection
        # --thread for parallel processing
        cmd = [
            "mafft",
            "--auto",
            "--thread", str(threads),
            str(input_fasta)
        ]
        
        with open(output_fasta, 'w') as outfile:
            result = subprocess.run(
                cmd,
                stdout=outfile,
                stderr=subprocess.PIPE,
                text=True,
                check=False
            )
        
        if result.returncode != 0:
            raise MLSTInspectError(f"MAFFT failed: {result.stderr}")
    
    except FileNotFoundError:
        raise MLSTInspectError("MAFFT not found. Please install mafft.")
    except Exception as e:
        raise MLSTInspectError(f"Error running MAFFT: {e}")


def run_fasttree(aligned_fasta: Path, output_tree: Path, is_nucleotide: bool = True) -> None:
    """
    Build a phylogenetic tree using FastTree (much faster than Clustal Omega).
    
    Args:
        aligned_fasta: Path to aligned FASTA file
        output_tree: Path to output Newick tree file
        is_nucleotide: Whether sequences are nucleotide (True) or protein (False)
    
    Raises:
        MLSTInspectError: If tree building fails
    """
    try:
        cmd = ["fasttree"]
        
        if is_nucleotide:
            cmd.append("-nt")
        
        cmd.append(str(aligned_fasta))
        
        with open(output_tree, 'w') as outfile:
            result = subprocess.run(
                cmd,
                stdout=outfile,
                stderr=subprocess.PIPE,
                text=True,
                check=False
            )
        
        if result.returncode != 0:
            raise MLSTInspectError(f"FastTree failed: {result.stderr}")
    
    except FileNotFoundError:
        raise MLSTInspectError("FastTree not found. Please install fasttree.")
    except Exception as e:
        raise MLSTInspectError(f"Error running FastTree: {e}")


def calculate_distance_matrix_from_alignment(aligned_fasta: Path, output_matrix: Path) -> None:
    """
    Calculate pairwise distance matrix from aligned sequences using BioPython.
    Much faster than running clustalo again.
    
    Args:
        aligned_fasta: Path to aligned FASTA file
        output_matrix: Path to output distance matrix file
    
    Raises:
        MLSTInspectError: If distance calculation fails
    """
    try:
        # Read alignment
        alignment = AlignIO.read(aligned_fasta, "fasta")
        
        # Calculate distances
        calculator = DistanceCalculator('identity')
        dm = calculator.get_distance(alignment)
        
        # Write distance matrix in tab-separated format
        with open(output_matrix, 'w') as f:
            # Write header with number of sequences
            f.write(f"{len(dm.names)}\n")
            
            # Write matrix
            for i, name in enumerate(dm.names):
                distances = [str(dm[name][other]) for other in dm.names]
                f.write(f"{name}\t" + "\t".join(distances) + "\n")
    
    except Exception as e:
        raise MLSTInspectError(f"Error calculating distance matrix: {e}")


def process_single_allele(
    allele_name: str,
    allele_file: Path,
    output_dir: Path,
    threads: int = 1
) -> Tuple[str, bool, Optional[str]]:
    """
    Process a single allele: align, build tree, and calculate distance matrix.
    
    Args:
        allele_name: Name of the allele
        allele_file: Path to allele FASTA file
        output_dir: Output directory for results
        threads: Number of threads for MAFFT
    
    Returns:
        Tuple of (allele_name, success, error_message)
    """
    try:
        # Read sequences
        sequences = read_allele_sequences(allele_file)
        
        if len(sequences) < 2:
            return (allele_name, False, "Not enough sequences for analysis")
        
        # Create temporary file for input
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as tmp:
            tmp_input = Path(tmp.name)
            for allele_id, seq in sequences.items():
                tmp.write(f">{allele_id}\n{seq}\n")
        
        # Create output files
        aligned_fasta = output_dir / f"{allele_name}_aligned.fasta"
        tree_file = output_dir / f"{allele_name}_tree.newick"
        dist_matrix = output_dir / f"{allele_name}_distance_matrix.tsv"
        
        # Run alignment with MAFFT
        run_mafft_alignment(tmp_input, aligned_fasta, threads=threads)
        
        # Build tree with FastTree
        run_fasttree(aligned_fasta, tree_file, is_nucleotide=True)
        
        # Calculate distance matrix
        calculate_distance_matrix_from_alignment(aligned_fasta, dist_matrix)
        
        # Clean up temp file
        tmp_input.unlink()
        
        return (allele_name, True, None)
    
    except Exception as e:
        return (allele_name, False, str(e))


def process_alleles(
    scheme_dir: Path,
    output_dir: Path,
    parallel: bool = False,
    max_workers: int = 4
) -> Dict[str, bool]:
    """
    Process all alleles in a scheme directory.
    
    Args:
        scheme_dir: Directory containing allele FASTA files
        output_dir: Output directory for results
        parallel: Whether to process alleles in parallel
        max_workers: Maximum number of parallel workers
    
    Returns:
        Dictionary mapping allele names to success status
    """
    from tqdm import tqdm
    
    # Find all .tfa files
    allele_files = list(scheme_dir.glob("*.tfa"))
    
    if not allele_files:
        raise MLSTInspectError(f"No .tfa files found in {scheme_dir}")
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    results = {}
    
    if parallel and len(allele_files) > 1:
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = {
                executor.submit(
                    process_single_allele,
                    allele_file.stem,
                    allele_file,
                    output_dir,
                    threads=1  # Each worker uses 1 thread
                ): allele_file.stem
                for allele_file in allele_files
            }
            
            with tqdm(total=len(allele_files), desc="Processing alleles") as pbar:
                for future in as_completed(futures):
                    allele_name, success, error = future.result()
                    results[allele_name] = success
                    if not success and error:
                        print(f"Warning: Failed to process {allele_name}: {error}")
                    pbar.update(1)
    else:
        # Sequential processing with more threads per alignment
        threads_per_job = max_workers if not parallel else 1
        for allele_file in tqdm(allele_files, desc="Processing alleles"):
            allele_name, success, error = process_single_allele(
                allele_file.stem,
                allele_file,
                output_dir,
                threads=threads_per_job
            )
            results[allele_name] = success
            if not success and error:
                print(f"Warning: Failed to process {allele_name}: {error}")
    
    return results


def extract_st_sequences(
    st: str,
    profile: Dict[str, str],
    allele_sequences: Dict[str, Dict[str, str]],
    allele_names: List[str]
) -> Dict[str, str]:
    """
    Extract sequences for a specific ST from allele files.
    
    Args:
        st: Sequence Type identifier
        profile: Dictionary mapping allele names to allele numbers for this ST
        allele_sequences: Nested dictionary {allele_name: {allele_id: sequence}}
        allele_names: List of allele names in order
    
    Returns:
        Dictionary mapping allele names to sequences for this ST
    """
    st_sequences = {}
    
    for allele_name in allele_names:
        allele_num = profile.get(allele_name)
        
        if allele_num is None:
            continue
        
        # Construct allele identifier (e.g., gapA_1)
        allele_id = f"{allele_name}_{allele_num}"
        
        # Get sequence from the allele sequences dictionary
        if allele_name in allele_sequences and allele_id in allele_sequences[allele_name]:
            st_sequences[allele_name] = allele_sequences[allele_name][allele_id]
    
    return st_sequences


def align_allele_sequences_parallel(
    allele_name: str, 
    sequences: Dict[str, str], 
    threads: int = 1
) -> Tuple[str, Dict[str, str]]:
    """
    Align sequences for a single allele using MAFFT (parallelisable).
    
    Args:
        allele_name: Name of the allele
        sequences: Dictionary mapping ST to sequence
        threads: Number of threads for MAFFT
    
    Returns:
        Tuple of (allele_name, aligned_sequences_dict)
    """
    if len(sequences) < 2:
        return (allele_name, sequences)
    
    # Create temporary input file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as tmp_in:
        tmp_input = Path(tmp_in.name)
        for st, seq in sequences.items():
            tmp_in.write(f">{st}\n{seq}\n")
    
    # Create temporary output file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as tmp_out:
        tmp_output = Path(tmp_out.name)
    
    try:
        # Run alignment with MAFFT
        run_mafft_alignment(tmp_input, tmp_output, threads=threads)
        
        # Read aligned sequences
        aligned = {}
        for record in SeqIO.parse(tmp_output, "fasta"):
            aligned[record.id] = str(record.seq)
        
        return (allele_name, aligned)
    
    finally:
        # Clean up temporary files
        tmp_input.unlink()
        tmp_output.unlink()


def concatenate_aligned_sequences(
    aligned_alleles: Dict[str, Dict[str, str]],
    allele_names: List[str],
    sts: List[str]
) -> Dict[str, str]:
    """
    Concatenate aligned allele sequences for each ST.
    
    Args:
        aligned_alleles: Nested dictionary {allele_name: {ST: aligned_sequence}}
        allele_names: List of allele names in order
        sts: List of ST identifiers
    
    Returns:
        Dictionary mapping ST to concatenated sequence
    """
    concatenated = {}
    
    for st in sts:
        concat_seq = ""
        for allele_name in allele_names:
            if allele_name in aligned_alleles and st in aligned_alleles[allele_name]:
                concat_seq += aligned_alleles[allele_name][st]
        
        if concat_seq:
            concatenated[st] = concat_seq
    
    return concatenated


def process_concatenated(
    scheme_dir: Path,
    profile_file: Path,
    output_dir: Path,
    scheme_name: str,
    parallel: bool = False,
    max_workers: int = 4
) -> bool:
    """
    Process concatenated allele sequences for all STs in a scheme.
    Uses MAFFT for fast alignment and FastTree for rapid tree building.
    
    Args:
        scheme_dir: Directory containing allele FASTA files
        profile_file: Path to profile file
        output_dir: Output directory
        scheme_name: Name of the scheme
        parallel: Whether to parallelise alignment steps
        max_workers: Maximum number of parallel workers
    
    Returns:
        bool: Success status
    """
    from tqdm import tqdm
    
    try:
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Read profile
        print("Loading profile...")
        allele_names, profiles = read_profile(profile_file)
        
        if not profiles:
            raise MLSTInspectError("No profiles found in profile file")
        
        print(f"Found {len(profiles)} STs with {len(allele_names)} alleles")
        
        # Load all allele sequences
        print("Loading allele sequences...")
        allele_sequences = {}
        for allele_name in allele_names:
            allele_file = scheme_dir / f"{allele_name}.tfa"
            if allele_file.exists():
                allele_sequences[allele_name] = read_allele_sequences(allele_file)
        
        # Extract sequences for each ST
        print("Extracting ST sequences...")
        st_allele_seqs = {}  # {allele_name: {ST: sequence}}
        
        for allele_name in allele_names:
            st_allele_seqs[allele_name] = {}
        
        for st, profile in profiles.items():
            st_seqs = extract_st_sequences(st, profile, allele_sequences, allele_names)
            for allele_name, seq in st_seqs.items():
                st_allele_seqs[allele_name][st] = seq
        
        # Align sequences for each allele IN PARALLEL
        print("Aligning allele sequences...")
        aligned_alleles = {}
        
        if parallel and len(allele_names) > 1:
            # TRUE PARALLELISATION: align different alleles in parallel
            with ProcessPoolExecutor(max_workers=max_workers) as executor:
                futures = {
                    executor.submit(
                        align_allele_sequences_parallel,
                        allele_name,
                        st_allele_seqs[allele_name],
                        threads=1  # Each worker uses 1 thread
                    ): allele_name
                    for allele_name in allele_names if allele_name in st_allele_seqs
                }
                
                with tqdm(total=len(futures), desc="Aligning alleles") as pbar:
                    for future in as_completed(futures):
                        allele_name, aligned = future.result()
                        aligned_alleles[allele_name] = aligned
                        pbar.update(1)
        else:
            # Sequential with more threads per alignment
            threads_per_job = max_workers if not parallel else 1
            for allele_name in tqdm(allele_names, desc="Aligning alleles"):
                if allele_name in st_allele_seqs:
                    _, aligned = align_allele_sequences_parallel(
                        allele_name,
                        st_allele_seqs[allele_name],
                        threads=threads_per_job
                    )
                    aligned_alleles[allele_name] = aligned
        
        # Concatenate aligned sequences
        print("Concatenating sequences...")
        sts = list(profiles.keys())
        concatenated = concatenate_aligned_sequences(aligned_alleles, allele_names, sts)
        
        # Write concatenated sequences to file
        concat_fasta = output_dir / f"{scheme_name}_concatenated.fasta"
        with open(concat_fasta, 'w') as f:
            for st, seq in concatenated.items():
                f.write(f">ST_{st}\n{seq}\n")
        
        print(f"Wrote {len(concatenated)} concatenated sequences")
        
        # Build tree using FastTree (MUCH faster)
        print("Building tree with FastTree...")
        tree_file = output_dir / f"{scheme_name}_concatenated_tree.newick"
        run_fasttree(concat_fasta, tree_file, is_nucleotide=True)
        
        # Calculate distance matrix using BioPython (faster)
        print("Calculating distance matrix...")
        dist_matrix = output_dir / f"{scheme_name}_concatenated_distance_matrix.tsv"
        calculate_distance_matrix_from_alignment(concat_fasta, dist_matrix)
        
        return True
    
    except Exception as e:
        raise MLSTInspectError(f"Error processing concatenated sequences: {e}")


def parse_distance_matrix(matrix_file: Path) -> Tuple[List[str], np.ndarray]:
    """
    Parse a distance matrix file.
    
    Args:
        matrix_file: Path to distance matrix file
    
    Returns:
        Tuple of (labels, distance_matrix)
    """
    with open(matrix_file, 'r') as f:
        lines = f.readlines()
    
    # Skip header lines (first line typically contains number of sequences)
    data_lines = [l.strip() for l in lines if l.strip() and not l.strip().isdigit()]
    
    labels = []
    matrix_data = []
    
    for line in data_lines:
        parts = line.split()
        if len(parts) > 1:
            labels.append(parts[0])
            # Extract numeric values
            distances = [float(x) for x in parts[1:]]
            matrix_data.append(distances)
    
    return labels, np.array(matrix_data)


def detect_outliers(
    matrix_file: Path,
    output_file: Path,
    std_threshold: float = 2.0
) -> int:
    """
    Detect outlier STs based on distance matrix analysis.
    
    Args:
        matrix_file: Path to distance matrix file
        output_file: Path to output outlier report file
        std_threshold: Number of standard deviations for outlier detection
    
    Returns:
        Number of outliers detected
    """
    try:
        labels, matrix = parse_distance_matrix(matrix_file)
        
        if matrix.size == 0:
            return 0
        
        # Calculate mean distance for each ST (row-wise mean excluding diagonal)
        n = len(labels)
        mean_distances = []
        
        for i in range(n):
            # Get all distances for this ST, excluding itself (diagonal)
            distances = [matrix[i][j] for j in range(n) if i != j]
            if distances:
                mean_distances.append(np.mean(distances))
            else:
                mean_distances.append(0.0)
        
        mean_distances = np.array(mean_distances)
        
        # Calculate overall statistics
        overall_mean = np.mean(mean_distances)
        overall_std = np.std(mean_distances)
        
        # Identify outliers
        threshold = overall_mean + (std_threshold * overall_std)
        outliers = []
        
        for i, (label, mean_dist) in enumerate(zip(labels, mean_distances)):
            if mean_dist > threshold:
                outliers.append({
                    'label': label,
                    'mean_distance': mean_dist,
                    'z_score': (mean_dist - overall_mean) / overall_std if overall_std > 0 else 0
                })
        
        # Write outlier report
        with open(output_file, 'w') as f:
            f.write("ST\tMean_Distance\tZ_Score\tThreshold\n")
            for outlier in outliers:
                f.write(f"{outlier['label']}\t{outlier['mean_distance']:.6f}\t"
                       f"{outlier['z_score']:.3f}\t{threshold:.6f}\n")
        
        return len(outliers)
    
    except Exception as e:
        raise MLSTInspectError(f"Error detecting outliers: {e}")


def validate_scheme_directory(scheme_dir: Path, scheme_name: str) -> Tuple[Path, Path]:
    """
    Validate that the scheme directory exists and contains required files.
    
    Args:
        scheme_dir: Path to scheme directory
        scheme_name: Name of the scheme
    
    Returns:
        Tuple of (scheme_dir, profile_file)
    
    Raises:
        MLSTInspectError: If validation fails
    """
    if not scheme_dir.exists():
        raise MLSTInspectError(f"Scheme directory not found: {scheme_dir}")
    
    if not scheme_dir.is_dir():
        raise MLSTInspectError(f"Not a directory: {scheme_dir}")
    
    # Look for profile file
    profile_file = scheme_dir / f"{scheme_name}.txt"
    
    if not profile_file.exists():
        raise MLSTInspectError(f"Profile file not found: {profile_file}")
    
    # Check for at least one .tfa file
    tfa_files = list(scheme_dir.glob("*.tfa"))
    if not tfa_files:
        raise MLSTInspectError(f"No .tfa files found in {scheme_dir}")
    
    return scheme_dir, profile_file
