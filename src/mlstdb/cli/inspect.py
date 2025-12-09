# SPDX-FileCopyrightText: 2025-present Himal Shrestha <stha.himal2007@gmail.com>
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""
CLI command for mlstdb inspect.
Performs phylogenetic analysis of MLST schemes using MAFFT and FastTree.
"""

import click
from pathlib import Path
import sys
from mlstdb.core.inspect import (
    check_tools_installed,
    validate_scheme_directory,
    process_alleles,
    process_concatenated,
    detect_outliers,
    MLSTInspectError
)
from mlstdb.utils import error, success, info


@click.command()
@click.help_option('-h', '--help')
@click.option(
    '--scheme', '-s',
    required=True,
    help='Name of the MLST scheme (e.g., klebsiella)'
)
@click.option(
    '--type', '-t',
    'analysis_type',
    type=click.Choice(['alleles', 'concatenated'], case_sensitive=False),
    required=True,
    help='Type of analysis: alleles (per-allele trees) or concatenated (concatenated ST tree)'
)
@click.option(
    '--pubmlst-dir', '-p',
    type=click.Path(exists=True, file_okay=False, dir_okay=True, path_type=Path),
    required=True,
    help='Path to pubmlst directory containing scheme data'
)
@click.option(
    '--outdir', '-o',
    type=click.Path(file_okay=False, dir_okay=True, path_type=Path),
    default=Path('./mlstdb_inspect_output'),
    help='Output directory for results (default: ./mlstdb_inspect_output)'
)
@click.option(
    '--outlier-threshold',
    type=float,
    default=2.0,
    help='Standard deviation threshold for outlier detection (default: 2.0)'
)
@click.option(
    '--parallel/--no-parallel',
    default=False,
    help='Enable parallel processing for alignment steps'
)
@click.option(
    '--max-workers',
    type=int,
    default=4,
    help='Maximum number of parallel workers (default: 4)'
)
@click.option(
    '--verbose', '-v',
    is_flag=True,
    help='Enable verbose output'
)
def inspect(
    scheme: str,
    analysis_type: str,
    pubmlst_dir: Path,
    outdir: Path,
    outlier_threshold: float,
    parallel: bool,
    max_workers: int,
    verbose: bool
):
    """
    Perform phylogenetic analysis of MLST schemes.
    
    This command analyses MLST scheme data to:
    - Build phylogenetic trees for individual alleles or concatenated STs
    - Calculate pairwise distance matrices
    - Detect outlier sequence types
    
    Uses MAFFT for fast multiple sequence alignment and FastTree for
    rapid phylogenetic tree construction.
    
    Examples:
    
        # Analyse individual alleles
        mlstdb inspect --scheme klebsiella --type alleles --pubmlst-dir ./pubmlst
        
        # Analyse concatenated STs with parallel processing
        mlstdb inspect --scheme klebsiella --type concatenated --pubmlst-dir ./pubmlst --parallel --max-workers 8
        
        # Use custom outlier threshold
        mlstdb inspect --scheme klebsiella --type concatenated --pubmlst-dir ./pubmlst --outlier-threshold 3.0
    """
    
    try:
        # Check if required tools are installed
        mafft_installed, fasttree_installed = check_tools_installed()
        
        if not mafft_installed:
            error("MAFFT is not installed or not in PATH.")
            error("Please install MAFFT: https://mafft.cbrc.jp/alignment/software/")
            error("  Ubuntu/Debian: sudo apt-get install mafft")
            error("  Conda: conda install -c bioconda mafft")
            sys.exit(1)
        
        if not fasttree_installed:
            error("FastTree is not installed or not in PATH.")
            error("Please install FastTree: http://www.microbesonline.org/fasttree/")
            error("  Ubuntu/Debian: sudo apt-get install fasttree")
            error("  Conda: conda install -c bioconda fasttree")
            sys.exit(1)
        
        if verbose:
            info("✓ MAFFT and FastTree are installed")
        
        # Validate scheme directory
        scheme_dir = pubmlst_dir / scheme
        
        if verbose:
            info(f"Validating scheme directory: {scheme_dir}")
        
        scheme_dir, profile_file = validate_scheme_directory(scheme_dir, scheme)
        
        # Create output directory
        output_dir = outdir / scheme
        output_dir.mkdir(parents=True, exist_ok=True)
        
        if verbose:
            info(f"Output directory: {output_dir}")
        
        # Process based on analysis type
        if analysis_type.lower() == 'alleles':
            info(f"Processing individual alleles for scheme: {scheme}")
            
            results = process_alleles(
                scheme_dir,
                output_dir,
                parallel=parallel,
                max_workers=max_workers
            )
            
            successful = sum(1 for v in results.values() if v)
            total = len(results)
            
            success(f"Processed {successful}/{total} alleles successfully")
            
            # Detect outliers for each allele
            if verbose:
                info("Detecting outliers for each allele...")
            
            total_outliers = 0
            for allele_name, succeeded in results.items():
                if succeeded:
                    dist_matrix = output_dir / f"{allele_name}_distance_matrix.tsv"
                    outlier_file = output_dir / f"{allele_name}_outliers.tsv"
                    
                    if dist_matrix.exists():
                        try:
                            n_outliers = detect_outliers(
                                dist_matrix,
                                outlier_file,
                                outlier_threshold
                            )
                            total_outliers += n_outliers
                            if verbose and n_outliers > 0:
                                info(f"  {allele_name}: {n_outliers} outliers detected")
                        except Exception as e:
                            if verbose:
                                error(f"  Failed to detect outliers for {allele_name}: {e}")
            
            if total_outliers > 0:
                info(f"Total outliers detected across all alleles: {total_outliers}")
            else:
                info("No outliers detected")
        
        elif analysis_type.lower() == 'concatenated':
            info(f"Processing concatenated sequences for scheme: {scheme}")
            
            success_status = process_concatenated(
                scheme_dir,
                profile_file,
                output_dir,
                scheme,
                parallel=parallel,
                max_workers=max_workers
            )
            
            if success_status:
                success("Concatenated analysis completed successfully")
                
                # Detect outliers
                dist_matrix = output_dir / f"{scheme}_concatenated_distance_matrix.tsv"
                outlier_file = output_dir / f"{scheme}_concatenated_outliers.tsv"
                
                if dist_matrix.exists():
                    info("Detecting outlier STs...")
                    n_outliers = detect_outliers(
                        dist_matrix,
                        outlier_file,
                        outlier_threshold
                    )
                    
                    if n_outliers > 0:
                        info(f"Detected {n_outliers} outlier STs")
                        info(f"Outlier report saved to: {outlier_file}")
                    else:
                        info("No outlier STs detected")
        
        success(f"\nAnalysis complete! Results saved to: {output_dir}")
    
    except MLSTInspectError as e:
        error(str(e))
        sys.exit(1)
    except KeyboardInterrupt:
        error("\nAnalysis interrupted by user")
        sys.exit(1)
    except Exception as e:
        error(f"Unexpected error: {e}")
        if verbose:
            import traceback
            error(traceback.format_exc())
        sys.exit(1)
