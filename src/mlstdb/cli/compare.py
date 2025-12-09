# SPDX-FileCopyrightText: 2025-present Himal Shrestha <stha.himal2007@gmail.com>
#
# SPDX-License-Identifier: GPL-3.0-or-later

import click
from pathlib import Path
import sys
from mlstdb.core.compare import (
    compare_databases,
    export_comparison_json,
    export_comparison_text,
    HashDatabase
)
from mlstdb.utils import error, success, info
from tqdm import tqdm


@click.command()
@click.help_option('-h', '--help')
@click.option('--old-db', '-o', required=True, type=click.Path(exists=True, path_type=Path),
              help='Path to the old database directory')
@click.option('--new-db', '-n', required=True, type=click.Path(exists=True, path_type=Path),
              help='Path to the new database directory')
@click.option('--outdir', '-d', required=True, type=click.Path(path_type=Path),
              help='Output directory for comparison results')
@click.option('--type', '-t', 'comparison_type', 
              type=click.Choice(['basic', 'full']), default='basic',
              help='Comparison type: basic (profiles only) or full (profiles + alleles)')
@click.option('--parallel/--no-parallel', default=True,
              help='Enable/disable parallel processing for multiple schemes')
@click.option('--max-workers', '-w', type=int, default=None,
              help='Maximum number of parallel workers (default: CPU count)')
@click.option('--scheme', '-s', multiple=True,
              help='Specific scheme(s) to compare (default: all schemes)')
@click.option('--create-hash-db', is_flag=True,
              help='Create SQLite database with allele hashes (for full comparison)')
def compare(old_db, new_db, outdir, comparison_type, parallel, max_workers, scheme, create_hash_db):
    """Compare two MLST databases to detect changes
    
    This tool compares two versions of an MLST database and reports:
    - Database version changes
    - Profile changes (STs added/removed/modified)
    - Allele sequence changes (with --type full)
    
    Examples:
    
      # Basic comparison (profiles only)
      mlstdb compare --old-db pubmlst_old/ --new-db pubmlst_new/ --outdir results/
    
      # Full comparison including allele sequences
      mlstdb compare --old-db pubmlst_old/ --new-db pubmlst_new/ --outdir results/ --type full
    
      # Compare specific schemes
      mlstdb compare --old-db pubmlst_old/ --new-db pubmlst_new/ --outdir results/ -s klebsiella -s listeria
    
      # Sequential processing (no parallelization)
      mlstdb compare --old-db pubmlst_old/ --new-db pubmlst_new/ --outdir results/ --no-parallel
    """
    
    try:
        # Validate inputs
        if not old_db.is_dir():
            error(f"Old database path is not a directory: {old_db}")
            sys.exit(1)
        
        if not new_db.is_dir():
            error(f"New database path is not a directory: {new_db}")
            sys.exit(1)
        
        # Create output directory
        outdir.mkdir(parents=True, exist_ok=True)
        info(f"Output directory: {outdir}")
        
        # Determine schemes to compare
        if scheme:
            # User specified specific schemes
            schemes_to_compare = list(scheme)
            info(f"Comparing {len(schemes_to_compare)} specific scheme(s): {', '.join(schemes_to_compare)}")
            
            # Validate that schemes exist
            for s in schemes_to_compare:
                old_exists = (old_db / s).exists()
                new_exists = (new_db / s).exists()
                if not old_exists and not new_exists:
                    error(f"Scheme '{s}' not found in either database")
                    sys.exit(1)
        else:
            # Compare all schemes
            old_schemes = {p.name for p in old_db.iterdir() if p.is_dir()}
            new_schemes = {p.name for p in new_db.iterdir() if p.is_dir()}
            schemes_to_compare = sorted(old_schemes | new_schemes)
            
            if not schemes_to_compare:
                error("No schemes found in either database")
                sys.exit(1)
            
            info(f"Found {len(schemes_to_compare)} scheme(s) to compare")
        
        # Show comparison settings
        info(f"Comparison type: {comparison_type}")
        if comparison_type == 'full':
            info("Will compare both profiles and allele sequences")
        else:
            info("Will compare profiles only (use --type full for allele sequences)")
        
        if parallel and len(schemes_to_compare) > 1:
            worker_info = f"{max_workers} workers" if max_workers else "CPU count workers"
            info(f"Parallel processing enabled ({worker_info})")
        else:
            info("Sequential processing")
        
        # Perform comparison
        click.echo()
        info("Starting comparison...")
        
        if scheme:
            # Compare specific schemes
            from mlstdb.core.compare import compare_scheme
            comparisons = []
            
            with tqdm(total=len(schemes_to_compare), desc="Comparing schemes", unit="scheme") as pbar:
                for scheme_name in schemes_to_compare:
                    try:
                        result = compare_scheme(old_db, new_db, scheme_name, comparison_type)
                        comparisons.append(result)
                    except Exception as e:
                        error(f"Error comparing scheme {scheme_name}: {e}")
                    pbar.update(1)
        else:
            # Compare all schemes
            comparisons = compare_databases(
                old_db, new_db, 
                comparison_type=comparison_type,
                parallel=parallel,
                max_workers=max_workers
            )
        
        if not comparisons:
            error("No comparisons completed successfully")
            sys.exit(1)
        
        # Create hash database if requested and full comparison
        if create_hash_db and comparison_type == 'full':
            info("Creating hash database...")
            hash_db_path = outdir / "allele_hashes.db"
            
            with HashDatabase(hash_db_path) as hash_db:
                from mlstdb.core.compare import hash_alleles_for_scheme
                
                # Hash alleles from new database
                for comp in tqdm(comparisons, desc="Hashing alleles", unit="scheme"):
                    new_scheme_path = new_db / comp.scheme_name
                    if new_scheme_path.exists():
                        hash_alleles_for_scheme(new_scheme_path, comp.scheme_name, hash_db)
                
                # Store comparison metadata
                hash_db.store_comparison_metadata(
                    str(old_db), str(new_db), comparison_type
                )
            
            success(f"Hash database created: {hash_db_path}")
        
        # Export results
        info("Exporting results...")
        
        # JSON output
        json_output = outdir / "comparison_results.json"
        export_comparison_json(comparisons, json_output)
        success(f"JSON report: {json_output}")
        
        # Text output
        text_output = outdir / "comparison_report.txt"
        export_comparison_text(comparisons, text_output)
        success(f"Text report: {text_output}")
        
        # Print summary to console
        click.echo()
        click.echo("=" * 80)
        click.echo("COMPARISON SUMMARY")
        click.echo("=" * 80)
        
        total_schemes = len(comparisons)
        schemes_with_changes = sum(1 for c in comparisons if c.has_changes)
        total_profiles_added = sum(c.profiles_added for c in comparisons)
        total_profiles_removed = sum(c.profiles_removed for c in comparisons)
        total_profiles_modified = sum(c.profiles_modified for c in comparisons)
        
        click.echo(f"Schemes compared: {total_schemes}")
        click.echo(f"Schemes with changes: {schemes_with_changes}")
        click.echo()
        click.echo("Profile Changes:")
        click.echo(f"  Added STs: {total_profiles_added}")
        click.echo(f"  Removed STs: {total_profiles_removed}")
        click.echo(f"  Modified STs: {total_profiles_modified}")
        
        if comparison_type == 'full':
            total_alleles_added = sum(c.alleles_added for c in comparisons)
            total_alleles_removed = sum(c.alleles_removed for c in comparisons)
            total_alleles_modified = sum(c.alleles_modified for c in comparisons)
            
            click.echo()
            click.echo("Allele Changes:")
            click.echo(f"  Added alleles: {total_alleles_added}")
            click.echo(f"  Removed alleles: {total_alleles_removed}")
            click.echo(f"  Modified alleles: {total_alleles_modified}")
            
            if total_alleles_modified > 0:
                click.echo()
                click.secho("WARNING: Some allele sequences were modified!", fg='yellow', bold=True)
                click.echo("Check the detailed report for affected alleles.")
        
        click.echo("=" * 80)
        click.echo()
        
        # Show schemes with changes
        if schemes_with_changes > 0:
            click.echo("Schemes with changes:")
            for comp in comparisons:
                if comp.has_changes:
                    changes = []
                    if comp.version_changed:
                        changes.append(f"version: {comp.old_version}→{comp.new_version}")
                    if comp.profiles_added:
                        changes.append(f"+{comp.profiles_added} STs")
                    if comp.profiles_removed:
                        changes.append(f"-{comp.profiles_removed} STs")
                    if comp.profiles_modified:
                        changes.append(f"~{comp.profiles_modified} STs")
                    if comparison_type == 'full':
                        if comp.alleles_added:
                            changes.append(f"+{comp.alleles_added} alleles")
                        if comp.alleles_removed:
                            changes.append(f"-{comp.alleles_removed} alleles")
                        if comp.alleles_modified:
                            changes.append(f"~{comp.alleles_modified} alleles")
                    
                    click.echo(f"  • {comp.scheme_name}: {', '.join(changes)}")
        else:
            info("No changes detected in any scheme")
        
        click.echo()
        success("Comparison completed successfully!")
        
    except Exception as e:
        error(f"Comparison failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
