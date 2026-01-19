#!/usr/bin/python3
"""
Helpers for running Amsample-related work in parallel (or in chunks).

This module does NOT do any multiprocessing / SLURM by itself.
Instead, it provides:
  - build_partial_amsample: run bam_to_am on a subset of chromosomes
  - merge_amsamples: merge several partial Amsamples into one full object

You can call these sequentially or from multiprocessing / job arrays.
"""

import copy
from typing import List, Sequence, Optional, Dict, Any
import multiprocessing as mp

import amsample as a


def build_partial_amsample(
    name: str,
    filename: str = "",
    filedir: Optional[str] = None,
    file_per_chrom: bool = False,
    library: Optional[str] = None,
    species: Optional[str] = None,
    chroms: Optional[Sequence[Any]] = None,
    trim_ends: bool = False,
    mapq_thresh: int = 20,
    qual_thresh: int = 20,
    gc_object: str = "",
) -> a.Amsample:
    """
    Convenience wrapper around Amsample.bam_to_am that restricts to a
    subset of chromosomes.

    This is intended to be called by different jobs/processes with different
    chroms subsets.

    Parameters mirror Amsample.bam_to_am, plus 'name'.

    Returns
    -------
    Amsample
        An Amsample object that has data only for the requested chroms.
    """
    if chroms is None:
        # default: behave like bam_to_am with its default chrom list
        chroms = list(range(23))

    ams = a.Amsample(name=name)
    ams.bam_to_am(
        filename=filename,
        filedir=filedir,
        file_per_chrom=file_per_chrom,
        library=library,
        species=species,
        chroms=list(chroms),  # ensure we have a concrete list
        trim_ends=trim_ends,
        mapq_thresh=mapq_thresh,
        qual_thresh=qual_thresh,
        gc_object=gc_object,
    )
    return ams


def merge_amsamples(
    ams_list: List[a.Amsample],
    chroms_full_order: Sequence[Any],
) -> a.Amsample:
    """
    Merge several partial Amsample objects into a single full Amsample.

    Each partial Amsample is assumed to:
      - represent the same biological sample (same name, species, library, etc.)
      - contain data for a subset of chromosomes (its own chr_names / no_* entries)

    Parameters
    ----------
    ams_list : list of Amsample
        Partially-filled Amsample objects (e.g., each from a different subset of chroms).
    chroms_full_order : sequence
        The desired global chromosome order, e.g. [1,2,...,22,"X"] or
        ["1","2",...,"22","X"]. This should match what you normally pass to
        bam_to_am as 'chroms'.

    Returns
    -------
    Amsample
        A new Amsample instance with all chromosomes filled in, in the requested order.
    """
    if not ams_list:
        raise ValueError("merge_amsamples: ams_list is empty")

    # Use first as a template for top-level metadata
    template = ams_list[0]

    n_chroms = len(chroms_full_order)

    # Create a fresh Amsample with the correct sized arrays
    merged = a.Amsample(
        name=template.name,
        species=template.species,
        reference=template.reference,
        library=template.library,
        chr_names=[None] * n_chroms,
        coord_per_position=template.coord_per_position,
        no_a=[None] * n_chroms,
        no_c=[None] * n_chroms,
        no_g=[None] * n_chroms,
        no_t=[None] * n_chroms,
        g_to_a=[None] * n_chroms,
        c_to_t=[None] * n_chroms,
        diagnostics=copy.deepcopy(template.diagnostics),
        p_filters=copy.deepcopy(template.p_filters),
        is_filtered=template.is_filtered,
        is_simulated=template.is_simulated,
        methylation=copy.deepcopy(template.methylation),
        d_rate=copy.deepcopy(template.d_rate),
        metadata=copy.deepcopy(template.metadata),
    )

    # Build mapping from chromosome "names" to positions in the final arrays.
    #
    # We don't assume how chroms_full_order is formatted; instead, we allow either
    # the raw value ("1","X",1) or a "chr" prefixed version ("chr1","chrX").
    chrom_to_pos: Dict[str, int] = {}
    for i, ch in enumerate(chroms_full_order):
        ch_str = str(ch)
        chrom_to_pos[ch_str] = i
        chrom_to_pos["chr" + ch_str.lstrip("chr")] = i

    # Merge each partial ams into the merged one
    for part in ams_list:
        # basic sanity checks: same sample identity
        if part.name != template.name:
            raise ValueError("merge_amsamples: mismatched sample names")
        if part.species != template.species:
            raise ValueError("merge_amsamples: mismatched species")
        if part.library != template.library:
            raise ValueError("merge_amsamples: mismatched library")

        for local_idx, name in enumerate(part.chr_names):
            if not name:
                continue  # None or empty, nothing filled here
            pos = chrom_to_pos.get(name)
            if pos is None:
                raise ValueError(
                    f"merge_amsamples: chromosome name '{name}' from partial "
                    "Amsample not found in chroms_full_order"
                )

            # Copy per-chrom arrays
            merged.no_a[pos] = part.no_a[local_idx]
            merged.no_c[pos] = part.no_c[local_idx]
            merged.no_g[pos] = part.no_g[local_idx]
            merged.no_t[pos] = part.no_t[local_idx]
            merged.c_to_t[pos] = part.c_to_t[local_idx]
            merged.g_to_a[pos] = part.g_to_a[local_idx]
            merged.chr_names[pos] = name

    merged.no_chrs = n_chroms
    return merged
    
def _split_into_chunks(seq, n_chunks):
    """
    Split a sequence into n_chunks in a round-robin way so that
    big chroms are more evenly spread (0, n, 2n, ...).
    """
    seq = list(seq)
    if n_chunks <= 1 or len(seq) <= 1:
        return [seq]

    chunks = [[] for _ in range(n_chunks)]
    for i, item in enumerate(seq):
        chunks[i % n_chunks].append(item)
    return chunks


def _worker_build_partial(kwargs):
    """Wrapper so we can use Pool.map; kwargs is a dict of args to build_partial_amsample."""
    return build_partial_amsample(**kwargs)


def build_amsample_parallel(
    name: str,
    filename: str = "",
    filedir: Optional[str] = None,
    file_per_chrom: bool = False,
    library: Optional[str] = None,
    species: Optional[str] = None,
    chroms_full: Optional[Sequence[Any]] = None,
    trim_ends: bool = False,
    mapq_thresh: int = 20,
    qual_thresh: int = 20,
    gc_object: str = "",
    n_workers: int = 2,
) -> a.Amsample:
    """
    High-level helper: run bam_to_am in parallel over subsets of chromosomes,
    then merge the resulting partial Amsamples.

    Parameters mirror bam_to_am, plus:
      - chroms_full : the global chrom list (same as you'd pass to bam_to_am)
      - n_workers   : how many worker processes to spawn

    Returns
    -------
    Amsample
        Fully merged Amsample covering chroms_full in the given order.
    """
    if chroms_full is None:
        chroms_full = list(range(23))
    chroms_full = list(chroms_full)

    # Decide worker chunks (round robin so big chroms are spread)
    chunks = _split_into_chunks(chroms_full, n_workers)

    # Build per-chunk argument dicts
    worker_args = []
    for chunk in chunks:
        if not chunk:
            continue
        worker_args.append(
            dict(
                name=name,
                filename=filename,
                filedir=filedir,
                file_per_chrom=file_per_chrom,
                library=library,
                species=species,
                chroms=chunk,
                trim_ends=trim_ends,
                mapq_thresh=mapq_thresh,
                qual_thresh=qual_thresh,
                gc_object=gc_object,
            )
        )

    if not worker_args:
        raise ValueError("No chromosome chunks to process")

    # If only one chunk or one worker, fall back to sequential execution
    if len(worker_args) == 1 or n_workers <= 1:
        partials = [build_partial_amsample(**worker_args[0])]
    else:
        with mp.Pool(processes=n_workers) as pool:
            partials = pool.map(_worker_build_partial, worker_args)

    # Merge all partial Amsamples into a full one
    merged = merge_amsamples(partials, chroms_full_order=chroms_full)
    return merged

