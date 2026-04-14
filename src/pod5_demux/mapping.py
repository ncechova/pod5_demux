# src/pod5_demux/mapping.py
import os
import re
import gzip
import glob
from pathlib import Path
from functools import partial
from multiprocessing import Pool
from typing import Tuple, Dict

try:
    import pysam
    PYSAM_AVAILABLE = True
except ImportError:
    PYSAM_AVAILABLE = False

from Bio import SeqIO
from pod5_demux.utils import detect_format

BARCODE_RE = re.compile(r"(barcode\d+)")
def parse_single_mapping_file(file_path: str, file_type: str, is_fmap: bool) -> Tuple[Dict, Dict]:
    """Vyparsuje jeden soubor a vrátí lokální slovníky."""
    mapping = {}
    filename_map = {}
    
    clean_filename = Path(file_path).name
    for ext in [".fastq.gz", ".fastq", f".{file_type}"]:
        clean_filename = clean_filename.replace(ext, "")

    if file_type in ["sam", "bam"]:
        if not PYSAM_AVAILABLE and file_type == "bam":
            return mapping, filename_map 
            
        elif not PYSAM_AVAILABLE and file_type == "sam":
            with open(file_path, 'r') as f:
                for line in f:
                    if line.startswith('@'): continue
                    parts = line.split('\t')
                    if len(parts) < 11: continue
                    read_id = parts[0]
                    found_bc = None
                    for field in parts[11:]:
                        if field.startswith("BC:Z:"):
                            found_bc = field.split(":")[2].strip()
                            break
                        elif field.startswith("RG:Z:"):
                            match = BARCODE_RE.search(field)
                            found_bc = match.group() if match else field
                    if found_bc:
                        mapping[read_id] = found_bc
                        if is_fmap: filename_map[read_id] = clean_filename
        else:
            mode = "rb" if file_type == "bam" else "r"
            try:
                with pysam.AlignmentFile(file_path, mode, check_sq=False) as samfile:
                    for read in samfile.fetch(until_eof=True):
                        read_id = read.query_name
                        try:
                            bc = read.get_tag("BC")
                        except KeyError:
                            try:
                                rg = read.get_tag("RG")
                                match = BARCODE_RE.search(str(rg))
                                bc = match.group() if match else str(rg)
                            except KeyError:
                                bc = None
                        if bc:
                            mapping[read_id] = bc
                            if is_fmap: filename_map[read_id] = clean_filename
            except Exception as e:
                print(f"!!! Chyba při čtení BAM/SAM {file_path}: {e}")

    elif file_type == "fastq":
        open_func = gzip.open if file_path.endswith('.gz') else open
        mode = "rt" if file_path.endswith('.gz') else "r"
        try:
            with open_func(file_path, mode) as handle:
                for record in SeqIO.parse(handle, "fastq"):
                    match = BARCODE_RE.search(record.description)
                    if match:
                        mapping[record.id] = match.group()
                        if is_fmap: filename_map[record.id] = clean_filename
        except Exception as e:
            print(f"!!! Chyba při čtení FASTQ {file_path}: {e}")

    return mapping, filename_map

def load_barcode_map_parallel(input_path: str, output_format: str, n_cores: int) -> Tuple[Dict, Dict]:
    """Najde soubory a paralelně z nich sestaví jeden velký slovník."""
    is_fmap = (output_format == "folder")
    input_format, path_type = detect_format(input_path)

    if not input_format:
        raise ValueError(f"Nelze detekovat formát vstupu pro {input_path}")

    if input_format == "fastq":
        files = glob.glob(os.path.join(input_path, "**", f"*.{input_format}*"), recursive=True) if path_type == "dir" else [input_path]
    else:
        files = glob.glob(os.path.join(input_path, "**", f"*.{input_format}"), recursive=True) if path_type == "dir" else [input_path]

    print(f"-> Nalezeno {len(files)} mapovacích souborů ({input_format.upper()}). Načítám paralelně...")

    barcode_map = {}
    filename_map = {}
    
    worker = partial(parse_single_mapping_file, file_type=input_format, is_fmap=is_fmap)
    
    with Pool(processes=n_cores) as pool:
        results = pool.map(worker, files)
        for m, fm in results:
            barcode_map.update(m)
            filename_map.update(fm)
            
    return barcode_map, filename_map