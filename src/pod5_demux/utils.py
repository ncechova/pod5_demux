# src/pod5_demux/utils.py
import os
from typing import Tuple

def ensure_unique_dir(path: str) -> str:
    """Zajistí unikátní název výstupní složky (přidáním čísla, pokud již existuje)."""
    base = path.rstrip("/\\")
    if not os.path.exists(base):
        return base
    i = 1
    while True:
        candidate = f"{base}({i})"
        if not os.path.exists(candidate):
            return candidate
        i += 1

def detect_format(input_path: str) -> Tuple[str, str]:
    """Detekuje formát vstupních dat (BAM, SAM, FASTQ) a zda jde o složku či soubor."""
    if not os.path.exists(input_path):
        return "", ""
        
    if os.path.isdir(input_path):
        for _, _, files in os.walk(input_path):
            for f_name in files:
                name_lower = f_name.lower()
                if name_lower.endswith(".bam"): return "bam", "dir"
                elif name_lower.endswith(".sam"): return "sam", "dir"
                elif name_lower.endswith((".fastq", ".fastq.gz")): return "fastq", "dir"
    elif os.path.isfile(input_path):
        name_lower = input_path.lower()
        if name_lower.endswith(".bam"): return "bam", "file"
        elif name_lower.endswith(".sam"): return "sam", "file"
        elif name_lower.endswith((".fastq", ".fastq.gz")): return "fastq", "file"
        
    return "", ""