# src/pod5_demux/mapping.py
import os
import re
import gzip
import io
import glob
from pathlib import Path
from functools import partial
from multiprocessing import Pool
from typing import Tuple, Dict, Optional
import struct
import zlib

try:
    import pysam
    PYSAM_AVAILABLE = True
except ImportError:
    PYSAM_AVAILABLE = False

from Bio import SeqIO
from pod5_demux.utils import detect_format

BARCODE_RE = re.compile(r"(barcode\d+)")

# BGZF metahlavička dle SAMv1 sekce 4.1, tabulka na str. 13:
# ID1 ID2 CM FLG(4B) + MTIME(4B) + XFL OS(2B) + XLEN(2B) + SI1 SI2 SLEN(2BH)
_BGZF_HEADER = struct.Struct('<4BI2BH2BH')
 
# CRC32 + ISIZE na konci každého BGZF bloku (SAMv1 sekce 4.1, str. 13)
_BGZF_TAIL = struct.Struct('<2I')
 
# Pevná část BAM záznamu dle SAMv1 sekce 4.2, str. 15:
# refID(4) pos(4) l_read_name(1) mapq(1) bin(2) n_cigar_op(2) flag(2) l_seq(4)
# next_refID(4) next_pos(4) tlen(4) = celkem 32 bytů
# Používáme pouze: l_read_name (byte 8), n_cigar_op (bytes 12-13), l_seq (bytes 16-19)
_BAM_RECORD_HEADER = struct.Struct('<4s4sBBHHHI4s4s4s')  # 32 bytů
 
 
# ---------------------------------------------------------------------------
# BGZF dekomprimace — SAMv1 sekce 4.1
# ---------------------------------------------------------------------------
 
def _read_bgzf_blocks(file_path: str) -> bytes:
    """
    Dekomprimuje BGZF soubor do jednoho bytearray.
 
    Dle SAMv1 sekce 4.1: BGZF soubor je posloupnost gzip bloků, každý
    obsahuje extra pole s BSIZE (celková velikost bloku - 1).
    Data jsou komprimována metodou raw deflate (zlib bez gzip hlavičky).
    EOF je signalizován prázdným blokem s ISIZE == 0.
    """
    result = bytearray()
    hsize = _BGZF_HEADER.size  # 18 bytů
 
    with open(file_path, 'rb') as f:
        while True:
            raw = f.read(hsize)
            if len(raw) < hsize:
                break
 
            ID1, ID2, CM, FLG, MTIME, XFL, OS, XLEN, SI1, SI2, SLEN = _BGZF_HEADER.unpack(raw)
 
            # Ověření BGZF signatury dle SAMv1 sekce 4.1, str. 13
            if not (ID1 == 31 and ID2 == 139 and CM == 8 and FLG == 4
                    and SI1 == 66 and SI2 == 67 and SLEN == 2):
                raise ValueError('Neplatná BGZF hlavička.')
 
            # BSIZE: celková velikost bloku - 1 (SAMv1 sekce 4.1)
            bsize = struct.unpack('<H', f.read(2))[0]
 
            # Velikost komprimovaných dat: BSIZE+1 - 18(hlavička) - 2(BSIZE pole) - 8(CRC32+ISIZE)
            # = BSIZE - XLEN - 19  (SAMv1 sekce 4.1, tabulka str. 13: CDATA = BSIZE-XLEN-19)
            cdata = f.read(bsize - XLEN - 19)
            crc32, isize = _BGZF_TAIL.unpack(f.read(8))
 
            # EOF marker: prázdný blok s ISIZE == 0 (SAMv1 sekce 4.1.2)
            if isize == 0:
                break
 
            # Raw deflate dekomprimace — zlib.decompressobj(-15) = bez gzip hlavičky
            # Stejný přístup jako htslib a bamnostic
            data = zlib.decompressobj(-15).decompress(cdata)
 
            # Ověření integrity dle SAMv1 sekce 4.1
            if zlib.crc32(data) & 0xFFFFFFFF != crc32 or len(data) != isize:
                raise ValueError('BGZF: chyba integrity (CRC32 nebo ISIZE nesouhlasí).')
 
            result.extend(data)
 
    return bytes(result)
 
 
# ---------------------------------------------------------------------------
# BAM parser — čte pouze read_id a tagy BC/RG
# ---------------------------------------------------------------------------
 
def _parse_bam_without_pysam(file_path: str, is_fmap: bool, clean_filename: str, known_bc: Optional[str]) -> Tuple[Dict, Dict]:
    """
    Čte BAM soubor bez pysam. Extrahuje pouze read_id a tagy BC/RG.
 
    Struktura BAM záznamu dle SAMv1 sekce 4.2, str. 15:
      32B  pevná hlavička (refID, pos, l_read_name, mapq, bin,
                           n_cigar_op, flag, l_seq, next_refID, next_pos, tlen)
      l_read_name B   read_name (NUL-terminated)
      n_cigar_op×4 B  cigar     ← přeskočíme
      (l_seq+1)//2 B  seq       ← přeskočíme
      l_seq B         qual      ← přeskočíme
      optional tagy   BC, RG    ← čteme
    """
    mapping = {}
    filename_map = {}
 
    try:
        buf = io.BytesIO(_read_bgzf_blocks(file_path))
 
        # --- BAM hlavička dle SAMv1 sekce 4.2, str. 15 ---
        if buf.read(4) != b'BAM\x01':
            raise ValueError('Neplatný BAM magic.')
 
        l_text = struct.unpack('<I', buf.read(4))[0]
        buf.read(l_text)  # plain-text SAM hlavička
 
        n_ref = struct.unpack('<I', buf.read(4))[0]
        for _ in range(n_ref):
            l_name = struct.unpack('<I', buf.read(4))[0]
            buf.read(l_name + 4)  # jméno reference + délka
 
        # --- Alignment záznamy dle SAMv1 sekce 4.2, str. 15 ---
        while True:
            bs = buf.read(4)
            if len(bs) < 4:
                break
            block_size = struct.unpack('<I', bs)[0]
            block = buf.read(block_size)
            if len(block) < block_size:
                break
 
            # Pevná část: 32 bytů
            # _BAM_RECORD_HEADER = refID(4s) pos(4s) l_read_name(B) mapq(B) bin(H)
            #                      n_cigar_op(H) flag(H) l_seq(I) next_refID(4s) next_pos(4s) tlen(4s)
            _, _, l_read_name, _, _, n_cigar_op, _, l_seq, _, _, _ = _BAM_RECORD_HEADER.unpack(block[:32])
 
            # read_name (UUID) — bez koncového NUL (SAMv1 sekce 4.2)
            read_id = block[32: 32 + l_read_name - 1].decode('ascii', errors='replace')

            if known_bc:
                found_bc = known_bc
            else:
                # Offset kde začínají optional tagy:
                # 32 (pevná část) + l_read_name + n_cigar_op×4 + (l_seq+1)//2 + l_seq
                tag_offset = 32 + l_read_name + n_cigar_op * 4 + (l_seq + 1) // 2 + l_seq
    
                # Parsování optional tagů dle SAMv1 sekce 4.2.4, str. 16–17
                found_bc = _extract_bc_tag(block, tag_offset)
 
            if found_bc:
                mapping[read_id] = found_bc
                if is_fmap:
                    filename_map[read_id] = clean_filename
 
    except Exception as e:
        print(f'!!! Chyba při čtení BAM (vlastní parser) {file_path}: {e}')
 
    return mapping, filename_map
 
 
def _extract_bc_tag(block: bytes, offset: int) -> str:
    """
    Prochází optional tagy BAM záznamu a vrátí hodnotu tagu BC nebo RG.
 
    Struktura optional tagu dle SAMv1 sekce 4.2.4, str. 16–17:
      2B  tag (např. b'BC')
      1B  val_type ('A','c','C','s','S','i','I','f','Z','H','B')
      xB  value (délka závisí na val_type)
 
    Typy pevné délky: A=1B, c/C=1B, s/S=2B, i/I/f=4B
    Typy proměnné délky: Z/H = NUL-terminated string, B = sub-typ + count + pole
    """
    # Velikosti hodnot číselných tag typů v bajtech 
    TAG_VALUE_SIZES = {b'A': 1, b'c': 1, b'C': 1, b's': 2, b'S': 2,
                       b'i': 4, b'I': 4, b'f': 4}
 
    rg_val = None
    i = offset
 
    while i < len(block) - 2:
        tag      = block[i:i+2]
        val_type = block[i+2:i+3]
        i += 3
 
        if val_type in TAG_VALUE_SIZES:
            size = TAG_VALUE_SIZES[val_type]
            val = block[i:i+size]
            i += size
 
            if tag in (b'BC', b'RG'):
                # Číselné typy — pro BC/RG nepravděpodobné, ale ošetříme
                decoded = val.decode('ascii', errors='replace').rstrip('\x00')
                if tag == b'BC':
                    return decoded
                elif tag == b'RG':
                    rg_val = decoded
 
        elif val_type in (b'Z', b'H'):
            # NUL-terminated string (SAMv1 sekce 4.2.4, str. 16)
            end = block.index(b'\x00', i)
            val = block[i:end].decode('ascii', errors='replace')
            i = end + 1
 
            if tag == b'BC':
                return val
            elif tag == b'RG' and rg_val is None:
                match = BARCODE_RE.search(val)
                rg_val = match.group() if match else val
 
        elif val_type == b'B':
            # Pole hodnot: sub-typ(1B) + count(4B) + count×size prvků
            # (SAMv1 sekce 4.2.4, str. 17)
            sub_type = block[i:i+1]
            count    = struct.unpack('<I', block[i+1:i+5])[0]
            size     = TAG_VALUE_SIZES.get(sub_type, 1)
            i += 5 + count * size
 
        else:
            break  # neznámý typ — přestaneme
 
    return rg_val  # BC nenalezeno, vrátíme RG (nebo None)

def parse_single_mapping_file(file_path: str, file_type: str, is_fmap: bool, known_bc: Optional[str]) -> Tuple[Dict, Dict]:
    """Vyparsuje jeden soubor a vrátí lokální slovníky."""
    mapping = {}
    filename_map = {}
    
    clean_filename = Path(file_path).name
    for ext in [".fastq.gz", ".fastq", f".{file_type}"]:
        clean_filename = clean_filename.replace(ext, "")

    if file_type in ["sam", "bam"]:
        if not PYSAM_AVAILABLE and file_type == "bam":
            mapping, filename_map = _parse_bam_without_pysam(file_path, is_fmap, clean_filename, known_bc)
            
        elif not PYSAM_AVAILABLE and file_type == "sam":
            with open(file_path, 'r') as f:
                for line in f:
                    if line.startswith('@'): continue
                    parts = line.split('\t')
                    if len(parts) < 11: continue
                    read_id = parts[0]
                    found_bc = None
                    if known_bc:
                        found_bc = known_bc
                    else:
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
                        if known_bc:
                            bc = known_bc
                        else:
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
                    if known_bc:
                        mapping[record.id] = known_bc
                    else:
                        match = BARCODE_RE.search(record.description)
                        if match:
                            mapping[record.id] = match.group()
                    if is_fmap: filename_map[record.id] = clean_filename
        except Exception as e:
            print(f"!!! Chyba při čtení FASTQ {file_path}: {e}")

    return mapping, filename_map

def load_barcode_map_parallel(input_path: str, output_format: str, known_bc: Optional[str], n_cores: int) -> Tuple[Dict, Dict]:
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
    
    worker = partial(parse_single_mapping_file, file_type=input_format, is_fmap=is_fmap, known_bc=known_bc)
    
    with Pool(processes=n_cores) as pool:
        results = pool.map(worker, files)
        for m, fm in results:
            barcode_map.update(m)
            filename_map.update(fm)
            
    return barcode_map, filename_map
