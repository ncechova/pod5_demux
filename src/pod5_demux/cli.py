# src/pod5_demux/cli.py
import os
import glob
import time
import shutil
from functools import partial
from multiprocessing import Pool
from typing import Annotated, Optional
import re
from enum import Enum

import typer

# Import funkcí
from pod5_demux.utils import ensure_unique_dir
from pod5_demux.mapping import load_barcode_map_parallel
from pod5_demux.pod5_core import split_one_pod5_by_barcode, merge_pod5_files

try:
    # Změna barvy pro [required]
    from typer import rich_utils as _ru
    _ru.STYLE_REQUIRED_SHORT = "bright_red"
    _ru.STYLE_REQUIRED_LONG  = "bright_red"
except (ImportError, AttributeError):
    pass

app = typer.Typer(
    help="Optimalizovaný nástroj pro demultiplexaci POD5 souborů podle mapy z BAM/SAM/FASTQ.",
    add_completion=False,
)
BARCODE_RE = re.compile(r"(barcode\d+)")

# Povolené hodnoty pro output_mode, aby uživatel nemohl zadat neplatnou možnost
class OutputMode(str, Enum):
    single_file = "single_file"
    folder = "folder"

def run_demultiplexing(input_seq: str, input_pod5: str, output_dir: str, known_bc: Optional[str], output_mode: OutputMode, n_cores: int):
    start_time = time.perf_counter()

    output_dir = ensure_unique_dir(output_dir)
    os.makedirs(output_dir, exist_ok=True)

    print("*" * 50)
    print(" POD5 DEMULTIPLEXER")
    print("*" * 50)
    print(f"Vstupní mapa:   {input_seq}")
    print(f"Vstupní POD5:   {input_pod5}")
    print(f"Výstupní složka: {output_dir}")
    print(f"Režim výstupu:  {output_mode}")
    print(f"Počet jader:    {n_cores}")
    print("=" * 50)

    # NAČTENÍ MAPY
    try:
        bc_map, f_map = load_barcode_map_parallel(input_seq, output_mode, known_bc, n_cores)
    except ValueError as e:
        print(f"Kritická chyba: {e}")
        return

    if not bc_map:
        print("!!! Proces ukončen: Nenalezeny žádné barkódy v mapovacích souborech.")
        return

    map_time = time.perf_counter()
    print(f"-> Mapa načtena: {len(bc_map)} ID (čas: {map_time - start_time:.2f} s)")

    # NALEZENÍ POD5
    pod5_files = glob.glob(os.path.join(input_pod5, "**", "*.pod5"), recursive=True) if os.path.isdir(input_pod5) else [input_pod5]
    if not pod5_files:
        print("!!! Proces ukončen: Nenalezeny žádné zdrojové POD5 soubory.")
        return

    # SPLIT FÁZE
    # vytvoření dočasné složky pro rozdělené části
    temp_dir = os.path.join(output_dir, "temp_parts")
    os.makedirs(temp_dir, exist_ok=True)

    print(f"-> Rozděluji {len(pod5_files)} POD5 souborů...")
    worker_pod5 = partial(split_one_pod5_by_barcode, mapping=bc_map, output_dir=output_dir, output_mode=output_mode, filename_map=f_map, known_bc=known_bc)

    with Pool(n_cores) as p:
        results = p.map(worker_pod5, pod5_files)

    total_reads = sum(r["processed"] for r in results)
    process_time = time.perf_counter()
    print(f"-> Dělení dokončeno: Zpracováno {total_reads} čtení (čas: {process_time - map_time:.2f} s)")

    # MERGE FÁZE
    print("-> Slučování rozdělených souborů (Merge)...")
    temp_files = glob.glob(os.path.join(temp_dir, "*.pod5"))

    groups = {}
    for f in temp_files:
        fname = os.path.basename(f)
        parts = fname.split("__")
        bc = parts[0]

        if output_mode == "folder":
            origin = parts[1] if len(parts) > 1 else "unknown"
            final_folder = os.path.join(output_dir, bc)
            final_name = f"{origin}.pod5"
        else:
            final_folder = output_dir
            final_name = f"{bc}.pod5"

        target_key = (final_folder, final_name)
        if target_key not in groups:
            groups[target_key] = []
        groups[target_key].append(f)

    # Vytvoření výstupní složky a příprava argumentů pro paralelní merge
    merge_args = []
    single_moves = []
    for (folder, name), parts in groups.items():
        os.makedirs(folder, exist_ok=True)
        final_path = os.path.join(folder, name)
        if len(parts) > 1:
            merge_args.append((parts, final_path))
        else:
            single_moves.append((parts[0], final_path))

    # Přesun jednoduchých (1 část) souborů bez nutnosti merge
    for src, dst in single_moves:
        shutil.move(src, dst)

    # Paralelní merge vícedílných skupin
    if merge_args:
        with Pool(n_cores) as p:
            p.map(merge_pod5_files, merge_args)

    shutil.rmtree(temp_dir, ignore_errors=True)

    end_time = time.perf_counter()

    # REPORT
    print("\n" + "="*50)
    print(" ZPRACOVÁNÍ DOKONČENO ")
    print("="*50)
    print(f"Celkový čas:          {end_time - start_time:.2f} sekund")
    print(f"Čas dělení (Split):   {process_time - map_time:.2f} sekund")
    print(f"Čas slučování (Merge): {end_time - process_time:.2f} sekund")

    with open(os.path.join(output_dir, "demux_stats.txt"), "w") as stats_file:
        stats_file.write(f"Zpracovaná čtení: {total_reads}\n")
        stats_file.write(f"Celkový čas: {end_time - start_time:.2f} s\n")
        stats_file.write(f"Režim: {output_mode}\n")
        stats_file.write(f"Unikátní barkódy: {len(set(bc_map.values()))}\n")
        stats_file.write(f"Vytvořeno souborů: {len(groups)}\n")


@app.command()
def main(
    seq: Annotated[str, typer.Option("-s", "--seq",
                                     exists=True,
                                     resolve_path=True,
                                     readable=True,
                                     help="Cesta k referenční mapě (BAM/SAM/FASTQ).")],
    pod5_path: Annotated[str, typer.Option("-p", "--pod5",
                                     exists=True,
                                     resolve_path=True,
                                     readable=True,
                                     help="Cesta k POD5 datům.")],
    output: Annotated[str, typer.Option("-o", "--output", help="Výstupní složka. (Výchozí: 'pod5_demux_output')")] = "pod5_demux_output",
    known_bc: Annotated[Optional[str], typer.Option("-b", "--bc", help="Název barkódu mapovací sekvence.")] = None,
    mode: Annotated[OutputMode, typer.Option("-m", "--mode", help="Režim: 'single_file' nebo 'folder'.")] = OutputMode.folder,
    threads: Annotated[Optional[int], typer.Option("-t", "--threads", help="Počet jader CPU použitých k výpočtům.")] = 8
):
    """
    Spustí proces demultiplexace.
    """

    run_demultiplexing(seq, pod5_path, output, known_bc, mode, threads)

