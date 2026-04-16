# src/pod5_demux/pod5_core.py
import os
import sys
import subprocess
from pathlib import Path
import pod5

def split_one_pod5_by_barcode(pod5_file: str, output_dir: str, mapping: dict, filename_map: dict, output_mode: str) -> dict:
    source_name = Path(pod5_file).stem
    writers = {}
    stats = {"processed": 0, "sorted": 0, "unclassified": 0}
    temp_dir = os.path.join(output_dir, "temp_parts")

    try:
        with pod5.Reader(pod5_file) as reader:
            for read_record in reader:
                stats["processed"] += 1
                read_id = str(read_record.read_id)
                bc_name = mapping.get(read_id, "unclassified")

                if output_mode == "folder" and filename_map:
                    origin_name = filename_map.get(read_id, "unknown")
                    writer_key = f"{bc_name}__{origin_name}"
                else:
                    writer_key = bc_name

                if writer_key not in writers:
                    out_name = f"{writer_key}__{source_name}.pod5"
                    out_path = os.path.join(temp_dir, out_name)
                    writers[writer_key] = pod5.Writer(out_path)
                    writers[writer_key]._add_run_info(read_record.run_info)

                try:
                    writers[writer_key].add_read(read_record.to_read())
                    if bc_name != "unclassified": stats["sorted"] += 1
                    else: stats["unclassified"] += 1
                except Exception as e:
                    print(f"!!! Chyba zápisu ID {read_id}: {e}")
    finally:
        for w in writers.values():
            w.close()
            
    return stats

def merge_pod5_files(args:tuple):
    """Pokusí se o rychlý subprocess merge, při chybě použije pomalejší writer."""
    input_files, output_file = args
    try:
        # 1. Pokus: Volání pod5 přímo
        subprocess.run(
            ["pod5", "merge"] + input_files + ["-o", output_file],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL
        )
    except (FileNotFoundError, subprocess.CalledProcessError):
        try:
            # 2. Pokus: Volání přes aktuální Python exe
            python_bin_dir = os.path.dirname(sys.executable)
            pod5_exe = os.path.join(python_bin_dir, "pod5")
            subprocess.run(
                [pod5_exe, "-m", "pod5", "merge"] + input_files + ["-o", output_file],
                check=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL
            )
        except (FileNotFoundError, subprocess.CalledProcessError):
            # Záchrana přes pod5 writer (nejpomalejší)
            print(f"  [Info] Subprocess merge selhal, používám Python API pro {Path(output_file).name}")
            with pod5.Writer(output_file) as writer:
                for file in input_files:
                    with pod5.Reader(file) as reader:
                        for i, read_rec in enumerate(reader.reads()):
                            if i == 0:
                                writer._add_run_info(read_rec.run_info)
                            writer.add_read(read_rec.to_read())
