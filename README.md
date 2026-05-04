# pod5-demux

Nástroj pro demultiplexaci surových nanopórových signálů z POD5 souborů na základě klasifikovaných sekvenačních dat (BAM/SAM/FASTQ).

---

## Obsah

1. [Popis nástroje](#popis-nástroje)
2. [Požadavky](#požadavky)
3. [Instalace](#instalace)
4. [Použití](#použití)
5. [Popis výstupů](#popis-výstupů)
6. [Ukázkový dataset](#ukázkový-dataset)
7. [Limitace](#limitace)

---

## Popis nástroje

`pod5-demux` rozdělí POD5 soubory podle barkódů identifikovaných při basecallingu. Vstupem jsou:

- **klasifikovaná sekvenační data** ve formátu BAM, SAM nebo FASTQ
- **surová signální data** ve formátu POD5

Nástroj vytvoří sadu výstupních POD5 souborů, kde každý soubor (nebo složka) obsahuje pouze čtení příslušející jednomu vzorku (barkódu).

Zpracování probíhá paralelně na více jádrech CPU pomocí strategie **split–merge**:
1. Každý POD5 soubor je rozdělen na dočasné části podle barkódu.
2. Části jsou sloučeny do finálních výstupních souborů.

---

## Požadavky

- Python >= 3.8
- Operační systém: **Linux** (doporučeno), Windows (s omezením — viz [Limitace](#limitace))

Závislosti se nainstalují automaticky při instalaci balíčku:

| Balíček    | Účel                                      |
|------------|-------------------------------------------|
| `pod5`     | Čtení a zápis POD5 souborů                |
| `biopython`| Parsování FASTQ souborů                   |
| `pysam`    | Čtení BAM/SAM souborů (pouze Linux/macOS) |
| `typer`    | Rozhraní příkazové řádky                  |

---

## Instalace

### Instalace z lokálního zdrojového kódu

'''bash
# 1. Klonujte nebo rozbalte zdrojové kódy do složky
cd pod5-demux

# 2. Nainstalujte balíček (včetně závislostí)
pip install.
'''

### Instalace z githubu

'''bash
pip install git+https://github.com/ncechova/pod5_demux.git
'''

Po instalaci bude dostupný příkaz `pod5_demux`:

'''bash
pod5_demux --help
'''

---

## Použití

### Základní syntaxe

'''bash
pod5_demux -s <cesta_k_mapě> -p <cesta_k_pod5> -o <výstupní_složka> [volby]
'''

### Parametry

| Parametr | Zkratka | Popis | Výchozí hodnota |
|---|---|---|---|
| `--seq` | `-s` | Cesta k mapovacímu souboru nebo složce (BAM/SAM/FASTQ) | -- (povinné) |
| `--pod5` | `-p` | Cesta k POD5 souboru nebo složce | -- (povinné) |
| `--output` | `-o` | Výstupní složka | `pod5_demux_output` |
| `--bc` | `-b` | Vynutí přiřazení konkrétního barkódu (např. `barcode05`) všem čtením z mapovacího souboru. | -- (volitelné) |
| `--mode` | `-m` | Režim výstupu: `folder` nebo `single_file` | `folder` |
| `--threads` | `-t` | Počet vláken CPU | 8 |

### Režimy výstupu

**'folder'** (výchozí) — zachová původní strukturu souborů, každý barkód dostane vlastní podsložku:

'''
vystup/
├── barcode01/
│   ├── run1.pod5
│   └── run2.pod5
├── barcode02/
│   └── run1.pod5
└── unclassified/
│   ├── run1.pod5
│   └── unclassified.pod5
'''

**`single_file`** — všechna čtení jednoho barkódu jsou sloučena do jednoho souboru:

'''
vystup/
├── barcode01.pod5
├── barcode02.pod5
└── unclassified.pod5
'''

### Filtrování a vynucení specifického barkódu (`--bc`)

Pokud je použit parametr `--bc <název_barkódu>`, nástroj nebude barkódy vyhledávat v metadatech mapovacích souborů, ale **všem čtením v mapě přiřadí zadaný barkód**. Zároveň dojde k vyfiltrování pouze tohoto vzorku – ostatní čtení v původních POD5 souborech budou přeskočena.

Tato funkce je ideální v situacích, kdy:
- Máte k dispozici vstupní soubor (např. FASTQ), který už obsahuje data pouze z jednoho vzorku, ale chybí v něm metadata o barkódech.
- Chcete z původního velkého POD5 datasetu efektivně a rychle vyextrahovat pouze data patřící jednomu konkrétnímu vzorku bez nutnosti zpracovávat zbytek.

### Příklady použití

**Demultiplexace s BAM mapou, výstup do složek:**

'''bash
pod5_demux \
  -s /data/basecalled/calls.bam \
  -p /data/raw/ \
  -o /data/demuxed/
'''

**Demultiplexace s FASTQ mapou, vše do jednoho souboru na barkód:**

'''bash
pod5_demux \
  -s /data/basecalled/fastq/ \
  -p /data/raw/ \
  -o /data/demuxed/ \
  -m single_file
'''

**Omezení počtu jader:**

'''bash
pod5_demux \
  -s /data/basecalled/calls.bam \
  -p /data/raw/ \
  -o /data/demuxed/ \
  -t 8
'''

**Extrakce a přiřazení jednoho specifického barkódu:**

'''bash
pod5_demux \
  -s /data/basecalled/sample_barcode05.fastq \
  -p /data/raw/ \
  -o /data/demuxed_single/ \
  -b barcode05
'''

---

## Popis výstupů

Po dokončení nástroj vytvoří ve výstupní složce:

- **POD5 soubory** rozdělené podle barkódu (ve struktuře dle zvoleného režimu)
- **`demux_stats.txt`** — souhrnný report obsahující:
  - celkový počet zpracovaných čtení
  - celkový čas zpracování
  - použitý režim výstupu
  - počet unikátních barkódů
  - počet vytvořených výstupních souborů

Příklad obsahu `demux_stats.txt`:

```
Zpracovaná čtení: 150432
Celkový čas: 47.83 s
Režim: folder
Unikátní barkódy: 12
Vytvořeno souborů: 24
```

---

## Ukázkový dataset


---

## Limitace

- **Windows:** Knihovna `pysam` není na Windows dostupná. Nástroj to automaticky detekuje a přepne na záložní vlastní parsery pro BAM i SAM (pomalejší než pysam). Formát FASTQ je plně podporován na všech platformách.
- **Nedemultiplexovaná čtení:** Čtení, jejichž UUID není nalezeno v mapovacím souboru, jsou zařazena do skupiny `unclassified`.
- **Formát vstupu:** Nástroj automaticky detekuje formát vstupního mapovacího souboru. Vstupní složka musí obsahovat soubory pouze jednoho formátu (nelze kombinovat BAM a FASTQ ve stejné složce).
