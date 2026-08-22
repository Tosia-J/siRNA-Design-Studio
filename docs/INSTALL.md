# Installation guide

Written for users who do not work with a terminal day to day. Every step is
spelled out.

---

## 1. What to download

All `.py` files into **one folder**. They import one another, so splitting
them into subfolders will break the imports.

```
sirna-design-studio/
    app.py              graphical interface
    main.py             command line
    validate.py         validation against measured efficacy
    thermo.py
    filters.py
    hosts.py
    design.py
    offtarget.py
    conservation.py
    scoring.py
    constructs.py
    oligos.py
    ncbi.py
    report.py
    blast_criteria.py
    requirements.txt
    tests/
        conftest.py
        test_thermo.py
        test_offtarget.py
        test_pipeline.py
```

The `tests` folder is the one exception — those files belong in a subfolder,
and `conftest.py` is what allows them to find the modules one level up.

---

## 2. Create an environment

An environment is an isolated set of packages. Installing something for this
project will not disturb anything else on the machine.

Open **Anaconda Prompt** on Windows, or **Terminal** on macOS and Linux:

```bash
conda create -n sirna python=3.11 -y
conda activate sirna
```

After the second command the prompt begins with `(sirna)`. The environment
must be activated **every time** a new terminal is opened.

---

## 3. Install the packages

```bash
pip install -r requirements.txt
```

If `ViennaRNA` fails to install — which happens on Windows when a C++ compiler
is absent — the pipeline still runs. Structural metrics are skipped, every
other metric is unaffected, and a warning is displayed. Install Microsoft
C++ Build Tools if the structural metrics are needed.

Check that it worked:

```bash
python -c "import RNA; print('ViennaRNA', RNA.__version__)"
```

---

## 4. Move into the folder

The terminal has to be inside the folder containing the code. The command is
`cd`.

**Windows**
```bash
cd C:\Users\YourName\sirna-design-studio
```

**macOS / Linux**
```bash
cd /Users/yourname/sirna-design-studio
```

A shortcut: type `cd ` with a trailing space, then **drag the folder from the
file manager onto the terminal window** — the path is filled in automatically.

Confirm you are in the right place:

```bash
dir        # Windows
ls         # macOS / Linux
```

The list should contain the `.py` files.

---

## 5. Run

### Graphical interface

```bash
streamlit run app.py
```

A browser opens at `http://localhost:8501`. To stop it, return to the
terminal and press **Ctrl+C**.

Note that `localhost` means *this computer*. The address is visible only to
you; nobody else can reach it, even knowing the address. The terminal **is**
the server — closing it stops the application.

### Command line

```bash
python main.py --target sequence.fasta --host plant --out results/
```

Three files appear in `results/`: `ranking.tsv`, `raport.txt` and
`wyniki.json`.

### Tests

```bash
python -m pytest tests/ -v
```

110 tests should pass. If they do not, a file is missing or is a different
version from the rest.

---

## 6. Where to obtain reference data

| Species | Source | File |
|---|---|---|
| *Arabidopsis thaliana* | [TAIR](https://www.arabidopsis.org/download/) | `TAIR10_cdna_*` |
| Tomato, potato, pepper | [Sol Genomics Network](https://solgenomics.net/ftp/) | `ITAG*_cDNA.fasta` |
| Other plants | [Ensembl Plants](https://plants.ensembl.org/info/data/ftp/index.html) | cDNA column |
| Human, mouse and others | [Ensembl](https://www.ensembl.org/info/data/ftp/index.html) | cDNA column |

**Use cDNA, not the genome.** A genome contains introns and intergenic
regions that are never transcribed, so an siRNA will never encounter them.
Hits there are biologically irrelevant and roughly double the analysis time.
For *Arabidopsis* the cDNA file is about 50 MB against 119 MB for the genome.

Files ending in `.gz` must be unpacked first — 7-Zip on Windows, a double
click on macOS.

For isolate sets, the **Data sources** tab in the application can retrieve
sequences directly from NCBI with quality filtering applied.

---

## 7. Common problems

**`ModuleNotFoundError: No module named 'RNA'`**
ViennaRNA is not installed, or the environment is not active. Check that the
prompt shows `(sirna)`.

**`ModuleNotFoundError: No module named 'thermo'`**
The terminal is in the wrong folder, or the files are scattered. Return to
step 4.

**`ERROR: file or directory not found: tests/`**
The `tests` subfolder does not exist. Create it and move the four test files
into it, or run `python -m pytest -v` without a path.

**`FileNotFoundError`**
The path to the FASTA file is wrong. Give the full path, or put the file in
the same folder as the scripts.

**Streamlit does not open a browser**
Go to `http://localhost:8501` manually.

**"No candidates passed the filters"**
The thresholds are too strict for that sequence. In advanced mode, widen the
GC range, lower the asymmetry requirement, or reduce the ORF-terminus
exclusions. For a short ORF, excluding 75 nt from the start and 50 from the
end may leave nothing.

**The application seems to hang while scanning**
Scanning a transcriptome takes about two minutes per 50 Mnt. A progress bar
shows the number of transcripts processed. If nothing moves at all, check
that the uploaded file is cDNA rather than a whole genome.
