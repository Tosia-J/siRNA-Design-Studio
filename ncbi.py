"""
=============================================================================
ncbi.py  --  NCBI E-UTILITIES INTEGRATION
=============================================================================

WHAT THIS MODULE IS FOR - AND WHAT IT IS NOT

E-utilities is the appropriate interface for retrieving a defined, bounded
set of records: a target gene, a set of viral isolates, a handful of
reference transcripts.

It is NOT an appropriate interface for whole transcriptomes. Retrieving tens
of thousands of cDNA records one request at a time is slow, places an
unreasonable load on a public service, and duplicates files that are
published as single downloads. Whole transcriptomes should be obtained from:

    Ensembl / Ensembl Plants   FTP, file *.cdna.all.fa.gz
    TAIR                       TAIR10_cdna_*
    Sol Genomics Network       ITAG*_cDNA.fasta

The distinction matters for the conservation analysis in particular. A set of
310 viral genomes is exactly the kind of bounded, query-defined collection
that E-utilities handles well, and assembling it by hand is tedious.

=============================================================================
USAGE LIMITS
=============================================================================

NCBI imposes:
    - 3 requests per second without an API key
    - 10 requests per second with a key (free, from an NCBI account)
    - identification by tool name and e-mail address

This module enforces the rate limit internally. Supplying an e-mail address
is a condition of use, not an optional courtesy: NCBI blocks unidentified
clients that generate heavy traffic.

Obtaining a key: https://www.ncbi.nlm.nih.gov/account/settings/

=============================================================================
TYPICAL WORKFLOW
=============================================================================

    from ncbi import NCBIClient

    cli = NCBIClient(email='you@example.org')

    # how many records match, before downloading anything
    n = cli.count('nuccore', 'tomato brown rugose fruit virus[Organism] '
                             'AND complete genome[Title]')

    # retrieve with quality filtering
    seqs = cli.fetch_sequences(
        'nuccore',
        'tomato brown rugose fruit virus[Organism] AND complete genome[Title]',
        max_records=310,
        min_length=6000,
        exclude_partial=True)

    cli.save_fasta(seqs, 'tobrfv_isolates.fasta')

Author: Antonina Jarecka
=============================================================================
"""

from dataclasses import dataclass, field
from typing import Callable, Dict, Iterator, List, Optional, Tuple
from urllib.parse import urlencode
from urllib.request import urlopen, Request
from urllib.error import HTTPError, URLError
import re
import time
import xml.etree.ElementTree as ET

BASE = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils'

# Terms in a record title that indicate an incomplete sequence
PARTIAL_MARKERS = ('partial', 'scaffold', 'contig', 'shotgun',
                   'unverified', 'UNVERIFIED')


@dataclass
class Record:
    accession: str
    title: str
    sequence: str
    length: int = 0

    def __post_init__(self):
        self.length = len(self.sequence)


@dataclass
class FetchReport:
    """Summary of what was retrieved and what was discarded, and why."""
    query: str
    found_total: int = 0
    requested: int = 0
    retrieved: int = 0
    rejected_partial: int = 0
    rejected_too_short: int = 0
    rejected_too_long: int = 0
    rejected_ambiguous: int = 0
    accepted: int = 0
    notes: List[str] = field(default_factory=list)

    def summary(self) -> str:
        return (f'Query returned {self.found_total} records; '
                f'{self.retrieved} downloaded, {self.accepted} accepted. '
                f'Rejected: {self.rejected_partial} partial, '
                f'{self.rejected_too_short} too short, '
                f'{self.rejected_too_long} too long, '
                f'{self.rejected_ambiguous} with excess ambiguous bases.')


class NCBIClient:
    """Client for NCBI E-utilities with enforced rate limiting."""

    def __init__(self, email: str, api_key: Optional[str] = None,
                 tool: str = 'siRNA-design-studio', timeout: int = 60):
        if not email or '@' not in email:
            raise ValueError(
                'A valid e-mail address is required. NCBI uses it to contact '
                'clients that generate excessive load and blocks '
                'unidentified traffic.')
        self.email = email
        self.api_key = api_key
        self.tool = tool
        self.timeout = timeout
        self.min_interval = 0.11 if api_key else 0.34   # margin over limit
        self._last_call = 0.0

    # ---------------- low level ----------------

    def _wait(self) -> None:
        elapsed = time.time() - self._last_call
        if elapsed < self.min_interval:
            time.sleep(self.min_interval - elapsed)
        self._last_call = time.time()

    def _get(self, endpoint: str, params: Dict[str, str],
             retries: int = 3) -> str:
        params = dict(params)
        params['tool'] = self.tool
        params['email'] = self.email
        if self.api_key:
            params['api_key'] = self.api_key

        url = f'{BASE}/{endpoint}.fcgi?{urlencode(params)}'
        last = None
        for attempt in range(retries):
            self._wait()
            try:
                req = Request(url, headers={'User-Agent': self.tool})
                with urlopen(req, timeout=self.timeout) as r:
                    return r.read().decode('utf-8', errors='replace')
            except HTTPError as e:
                last = e
                if e.code in (429, 500, 502, 503):
                    time.sleep(2 ** attempt)     # exponential back-off
                    continue
                raise
            except URLError as e:
                last = e
                time.sleep(2 ** attempt)
        raise RuntimeError(f'NCBI request failed after {retries} attempts: '
                           f'{last}')

    # ---------------- search ----------------

    def count(self, db: str, term: str) -> int:
        """Number of matching records, without downloading them."""
        xml = self._get('esearch', {'db': db, 'term': term, 'retmax': '0'})
        m = re.search(r'<Count>(\d+)</Count>', xml)
        return int(m.group(1)) if m else 0

    def search_ids(self, db: str, term: str, max_records: int = 500,
                   progress: Optional[Callable[[str], None]] = None
                   ) -> Tuple[List[str], int]:
        """
        Retrieves record identifiers. Returns (ids, total_found).

        Identifiers are fetched in batches of 500, the maximum E-utilities
        will return in a single esearch call.
        """
        ids: List[str] = []
        total = self.count(db, term)
        if progress:
            progress(f'{total} records match the query')

        retstart = 0
        while len(ids) < min(max_records, total):
            batch = min(500, max_records - len(ids))
            xml = self._get('esearch', {
                'db': db, 'term': term,
                'retstart': str(retstart), 'retmax': str(batch)})
            found = re.findall(r'<Id>(\d+)</Id>', xml)
            if not found:
                break
            ids.extend(found)
            retstart += len(found)
            if progress:
                progress(f'{len(ids)} identifiers collected')
        return ids[:max_records], total

    # ---------------- fetch ----------------

    def fetch_fasta(self, db: str, ids: List[str], batch_size: int = 200,
                    progress: Optional[Callable[[str], None]] = None
                    ) -> List[Record]:
        """
        Downloads sequences in FASTA format, in batches.

        Batching is essential: one request per record would take hours for
        several hundred sequences and would violate the usage limits.
        """
        records: List[Record] = []
        for i in range(0, len(ids), batch_size):
            chunk = ids[i:i + batch_size]
            txt = self._get('efetch', {
                'db': db, 'id': ','.join(chunk),
                'rettype': 'fasta', 'retmode': 'text'})
            records.extend(_parse_fasta(txt))
            if progress:
                progress(f'{len(records)} / {len(ids)} sequences downloaded')
        return records

    # ---------------- high level ----------------

    def fetch_sequences(self, db: str, term: str,
                        max_records: int = 500,
                        min_length: Optional[int] = None,
                        max_length: Optional[int] = None,
                        exclude_partial: bool = True,
                        max_ambiguous_fraction: float = 0.01,
                        progress: Optional[Callable[[str], None]] = None
                        ) -> Tuple[List[Record], FetchReport]:
        """
        Search, download and filter in a single call.

        Filtering criteria and their rationale:

        exclude_partial
            Records whose title contains 'partial', 'scaffold', 'contig',
            'shotgun' or 'UNVERIFIED'. A partial sequence lowers apparent
            conservation at its ends purely because the region is absent,
            not because it is variable.

        min_length / max_length
            Removes fragments and mis-assigned records. For a viral genome,
            a sensible window is roughly 80 to 120 per cent of the reference
            length.

        max_ambiguous_fraction
            Records with an excess of N. Ambiguous positions cannot be
            evaluated for conservation and effectively lower coverage.
        """
        report = FetchReport(query=term, requested=max_records)

        ids, total = self.search_ids(db, term, max_records, progress)
        report.found_total = total
        if not ids:
            report.notes.append('Query returned no records.')
            return [], report

        raw = self.fetch_fasta(db, ids, progress=progress)
        report.retrieved = len(raw)

        accepted: List[Record] = []
        for r in raw:
            if exclude_partial and any(m.lower() in r.title.lower()
                                       for m in PARTIAL_MARKERS):
                report.rejected_partial += 1
                continue
            if min_length and r.length < min_length:
                report.rejected_too_short += 1
                continue
            if max_length and r.length > max_length:
                report.rejected_too_long += 1
                continue
            n_amb = sum(1 for c in r.sequence if c not in 'ACGT')
            if r.length and n_amb / r.length > max_ambiguous_fraction:
                report.rejected_ambiguous += 1
                continue
            accepted.append(r)

        report.accepted = len(accepted)

        if accepted:
            lengths = [r.length for r in accepted]
            report.notes.append(
                f'Length range of accepted records: '
                f'{min(lengths)}-{max(lengths)} nt '
                f'(median {sorted(lengths)[len(lengths) // 2]}).')
        if report.accepted < report.retrieved * 0.5:
            report.notes.append(
                'More than half of the downloaded records were rejected. '
                'Consider relaxing the filters or refining the query.')

        return accepted, report

    @staticmethod
    def save_fasta(records: List[Record], path: str,
                   line_width: int = 70) -> None:
        with open(path, 'w') as fh:
            for r in records:
                fh.write(f'>{r.accession} {r.title}\n')
                for i in range(0, len(r.sequence), line_width):
                    fh.write(r.sequence[i:i + line_width] + '\n')


def _parse_fasta(text: str) -> List[Record]:
    recs, acc, title, buf = [], None, '', []
    for line in text.splitlines():
        line = line.rstrip()
        if line.startswith('>'):
            if acc is not None:
                recs.append(Record(acc, title, ''.join(buf).upper()))
            head = line[1:]
            parts = head.split(None, 1)
            acc = parts[0]
            title = parts[1] if len(parts) > 1 else ''
            buf = []
        elif acc is not None:
            buf.append(line)
    if acc is not None:
        recs.append(Record(acc, title, ''.join(buf).upper()))
    return recs


# ============================================================================
# QUERY TEMPLATES
# ============================================================================

QUERY_TEMPLATES = {
    'viral_genomes': {
        'template': '{organism}[Organism] AND complete genome[Title]',
        'db': 'nuccore',
        'description': 'Complete genomes of a given virus - the standard '
                       'starting point for conservation analysis',
        'example': 'tomato brown rugose fruit virus',
        'suggested_filters': {'exclude_partial': True,
                              'max_ambiguous_fraction': 0.01},
    },
    'viral_recent': {
        'template': ('{organism}[Organism] AND complete genome[Title] '
                     'AND ("{year_from}"[PDAT] : "3000"[PDAT])'),
        'db': 'nuccore',
        'description': 'Genomes deposited from a given year onwards - '
                       'useful for tracking a currently circulating '
                       'population',
        'example': 'tomato brown rugose fruit virus, 2020',
    },
    'gene_by_name': {
        'template': '{gene}[Gene Name] AND {organism}[Organism]',
        'db': 'nuccore',
        'description': 'A specific gene in a specific organism',
        'example': 'PDS, Solanum lycopersicum',
    },
    'refseq_mrna': {
        'template': ('{organism}[Organism] AND biomol_mrna[PROP] '
                     'AND refseq[filter]'),
        'db': 'nuccore',
        'description': 'RefSeq mRNA records. For a whole transcriptome '
                       'use an FTP download instead - see module docstring.',
    },
    'accession_list': {
        'template': '{accessions}',
        'db': 'nuccore',
        'description': 'Explicit list of accession numbers, separated by '
                       'OR - for reproducing a published dataset',
        'example': 'MN882030 OR MT018320 OR MK648157',
    },
}


def build_query(template_key: str, **kwargs) -> str:
    if template_key not in QUERY_TEMPLATES:
        raise ValueError(f'Unknown template: {template_key}. '
                         f'Available: {sorted(QUERY_TEMPLATES)}')
    return QUERY_TEMPLATES[template_key]['template'].format(**kwargs)


if __name__ == '__main__':
    print(__doc__.split('=' * 77)[1])
    print('AVAILABLE QUERY TEMPLATES')
    print('=' * 74)
    for k, v in QUERY_TEMPLATES.items():
        print(f'\n{k}')
        print(f'  {v["description"]}')
        print(f'  template: {v["template"]}')
        if 'example' in v:
            print(f'  example:  {v["example"]}')

    print('\n' + '=' * 74)
    print('Offline test of the FASTA parser and filters')
    print('=' * 74)
    demo = ('>NC_000001.1 Test virus, complete genome\n'
            'ACGTACGTACGTACGTACGT\n'
            '>XX_000002.1 Test virus, partial cds\n'
            'ACGTACGT\n'
            '>XX_000003.1 Test virus, complete genome\n'
            'ACGTNNNNNNNNNNNNACGT\n')
    recs = _parse_fasta(demo)
    print(f'Parsed {len(recs)} records:')
    for r in recs:
        n_amb = sum(1 for c in r.sequence if c not in 'ACGT')
        partial = any(m.lower() in r.title.lower() for m in PARTIAL_MARKERS)
        print(f'  {r.accession:<14} {r.length:>3} nt  '
              f'ambiguous {n_amb:>2}  '
              f'{"PARTIAL" if partial else "complete"}')
