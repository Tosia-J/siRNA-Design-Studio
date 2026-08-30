"""
=============================================================================
design.py  --  GENEROWANIE KANDYDATOW I ANALIZA STRUKTURY
=============================================================================

STATUS: HYBRYDA -- wlasna logika + ViennaRNA jako biblioteka

=============================================================================
DLACZEGO NIE ZASTEPUJEMY ViennaRNA (a BLAST zastepujemy)
=============================================================================

Pytalas o zastapienie zewnetrznych programow wlasnymi obliczeniami. Odpowiedz
jest rozna dla roznych programow i warto rozumiec, na czym polega roznica.

siRNA Wizard  ->  ZASTEPUJEMY W CALOSCI
    To zbior regul heurystycznych (GC, homopolimery, asymetria, preferencja
    5-konca). Reguly sa opublikowane, wiec da sie je zaimplementowac
    samodzielnie i miec nad nimi pelna kontrole. Dodatkowo Wizard stosuje
    reguly SSACZE -- filtr motywow immunostymulacyjnych, ktory w roslinach
    nie ma sensu (patrz filters.py). Wlasna implementacja jest tu nie tylko
    mozliwa, ale LEPSZA, bo dostosowana do gospodarza.

BLAST  ->  ZASTEPUJEMY (offtarget.py)
    Nie dlatego, ze BLAST jest zly, tylko dlatego, ze jest narzedziem do
    innego zadania. Dla zapytan 21-nt jego word_size i filtr E-value
    systematycznie gubia trafienia istotne biologicznie. Wlasny indeks
    k-merowy jest tu szybszy, dokladniejszy i pozwala na wagi pozycyjne,
    ktorych BLAST nie ma. To jest realne ulepszenie, nie tylko niezaleznosc.

ViennaRNA / RNAfold  ->  NIE ZASTEPUJEMY, uzywamy jako biblioteki
    I to jest wazna decyzja, ktora trzeba umiec obronic.

    RNAfold rozwiazuje problem MINIMALIZACJI ENERGII SWOBODNEJ po przestrzeni
    wszystkich mozliwych struktur drugorzednych. Dla sekwencji dlugosci n
    liczba struktur rosnie wykladniczo, wiec robi sie to programowaniem
    dynamicznym (algorytm Zukera, O(n^3) czasu, O(n^2) pamieci), a rozklad
    Boltzmanna po strukturach -- algorytmem McCaskilla.

    Napisanie tego samodzielnie oznacza zaimplementowanie:
      - pelnego modelu energetycznego Turnera 2004 (kilkaset parametrow dla
        petli wewnetrznych, wybrzuszen, petli spinkowych, koaksjalnego
        stackingu, terminalnych niesparowan)
      - rekurencji Zukera z obsluga wszystkich typow petli
      - algorytmu McCaskilla, jesli chcemy prawdopodobienstwa parowania

    To jest praca na miesiace, a wynik bylby w najlepszym razie identyczny
    z ViennaRNA, a realnie -- gorszy, bo parametry Turnera sa mozolnie
    dopasowane do danych eksperymentalnych i latwo je zle przepisac.

    ZASADA OGOLNA: implementuj samodzielnie tam, gdzie masz przewage
    (znajomosc problemu, mozliwosc dostosowania, prostota zadania). Uzywaj
    bibliotek tam, gdzie standardowe narzedzie jest referencyjna
    implementacja modelu naukowego. ViennaRNA to drugi przypadek.

    Co ZA TO implementujemy samodzielnie w warstwie termodynamicznej:
    dupleks o znanej, pelnej komplementarnosci (thermo.py). Tam nie ma
    przeszukiwania, wiec zadanie jest proste i wlasna implementacja daje
    pelna kontrole nad asymetria -- najwazniejsza metryka w calym pipelinie.

    W metodyce zapisz to tak: "Faldowanie struktur drugorzednych wykonano
    przy uzyciu ViennaRNA 2.7 (Lorenz i wsp. 2011) jako referencyjnej
    implementacji modelu Turnera. Termodynamike dupleksu, analize
    off-target i punktacje zaimplementowano samodzielnie."

Autor: Antonina Jarecka
=============================================================================
"""

from typing import Dict, List, Optional, Tuple
import math
import warnings

# --- ViennaRNA: uzywane swiadomie (patrz uzasadnienie wyzej), ale opcjonalne
# --- Na Windows instalacja bywa klopotliwa (bioconda nie ma buildow win-64).
# --- Bez ViennaRNA pipeline dziala w trybie ograniczonym: pomija metryki
# --- strukturalne, reszta liczona normalnie.
try:
    import RNA
    VIENNA_DOSTEPNE = True
except ImportError:
    RNA = None
    VIENNA_DOSTEPNE = False
    warnings.warn(
        'ViennaRNA niedostepne. Pipeline dziala w TRYBIE OGRANICZONYM:\n'
        '  - MFE nici prowadzacej: nieliczone (wartosc neutralna 0.0)\n'
        '  - dostepnosc miejsca docelowego: nieliczona (wartosc neutralna)\n'
        '  - pozostale metryki (asymetria, GC, off-target): BEZ ZMIAN\n'
        'Ranking pozostaje uzyteczny, ale dwie z szesciu metryk sa wylaczone.\n'
        'Instalacja: pip install ViennaRNA',
        RuntimeWarning)

from thermo import (to_rna, revcomp_rna, gc_content, has_homopolymer,
                    asymmetry, duplex_dg, tm)
from filters import physicochemical_filters, check_mammalian_motifs

# Klasy dlugosci istotne w roslinach
DLUGOSCI_ROSLINNE = {
    21: 'DCL4 -- PTGS, glowna sciezka przeciwwirusowa, ciecie mRNA',
    22: 'DCL2 -- wyzwala transitivity (RDR6), wtorne siRNA, ruch systemiczny',
    24: 'DCL3 -- RdDM, metylacja DNA, wyciszanie transkrypcyjne',
}


# ============================================================================
# STRUKTURA DRUGORZEDOWA (ViennaRNA)
# ============================================================================

def fold_guide(guide_rna: str) -> Tuple[str, float]:
    """
    Struktura wlasna nici prowadzacej.

    Silnie sfaldowana nic prowadzaca nie zwiaze celu, bo jest zajeta sama
    soba. Prog roboczy: MFE > -3 kcal/mol.
    """
    if not VIENNA_DOSTEPNE:
        return '.' * len(guide_rna), 0.0    # wartosc neutralna
    ss, mfe = RNA.fold(to_rna(guide_rna))
    return ss, mfe


def dostepnosc_celu(mrna: str, start_0based: int, dlugosc: int,
                    okno: int = 25) -> float:
    """
    Dostepnosc miejsca docelowego, mierzona jako MFE lokalnego okna mRNA
    znormalizowane na nukleotyd.

    Logika: jesli okolica miejsca docelowego jest silnie ustrukturyzowana
    (bardzo ujemne MFE), RISC bedzie mial trudniej dotrzec do celu.
    Wartosc MNIEJ ujemna = miejsce bardziej dostepne = lepiej.

    UWAGA METODYCZNA: to jest przyblizenie. Scisle rzecz biorac nalezaloby
    liczyc prawdopodobienstwo niesparowania kazdej pozycji z rozkladu
    Boltzmanna (RNA.pfl_fold), co jest dokladniejsze, ale wolniejsze.
    Funkcja dostepnosc_celu_boltzmann() ponizej robi to dokladnie -- warto
    jej uzyc w wersji finalnej.
    """
    if not VIENNA_DOSTEPNE:
        return 0.0                          # wartosc neutralna
    lo = max(0, start_0based - okno)
    hi = min(len(mrna), start_0based + dlugosc + okno)
    win = to_rna(mrna[lo:hi])
    _, mfe = RNA.fold(win)
    return mfe / len(win)


def dostepnosc_celu_boltzmann(mrna: str, start_0based: int,
                              dlugosc: int) -> float:
    """
    Dokladniejsza dostepnosc: srednie prawdopodobienstwo NIESPAROWANIA
    pozycji w regionie docelowym, liczone z rozkladu Boltzmanna.

    Zwraca wartosc w [0,1]. Wieksza = bardziej dostepne.
    """
    if not VIENNA_DOSTEPNE:
        return 0.5                          # wartosc neutralna
    seq = to_rna(mrna)
    fc = RNA.fold_compound(seq)
    fc.pf()                       # funkcja podzialu
    bpp = fc.bpp()                # macierz prawdopodobienstw parowania
    n = len(seq)

    sumy = 0.0
    for i in range(start_0based + 1, start_0based + dlugosc + 1):  # 1-based
        p_sparowana = 0.0
        for j in range(1, n + 1):
            if i < j:
                p_sparowana += bpp[i][j]
            elif j < i:
                p_sparowana += bpp[j][i]
        sumy += max(0.0, 1.0 - p_sparowana)
    return sumy / dlugosc


# ============================================================================
# GENEROWANIE KANDYDATOW
# ============================================================================

def generuj_kandydatow(mrna: str,
                       dlugosci: Tuple[int, ...] = (21, 22, 24),
                       wyklucz_5prim_nt: int = 75,
                       wyklucz_3prim_nt: int = 50,
                       gc_min: float = 30.0,
                       gc_max: float = 52.0,
                       min_asymetria: float = 0.3,
                       max_mfe_guide: float = -3.0,
                       uzyj_boltzmann: bool = False,
                       verbose: bool = True) -> List[Dict]:
    """
    Generuje i filtruje kandydatow siRNA dla podanego mRNA.

    Zwraca liste slownikow z metrykami gotowymi do modulu scoringu.
    """
    mrna = mrna.upper().replace('U', 'T')
    wyniki: List[Dict] = []
    statystyki = {d: {'wszystkie': 0, 'po_filtrach': 0} for d in dlugosci}

    for L in dlugosci:
        lo_bound = wyklucz_5prim_nt
        hi_bound = len(mrna) - L - wyklucz_3prim_nt

        if hi_bound <= lo_bound:
            if verbose:
                print(f'[!] Dlugosc {L}: po wykluczeniu koncow nie zostaje '
                      f'zadne okno. ORF ma {len(mrna)} nt -- poluzuj '
                      f'wyklucz_5prim_nt/wyklucz_3prim_nt.')
            continue

        for i in range(lo_bound, hi_bound + 1):
            statystyki[L]['wszystkie'] += 1
            sense = mrna[i:i + L]
            if set(sense) - set('ACGT'):
                continue

            ok, powody = physicochemical_filters(
                sense, gc_min=gc_min, gc_max=gc_max,
                min_asymmetry=min_asymetria)
            if not ok:
                continue

            guide = revcomp_rna(sense)
            ss_g, mfe_g = fold_guide(guide)
            if VIENNA_DOSTEPNE and mfe_g < max_mfe_guide:
                continue

            if uzyj_boltzmann:
                dost = dostepnosc_celu_boltzmann(mrna, i, L)
            else:
                dost = dostepnosc_celu(mrna, i, L)

            dg_g, dg_p, asym = asymmetry(guide)
            g_pct = gc_content(sense)

            statystyki[L]['po_filtrach'] += 1
            wyniki.append({
                'nazwa': f'{L}nt_poz{i + 1}',
                'dlugosc': L,
                'klasa_DCL': DLUGOSCI_ROSLINNE.get(L, 'undefined'),
                'poz_start_1based': i + 1,
                'poz_end_1based': i + L,
                'sense_dna': sense,
                'guide_rna': guide,
                'passenger_rna': to_rna(sense),
                'seed_2_8': guide[1:8],
                'gc_proc': round(g_pct, 1),
                'dG_dupleksu': round(duplex_dg(guide), 2),
                'Tm': round(tm(guide), 1),
                'dG5_guide': round(dg_g, 2),
                'dG5_passenger': round(dg_p, 2),

                # --- metryki do scoringu ---
                'asymetria': round(asym, 3),
                'dostepnosc_celu': round(dost, 4),
                'mfe_guide': round(mfe_g, 2),
                'gc_optymalnosc': round(1.0 - abs(g_pct - 42.0) / 42.0, 4),
                'koszt_offtarget': None,    # uzupelnia offtarget.py
                'seed_czystosc': None,      # uzupelnia offtarget.py

                'struktura_guide': ss_g,
                'motywy_ssacze': check_mammalian_motifs(guide),
            })

    if verbose:
        if not VIENNA_DOSTEPNE:
            print('  [TRYB OGRANICZONY] ViennaRNA niedostepne - metryki '
                  'strukturalne pominiete')
        print('GENERATION STATISTICS')
        for L in dlugosci:
            s = statystyki[L]
            pct = (100 * s['po_filtrach'] / s['wszystkie']
                   if s['wszystkie'] else 0)
            print(f'  {L} nt: {s["po_filtrach"]:>4} / {s["wszystkie"]:>4} '
                  f'przeszlo filtry ({pct:.1f} %)  --  {DLUGOSCI_ROSLINNE.get(L, "")}')

    return wyniki


def serie_zagniezdzone(kandydaci: List[Dict],
                       tolerancja_poz: int = 3) -> List[Dict]:
    """
    Znajduje zestawy kandydatow o roznych dlugosciach celujace w TO SAMO
    miejsce -- czyli serie nadajace sie do porownania dlugosci.

    Bez tego porownanie 21 vs 22 vs 24 miesza dwie zmienne: dlugosc i
    kontekst sekwencyjny. Serie zagniezdzone eliminuja te druga.

    Zwraca liste serii, kazda jako slownik {dlugosc: kandydat}.
    """
    wg_dlugosci: Dict[int, List[Dict]] = {}
    for k in kandydaci:
        wg_dlugosci.setdefault(k['dlugosc'], []).append(k)

    dlugosci = sorted(wg_dlugosci)
    if len(dlugosci) < 2:
        return []

    baza = wg_dlugosci[dlugosci[0]]
    serie = []
    for k0 in baza:
        seria = {dlugosci[0]: k0}
        kompletna = True
        for L in dlugosci[1:]:
            dopasowany = None
            for k in wg_dlugosci[L]:
                if abs(k['poz_start_1based'] - k0['poz_start_1based']) <= tolerancja_poz:
                    dopasowany = k
                    break
            if dopasowany is None:
                kompletna = False
                break
            seria[L] = dopasowany
        if kompletna:
            serie.append(seria)
    return serie


if __name__ == '__main__':
    GFP = ("ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAG"
           "TTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTG"
           "CCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAG"
           "CACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAG"
           "ACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGC"
           "AACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATC"
           "AAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATC"
           "GGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGC"
           "GATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAA")

    kand = generuj_kandydatow(GFP)
    print(f'\nLacznie {len(kand)} kandydatow.\n')

    serie = serie_zagniezdzone(kand)
    print(f'Znaleziono {len(serie)} serii zagniezdzonych '
          f'(to samo miejsce, rozne dlugosci):')
    for s in serie:
        poz = s[21]['poz_start_1based']
        print(f'  poz ~{poz}: ' + ', '.join(
            f'{L}nt(GC {s[L]["gc_proc"]}%, asym {s[L]["asymetria"]:+.2f})'
            for L in sorted(s)))
