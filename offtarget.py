"""
=============================================================================
offtarget.py  --  ANALIZA OFF-TARGET (skanowanie strumieniowe)
=============================================================================

=============================================================================
UWAGA O WERSJI POPRZEDNIEJ
=============================================================================

Pierwsza implementacja budowala pelny indeks k-merowy transkryptomu -
slownik odwzorowujacy kazdy 8-mer na liste wszystkich jego pozycji.
Dla transkryptomu Arabidopsis (okolo 50 mln nukleotydow) oznacza to
50 milionow krotek przechowywanych w pamieci, czyli kilka gigabajtow.
Na typowym komputerze prowadzi to do zrzucania pamieci na dysk i pozornego
zawieszenia programu.

Bledem projektowym bylo indeksowanie CALEGO transkryptomu, podczas gdy
interesuje nas jedynie kilkadziesiat sekwencji zapytan.

=============================================================================
ROZWIAZANIE: ODWROCENIE KIERUNKU INDEKSOWANIA
=============================================================================

Zamiast indeksowac transkryptom i odpytywac go zapytaniami, indeksujemy
ZAPYTANIA i przepuszczamy przez nie transkryptom strumieniowo.

    poprzednio:  indeks(transkryptom)  ->  pamiec O(dlugosc transkryptomu)
    obecnie:     indeks(zapytania)     ->  pamiec O(liczba zapytan)

Dla 100 kandydatow o dlugosci 21 nt liczba sond wynosi:
    100 zapytan * 14 osiem-merow = 1400 sond
    + sondy regionu seed
    -> okolo 2000 pozycji w slowniku zamiast 50 milionow.

Zuzycie pamieci spada z kilku gigabajtow do kilku megabajtow, a zlozonosc
czasowa pozostaje liniowa wzgledem dlugosci transkryptomu: przesuwamy okno
o dlugosci 8 i sprawdzamy przynaleznosc do zbioru, co jest operacja o stalym
koszcie.

Transkryptom jest czytany rekord po rekordzie i nigdy nie znajduje sie
w pamieci w calosci.

=============================================================================
CO WYKRYWAMY
=============================================================================

A) TRAFIENIA REGIONU SEED - wiazanie typu miRNA

   8mer     : komplementarnosc pozycji 2-8 nici prowadzacej + A naprzeciw
              pozycji 1  (klasa najsilniejsza)
   7mer-m8  : komplementarnosc pozycji 2-8
   7mer-A1  : komplementarnosc pozycji 2-7 + A naprzeciw pozycji 1

   Klasyfikacja wedlug Bartel 2009 (Cell 136:215).

B) TRAFIENIA PELNEJ DLUGOSCI - ryzyko ciecia

   Zasiew osiem-merem, rozszerzenie w obie strony, punktacja niedopasowan
   z WAGAMI POZYCYJNYMI:

       pozycja 1        waga 0.5   nie paruje z celem, siedzi w kieszeni AGO
       pozycje 2-8      waga 3.0   region seed, decyduje o rozpoznaniu
       pozycje 9-11     waga 2.0   miejsce ciecia
       pozycje 12+      waga 1.0   region 3-prim, wspomaga ale nie decyduje

   Wagi pozycyjne sa glowna przewaga nad blastn, ktory traktuje wszystkie
   pozycje rownowaznie, a przy domyslnym word_size = 11 w ogole nie zasiewa
   dopasowan regionu seed.

Autor: Antonina Jarecka
=============================================================================
"""

from collections import defaultdict
from typing import Callable, Dict, Iterator, List, Optional, Tuple
import os

COMP = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'U': 'A', 'N': 'N'}

DLUGOSC_SONDY = 8      # dlugosc k-meru zasiewajacego


def position_weight(pos_1based: int) -> float:
    if pos_1based == 1:
        return 0.5
    if 2 <= pos_1based <= 8:
        return 3.0
    if 9 <= pos_1based <= 11:
        return 2.0
    return 1.0


def to_dna(s: str) -> str:
    return s.upper().replace('U', 'T')


def revcomp(s: str) -> str:
    return ''.join(COMP.get(c, 'N') for c in reversed(s.upper()))


def read_fasta(path: str) -> Iterator[Tuple[str, str]]:
    """Czyta FASTA rekord po rekordzie. Nie trzyma pliku w pamieci."""
    nazwa, buf = None, []
    with open(path) as fh:
        for linia in fh:
            linia = linia.rstrip()
            if linia.startswith('>'):
                if nazwa is not None:
                    yield nazwa, ''.join(buf).upper()
                nazwa, buf = linia[1:].split()[0], []
            else:
                buf.append(linia)
    if nazwa is not None:
        yield nazwa, ''.join(buf).upper()


def rozmiar_pliku_mb(path: str) -> float:
    try:
        return os.path.getsize(path) / (1024 * 1024)
    except OSError:
        return 0.0


# ============================================================================
# SKANER
# ============================================================================

class OffTargetScanner:
    """
    Skaner off-target o stalym zuzyciu pamieci.

    Buduje slownik sond WYLACZNIE z sekwencji zapytan, a nastepnie przepuszcza
    transkryptom strumieniowo.
    """

    def __init__(self, guides: Dict[str, str], max_koszt: float = 12.0,
                 mrna_celu: Optional[str] = None,
                 kontekst_celu: int = 45):
        """
        guides       : {nazwa: nic_prowadzaca}  (RNA lub DNA, 5'->3')
        mrna_celu    : sekwencja mRNA, przeciwko ktorej projektowano siRNA.
                       Sluzy do ODROZNIENIA trafienia zamierzonego od
                       off-target - patrz nizej.
        kontekst_celu: dlugosc otoczenia sprawdzanego przy rozpoznawaniu
                       trafienia zamierzonego

        ROZPOZNAWANIE TRAFIENIA ZAMIERZONEGO

        Jesli gen docelowy wystepuje w transkryptomie gospodarza - a tak jest
        zawsze, gdy celem jest gen endogenny - kazdy kandydat trafia we
        wlasny cel z kosztem zero. Bez rozroznienia trafienia zamierzonego
        od niezamierzonego wszystkie kandydaty zostalyby odrzucone.

        Rozpoznanie: dla kazdego trafienia sprawdzamy, czy jego otoczenie
        w transkrypcie wystepuje rowniez w mRNA celu. Jesli tak, jest to
        ten sam gen (albo jego bliski paralog) i trafienie klasyfikujemy
        jako ZAMIERZONE.
        """
        self.guides = {n: to_dna(g) for n, g in guides.items()}
        self.max_koszt = max_koszt
        self.kontekst_celu = kontekst_celu
        self.mrna_celu = to_dna(mrna_celu) if mrna_celu else None
        self._kmery_celu = None
        if self.mrna_celu:
            kk = 24
            self._kmery_celu = {self.mrna_celu[i:i + kk]
                                for i in range(len(self.mrna_celu) - kk + 1)}
            self._dl_kmeru_celu = kk

        # sonda -> lista (nazwa_zapytania, przesuniecie_w_celu)
        self.sondy_pelne: Dict[str, List[Tuple[str, int]]] = defaultdict(list)
        # sonda -> lista (nazwa_zapytania, klasa_seed)
        self.sondy_seed: Dict[str, List[Tuple[str, str]]] = defaultdict(list)
        # nazwa -> idealna sekwencja mRNA
        self.cele: Dict[str, str] = {}

        self._zbuduj_sondy()

    def _zbuduj_sondy(self) -> None:
        for nazwa, g in self.guides.items():
            cel = revcomp(g)              # idealne miejsce docelowe w mRNA
            self.cele[nazwa] = cel

            # --- sondy do zasiewu trafien pelnej dlugosci ---
            for off in range(len(cel) - DLUGOSC_SONDY + 1):
                self.sondy_pelne[cel[off:off + DLUGOSC_SONDY]].append(
                    (nazwa, off))

            # --- sondy regionu seed ---
            seed_2_8 = g[1:8]
            seed_2_7 = g[1:7]
            # Wszystkie sondy sprowadzamy do dlugosci 8, zeby w petli
            # skanujacej wystarczylo JEDNO wyszukanie w slowniku na pozycje.
            # Sondy 7-nt rozwijamy na cztery warianty przez dopisanie
            # kazdego mozliwego nukleotydu.
            self.sondy_seed[revcomp(seed_2_8) + 'A'].append((nazwa, '8mer'))
            for b in 'ACGT':
                self.sondy_seed[revcomp(seed_2_8) + b].append(
                    (nazwa, '7mer_m8'))
                self.sondy_seed[revcomp(seed_2_7) + 'A' + b].append(
                    (nazwa, '7mer_A1'))

    @property
    def liczba_sond(self) -> int:
        return len(self.sondy_pelne) + len(self.sondy_seed)

    def _oceń(self, nazwa: str, okno: str) -> Tuple[float, int, List[int]]:
        """Wazony koszt niedopasowania miedzy zapytaniem a oknem mRNA."""
        cel = self.cele[nazwa]
        L = len(cel)
        koszt, n_mm, pozycje = 0.0, 0, []
        for j in range(L):
            if okno[j] != cel[j]:
                poz_guide = L - j     # antyrownolegle
                koszt += position_weight(poz_guide)
                n_mm += 1
                pozycje.append(poz_guide)
        return koszt, n_mm, sorted(pozycje)

    def _czy_cel_zamierzony(self, seq: str, start: int, dl: int) -> bool:
        """
        Sprawdza, czy trafienie pochodzi z genu docelowego.

        Pobiera otoczenie trafienia i sprawdza, czy ktorykolwiek jego
        24-mer wystepuje w mRNA celu. Uzycie 24-merow zamiast calego
        okna daje odpornosc na roznice w wariantach transkryptu
        i na granice trafienia.
        """
        if self._kmery_celu is None:
            return False
        kk = self._dl_kmeru_celu
        lo = max(0, start - self.kontekst_celu)
        hi = min(len(seq), start + dl + self.kontekst_celu)
        otoczenie = seq[lo:hi]
        for i in range(len(otoczenie) - kk + 1):
            if otoczenie[i:i + kk] in self._kmery_celu:
                return True
        return False

    def scan(self, fasta_path: str,
             progress: Optional[Callable[[int, str], None]] = None,
             limit_trafien_na_zapytanie: int = 500) -> Dict[str, Dict]:
        """
        Jednokrotne przejscie przez transkryptom.

        progress : funkcja wywolywana co 2000 transkryptow, otrzymuje
                   (liczba_przetworzonych, nazwa_ostatniego)
        limit_trafien_na_zapytanie : gorne ograniczenie liczby zapamietanych
                   trafien pelnej dlugosci; chroni pamiec przy sekwencjach
                   powtarzalnych. Licznik zlicza wszystkie trafienia
                   niezaleznie od limitu.
        """
        wyniki = {n: {'guide': g, 'seed_2_8': g[1:8],
                      'n_8mer': 0, 'n_7mer_m8': 0, 'n_7mer_A1': 0,
                      'geny_8mer': [], 'n_trafien_pelnych': 0,
                      'n_trafien_zamierzonych': 0,
                      'n_krytycznych': 0, 'geny_krytyczne': [],
                      'min_koszt': None, 'najgorszy_transkrypt': None,
                      'trafienia': []}
                   for n, g in self.guides.items()}

        widziane = {n: set() for n in self.guides}
        k = DLUGOSC_SONDY
        n_tx = 0
        laczna_dlugosc = 0

        for naglowek, seq in read_fasta(fasta_path):
            n_tx += 1
            if progress and n_tx % 2000 == 0:
                progress(n_tx, naglowek)

            L = len(seq)
            laczna_dlugosc += L
            if L < k:
                continue

            # Odwolania lokalne - w Pythonie dostep do zmiennej lokalnej
            # jest istotnie szybszy niz do atrybutu obiektu, a petla wykonuje
            # sie tyle razy, ile transkryptom ma nukleotydow.
            get_seed = self.sondy_seed.get
            get_pelne = self.sondy_pelne.get

            for i in range(L - k + 1):
                okno_k = seq[i:i + k]

                # --- trafienia regionu seed ---
                dopasowania_seed = get_seed(okno_k)
                if dopasowania_seed is not None:
                    if self._czy_cel_zamierzony(seq, i, k):
                        pass          # trafienie w gen docelowy - pomijamy
                    else:
                      for nazwa, klasa in dopasowania_seed:
                        w = wyniki[nazwa]
                        w[f'n_{klasa}'] += 1
                        if klasa == '8mer' and len(w['geny_8mer']) < 20:
                            if naglowek not in w['geny_8mer']:
                                w['geny_8mer'].append(naglowek)

                # --- trafienia pelnej dlugosci ---
                dopasowania = get_pelne(okno_k)
                if dopasowania is None:
                    continue
                for nazwa, off in dopasowania:
                    start = i - off
                    dl = len(self.cele[nazwa])
                    if start < 0 or start + dl > L:
                        continue
                    klucz = (n_tx, start)
                    if klucz in widziane[nazwa]:
                        continue
                    widziane[nazwa].add(klucz)

                    okno = seq[start:start + dl]
                    koszt, n_mm, poz = self._oceń(nazwa, okno)
                    if koszt > self.max_koszt:
                        continue

                    w = wyniki[nazwa]

                    # trafienie w gen docelowy nie jest off-target
                    if self._czy_cel_zamierzony(seq, start, dl):
                        w['n_trafien_zamierzonych'] += 1
                        continue

                    w['n_trafien_pelnych'] += 1

                    # Trafienie KRYTYCZNE: najwyzej 3 niedopasowania i zaden
                    # w regionie seed. Tylko takie niosa realne ryzyko ciecia.
                    # Prawdopodobienstwo przypadkowego wystapienia takiego
                    # dopasowania nawet w genomie o rozmiarze 10^8 jest
                    # pomijalne, wiec kryterium nie wymaga poprawki na
                    # wielokrotne testowanie.
                    if n_mm <= 3 and not any(2 <= p <= 8 for p in poz):
                        w['n_krytycznych'] += 1
                        if len(w['geny_krytyczne']) < 20:
                            w['geny_krytyczne'].append(naglowek)

                    if w['min_koszt'] is None or koszt < w['min_koszt']:
                        w['min_koszt'] = round(koszt, 2)
                        w['najgorszy_transkrypt'] = naglowek
                    if len(w['trafienia']) < limit_trafien_na_zapytanie:
                        w['trafienia'].append({
                            'transkrypt': naglowek, 'pozycja': start,
                            'koszt': round(koszt, 2), 'n_niedopasowan': n_mm,
                            'pozycje_niedopasowan': poz})

        # oczekiwana liczba trafien przypadkowych dla sondy 8-nt
        oczek_8mer = laczna_dlugosc / (4 ** 8)

        for n, w in wyniki.items():
            w['n_transkryptow'] = n_tx
            w['laczna_dlugosc'] = laczna_dlugosc
            w['oczek_8mer_losowo'] = round(oczek_8mer, 1)
            w['wzbogacenie_8mer'] = (round(w['n_8mer'] / oczek_8mer, 2)
                                     if oczek_8mer > 0 else None)
            w['flaga'] = _flaga(w['n_krytycznych'], w['n_8mer'], oczek_8mer)
            w['trafienia'].sort(key=lambda t: t['koszt'])

        return wyniki


def _flaga(n_krytycznych: int, n_8mer: int, oczek_8mer: float) -> str:
    """
    Klasyfikacja z poprawka na rozmiar przeszukiwanej bazy.

    UZASADNIENIE

    Wczesniejsza wersja odrzucala kandydata przy jakimkolwiek trafieniu
    klasy 8mer oraz przy minimalnym koszcie ponizej progu. Oba kryteria
    zawodza przy duzych bazach:

      - dokladne dopasowanie 8 nt wystepuje losowo raz na 4^8 = 65 536
        pozycji, czyli okolo 1800 razy w bazie o dlugosci 1.2 * 10^8;
      - minimum kosztu wyznaczane po dziesiatkach tysiecy pozycji losowych
        jest z koniecznosci niskie - jest to problem wielokrotnego
        testowania.

    Obecne kryterium opiera sie na trafieniach KRYTYCZNYCH: najwyzej trzy
    niedopasowania i zaden w regionie seed. Oczekiwana liczba takich
    dopasowan w bazie o dlugosci 10^8 jest pomijalna, wiec ich occurrence
    is a real signal.

    The number of 8mer hits is instead reported as ENRICHMENT relative to
    chance expectation, rather than as an absolute value.
    """
    if n_krytycznych > 0:
        return 'REJECT (off-target cleavage risk)'
    if oczek_8mer > 0 and n_8mer > max(3 * oczek_8mer, 5):
        return 'REVIEW (excess seed hits)'
    return 'OK'


# ============================================================================
# INTERFEJS WYSOKIEGO POZIOMU
# ============================================================================

def analyze_guides(guides: Dict[str, str], fasta_path: str,
                   max_koszt: float = 12.0,
                   mrna_celu: Optional[str] = None,
                   progress: Optional[Callable[[int, str], None]] = None
                   ) -> Dict[str, Dict]:
    """
    Pelna analiza off-target dla slownika {nazwa: nic_prowadzaca}.

    Jedno przejscie przez plik, stala pamiec.

    mrna_celu : jesli podane, trafienia w gen docelowy sa rozpoznawane
                i wylaczane z oceny off-target.
    """
    return OffTargetScanner(guides, max_koszt=max_koszt,
                            mrna_celu=mrna_celu).scan(
        fasta_path, progress=progress)


def szacuj_czas(fasta_path: str) -> Dict[str, float]:
    """
    Szacuje czas skanowania na podstawie rozmiaru pliku.

    Podstawa: pomiar rzeczywisty, okolo 0.4-0.7 mln nukleotydow na sekunde
    w czystym Pythonie (bez bibliotek kompilowanych).
    """
    mb = rozmiar_pliku_mb(fasta_path)
    nt_mln = mb * 0.97          # FASTA to niemal same nukleotydy
    return {'rozmiar_mb': round(mb, 1),
            'szac_nt_mln': round(nt_mln, 1),
            'szac_sekund_min': round(nt_mln / 0.7, 0),
            'szac_sekund_max': round(nt_mln / 0.4, 0),
            'szac_minut': round(nt_mln / 0.55 / 60, 1)}


if __name__ == '__main__':
    import time

    # mini-transkryptom demonstracyjny
    with open('/tmp/demo_tx.fasta', 'w') as fh:
        fh.write('>TX_bez_trafien\n')
        fh.write('ATGCGCGCGATATATCGCGCGATATATGCGCGCATATATCGCGCGATATAT\n')
        fh.write('>TX_z_seedem\n')
        fh.write('GGGGCCCCAAACTTCAAGGGGCCCCTTTTAAAACCCCGGGG\n')
        fh.write('>TX_z_pelnym_celem\n')
        fh.write('AAAAAACGGCATCAAGGTGAACTTCAAGGGGGG\n')

    guides = {'A-21': 'UUGAAGUUCACCUUGAUGCCG',
              'B-21': 'UAGUUGUACUCCAGCUUGUGC'}

    sk = OffTargetScanner(guides)
    print(f'Zapytan: {len(guides)}, sond w slowniku: {sk.liczba_sond}')
    print('(dla porownania: indeks pelnego transkryptomu Arabidopsis '
          'wymagalby okolo 50 mln pozycji)')

    t0 = time.time()
    rap = sk.scan('/tmp/demo_tx.fasta')
    print(f'Czas: {time.time() - t0:.3f} s\n')

    for n, r in rap.items():
        print(f'{n}: seed={r["seed_2_8"]}  8mer={r["n_8mer"]}  '
              f'7mer_m8={r["n_7mer_m8"]}  pelne={r["n_trafien_pelnych"]}  '
              f'min_koszt={r["min_koszt"]}')
        print(f'   geny z 8mer: {r["geny_8mer"] or "brak"}')
        print(f'   FLAGA: {r["flaga"]}')
