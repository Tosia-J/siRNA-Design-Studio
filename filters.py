"""
=============================================================================
filters.py  --  FILTRY SEKWENCYJNE
=============================================================================

STATUS: WLASNA IMPLEMENTACJA -- zastepuje reguly siRNA Wizard

=============================================================================
KOREKTA MERYTORYCZNA: "SEKWENCJE IMMUNOGENNE" W ROSLINACH
=============================================================================

W planie mialas krok "wykluczenie sekwencji immunogennych". Ten filtr, w
postaci w jakiej wystepuje w narzedziach typu siRNA Wizard czy Dharmacon,
JEST FILTREM SSACZYM I NIE MA ZASTOSOWANIA W ROSLINACH.

Skad sie wzial: u ssakow jednoniciowe RNA jest rozpoznawane przez receptory
TLR7 i TLR8 w endosomach, a dwuniciowe przez TLR3, PKR i RIG-I. Uruchamia to
odpowiedz interferonowa. Zidentyfikowano konkretne motywy nasilajace to
rozpoznanie:
    5'-UGUGU-3'        (Judge i wsp. 2005, Nat Biotechnol 23:457)
    5'-GUCCUUCAA-3'    (Hornung i wsp. 2005, Nat Med 11:263)
    sekwencje bogate w UGU
Odfiltrowanie ich ma sens, gdy projektujesz terapeutyk dla czlowieka.

Dlaczego to nie dotyczy roslin:
    - Rosliny NIE MAJA receptorow Toll-like. Rodzina TLR jest zwierzeca.
    - Rosliny NIE MAJA interferonow ani szlaku interferonowego.
    - Rosliny NIE MAJA PKR ani RIG-I.

Zostawienie tego filtru w pipelinie roslinnym oznacza odrzucanie dobrych
kandydatow z powodu, ktory w tym organizmie nie istnieje. Recenzent, ktory
to zauwazy, slusznie zapyta, czy rozumiesz, co robi Twoj wlasny kod.

=============================================================================
CO JEST ODPOWIEDNIKIEM W ROSLINACH
=============================================================================

Rosliny maja wlasny mechanizm rozpoznawania dsRNA, ale dziala inaczej:

1. dsRNA-PTI -- roslina rozpoznaje dwuniciowe RNA jako wzorzec zagrozenia i
   uruchamia odpornosc wyzwalana wzorcami. KLUCZOWE: ten mechanizm jest
   NIEZALEZNY OD SEKWENCJI. Nie da sie go odfiltrowac na poziomie sekwencji,
   bo reaguje na sama dwuniciowosc.
   -> Kontrola eksperymentalna (dsRNA o obojetnej sekwencji), nie filtr.

2. KONKURENCJA O AGO -- to jest realne ryzyko sekwencyjne i JEST filtrowalne.
   Jesli Twoja nic prowadzaca ma seed identyczny z endogennym miRNA
   gospodarza, moze:
     a) konkurowac z tym miRNA o zaladowanie do AGO1 (efekt: rozregulowanie
        naturalnej regulacji rozwojowej rosliny),
     b) wyciszac naturalne cele tego miRNA (efekt: fenotyp niezwiazany
        z Twoim celem).
   Objawy: karlowatosc, deformacje lisci, zaburzenia kwitnienia. Latwo
   pomylic je z toksycznoscia nosnika.

3. KOLIZJA Z PREKURSORAMI miRNA -- jesli konstrukt jest komplementarny do
   pre-miRNA, moze zaburzyc jego dojrzewanie.

TEN MODUL FILTRUJE PUNKTY 2 I 3, a punkt 1 zostawia kontroli eksperymentalnej.

=============================================================================
POZOSTALE FILTRY I ICH UZASADNIENIE
=============================================================================

GC 30-52 %
    Ponizej 30 %: dupleks za slaby, siRNA nie utrzyma sie na celu.
    Powyzej 52 %: dupleks za silny, RISC nie rozplecie nici; dodatkowo rosnie
    ryzyko struktur drugorzednych w samej nici prowadzacej.
    Zakres wezszy niz w narzedziach ssaczych (tam zwykle 30-64 %), bo AGO1
    roslinne jest bardziej wrazliwe na stabilnosc dupleksu.

5'-koniec nici prowadzacej = U lub A
    AGO1 roslinne ma kieszen MID selektywna wobec 5'-U. AGO2 preferuje 5'-A,
    AGO5 preferuje 5'-C. Poniewaz szlak przeciwwirusowy w roslinach idzie
    glownie przez AGO1 i AGO2, wymuszamy U albo A.
    Zrodlo: Mi i wsp. 2008 (Cell 133:116), Takeda i wsp. 2008 (Plant Cell
    Physiol 49:493).

5'-koniec nici PASAZERSKIEJ = G lub C
    Odwrotnosc powyzszego -- wzmacnia asymetrie w pozadanym kierunku.

Brak ciagow >=4 identycznych nukleotydow
    Utrudniaja synteze chemiczna (bledy przy sprzeganiu), a ciagi G moga
    tworzyc struktury G-kwadrupleksowe.

Wykluczenie pierwszych 75 nt i ostatnich 50 nt ORF
    Poczatek ORF jest oslaniany przez kompleks inicjacyjny translacji, koniec
    przez czynniki terminacyjne. Dostepnosc dla RISC jest tam nizsza.
    Uwaga: to reguła empiryczna, nie prawo. Przy krotkich ORF trzeba ja
    poluzowac, bo inaczej nie zostanie nic do wyboru.

Autor: Antonina Jarecka
=============================================================================
"""

from typing import Dict, List, Optional, Set, Tuple
from thermo import to_rna, revcomp_rna, gc_content, has_homopolymer, asymmetry

# --- Motywy immunostymulacyjne SSACZE -- trzymamy WYLACZNIE po to, zeby
# --- moc jawnie zaraportowac, ze ich NIE stosujemy w trybie roslinnym.
MAMMALIAN_IMMUNO_MOTIFS = ('UGUGU', 'GUCCUUCAA', 'UGGC')


def check_mammalian_motifs(guide: str) -> List[str]:
    """
    Zwraca liste znalezionych motywow ssaczych.

    NIE UZYWAC jako filtru w trybie roslinnym. Funkcja istnieje po to, zeby
    raport mogl zawierac zdanie: "sprawdzono, motywy obecne/nieobecne,
    nieistotne w tym gospodarzu" -- co jest mocniejsze niz przemilczenie.
    """
    g = to_rna(guide)
    return [m for m in MAMMALIAN_IMMUNO_MOTIFS if m in g]


class MirnaGuard:
    """
    Sprawdza kolizje z endogennymi miRNA gospodarza.

    Wejscie: plik FASTA z dojrzalymi miRNA gospodarza.
    Zrodlo danych: miRBase (mirbase.org) -> mature.fa, przefiltrowane po
    prefiksie gatunku ('ath-' dla Arabidopsis, 'sly-' dla pomidora).

    Przyklad przygotowania pliku:
        grep -A1 '^>ath-' mature.fa | grep -v '^--' > ath_mature.fa
    """

    def __init__(self):
        self.mirnas: Dict[str, str] = {}
        self.seed_map: Dict[str, List[str]] = {}

    def load(self, fasta_path: str, species_prefix: Optional[str] = None
             ) -> 'MirnaGuard':
        name, buf = None, []

        def _flush():
            if name is None:
                return
            seq = to_rna(''.join(buf))
            if species_prefix and not name.startswith(species_prefix):
                return
            self.mirnas[name] = seq
            if len(seq) >= 8:
                self.seed_map.setdefault(seq[1:8], []).append(name)

        with open(fasta_path) as fh:
            for line in fh:
                line = line.rstrip()
                if line.startswith('>'):
                    _flush()
                    name, buf = line[1:].split()[0], []
                else:
                    buf.append(line)
        _flush()
        return self

    def check(self, guide: str) -> Dict[str, object]:
        """
        Sprawdza, czy nic prowadzaca dzieli seed z ktoryms miRNA gospodarza.
        """
        g = to_rna(guide)
        seed = g[1:8]
        kolizje = self.seed_map.get(seed, [])

        # dodatkowo: czy cala nic jest bardzo podobna do ktoregos miRNA
        podobne = []
        for nazwa, mir in self.mirnas.items():
            n = min(len(g), len(mir))
            ident = sum(1 for i in range(n) if g[i] == mir[i])
            if ident / n >= 0.85:
                podobne.append((nazwa, round(100 * ident / n, 1)))

        return {
            'seed': seed,
            'kolizje_seed_z_miRNA': kolizje,
            'n_kolizji': len(kolizje),
            'miRNA_podobne_85proc': podobne,
            'flaga': 'ODRZUC' if kolizje else ('SPRAWDZ' if podobne else 'OK'),
        }


def physicochemical_filters(sense_dna: str,
                            gc_min: float = 30.0,
                            gc_max: float = 52.0,
                            max_homopolymer: int = 4,
                            require_guide_5_AU: bool = True,
                            require_passenger_5_GC: bool = True,
                            min_asymmetry: float = 0.3
                            ) -> Tuple[bool, List[str]]:
    """
    Komplet filtrow fizykochemicznych. Zwraca (czy_przeszlo, lista_powodow).

    `sense_dna` to fragment mRNA (nic sensowna) -- z niego wyliczamy nic
    prowadzaca jako odwrotna komplementarnosc.
    """
    powody: List[str] = []
    sense_rna = to_rna(sense_dna)
    guide = revcomp_rna(sense_dna)

    g = gc_content(sense_rna)
    if not (gc_min <= g <= gc_max):
        powody.append(f'GC {g:.1f}% poza zakresem {gc_min}-{gc_max}%')

    if has_homopolymer(sense_rna, max_homopolymer):
        powody.append(f'ciag >={max_homopolymer} identycznych nt')

    if require_guide_5_AU and guide[0] not in 'AU':
        powody.append(f"5'-koniec nici prowadzacej = {guide[0]}, wymagane U/A")

    if require_passenger_5_GC and sense_rna[0] not in 'GC':
        powody.append(f"5'-koniec nici pasazerskiej = {sense_rna[0]}, "
                      f'wymagane G/C')

    _, _, asym = asymmetry(guide)
    if asym <= min_asymmetry:
        powody.append(f'asymetria {asym:+.2f} <= progu {min_asymmetry}')

    return (len(powody) == 0), powody


if __name__ == '__main__':
    print('--- Test filtrow na kandydacie A-21 ---')
    sense = 'CGGCATCAAGGTGAACTTCAA'
    ok, powody = physicochemical_filters(sense)
    print(f'sense    : {sense}')
    print(f'guide    : {revcomp_rna(sense)}')
    print(f'przeszlo : {ok}')
    if powody:
        for p in powody:
            print(f'   - {p}')

    print('\n--- Motywy ssacze (raportowane, NIE filtrowane) ---')
    m = check_mammalian_motifs(revcomp_rna(sense))
    print(f'znalezione: {m if m else "brak"}  '
          f'(nieistotne w roslinach -- brak TLR7/8)')
