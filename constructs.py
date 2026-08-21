"""
=============================================================================
constructs.py  --  SKLADANIE shRNA I KASET EKSPRESYJNYCH
=============================================================================

STATUS: WLASNA IMPLEMENTACJA

=============================================================================
ODPOWIEDZ NA PYTANIE: CZY ROZNE PLAZMIDY WYMAGAJA INNYCH SEKWENCJI
NA KONCACH KASETY
=============================================================================

Tak, i to na dwoch niezaleznych poziomach. Trzeba je rozroznic, bo mieszanie
ich jest zrodlem wiekszosci nieudanych klonowan.

POZIOM 1 -- WYMAGANIA POLIMERAZY (biologia)
--------------------------------------------
Zalezy od PROMOTORA, nie od plazmidu. Decyduje o tym, czy transkrypt w ogole
powstanie i jak bedzie wygladal.

  Pol III (U6, U3, 7SL)
    - transkrypt MUSI zaczynac sie od okreslonego nukleotydu:
        U6 -> G
        U3 -> A
      Jesli Twoja nic prowadzaca zaczyna sie inaczej, trzeba DODAC ten
      nukleotyd na 5'. Powstaje wtedy jeden nukleotyd nadmiarowy, ktory
      DCL zwykle toleruje, ale ktory nalezy zaraportowac w metodyce.
    - terminacja: ciag T (zwykle 6, minimum 4-5) bezposrednio za insertem
    - transkrypt NIE dostaje czapeczki 5' ani ogona poliA
    - limit dlugosci praktyczny: do okolo 250 nt
    - ZALETA dla shRNA: transkrypt ma zdefiniowane konce, wiec spinka
      powstaje czysta

  Pol II (35S, UBQ10, nos)
    - brak wymogu co do pierwszego nukleotydu
    - terminator wymagany jako osobny element (NOS-ter, 35S-ter, HSP18.2)
    - transkrypt dostaje czapeczke i ogon poliA
    - brak limitu dlugosci -- mozna dawac dlugie hairpiny z intronem
    - WADA dla krotkich shRNA: czapeczka i ogon poliA moga utrudniac
      dostep DCL do konca dupleksu

  WNIOSEK PRAKTYCZNY dla Twojego projektu: dla krotkich shRNA (54-60 nt)
  Pol III jest naturalniejszy. Ale wtedy konstrukty A-21, A-22, A-24, B-21,
  B-22, B-24 i SCRAMBLED WSZYSTKIE wymagaja dodania G na 5', bo zaczynaja
  sie od U albo A.

POZIOM 2 -- WYMAGANIA METODY KLONOWANIA (technika)
---------------------------------------------------
Zalezy od tego, JAK wkladasz kasete do plazmidu. Nie ma zwiazku z biologia.

  KLONOWANIE RESTRYKCYJNE
    - potrzebne miejsca restrykcyjne na obu koncach insertu
    - miejsca MUSZA byc unikatowe w calym plazmidzie, nie tylko w insercie
    - dodaje sie 4-6 nt "wypelniacza" przed miejscem, bo enzymy tna slabo
      przy samym koncu fragmentu DNA
    - te dodatkowe nukleotydy ZOSTAJA w gotowym konstrukcie

  GOLDEN GATE (BsaI, BpiI/BbsI, Esp3I)
    - potrzebne miejsce rozpoznania enzymu + 4-nt miejsce fuzji
    - miejsce fuzji definiuje kolejnosc skladania -- musi byc zgodne ze
      standardem (np. MoClo), inaczej elementy zloza sie w zlej kolejnosci
    - KRYTYCZNE: insert NIE MOZE zawierac wewnetrznych miejsc uzywanego
      enzymu. Jesli zawiera, trzeba je usunac cichymi mutacjami -- co przy
      shRNA jest niemozliwe, bo sekwencja jest zdeterminowana przez cel.
      Wtedy zmieniamy enzym.
    - w gotowym konstrukcie zostaje TYLKO 4-nt miejsce fuzji (szew)

  GATEWAY (attB/attP/attL/attR)
    - potrzebne miejsca attB1 i attB2 flankujace insert
    - w gotowym konstrukcie zostaje attB (okolo 25 nt po kazdej stronie)
    - te 25 nt to sporo przy 54-nt shRNA -- warto to policzyc
    - system wygodny przy pHellsgate/pHANNIBAL (klasyka hairpin RNAi
      w roslinach)

  GIBSON / NEBuilder
    - potrzebne ramiona homologii 20-40 nt zgodne z miejscem wstawienia
    - w gotowym konstrukcie NIE zostaje nic dodatkowego (szew bezszwowy)
    - najczystszy wynik, ale wymaga znajomosci dokladnej sekwencji wektora

=============================================================================
UWAGA O SEKWENCJACH WEKTOROW
=============================================================================

W tym module NIE MA wpisanych sekwencji konkretnych wektorow. To celowe.

Sekwencje pCAMBIA roznia sie miedzy wariantami i miedzy laboratoriami
(wersje pochodne, modyfikacje MCS). Wpisanie tu sekwencji "z pamieci"
oznaczaloby zaprojektowanie konstruktu, ktory nie pasuje do plazmidu
lezacego w Twojej zamrazarce. Zamiast tego rejestr koduje WYMAGANIA, a
konkretne sekwencje wklejasz z mapy swojego wektora.

To samo dotyczy motywu TLS2. Sekwencja pochodzi z konkretnej publikacji i
musi byc przepisana ze zrodla, nie odtworzona z pamieci.

Autor: Antonina Jarecka
=============================================================================
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple
from thermo import to_rna, revcomp_rna

COMP = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'U': 'A'}


def revcomp_dna(s: str) -> str:
    return ''.join(COMP[c] for c in reversed(s.upper().replace('U', 'T')))


def to_dna(s: str) -> str:
    return s.upper().replace('U', 'T')


# ============================================================================
# PETLE
# ============================================================================

PETLE = {
    'standard': {
        'seq': 'GTTCAAGAGA',
        'opis': 'klasyczna petla shRNA, 10 nt, stabilna, szeroko stosowana',
    },
    'standard_12': {
        'seq': 'GTTCAAGAGAAC',
        'opis': '12-nt wariant, uzyty w konstruktach anty-GFP',
    },
    'mir_like': {
        'seq': 'CTGCTGCTGCT',
        'opis': 'petla o niskiej strukturze wlasnej, stosowana w wektorach '
                'roslinnych',
    },
}


# ============================================================================
# REJESTR PROMOTOROW
# ============================================================================

@dataclass
class Promotor:
    nazwa: str
    polimeraza: str                 # 'PolII' albo 'PolIII'
    wymagany_pierwszy_nt: Optional[str]   # None jesli brak wymogu
    terminator: str                 # opis wymaganego terminatora
    max_dlugosc_nt: Optional[int]
    uwagi: str


PROMOTORY: Dict[str, Promotor] = {
    'AtU6-1': Promotor(
        nazwa='AtU6-1', polimeraza='PolIII', wymagany_pierwszy_nt='G',
        terminator='TTTTTT', max_dlugosc_nt=250,
        uwagi='Transkrypt musi zaczynac sie od G. Konstytutywny. '
              'Naturalny wybor dla krotkich shRNA.'),
    'AtU3': Promotor(
        nazwa='AtU3', polimeraza='PolIII', wymagany_pierwszy_nt='A',
        terminator='TTTTTT', max_dlugosc_nt=250,
        uwagi='Transkrypt musi zaczynac sie od A. Alternatywa dla U6, '
              'gdy nic prowadzaca zaczyna sie od A -- wtedy zaden dodatkowy '
              'nukleotyd nie jest potrzebny.'),
    'CaMV35S': Promotor(
        nazwa='CaMV35S', polimeraza='PolII', wymagany_pierwszy_nt=None,
        terminator='NOS-ter albo 35S-ter (element osobny)', max_dlugosc_nt=None,
        uwagi='Silny, konstytutywny. Transkrypt dostaje czapeczke i poliA, '
              'co moze utrudniac dostep DCL do krotkiej spinki.'),
    'AtUBQ10': Promotor(
        nazwa='AtUBQ10', polimeraza='PolII', wymagany_pierwszy_nt=None,
        terminator='NOS-ter', max_dlugosc_nt=None,
        uwagi='Konstytutywny, bardziej rownomierny niz 35S w roznych tkankach.'),
}


# ============================================================================
# REJESTR METOD KLONOWANIA
# ============================================================================

@dataclass
class MetodaKlonowania:
    nazwa: str
    flanka_5: str                   # co dodajemy przed insertem
    flanka_3: str                   # co dodajemy za insertem
    zabronione_motywy: List[str]    # nie moga wystapic wewnatrz insertu
    zostaje_w_konstrukcie: str
    uwagi: str


KLONOWANIE: Dict[str, MetodaKlonowania] = {
    'restrykcyjne_BamHI_EcoRI': MetodaKlonowania(
        nazwa='restrykcyjne BamHI/EcoRI',
        flanka_5='AAAAGGATCC',       # wypelniacz + BamHI
        flanka_3='GAATTCAAAA',       # EcoRI + wypelniacz
        zabronione_motywy=['GGATCC', 'GAATTC'],
        zostaje_w_konstrukcie='pelne miejsca restrykcyjne + wypelniacze',
        uwagi='Sprawdz, czy oba miejsca sa UNIKATOWE w calym plazmidzie, '
              'a nie tylko nieobecne w insercie.'),
    'goldengate_BsaI': MetodaKlonowania(
        nazwa='Golden Gate (BsaI)',
        flanka_5='GGTCTCN' + 'AATG',   # BsaI + miejsce fuzji (przyklad MoClo)
        flanka_3='GCTT' + 'NGAGACC',
        zabronione_motywy=['GGTCTC', 'GAGACC'],
        zostaje_w_konstrukcie='tylko 4-nt miejsce fuzji (szew)',
        uwagi='Miejsca fuzji MUSZA byc zgodne ze standardem MoClo dla Twojego '
              'zestawu. Podane tu sa przykladowe -- sprawdz w dokumentacji '
              'zestawu.'),
    'gateway_attB': MetodaKlonowania(
        nazwa='Gateway (attB1/attB2)',
        flanka_5='[attB1 -- wklej z dokumentacji zestawu, ok. 25 nt]',
        flanka_3='[attB2 -- wklej z dokumentacji zestawu, ok. 25 nt]',
        zabronione_motywy=[],
        zostaje_w_konstrukcie='attB1 i attB2, po ok. 25 nt z kazdej strony',
        uwagi='Przy 54-nt shRNA flanki attB sa dluzsze niz sam insert. '
              'Policz, czy to nie zaburzy faldowania spinki.'),
    'gibson': MetodaKlonowania(
        nazwa='Gibson / NEBuilder',
        flanka_5='[20-40 nt homologii do wektora -- wklej z mapy]',
        flanka_3='[20-40 nt homologii do wektora -- wklej z mapy]',
        zabronione_motywy=[],
        zostaje_w_konstrukcie='nic -- szew bezszwowy',
        uwagi='Najczystszy wynik. Wymaga dokladnej sekwencji miejsca '
              'wstawienia w wektorze.'),
}


# ============================================================================
# MOTYWY DODATKOWE
# ============================================================================

MOTYWY_DODATKOWE = {
    'TLS2': {
        'seq': None,   # CELOWO PUSTE
        'opis': 'Skrocony motyw tRNA-podobny (bez petli D i T). Nadaje '
                'transkryptowi kompetencje do transportu dalekodystansowego. '
                'Wersja skrocona dziala tak samo jak pelna TLS1, wiec '
                'preferujemy ja jako krotsza.',
        'status': 'SEKWENCJA DO UZUPELNIENIA ZE ZRODLA -- nie odtwarzac '
                  'z pamieci. Wklej z publikacji zrodlowej.',
        'ostrzezenie': 'Prace nad TLS dotycza mRNA, nie shRNA. Nie wiadomo, '
                       'czy krotka spinka korzysta z tego samego mechanizmu. '
                       'To jest hipoteza do przetestowania, nie fakt.',
    },
}


# ============================================================================
# SKLADANIE KONSTRUKTU
# ============================================================================

@dataclass
class Konstrukt:
    nazwa: str
    guide: str
    passenger: str
    petla: str
    shrna_dna: str
    kaseta_dna: str
    promotor: str
    metoda_klonowania: str
    dodany_nt_5: Optional[str] = None
    ostrzezenia: List[str] = field(default_factory=list)


def zbuduj_shrna(guide_rna: str,
                 petla_klucz: str = 'standard_12',
                 kolejnosc: str = 'guide_first') -> str:
    """
    Sklada spinke shRNA w DNA.

    kolejnosc:
      'guide_first'     -- nic prowadzaca na 5' (zalecane: lepsze zaladowanie
                           do AGO, bo koniec 5' spinki jest przetwarzany
                           najpierw)
      'passenger_first' -- nic pasazerska na 5'
    """
    g = to_dna(guide_rna)
    p = revcomp_dna(g)
    loop = PETLE[petla_klucz]['seq']
    if kolejnosc == 'guide_first':
        return g + loop + p
    return p + loop + g


def zbuduj_kasete(nazwa: str,
                  guide_rna: str,
                  promotor_klucz: str,
                  metoda_klonowania_klucz: str,
                  petla_klucz: str = 'standard_12',
                  motyw_dodatkowy: Optional[str] = None) -> Konstrukt:
    """
    Sklada kompletna kasete gotowa do zamowienia syntezy.

    Automatycznie:
      - dodaje nukleotyd wymagany przez polimeraze, jesli trzeba
      - dodaje terminator dla Pol III
      - dodaje flanki metody klonowania
      - sprawdza kolizje z zabronionymi motywami
      - zbiera ostrzezenia
    """
    prom = PROMOTORY[promotor_klucz]
    klon = KLONOWANIE[metoda_klonowania_klucz]
    ostrzezenia: List[str] = []

    shrna = zbuduj_shrna(guide_rna, petla_klucz)
    dodany_nt = None

    # --- wymog pierwszego nukleotydu (Pol III) ---
    if prom.wymagany_pierwszy_nt:
        if shrna[0] != prom.wymagany_pierwszy_nt:
            dodany_nt = prom.wymagany_pierwszy_nt
            shrna = dodany_nt + shrna
            ostrzezenia.append(
                f'Dodano {dodany_nt} na 5-koncu, bo {prom.nazwa} '
                f'({prom.polimeraza}) wymaga transkryptu zaczynajacego sie od '
                f'{prom.wymagany_pierwszy_nt}. Ten nukleotyd bedzie obecny '
                f'w transkrypcie -- zaraportuj to w metodyce.')

    # --- motyw dodatkowy (np. TLS2) ---
    if motyw_dodatkowy:
        m = MOTYWY_DODATKOWE.get(motyw_dodatkowy)
        if m is None:
            raise ValueError(f'Nieznany motyw: {motyw_dodatkowy}')
        if m['seq'] is None:
            ostrzezenia.append(
                f'Motyw {motyw_dodatkowy} zazadany, ale jego sekwencja nie '
                f'jest wpisana w rejestrze. {m["status"]}')
        else:
            shrna = shrna + m['seq']
        if 'ostrzezenie' in m:
            ostrzezenia.append(f'{motyw_dodatkowy}: {m["ostrzezenie"]}')

    # --- terminator ---
    rdzen = shrna
    if prom.polimeraza == 'PolIII':
        rdzen = rdzen + prom.terminator
    else:
        ostrzezenia.append(
            f'{prom.nazwa} to {prom.polimeraza} -- terminator '
            f'({prom.terminator}) musi byc osobnym elementem wektora, '
            f'nie czescia insertu.')

    # --- limit dlugosci ---
    if prom.max_dlugosc_nt and len(rdzen) > prom.max_dlugosc_nt:
        ostrzezenia.append(
            f'Dlugosc transkryptu {len(rdzen)} nt przekracza praktyczny limit '
            f'{prom.max_dlugosc_nt} nt dla {prom.polimeraza}.')

    # --- kolizje z motywami zabronionymi ---
    for motyw in klon.zabronione_motywy:
        if motyw in rdzen or motyw in revcomp_dna(rdzen):
            ostrzezenia.append(
                f'KOLIZJA: insert zawiera {motyw}, wymagany przez metode '
                f'{klon.nazwa}. Zmien enzym albo metode -- sekwencji shRNA '
                f'nie da sie zmienic bez utraty celu.')

    kaseta = klon.flanka_5 + rdzen + klon.flanka_3

    if '[' in kaseta:
        ostrzezenia.append(
            'Kaseta zawiera placeholder w nawiasach kwadratowych. '
            'Uzupelnij sekwencja z mapy swojego wektora przed zamowieniem.')

    return Konstrukt(
        nazwa=nazwa, guide=to_rna(guide_rna),
        passenger=to_rna(revcomp_dna(guide_rna)),
        petla=PETLE[petla_klucz]['seq'],
        shrna_dna=rdzen, kaseta_dna=kaseta,
        promotor=prom.nazwa, metoda_klonowania=klon.nazwa,
        dodany_nt_5=dodany_nt, ostrzezenia=ostrzezenia)


if __name__ == '__main__':
    print('=' * 78)
    print('TEST: ten sam guide w trzech roznych systemach')
    print('=' * 78)

    guide = 'UUGAAGUUCACCUUGAUGCCG'   # A-21 anty-GFP

    for prom, klon in (('AtU6-1', 'restrykcyjne_BamHI_EcoRI'),
                       ('AtU6-1', 'goldengate_BsaI'),
                       ('CaMV35S', 'restrykcyjne_BamHI_EcoRI')):
        k = zbuduj_kasete(f'A-21_{prom}', guide, prom, klon)
        print(f'\n--- {prom} + {klon} ---')
        print(f'  dodany nt na 5 : {k.dodany_nt_5 or "brak"}')
        print(f'  rdzen ({len(k.shrna_dna)} nt): {k.shrna_dna}')
        print(f'  kaseta ({len(k.kaseta_dna)} nt): {k.kaseta_dna}')
        for o in k.ostrzezenia:
            print(f'  [!] {o}')
