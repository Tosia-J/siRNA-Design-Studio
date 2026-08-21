#!/usr/bin/env python3
"""
=============================================================================
main.py  --  ORKIESTRATOR PIPELINE'U
=============================================================================

Jeden punkt wejscia. Uruchamia wszystkie moduly po kolei i przekazuje dane
miedzy nimi automatycznie - nie ma potrzeby recznego kopiowania sekwencji.

PRZEPLYW

    plik FASTA z mRNA celu
            |
            v
    [design]      generowanie okien 21/22/24 nt
            |
            v
    [filters]     filtry fizykochemiczne wg profilu gospodarza
            |
            v
    [thermo]      asymetria, dG dupleksu, Tm
            |
            v
    [design]      MFE nici prowadzacej, dostepnosc celu (ViennaRNA)
            |
            v
    [offtarget]   indeks k-merowy transkryptomu, seed + seed-and-extend
            |     (pomijane, jesli nie podano transkryptomu)
            v
    [scoring]     normalizacja, ranking, analiza wrazliwosci, front Pareto
            |
            v
    [constructs]  shRNA i kaseta pod wybrany plazmid (opcjonalne)
            |
            v
    pliki wynikowe: TSV, raport tekstowy

UZYCIE

    # minimum - sam projekt, bez off-target
    python3 main.py --target gfp.fasta --host plant --out wyniki/

    # z analiza off-target
    python3 main.py --target gfp.fasta --host plant \\
                    --transcriptome TAIR10_cdna.fasta --out wyniki/

    # cel ssaczy, jedna dlugosc
    python3 main.py --target mrna.fasta --host mammal --lengths 21 --out wyniki/

    # z budowa konstruktow
    python3 main.py --target gfp.fasta --host plant \\
                    --promoter AtU6-1 --cloning goldengate_BsaI --out wyniki/

    # indeks transkryptomu buduje sie raz, potem wczytuje
    python3 main.py --build-index TAIR10_cdna.fasta --index-out tair10.idx
    python3 main.py --target gfp.fasta --index tair10.idx --out wyniki/

WYMAGANIA

    Wszystkie moduly (thermo.py, filters.py, hosts.py, design.py,
    offtarget.py, scoring.py, constructs.py) musza znajdowac sie
    W TYM SAMYM KATALOGU co main.py.

    pip install ViennaRNA

Autor: Antonina Jarecka
=============================================================================
"""

import argparse
import csv
import json
import os
import sys
from typing import Dict, List, Optional

# --- moduly wlasne ---
import hosts
import thermo
import filters
import design
import scoring

try:
    import offtarget
    OFFTARGET_DOSTEPNY = True
except ImportError:
    OFFTARGET_DOSTEPNY = False

try:
    import constructs
    CONSTRUCTS_DOSTEPNY = True
except ImportError:
    CONSTRUCTS_DOSTEPNY = False


# ============================================================================
# WEJSCIE
# ============================================================================

def wczytaj_fasta_pojedynczy(sciezka: str) -> tuple:
    """Wczytuje pierwszy rekord z pliku FASTA. Zwraca (naglowek, sekwencja)."""
    naglowek, buf = None, []
    with open(sciezka) as fh:
        for linia in fh:
            linia = linia.rstrip()
            if linia.startswith('>'):
                if naglowek is not None:
                    break
                naglowek = linia[1:].split()[0]
            elif naglowek is not None:
                buf.append(linia)
    if naglowek is None:
        raise ValueError(f'Brak rekordu FASTA w pliku {sciezka}')
    seq = ''.join(buf).upper().replace('U', 'T')
    niedozwolone = set(seq) - set('ACGTN')
    if niedozwolone:
        raise ValueError(f'Niedozwolone znaki w sekwencji: {sorted(niedozwolone)}')
    return naglowek, seq


# ============================================================================
# ETAPY
# ============================================================================

def etap_generowanie(mrna: str, profil: hosts.HostProfile,
                     dlugosci: Optional[tuple], verbose: bool) -> List[Dict]:
    if verbose:
        print('\n[1/5] GENEROWANIE I FILTRY')
        print('-' * 70)
        print(f'  profil gospodarza : {profil.nazwa} ({profil.opis})')
        print(f'  zakres GC         : {profil.gc_min}-{profil.gc_max}%')
        print(f'  5-koniec guide    : '
              f'{profil.guide_5_dozwolone or "brak wymogu"}')
        print(f'  motywy filtrowane : '
              f'{", ".join(profil.motywy_zabronione) or "brak"}')

    L = dlugosci if dlugosci else profil.dlugosci

    kandydaci = design.generuj_kandydatow(
        mrna,
        dlugosci=L,
        wyklucz_5prim_nt=profil.wyklucz_5prim_nt,
        wyklucz_3prim_nt=profil.wyklucz_3prim_nt,
        gc_min=profil.gc_min,
        gc_max=profil.gc_max,
        min_asymetria=profil.min_asymetria,
        max_mfe_guide=profil.max_mfe_guide,
        verbose=verbose,
    )

    # filtr motywow zaleznych od gospodarza
    if profil.motywy_zabronione:
        przed = len(kandydaci)
        kandydaci = [k for k in kandydaci
                     if not any(m in k['guide_rna']
                                for m in profil.motywy_zabronione)]
        if verbose:
            print(f'  filtr motywow: odrzucono {przed - len(kandydaci)} '
                  f'kandydatow')

    return kandydaci


def etap_offtarget(kandydaci: List[Dict], indeks, verbose: bool) -> List[Dict]:
    if verbose:
        print('\n[2/5] ANALIZA OFF-TARGET')
        print('-' * 70)

    if indeks is None:
        if verbose:
            print('  POMINIETO - nie podano transkryptomu.')
            print('  UWAGA: bez tej analizy ranking nie uwzglednia '
                  'bezpieczenstwa.')
        for k in kandydaci:
            k['koszt_offtarget'] = 20.0     # wartosc neutralna
            k['seed_czystosc'] = 0
            k['offtarget_flaga'] = 'NIESPRAWDZONE'
        return kandydaci

    slownik = {k['nazwa']: k['guide_rna'] for k in kandydaci}
    raport = offtarget.analyze_guides(slownik, indeks)

    for k in kandydaci:
        r = raport[k['nazwa']]
        k['koszt_offtarget'] = (r['min_koszt'] if r['min_koszt'] is not None
                                else 30.0)
        k['seed_czystosc'] = r['n_8mer']
        k['offtarget_flaga'] = r['flaga']
        k['offtarget_geny'] = ';'.join(r['geny_8mer'][:5])

    if verbose:
        odrzucone = sum(1 for k in kandydaci
                        if k['offtarget_flaga'].startswith('ODRZUC'))
        print(f'  sprawdzono {len(kandydaci)} kandydatow')
        print(f'  oznaczonych do odrzucenia: {odrzucone}')

    return kandydaci


def etap_scoring(kandydaci: List[Dict], verbose: bool) -> Dict:
    if verbose:
        print('\n[3/5] PUNKTACJA I RANKING')
        print('-' * 70)

    ranking = scoring.normalize_and_score(kandydaci)
    wrazliwosc = scoring.sensitivity_analysis(kandydaci, top_n=5)
    pareto = scoring.pareto_front(kandydaci)
    nazwy_pareto = {c['nazwa'] for c in pareto}

    for c in ranking:
        c['na_froncie_pareto'] = c['nazwa'] in nazwy_pareto

    if verbose:
        print(f'  kandydatow w rankingu     : {len(ranking)}')
        print(f'  na froncie Pareto         : {len(pareto)}')
        print(f'  stabilnosc rankingu       : {wrazliwosc["stabilnosc"]}')
        print(f'  {wrazliwosc["interpretacja"]}')

    return {'ranking': ranking, 'wrazliwosc': wrazliwosc,
            'pareto': nazwy_pareto}


def etap_rekomendacja(ranking: List[Dict], wrazliwosc: Dict,
                      n: int, verbose: bool) -> List[Dict]:
    """
    Laczy ranking wazony z frontem Pareto.

    Kandydat rekomendowany musi spelniac oba kryteria:
      - wysoka pozycja w rankingu wazonym
      - obecnosc na froncie Pareto (nie jest zdominowany przez zadnego innego)

    Kandydat wysoko oceniony, ale poza frontem, oznacza sytuacje, w ktorej
    istnieje inny kandydat lepszy we WSZYSTKICH kryteriach - wysoka pozycja
    wynika wtedy z przypadkowego doboru wag, nie z jakosci sekwencji.
    """
    if verbose:
        print('\n[4/5] REKOMENDACJA')
        print('-' * 70)

    rdzen = set(wrazliwosc['stabilny_rdzen'])
    rekomendowane = []
    for c in ranking:
        if not c['na_froncie_pareto']:
            continue
        if c.get('offtarget_flaga', '').startswith('ODRZUC'):
            continue
        c['stabilny_przy_zmianie_wag'] = c['nazwa'] in rdzen
        rekomendowane.append(c)
        if len(rekomendowane) >= n:
            break

    if verbose:
        for c in rekomendowane:
            gwiazdka = ' *' if c.get('stabilny_przy_zmianie_wag') else '  '
            print(f'  {c["ranga"]:>3}.{gwiazdka} {c["nazwa"]:<18} '
                  f'wynik {c["wynik"]:.4f}  '
                  f'{c["dlugosc"]} nt  poz {c["poz_start_1based"]}')
        print('  * = pozycja stabilna przy zmianie wag o +/-25%')

    return rekomendowane


def etap_konstrukty(rekomendowane: List[Dict], promotor: Optional[str],
                    klonowanie: Optional[str], verbose: bool) -> List:
    if verbose:
        print('\n[5/5] KONSTRUKTY')
        print('-' * 70)

    if not (promotor and klonowanie and CONSTRUCTS_DOSTEPNY):
        if verbose:
            print('  POMINIETO - nie podano promotora i metody klonowania.')
        return []

    wynik = []
    for c in rekomendowane:
        k = constructs.zbuduj_kasete(
            nazwa=c['nazwa'], guide_rna=c['guide_rna'],
            promotor_klucz=promotor, metoda_klonowania_klucz=klonowanie)
        wynik.append(k)
        if verbose:
            print(f'  {k.nazwa}: shRNA {len(k.shrna_dna)} nt, '
                  f'kaseta {len(k.kaseta_dna)} nt')
            for o in k.ostrzezenia:
                print(f'      [!] {o}')
    return wynik


# ============================================================================
# WYJSCIE
# ============================================================================

def zapisz_tsv(ranking: List[Dict], sciezka: str) -> None:
    if not ranking:
        return
    kolumny = ['ranga', 'nazwa', 'wynik', 'na_froncie_pareto', 'dlugosc',
               'klasa_DCL', 'poz_start_1based', 'poz_end_1based',
               'guide_rna', 'passenger_rna', 'seed_2_8', 'gc_proc',
               'asymetria', 'dG_dupleksu', 'Tm', 'mfe_guide',
               'dostepnosc_celu', 'koszt_offtarget', 'seed_czystosc',
               'offtarget_flaga', 'struktura_guide']
    with open(sciezka, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=kolumny, delimiter='\t',
                           extrasaction='ignore')
        w.writeheader()
        for c in ranking:
            w.writerow(c)


def zapisz_raport(sciezka: str, naglowek_celu: str, mrna: str,
                  profil: hosts.HostProfile, wynik: Dict,
                  rekomendowane: List[Dict], transkryptom: Optional[str]) -> None:
    L = []
    L.append('RAPORT PROJEKTOWANIA siRNA')
    L.append('=' * 78)
    L.append(f'Cel              : {naglowek_celu} ({len(mrna)} nt)')
    L.append(f'Profil gospodarza: {profil.nazwa} - {profil.opis}')
    L.append(f'Transkryptom     : {transkryptom or "NIE PODANO"}')
    L.append('')
    L.append('PARAMETRY PROFILU')
    L.append('-' * 78)
    L.append(f'  GC                    : {profil.gc_min}-{profil.gc_max}%')
    L.append(f'  5-koniec guide        : {profil.guide_5_dozwolone or "brak wymogu"}')
    L.append(f'  5-koniec passenger    : {profil.passenger_5_dozwolone or "brak wymogu"}')
    L.append(f'  min. asymetria        : {profil.min_asymetria} kcal/mol')
    L.append(f'  motywy filtrowane     : {", ".join(profil.motywy_zabronione) or "brak"}')
    L.append(f'  motywy raportowane    : {", ".join(profil.motywy_raportowane) or "brak"}')
    L.append('')
    for u in profil.uwagi:
        L.append(f'  UWAGA: {u}')
    L.append('')
    L.append(scoring.raport_tekstowy(wynik['ranking'], top_n=15))
    L.append('')
    L.append('ANALIZA WRAZLIWOSCI NA WAGI')
    L.append('-' * 78)
    w = wynik['wrazliwosc']
    L.append(f'  Top-5 bazowy       : {", ".join(w["top_bazowy"])}')
    L.append(f'  Stabilny rdzen     : {", ".join(w["stabilny_rdzen"]) or "pusty"}')
    L.append(f'  Stabilnosc         : {w["stabilnosc"]}')
    L.append(f'  Interpretacja      : {w["interpretacja"]}')
    L.append('')
    L.append('  Znaczenie: jesli pierwsza piatka zmienia sie przy zmianie wag')
    L.append('  o +/-25%, oznacza to, ze ranking odzwierciedla arbitralny dobor')
    L.append('  wag, a nie wlasciwosci sekwencji. Stabilnosc >= 0.8 pozwala')
    L.append('  zaraportowac ranking jako niezalezny od tego doboru.')
    L.append('')
    L.append('KANDYDACI REKOMENDOWANI')
    L.append('-' * 78)
    L.append('  Kryteria: obecnosc na froncie Pareto ORAZ brak flagi odrzucenia')
    L.append('  w analizie off-target.')
    L.append('')
    for c in rekomendowane:
        L.append(f'  {c["nazwa"]}  ({c["dlugosc"]} nt, poz. '
                 f'{c["poz_start_1based"]}-{c["poz_end_1based"]})')
        L.append(f'      guide     : 5-{c["guide_rna"]}-3')
        L.append(f'      passenger : 5-{c["passenger_rna"]}-3')
        L.append(f'      GC {c["gc_proc"]}%, asymetria {c["asymetria"]:+.2f}, '
                 f'Tm {c["Tm"]} C')
        L.append(f'      off-target: {c.get("offtarget_flaga", "n/d")}')
        L.append('')
    with open(sciezka, 'w') as fh:
        fh.write('\n'.join(L))


# ============================================================================
# CLI
# ============================================================================

def main() -> int:
    p = argparse.ArgumentParser(
        description='Pipeline projektowania siRNA',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__)

    p.add_argument('--target', help='plik FASTA z mRNA celu')
    p.add_argument('--host', default='generic',
                   choices=sorted(hosts.PROFILE),
                   help='profil gospodarza (domyslnie: generic)')
    p.add_argument('--lengths', type=int, nargs='+',
                   help='dlugosci siRNA; domyslnie wg profilu')
    p.add_argument('--transcriptome',
                   help='FASTA transkryptomu gospodarza do analizy off-target')
    p.add_argument('--index', help='gotowy indeks transkryptomu (.idx)')
    p.add_argument('--build-index',
                   help='zbuduj indeks z podanego FASTA i zakoncz')
    p.add_argument('--index-out', default='transkryptom.idx',
                   help='sciezka zapisu indeksu')
    p.add_argument('--promoter', help='promotor do budowy kasety')
    p.add_argument('--cloning', help='metoda klonowania')
    p.add_argument('--top', type=int, default=5,
                   help='liczba kandydatow rekomendowanych')
    p.add_argument('--out', default='wyniki',
                   help='katalog wynikowy')
    p.add_argument('--quiet', action='store_true')
    p.add_argument('--list-profiles', action='store_true')

    a = p.parse_args()
    verbose = not a.quiet

    if a.list_profiles:
        print(hosts.lista_profili())
        return 0

    # budowa indeksu i wyjscie
    if a.build_index:
        if not OFFTARGET_DOSTEPNY:
            print('BLAD: modul offtarget niedostepny.', file=sys.stderr)
            return 1
        print(f'Buduje indeks z {a.build_index} ...')
        idx = offtarget.TranscriptomeIndex(k=8).build(a.build_index,
                                                      verbose=verbose)
        idx.save(a.index_out)
        print(f'Zapisano: {a.index_out}')
        return 0

    if not a.target:
        p.error('wymagane --target (albo --build-index / --list-profiles)')

    os.makedirs(a.out, exist_ok=True)
    profil = hosts.get_profile(a.host)
    naglowek, mrna = wczytaj_fasta_pojedynczy(a.target)

    if verbose:
        print('=' * 70)
        print('PIPELINE PROJEKTOWANIA siRNA')
        print('=' * 70)
        print(f'  cel: {naglowek}, {len(mrna)} nt')

    # indeks off-target
    indeks = None
    if a.index:
        if not OFFTARGET_DOSTEPNY:
            print('BLAD: modul offtarget niedostepny.', file=sys.stderr)
            return 1
        indeks = offtarget.TranscriptomeIndex.load(a.index)
    elif a.transcriptome:
        if not OFFTARGET_DOSTEPNY:
            print('BLAD: modul offtarget niedostepny.', file=sys.stderr)
            return 1
        indeks = offtarget.TranscriptomeIndex(k=8).build(a.transcriptome,
                                                         verbose=verbose)

    # etapy
    kandydaci = etap_generowanie(mrna, profil,
                                 tuple(a.lengths) if a.lengths else None,
                                 verbose)
    if not kandydaci:
        print('\nBRAK KANDYDATOW po filtrach. Rozwaz poluzowanie progow '
              '(zakres GC, asymetria) albo wykluczen koncow ORF.',
              file=sys.stderr)
        return 2

    kandydaci = etap_offtarget(kandydaci, indeks, verbose)
    wynik = etap_scoring(kandydaci, verbose)
    rekomendowane = etap_rekomendacja(wynik['ranking'], wynik['wrazliwosc'],
                                      a.top, verbose)
    etap_konstrukty(rekomendowane, a.promoter, a.cloning, verbose)

    # zapis
    tsv = os.path.join(a.out, 'ranking.tsv')
    raport = os.path.join(a.out, 'raport.txt')
    js = os.path.join(a.out, 'wyniki.json')
    zapisz_tsv(wynik['ranking'], tsv)
    zapisz_raport(raport, naglowek, mrna, profil, wynik, rekomendowane,
                  a.transcriptome or a.index)
    with open(js, 'w') as fh:
        json.dump({'cel': naglowek, 'profil': profil.nazwa,
                   'ranking': wynik['ranking'],
                   'wrazliwosc': wynik['wrazliwosc']}, fh, indent=1)

    if verbose:
        print(f'\nZapisano:\n  {tsv}\n  {raport}\n  {js}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
