# Varmesimulering

Interaktiv simulering av temperaturen i isolert metall, enten en stav på 1 meter (1D) eller en plate på 1 × 1 meter (2D).

## Start

Krever Python 3.10 eller nyere.

```bash
python -m venv .venv
.venv\Scripts\python.exe -m pip install -r requirements.txt
.venv\Scripts\python.exe simulering.py
```

På Windows kan du etterpå dobbeltklikke `start.bat`. På macOS/Linux er stien `.venv/bin/python` i stedet for `.venv\Scripts\python.exe`.

## Stav eller plate

Simuleringen starter med staven. Knappen *Bytt til 2D* / *Bytt til 1D* veksler mellom stav og plate. Materialet og innstillingene for varmekilder følger med.

- **Stav (1D):** den heltrukne kurven er den numeriske løsningen, den stiplede er Fourier-rekka. Fargestripen under viser staven.
- **Plate (2D):** det venstre fargekartet er den numeriske løsningen, det høyre er Fourier-rekka.

## Rediger temperaturen

Gjøres direkte på grafen, staven eller fargekartene:

| Handling | Stav | Plate |
|---|---|---|
| Venstreklikk | Legger til en varmekilde | Legger til en kvadratisk varmekilde |
| Venstre + dra | Tegner temperaturkurven fritt | Maler varme |
| Høyreklikk / høyre + dra | Fjerner varme (tilbake til 0 °C) | Fjerner varme (tilbake til 0 °C) |
| *Nullstill* | Fjerner alt | Fjerner alt |

Varmekildene får temperatur og bredde fra feltene *Kilde* og *Bredde*. Redigerer du midt i en simulering, blir tilstanden der den nye starten (t = 0).

## Tid

| Kontroll | Effekt |
|---|---|
| *Spill av* / mellomrom | Starter og stopper avspillingen |
| Tidsbryteren | Flytter tiden; starten av simuleringen har fått mest plass |
| ← → / musehjul | Små steg i tid (Shift gir 10 ganger større steg) |
| *t [s]* | Skriv inn et eksakt tidspunkt |
| *Slutt [s]* | Hvor langt tidsbryteren går; mindre verdi gir finere kontroll |
| *Fart* | Avspillingsfart fra ×0,05 til ×8 |

## Andre innstillinger

| Felt | Effekt |
|---|---|
| *k* | Varmeledningsevne [W/(m·K)] |
| *ρ* | Tetthet [kg/m³] |
| *c<sub>p</sub>* | Spesifikk varmekapasitet [J/(kg·K)] |
| *Kilde [°C]* | Temperaturen til nye varmekilder |
| *Bredde [cm]* | Bredden til nye varmekilder (sidekanten på kvadratet for platen) |
| *Fourier-ledd* | Antall ledd i Fourier-rekka (for platen: antall i hver retning) |
| *Farger* (bare plate) | *Lineære*: fargen følger temperaturen jevnt. *Forsterket*: lave temperaturer er strukket ut, så svak varme synes bedre |

Standardverdiene for materialet er for kobber. Trykk Enter etter å ha skrevet i et felt.
