# Varmesimulering

Interaktiv simulering av temperaturen i en isolert metallstav på 1 meter.

## Start

Krever Python 3.10 eller nyere.

```bash
python -m venv .venv
.venv\Scripts\python.exe -m pip install -r requirements.txt
.venv\Scripts\python.exe simulering.py
```

På Windows kan du etterpå dobbeltklikke `start.bat`. På macOS/Linux er stien `.venv/bin/python` i stedet for `.venv\Scripts\python.exe`.

## Rediger temperaturen

Gjøres direkte på grafen eller på staven:

| Handling | Effekt |
|---|---|
| Venstreklikk | Legger til en varmekilde (temperatur og bredde fra feltene *Kilde* og *Bredde*) |
| Venstre + dra | Tegner temperaturkurven fritt |
| Høyreklikk / høyre + dra | Fjerner varme (tilbake til 0 °C) |
| *Nullstill* | Fjerner alt |

Redigerer du midt i en simulering, blir tilstanden der den nye starten (t = 0).

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
| *Bredde [cm]* | Bredden til nye varmekilder |
| *Fourier-ledd* | Antall ledd i Fourier-rekka (stiplet linje) |

Standardverdiene for materialet er for kobber. Trykk Enter etter å ha skrevet i et felt.
