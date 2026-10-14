# Varmesimulering

Interaktiv simulering av temperaturen i isolert metall, enten en stav på 1 meter (1D) eller en plate på 1 × 1 meter (2D).

## Start

Krever Python 3.10 eller nyere, fra [python.org](https://www.python.org).

**Windows:** dobbeltklikk `start.bat`. Første gang installeres pakkene simuleringen trenger, og det tar et par minutter. Etter det starter simuleringen med en gang.

**Manuelt** (også på macOS/Linux, der stien er `.venv/bin/python` i stedet for `.venv\Scripts\python.exe`):

```bash
python -m venv .venv
.venv\Scripts\python.exe -m pip install -r requirements.txt
.venv\Scripts\python.exe simulering.py
```

## Stav eller plate

Simuleringen starter med staven. Knappen *Bytt til 2D* / *Bytt til 1D* veksler mellom stav og plate. Materialet og innstillingene for varmekilder og rand følger med.

- **Stav (1D):** den heltrukne kurven er den numeriske løsningen, den stiplede er Fourier-rekka. Fargestripen under viser staven.
- **Plate (2D):** det venstre fargekartet er den numeriske løsningen, det høyre er Fourier-rekka.

## Rand

Knappen *Rand* velger hva som skjer i kantene:

- **Rand: isolert** – ingen varme slipper ut. Temperaturen jevner seg ut til snittet.
- **Rand: fast temp.** – kantene holdes på en fast temperatur, som du skriver inn i feltene som dukker opp. For staven kan *Venstre* og *Høyre* ende ha ulik temperatur; for platen har alle fire kantene temperaturen *Kant*.

Bytte av rand, eller ny randtemperatur, starter simuleringen på nytt fra samme starttemperatur.

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

## Les av temperaturen

Skriv inn en posisjon i *Les av ved x [cm]* (og *y [cm]* for platen). Temperaturen i det punktet ved gjeldende tid vises til høyre, både for den numeriske løsningen og for Fourier-rekka. Et merke i figuren viser hvor det leses av.

For å finne temperaturen i et bestemt punkt ved et bestemt tidspunkt: skriv inn tiden i *t [s]* og posisjonen i avlesningsfeltene.

Du kan også holde musepekeren over grafen, staven eller fargekartene. Da vises posisjonen og temperaturen akkurat der pekeren er, i en liten etikett ved siden av den.

## Kopier grafen

Tiden står i selve grafen (øverst til venstre for staven, i titlene for platen), så den kommer med på skjermbilder.

Knappen *Kopier grafen* legger grafen på utklippstavlen, slik at den kan limes rett inn i et dokument. For platen er det én knapp for hvert fargekart, *Kopier numerisk* og *Kopier Fourier*, der fargesøylen følger med, og *Kopier begge*, som kopierer de to kartene ved siden av hverandre slik de vises i vinduet. På andre systemer enn Windows lagres bildet i stedet som en PNG-fil.

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
