# Varmesimulering

Interaktiv simulering av varmelikningen i en isolert metallstav. Den numeriske løsningen (endelige differanser) vises sammen med den analytiske løsningen som Fourier-cosinusrekke, slik at de to kan sammenlignes direkte.

## Modellen

Temperaturen $u(x,t)$ i en stav med lengde $L$ følger varmelikningen

$$\frac{\partial u}{\partial t} = \alpha \frac{\partial^2 u}{\partial x^2}, \qquad \alpha = \frac{k}{\rho\, c_p}$$

med isolerte ender (homogene Neumann-betingelser): $\dfrac{\partial u}{\partial x} = 0$ i $x = 0$ og $x = L$.

Løsningen kan skrives som en Fourier-cosinusrekke:

$$u(x,t) = a_0 + \sum_{n=1}^{\infty} a_n \cos\left(\frac{n\pi x}{L}\right) e^{-\alpha (n\pi/L)^2 t}$$

## Kom i gang

Krever Python 3.10 eller nyere.

```bash
python -m venv .venv
.venv\Scripts\python.exe -m pip install -r requirements.txt
.venv\Scripts\python.exe simulering.py
```

På Windows kan du etterpå dobbeltklikke `start.bat`. På macOS/Linux er stien `.venv/bin/python` i stedet for `.venv\Scripts\python.exe`.

## Bruk

**Rediger temperaturen** direkte på grafen eller på staven:

| Handling | Effekt |
|---|---|
| Venstreklikk | Legger til en varmekilde (temperatur og bredde fra feltene *Kilde* og *Bredde*) |
| Venstre + dra | Tegner temperaturkurven fritt |
| Høyreklikk / høyre + dra | Fjerner varme (tilbake til 20 °C) |
| *Nullstill* | Fjerner alt |

Redigerer du midt i en simulering, blir tilstanden der den nye starten ($t = 0$).

**Tid:**

| Kontroll | Effekt |
|---|---|
| *Spill av* / mellomrom | Starter og stopper avspillingen |
| Tidsbryteren | Skalert i $\sqrt{t}$, så starten får mest plass |
| ← → / musehjul | Små steg i tid (Shift gir 10 ganger større steg) |
| *t [s]* | Skriv inn et eksakt tidspunkt |
| *Slutt [s]* | Hvor langt tidsbryteren går; mindre verdi gir finere kontroll |
| *Fart* | Avspillingsfart fra ×0,05 til ×8 |

**Materiale:** Skriv inn varmeledningsevne $k$ [W/(m·K)], tetthet $\rho$ [kg/m³] og spesifikk varmekapasitet $c_p$ [J/(kg·K)]. Standardverdiene er for kobber.

**Fourier-ledd:** Antall ledd i rekka (stiplet linje). Få ledd ved $t = 0$ viser Gibbs-fenomenet ved skarpe kanter.

## Oppbygging

| Fil | Innhold |
|---|---|
| `simulering.py` | Åpner simuleringsvinduet |
| `varme/materialer.py` | Materialdata for noen metaller |
| `varme/rutenett.py` | Rutenett i 1, 2 eller 3 dimensjoner |
| `varme/startbetingelser.py` | Starttemperatur med varmt segment |
| `varme/losere.py` | Numerisk løser (FTCS) og Fourier-cosinusrekke |
| `varme/visualisering.py` | Det interaktive vinduet |
