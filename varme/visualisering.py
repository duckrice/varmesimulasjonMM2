"""Interaktivt simuleringsvindu for temperaturen langs en stav.

Alle kontroller er matplotlib-widgets inne i figuren, så samme vindu virker
både som eget vindu (`python simulering.py`) og i Jupyter med `%matplotlib widget`.

Redigering av starttemperaturen (gjøres på grafen eller staven):
* Venstreklikk:        legg til en varmekilde (temperatur og bredde fra feltene nederst).
* Venstre + dra:       tegn temperaturkurven fritt med musa.
* Høyreklikk / -dra:   fjern varme (sett tilbake til grunntemperaturen).
Redigerer du midt i en simulering, blir tilstanden der den nye starten (t = 0).

Tid:
* Tidsbryteren er skalert i √t: halvveis på bryteren er en firedel av sluttiden.
  Da får starten, der det skjer mest, mest plass.
* Piltaster ← → eller musehjulet flytter tiden et lite steg (hold Shift for 10 steg).
* Mellomrom starter og stopper avspillingen.
* «Slutt [s]» bestemmer hvor langt bryteren går. Mindre sluttid gir finere kontroll.
"""

import time

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.widgets import Button, Slider, TextBox

from .losere import FourierCosinusrekke, simuler_eksplisitt
from .materialer import Materiale
from .rutenett import Rutenett
from .startbetingelser import varmt_segment

VENSTRE, HOYRE = 1, 3

# Avspillingsfarter. Ved ×1 tar hele tidsbryteren 40 sekunder.
FARTER = (0.05, 0.1, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0)
SEKUNDER_VED_FART_1 = 40.0


def formater_tall(verdi: float) -> str:
    return f"{verdi:g}"


def les_tall(tekst: str) -> float:
    """Tolker tekst som tall, og godtar både komma og punktum som desimaltegn."""
    return float(tekst.strip().replace(",", "."))


class InteraktivStav:
    def __init__(
        self,
        materiale: Materiale,
        lengde: float = 1.0,  # [m]
        antall_celler: int = 200,
        T_omgivelse: float = 20.0,  # grunntemperatur [°C]
        T_kilde: float = 300.0,  # temperatur for nye varmekilder [°C]
        bredde: float = 0.05,  # bredde for nye varmekilder [m]
        T_akse: tuple[float, float] = (0.0, 400.0),  # temperaturområde på aksen [°C]
        antall_bilder: int = 1000,  # lagrede tidspunkt; tider imellom regnes ut ved behov
        antall_steg: int = 1000,  # piltast/musehjul flytter 1/antall_steg av bryteren
        fourier_ledd: int = 20,
        maks_fourier_ledd: int = 100,
    ):
        self.rutenett = Rutenett.stav(lengde, antall_celler)
        self.L = lengde
        self.x = self.rutenett.sentre[0]
        # Plotting: numerisk løsning i cellesentrene, forlenget flatt ut til endeflatene
        # (ingen varmestrøm over kanten), og Fourier-rekka på et fint rutenett fra 0 til L.
        self.x_numerisk = np.concatenate([[0.0], self.x, [lengde]])
        self.x_fourier = np.linspace(0.0, lengde, 1001)
        self.k = materiale.varmeledningsevne
        self.rho = materiale.tetthet
        self.c_p = materiale.spesifikk_varme
        self.T_omgivelse = T_omgivelse
        self.T_kilde = T_kilde
        self.bredde = bredde
        self.T_akse = T_akse
        self.antall_bilder = antall_bilder
        self.antall_steg = antall_steg
        self.maks_fourier_ledd = maks_fourier_ledd

        self.u0 = varmt_segment(self.rutenett, lengde / 2, bredde, T_omgivelse, T_kilde)
        self.t, self.u = 0.0, self.u0
        self._endrer_kontroller = False
        self._redigering = None
        self._spiller = False

        self._bygg_figur()
        self._bygg_kontroller(fourier_ledd)
        self._nytt_materiale()

    @property
    def alfa(self) -> float:
        """Termisk diffusivitet alfa = k / (rho · c_p)  [m²/s]."""
        return self.k / (self.rho * self.c_p)

    # --- Beregning -----------------------------------------------------------

    def _nytt_materiale(self) -> None:
        """Ny alfa gir ny tidsskala: etter t ≈ 0.3·L²/alfa er staven nesten jevnt varm."""
        self.t_slutt = float(f"{0.3 * self.L**2 / self.alfa:.3g}")
        self._sett_felt(self.slutt_felt, formater_tall(self.t_slutt))
        self.alfa_tekst.set_text(f"α = k/(ρ·c$_p$) = {self.alfa:.3g} m²/s")
        self._ny_tidsakse()

    def _ny_tidsakse(self, behold_tid: bool = False) -> None:
        # Lagrede bilder jevnt fordelt i √t, som tidsbryteren
        self.tider = self.t_slutt * np.linspace(0.0, 1.0, self.antall_bilder + 1) ** 2
        self.simuler(t=min(self.t, self.t_slutt) if behold_tid else 0.0)

    def tid_fra_bryter(self, posisjon: float) -> float:
        return self.t_slutt * posisjon**2

    def bryter_fra_tid(self, t: float) -> float:
        return float(np.sqrt(t / self.t_slutt))

    def simuler(self, t: float = 0.0) -> None:
        """Regner ut begge løsningene fra gjeldende starttemperatur og viser tid t."""
        self.bilder = simuler_eksplisitt(self.u0, self.rutenett, self.alfa, self.tider)
        self.fourier = FourierCosinusrekke(self.u0, self.rutenett, self.alfa, self.maks_fourier_ledd)
        snitt = self.fourier.a[0]
        self.likevekt_linje.set_ydata([snitt, snitt])
        self.forklaring.get_texts()[2].set_text(f"Likevekt: {snitt:.1f} °C")
        self.fourier_linje.set_visible(True)
        self._juster_temperaturakse()
        self.vis_tid(t)

    def vis_tid(self, t: float) -> None:
        """Viser temperaturen ved et vilkårlig tidspunkt t [s], uten avrunding.

        Starter fra nærmeste lagrede bilde før t og simulerer resten av veien.
        """
        t = float(np.clip(t, 0.0, self.t_slutt))
        i = int(np.searchsorted(self.tider, t, side="right")) - 1
        u = self.bilder[i]
        if t > self.tider[i]:
            u = simuler_eksplisitt(u, self.rutenett, self.alfa, [t - self.tider[i]])[0]
        self._sett_tilstand(t, u)

    def steg(self, antall: int) -> None:
        """Flytter tidsbryteren `antall` små steg (negativt = bakover)."""
        self._stopp()
        posisjon = np.clip(self.tid_bryter.val + antall / self.antall_steg, 0.0, 1.0)
        self.vis_tid(self.tid_fra_bryter(posisjon))

    def _sett_tilstand(self, t: float, u: np.ndarray) -> None:
        self.t, self.u = t, u
        # Oppdater kontrollene uten at de utløser nye beregninger
        self._endrer_kontroller = True
        try:
            self.tid_bryter.set_val(self.bryter_fra_tid(t))
            # Direkte på teksten: TextBox.set_val tegner hele figuren på nytt, som gjør avspillingen treg
            self.tid_felt.text_disp.set_text(f"{t:.6g}")
        finally:
            self._endrer_kontroller = False
        self.tegn()

    def _sett_felt(self, felt: TextBox, tekst: str) -> None:
        self._endrer_kontroller = True
        try:
            felt.set_val(tekst)
        finally:
            self._endrer_kontroller = False

    # --- Redigering med mus --------------------------------------------------

    def _mus_til_data(self, hendelse) -> tuple[float, float]:
        """Musposisjon i (x [m], T [°C]), også når musa er utenfor grafen."""
        x_cm, T = self.ax_temp.transData.inverted().transform((hendelse.x, hendelse.y))
        T_min, T_maks = self.ax_temp.get_ylim()
        return float(np.clip(x_cm / 100, 0.0, self.L)), float(np.clip(T, T_min, T_maks))

    def _trykk(self, hendelse) -> None:
        if hendelse.inaxes not in (self.ax_temp, self.ax_stav) or hendelse.button not in (VENSTRE, HOYRE):
            return
        self._stopp()
        # Redigering midt i simuleringen: tilstanden nå blir den nye starten
        self.u0 = np.array(self.u, dtype=float)
        x, T = self._mus_til_data(hendelse)
        self._redigering = {
            "knapp": hendelse.button,
            "på_grafen": hendelse.inaxes is self.ax_temp,
            "start_piksel": (hendelse.x, hendelse.y),
            "forrige": (x, T),
            "dratt": False,
        }
        if hendelse.button == HOYRE:
            self._fjern_varme(x, x)
            self._vis_redigering()

    def _bevegelse(self, hendelse) -> None:
        r = self._redigering
        if r is None:
            return
        x, T = self._mus_til_data(hendelse)
        x0, T0 = r["forrige"]
        if r["knapp"] == HOYRE:
            self._fjern_varme(x0, x)
        elif r["på_grafen"]:
            if not r["dratt"]:
                px, py = r["start_piksel"]
                if np.hypot(hendelse.x - px, hendelse.y - py) < 4:
                    return  # for liten bevegelse: fortsatt et klikk
                r["dratt"] = True
            self._tegn_strek(x0, T0, x, T)
        r["forrige"] = (x, T)
        self._vis_redigering()

    def _slipp(self, hendelse) -> None:
        r = self._redigering
        if r is None:
            return
        self._redigering = None
        if r["knapp"] == VENSTRE and not r["dratt"]:
            x, _ = r["forrige"]
            self.u0[np.abs(self.x - x) <= self.bredde / 2] = self.T_kilde
            self.u0[np.argmin(np.abs(self.x - x))] = self.T_kilde  # minst én celle
        self.simuler()

    def _tegn_strek(self, x0: float, T0: float, x1: float, T1: float) -> None:
        """Setter cellene mellom x0 og x1 til en rett linje fra T0 til T1."""
        if x1 < x0:
            x0, T0, x1, T1 = x1, T1, x0, T0
        inni = (self.x >= x0) & (self.x <= x1)
        if x1 > x0:
            self.u0[inni] = T0 + (T1 - T0) * (self.x[inni] - x0) / (x1 - x0)
        self.u0[np.argmin(np.abs(self.x - x1))] = T1

    def _fjern_varme(self, x0: float, x1: float) -> None:
        """Setter cellene innenfor en bredde rundt strekningen x0–x1 til grunntemperaturen."""
        lav, hoy = min(x0, x1) - self.bredde / 2, max(x0, x1) + self.bredde / 2
        self.u0[(self.x >= lav) & (self.x <= hoy)] = self.T_omgivelse

    def _vis_redigering(self) -> None:
        """Viser starttemperaturen mens den redigeres (simuleringen kjøres når musa slippes)."""
        self.fourier_linje.set_visible(False)
        self._sett_tilstand(0.0, self.u0)

    # --- Avspilling ----------------------------------------------------------

    def _animerte_artister(self) -> list:
        """Det som endrer seg under avspilling. Bare dette tegnes på nytt per bilde (blitting)."""
        b = self.tid_bryter
        return [self.numerisk_linje, self.fourier_linje, self.stav_bilde, b.poly, b._handle, self.tid_felt.text_disp]

    def _spill_pause(self, _=None) -> None:
        if self._spiller:
            self._stopp()
            return
        if self.t >= self.t_slutt:
            self.vis_tid(0.0)
        self._spiller = True
        self.spill_knapp.label.set_text("Pause")
        # Tegn alt unntatt de animerte delene én gang, og ta vare på bakgrunnen
        self.tid_bryter.drawon = False
        for artist in self._animerte_artister():
            artist.set_animated(True)
        self.fig.canvas.draw()
        self._bakgrunn = self.fig.canvas.copy_from_bbox(self.fig.bbox)
        self._posisjon = self.tid_bryter.val
        self._forrige_tikk = time.perf_counter()
        self.timer.start()

    def _stopp(self) -> None:
        if not self._spiller:
            return
        self._spiller = False
        self.timer.stop()
        self.tid_bryter.drawon = True
        for artist in self._animerte_artister():
            artist.set_animated(False)
        self.spill_knapp.label.set_text("Spill av")
        self.fig.canvas.draw_idle()

    @property
    def fart(self) -> float:
        return FARTER[int(self.fart_bryter.val)]

    def _tikk(self) -> None:
        """Går fram etter klokka, så farten er lik uansett hvor fort det tegnes.
        Hvert bilde regnes ut ved eksakt tid, så også sakte avspilling er jevn."""
        if not self._spiller:
            return
        naa = time.perf_counter()
        self._posisjon = min(self._posisjon + (naa - self._forrige_tikk) * self.fart / SEKUNDER_VED_FART_1, 1.0)
        self._forrige_tikk = naa
        self.vis_tid(self.tid_fra_bryter(self._posisjon))
        if self._posisjon >= 1.0:
            self._stopp()

    # --- Tekstfelt -----------------------------------------------------------

    def _tallfelt(self, felt: TextBox, attributt: str, skala: float = 1.0, positiv: bool = True, etter=None):
        """Kobler et tekstfelt til en egenskap. Ugyldige tall settes tilbake."""

        def ved_innsending(tekst: str) -> None:
            if self._endrer_kontroller:
                return
            gammel = getattr(self, attributt)
            try:
                verdi = les_tall(tekst) * skala
                if positiv and verdi <= 0:
                    raise ValueError
            except ValueError:
                verdi = gammel
            self._sett_felt(felt, formater_tall(verdi / skala))
            if verdi != gammel:
                setattr(self, attributt, verdi)
                if etter is not None:
                    etter()

        felt.on_submit(ved_innsending)

    def _tid_skrevet(self, tekst: str) -> None:
        if self._endrer_kontroller:
            return
        try:
            t = les_tall(tekst)
        except ValueError:
            t = self.t
        self._stopp()
        self.vis_tid(t)

    # --- Figur ---------------------------------------------------------------

    def _bygg_figur(self) -> None:
        L_cm = self.L * 100

        self.fig = plt.figure(figsize=(10, 6.5))
        lerret = self.fig.canvas
        if lerret.manager is not None:
            lerret.manager.set_window_title("Varmesimulering")
            # Slå av matplotlibs hurtigtaster (f.eks. «l» = logaritmisk akse)
            lerret.mpl_disconnect(lerret.manager.key_press_handler_id)
        for egenskap in ("header_visible", "footer_visible", "toolbar_visible", "resizable"):
            if hasattr(lerret, egenskap):  # ipympl i Jupyter
                setattr(lerret, egenskap, False)

        self.ax_temp = self.fig.add_axes((0.07, 0.40, 0.91, 0.57))
        self.ax_stav = self.fig.add_axes((0.07, 0.335, 0.91, 0.045), sharex=self.ax_temp)

        (self.numerisk_linje,) = self.ax_temp.plot(
            self.x_numerisk * 100, self._med_ender(self.u0), color="tab:red", lw=2, label="Numerisk")
        (self.fourier_linje,) = self.ax_temp.plot(
            self.x_fourier * 100, np.zeros_like(self.x_fourier), color="tab:blue", lw=1.5, ls="--", label="Fourier"
        )
        (self.likevekt_linje,) = self.ax_temp.plot([0, L_cm], [0, 0], color="gray", lw=1, ls=":", label="Likevekt")
        self.forklaring = self.ax_temp.legend(loc="upper right", fontsize=8, framealpha=0.8)
        self.ax_temp.set_xlim(0, L_cm)
        self.ax_temp.set_ylabel("Temperatur [°C]")
        self.ax_temp.grid(alpha=0.3)
        self.ax_temp.tick_params(labelbottom=False)

        self.stav_bilde = self.ax_stav.imshow(
            self.u0[np.newaxis, :], aspect="auto", extent=(0, L_cm, 0, 1), cmap="inferno"
        )
        self.ax_stav.set_yticks([])
        self.ax_stav.set_xlabel("x [cm]", labelpad=1)

        lerret.mpl_connect("button_press_event", self._trykk)
        lerret.mpl_connect("motion_notify_event", self._bevegelse)
        lerret.mpl_connect("button_release_event", self._slipp)
        lerret.mpl_connect("draw_event", self._full_tegning)
        lerret.mpl_connect("key_press_event", self._tast)
        lerret.mpl_connect("scroll_event", self._musehjul)

    def _tast(self, hendelse) -> None:
        if any(f.capturekeystrokes for f in self._tekstfelt):
            return  # det skrives i et tekstfelt
        if hendelse.key in ("left", "right", "shift+left", "shift+right"):
            retning = 1 if hendelse.key.endswith("right") else -1
            self.steg(retning * (10 if hendelse.key.startswith("shift") else 1))
        elif hendelse.key == " ":
            self._spill_pause()

    def _musehjul(self, hendelse) -> None:
        retning = 1 if hendelse.step > 0 else -1
        self.steg(retning * (10 if hendelse.key == "shift" else 1))

    def _full_tegning(self, _) -> None:
        """Hvis hele figuren tegnes på nytt under avspilling (f.eks. ved endret
        vindusstørrelse), må bakgrunnen for blittingen tas vare på på nytt."""
        if self._spiller:
            self._bakgrunn = self.fig.canvas.copy_from_bbox(self.fig.bbox)

    def _bygg_kontroller(self, fourier_ledd: int) -> None:
        def akse(x, y, b, h=0.045):
            return self.fig.add_axes((x, y, b, h))

        # Rad 1: avspilling og tid
        self.spill_knapp = Button(akse(0.07, 0.19, 0.09), "Spill av")
        self.spill_knapp.on_clicked(self._spill_pause)
        self.tid_bryter = Slider(akse(0.19, 0.19, 0.40), "", 0.0, 1.0, valinit=0.0)
        self.tid_bryter.valtext.set_visible(False)
        self.tid_bryter.on_changed(
            lambda p: None if self._endrer_kontroller else self.vis_tid(self.tid_fra_bryter(p))
        )
        self.tid_felt = TextBox(akse(0.66, 0.19, 0.10), "t [s] ", textalignment="right")
        self.tid_felt.on_submit(self._tid_skrevet)
        self.slutt_felt = TextBox(akse(0.88, 0.19, 0.10), "Slutt [s] ", textalignment="right")
        self._tallfelt(self.slutt_felt, "t_slutt", etter=lambda: self._ny_tidsakse(behold_tid=True))

        # Rad 2: materiale
        self.k_felt = TextBox(akse(0.10, 0.11, 0.08), "k ", formater_tall(self.k), textalignment="right")
        self.rho_felt = TextBox(akse(0.23, 0.11, 0.08), "ρ ", formater_tall(self.rho), textalignment="right")
        self.cp_felt = TextBox(akse(0.37, 0.11, 0.08), "c$_p$ ", formater_tall(self.c_p), textalignment="right")
        for felt, navn in ((self.k_felt, "k"), (self.rho_felt, "rho"), (self.cp_felt, "c_p")):
            self._tallfelt(felt, navn, etter=self._nytt_materiale)
        self.alfa_tekst = self.fig.text(0.47, 0.13, "", va="center", fontsize=9)

        # Rad 3: varmekilder, Fourier og fart
        self.kilde_felt = TextBox(akse(0.15, 0.03, 0.07), "Kilde [°C] ", formater_tall(self.T_kilde), textalignment="right")
        self.bredde_felt = TextBox(akse(0.33, 0.03, 0.06), "Bredde [cm] ", formater_tall(self.bredde * 100), textalignment="right")
        self._tallfelt(self.kilde_felt, "T_kilde", positiv=False, etter=self._juster_temperaturakse)
        self._tallfelt(self.bredde_felt, "bredde", skala=0.01)
        self.ledd_bryter = Slider(
            akse(0.50, 0.03, 0.15), "Fourier-ledd", 1, self.maks_fourier_ledd, valinit=fourier_ledd, valstep=1
        )
        self.ledd_bryter.on_changed(lambda _: self.tegn())
        self.fart_bryter = Slider(
            akse(0.76, 0.03, 0.08), "Fart", 0, len(FARTER) - 1, valinit=FARTER.index(1.0), valstep=1
        )
        self.fart_bryter.valtext.set_text("×1")
        self.fart_bryter.on_changed(lambda _: self.fart_bryter.valtext.set_text(f"×{self.fart:g}"))
        self.nullstill_knapp = Button(akse(0.88, 0.03, 0.10), "Nullstill")
        self.nullstill_knapp.on_clicked(self._nullstill)

        self._tekstfelt = [self.tid_felt, self.slutt_felt, self.k_felt, self.rho_felt, self.cp_felt,
                           self.kilde_felt, self.bredde_felt]
        self.timer = self.fig.canvas.new_timer(interval=15)
        self.timer.add_callback(self._tikk)

    def _nullstill(self, _=None) -> None:
        self._stopp()
        self.u0 = np.full_like(self.u0, self.T_omgivelse)
        self.simuler()

    def _juster_temperaturakse(self) -> None:
        """Temperaturaksen dekker minst T_akse, og utvides hvis noe er varmere/kaldere."""
        T_min = min(self.T_akse[0], self.u0.min(), self.T_kilde)
        T_maks = max(self.T_akse[1], self.u0.max(), self.T_kilde)
        margin = 0.03 * (T_maks - T_min)
        self.ax_temp.set_ylim(T_min - margin, T_maks + margin)
        self.stav_bilde.set_clim(T_min, T_maks)
        self.fig.canvas.draw_idle()

    @staticmethod
    def _med_ender(u: np.ndarray) -> np.ndarray:
        return np.concatenate([[u[0]], u, [u[-1]]])

    def tegn(self) -> None:
        ledd = int(self.ledd_bryter.val)
        self.numerisk_linje.set_ydata(self._med_ender(self.u))
        if self.fourier_linje.get_visible():
            self.fourier_linje.set_ydata(self.fourier.evaluer(self.x_fourier, self.t, ledd))
        self.forklaring.get_texts()[1].set_text(f"Fourier ({ledd} ledd)")
        self.stav_bilde.set_data(self.u[np.newaxis, :])
        if self._spiller:
            lerret = self.fig.canvas
            lerret.restore_region(self._bakgrunn)
            for artist in self._animerte_artister():
                self.fig.draw_artist(artist)
            lerret.blit(self.fig.bbox)
        else:
            self.fig.canvas.draw_idle()
