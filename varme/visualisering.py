"""Interaktivt simuleringsvindu for varmeledning i en stav (1D) eller en plate (2D).

Knappen «Bytt til 2D» / «Bytt til 1D» veksler mellom de to. Materialverdiene og
innstillingene for varmekilder og rand følger med.

Rand: knappen «Rand: isolert» / «Rand: fast temp.» velger mellom isolerte kanter
(Neumann) og kanter som holdes på en fast temperatur (Dirichlet). For staven kan
de to endene ha ulik temperatur; for platen har alle fire kantene samme.

Redigering av starttemperaturen:
* Stav:  venstreklikk legger til en varmekilde, venstre + dra tegner kurven fritt.
* Plate: venstreklikk legger til en kvadratisk varmekilde, venstre + dra maler varme.
* Høyreklikk / høyre + dra fjerner varme (setter tilbake til grunntemperaturen).
Redigerer du midt i en simulering, blir tilstanden der den nye starten (t = 0).

Avlesning: skriv inn en posisjon i «Les av ved x [cm]» (og «y [cm]» for platen),
så vises temperaturen der ved gjeldende tid, for både den numeriske løsningen og
Fourier-rekka. Sammen med feltet «t [s]» gir det eksakte verdier i et valgt punkt
ved et valgt tidspunkt. Et merke i figuren viser hvor det leses av.

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
from matplotlib.colors import Normalize, PowerNorm
from matplotlib.widgets import Button, Slider, TextBox
from scipy.interpolate import RegularGridInterpolator

from .losere import (
    FourierCosinusrekke,
    FourierCosinusrekke2D,
    FourierSinusrekke,
    FourierSinusrekke2D,
    simuler_eksplisitt,
)
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


class _Visning:
    """Det som er felles for stav og plate: tid, avspilling, materialfelt og
    lagring av beregnede bilder. Underklassene tegner og redigerer temperaturen.

    Tegning: Å tegne hele figuren tar lang tid, mest på grunn av all teksten.
    Derfor tas det vare på et bilde av alt som står stille (bakgrunnen), og ved
    hver oppdatering tegnes bare det som endrer seg oppå det (blitting). Hele
    figuren tegnes bare på nytt når noe i bakgrunnen endres, f.eks. en akse.
    """

    ANTALL_BILDER = 1000  # lagrede tidspunkt; tider imellom regnes ut ved behov
    ANTALL_STEG = 1000  # piltast/musehjul flytter 1/ANTALL_STEG av tidsbryteren
    MAKS_FOURIER_LEDD = 100
    BYTT_TIL = ""  # teksten på knappen som bytter til den andre visningen

    def __init__(self, fig, rutenett: Rutenett, innstillinger: dict, bytt_modus):
        self.fig = fig
        self.rutenett = rutenett
        self.L = rutenett.lengder[0]
        self.x = rutenett.sentre[0]
        self.k = innstillinger["k"]
        self.rho = innstillinger["rho"]
        self.c_p = innstillinger["c_p"]
        self.T_omgivelse = innstillinger["T_omgivelse"]
        self.T_kilde = innstillinger["T_kilde"]
        self.bredde = innstillinger["bredde"]
        self.fast_rand = innstillinger["fast_rand"]  # False: isolerte kanter, True: fast temperatur
        self.T_venstre = innstillinger["T_venstre"]  # stavens ender
        self.T_hoyre = innstillinger["T_hoyre"]
        self.T_kant = innstillinger["T_kant"]  # platens kanter
        self.les_x = innstillinger["les_x"]  # punktet der temperaturen leses av [m]
        self.les_y = innstillinger["les_y"]
        self._bytt_modus = bytt_modus

        self.u0 = self._startbetingelse()
        self.t, self.u = 0.0, self.u0
        self.bilder = np.empty((self.ANTALL_BILDER + 1, *self.u0.shape))
        self._endrer_kontroller = False
        self._redigering = None
        self._spiller = False
        self._bakgrunn = None
        self._widgets = []
        self._koblinger = []

        self._bygg_figur()
        self._bygg_kontroller(innstillinger["fourier_ledd"])
        for hendelse, funksjon in (
            ("button_press_event", self._trykk),
            ("motion_notify_event", self._bevegelse),
            ("button_release_event", self._slipp),
            ("draw_event", self._full_tegning),
            ("key_press_event", self._tast),
            ("scroll_event", self._musehjul),
        ):
            self._koblinger.append(fig.canvas.mpl_connect(hendelse, funksjon))
        for artist in self._animerte_artister():
            artist.set_animated(True)  # tegnes av oss, ikke sammen med bakgrunnen
        self._nytt_materiale()

    def innstillinger(self) -> dict:
        """Verdiene som følger med når man bytter mellom stav og plate."""
        return {
            "k": self.k, "rho": self.rho, "c_p": self.c_p,
            "T_omgivelse": self.T_omgivelse, "T_kilde": self.T_kilde, "bredde": self.bredde,
            "fast_rand": self.fast_rand, "T_venstre": self.T_venstre, "T_hoyre": self.T_hoyre,
            "T_kant": self.T_kant, "les_x": self.les_x, "les_y": self.les_y,
            "fourier_ledd": int(self.ledd_bryter.val),
        }

    def lukk(self) -> None:
        """Kobler fra alt før figuren tømmes og den andre visningen bygges."""
        self._stopp()
        self.timer.stop()
        for kobling in self._koblinger:
            self.fig.canvas.mpl_disconnect(kobling)
        for widget in self._widgets:
            widget.disconnect_events()

    @property
    def alfa(self) -> float:
        """Termisk diffusivitet alfa = k / (rho · c_p)  [m²/s]."""
        return self.k / (self.rho * self.c_p)

    # --- Underklassene fyller inn ----------------------------------------------

    def _startbetingelse(self) -> np.ndarray:
        raise NotImplementedError

    def _bygg_figur(self) -> None:
        raise NotImplementedError

    def _rand(self):
        """Faste randtemperaturer slik løseren vil ha dem, eller None når kantene er isolert."""
        raise NotImplementedError

    def _bygg_randfelt(self, y: float) -> list:
        """Lager tekstfeltene for randtemperaturene (på raden med høyde y) og returnerer dem."""
        raise NotImplementedError

    def _ny_start(self) -> None:
        """Kalles når starttemperaturen eller randen er endret: lager ny Fourier-rekke."""
        raise NotImplementedError

    def _likevekt(self) -> tuple:
        """(tekstobjektet som viser likevektstemperaturen, teksten det skal ha)."""
        raise NotImplementedError

    def _data_artister(self) -> list:
        """Det i figuren som endrer seg med temperaturen, i den rekkefølgen det skal tegnes."""
        raise NotImplementedError

    def _avlesning(self) -> tuple[float, float]:
        """(numerisk, Fourier): temperaturen i avlesningspunktet ved gjeldende tid."""
        raise NotImplementedError

    def _flytt_merke(self) -> None:
        """Flytter merket i figuren til avlesningspunktet."""
        raise NotImplementedError

    def _oppdater_data(self) -> None:
        """Setter gjeldende temperatur (self.u ved tid self.t) inn i figuren."""
        raise NotImplementedError

    def _trykk(self, hendelse) -> None:
        raise NotImplementedError

    def _bevegelse(self, hendelse) -> None:
        raise NotImplementedError

    def _slipp(self, hendelse) -> None:
        if self._redigering is not None:
            self._redigering = None
            self._vis_likevekt()

    def _juster_temperaturskala(self) -> None:
        raise NotImplementedError

    def _vis_likevekt(self) -> None:
        """Oppdaterer teksten med likevektstemperaturen. Den er en del av
        bakgrunnen (tekst er tregt å tegne), så mens man drar med musa venter
        vi til knappen slippes."""
        tekstobjekt, tekst = self._likevekt()
        if self._redigering is None and tekstobjekt.get_text() != tekst:
            tekstobjekt.set_text(tekst)
            self.fig.canvas.draw_idle()

    # --- Beregning -----------------------------------------------------------

    def _nytt_materiale(self) -> None:
        """Ny alfa gir ny tidsskala: etter t ≈ 0.3·L²/alfa er temperaturen nesten jevn."""
        self.t_slutt = float(f"{0.3 * self.L**2 / self.alfa:.3g}")
        self._sett_felt(self.slutt_felt, formater_tall(self.t_slutt))
        self.alfa_tekst.set_text(f"α = k/(ρ·c$_p$) = {self.alfa:.3g} m²/s")
        self._ny_tidsakse()

    def _ny_tidsakse(self, behold_tid: bool = False) -> None:
        # Lagrede bilder jevnt fordelt i √t, som tidsbryteren
        self.tider = self.t_slutt * np.linspace(0.0, 1.0, self.ANTALL_BILDER + 1) ** 2
        self.simuler(t=min(self.t, self.t_slutt) if behold_tid else 0.0)

    def tid_fra_bryter(self, posisjon: float) -> float:
        return self.t_slutt * posisjon**2

    def bryter_fra_tid(self, t: float) -> float:
        return float(np.sqrt(t / self.t_slutt))

    def simuler(self, t: float = 0.0) -> None:
        """Starter på nytt fra gjeldende starttemperatur og viser tid t.

        Bildene regnes ikke ut på forhånd, men etter hvert som de trengs
        (se `_beregn_til`). Derfor går redigering raskt, også for platen.
        """
        self.bilder[0] = self.u0
        self._antall_beregnet = 1
        self._forrige_viste = None
        self._ny_start()
        self._vis_likevekt()
        self._juster_temperaturskala()
        self.vis_tid(t)

    def _beregn_til(self, i: int) -> None:
        """Sørger for at de lagrede bildene 0 … i er regnet ut."""
        siste = self._antall_beregnet - 1
        if i > siste:
            tider = self.tider[siste + 1 : i + 1] - self.tider[siste]
            self.bilder[siste + 1 : i + 1] = simuler_eksplisitt(
                self.bilder[siste], self.rutenett, self.alfa, tider, rand=self._rand()
            )
            self._antall_beregnet = i + 1

    def vis_tid(self, t: float, fortsett: bool = False) -> None:
        """Viser temperaturen ved et vilkårlig tidspunkt t [s], uten avrunding.

        Starter fra nærmeste lagrede bilde før t og simulerer resten av veien.
        Samme t gir dermed alltid nøyaktig samme tall. Under avspilling
        (`fortsett=True`) fortsetter vi i stedet fra forrige viste tilstand når
        den er nærmere, for det går raskere.
        """
        t = float(np.clip(t, 0.0, self.t_slutt))
        i = int(np.searchsorted(self.tider, t, side="right")) - 1
        self._beregn_til(i)
        t_fra, u = self.tider[i], self.bilder[i]
        if fortsett and self._forrige_viste is not None and t_fra <= self._forrige_viste[0] <= t:
            t_fra, u = self._forrige_viste
        if t > t_fra:
            u = simuler_eksplisitt(u, self.rutenett, self.alfa, [t - t_fra], rand=self._rand())[0]
        self._forrige_viste = (t, u)
        self._sett_tilstand(t, u)

    def steg(self, antall: int) -> None:
        """Flytter tidsbryteren `antall` små steg (negativt = bakover)."""
        self._stopp()
        posisjon = np.clip(self.tid_bryter.val + antall / self.ANTALL_STEG, 0.0, 1.0)
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

    def _start_redigering(self) -> None:
        """Redigering midt i simuleringen: tilstanden nå blir den nye starten."""
        self._stopp()
        self.u0 = np.array(self.u, dtype=float)

    def _bytt_rand(self, _=None) -> None:
        """Veksler mellom isolerte kanter og fast temperatur på kantene, og
        starter simuleringen på nytt fra samme starttemperatur."""
        self._stopp()
        self.fast_rand = not self.fast_rand
        self._vis_randvalg()
        self.simuler()

    def _vis_randvalg(self) -> None:
        """Feltene for randtemperatur vises bare når kantene har fast temperatur."""
        self.rand_knapp.label.set_text("Rand: fast temp." if self.fast_rand else "Rand: isolert")
        for felt in self._randfelt:
            felt.ax.set_visible(self.fast_rand)
            felt.active = self.fast_rand
        self.fig.canvas.draw_idle()

    def _nullstill(self, _=None) -> None:
        self._stopp()
        self.u0 = np.full_like(self.u0, self.T_omgivelse)
        self.simuler()

    # --- Avspilling ----------------------------------------------------------

    def _animerte_artister(self) -> list:
        """Alt som tegnes oppå bakgrunnen: temperaturen, bryterne og tiden."""
        brytere = []
        for bryter in (self.tid_bryter, self.ledd_bryter, self.fart_bryter):
            brytere += [bryter.poly, bryter._handle, bryter.valtext]
        avlesning = [self.les_numerisk_tekst, self.les_fourier_tekst]
        return [*self._data_artister(), *brytere, self.tid_felt.text_disp, *avlesning]

    def _spill_pause(self, _=None) -> None:
        if self._spiller:
            self._stopp()
            return
        if self.t >= self.t_slutt:
            self.vis_tid(0.0)
        self._spiller = True
        self.spill_knapp.label.set_text("Pause")
        self.fig.canvas.draw_idle()
        self._posisjon = self.tid_bryter.val
        self._forrige_tikk = time.perf_counter()
        self.timer.start()

    def _stopp(self) -> None:
        if not self._spiller:
            return
        self._spiller = False
        self.timer.stop()
        self.spill_knapp.label.set_text("Spill av")
        self.fig.canvas.draw_idle()

    @property
    def fart(self) -> float:
        return FARTER[int(self.fart_bryter.val)]

    def _fart_endret(self, _) -> None:
        self.fart_bryter.valtext.set_text(f"×{self.fart:g}")
        self._blit()

    def _tikk(self) -> None:
        """Går fram etter klokka, så farten er lik uansett hvor fort det tegnes.
        Hvert bilde regnes ut ved eksakt tid, så også sakte avspilling er jevn."""
        if not self._spiller:
            return
        naa = time.perf_counter()
        self._posisjon = min(self._posisjon + (naa - self._forrige_tikk) * self.fart / SEKUNDER_VED_FART_1, 1.0)
        self._forrige_tikk = naa
        self.vis_tid(self.tid_fra_bryter(self._posisjon), fortsett=True)
        if self._posisjon >= 1.0:
            self._stopp()

    def _full_tegning(self, _) -> None:
        """Hele figuren er nettopp tegnet (uten de animerte delene): ta vare på
        den som bakgrunn, og tegn de animerte delene oppå."""
        if self.fig.canvas.is_saving():
            return  # ved lagring til fil tegner matplotlib alt selv
        self._bakgrunn = self.fig.canvas.copy_from_bbox(self.fig.bbox)
        for artist in self._animerte_artister():
            self.fig.draw_artist(artist)

    def _blit(self) -> None:
        """Viser bakgrunnen med de animerte delene tegnet oppå."""
        lerret = self.fig.canvas
        if self._bakgrunn is None:  # figuren er ikke tegnet ennå
            lerret.draw_idle()
            return
        lerret.restore_region(self._bakgrunn)
        for artist in self._animerte_artister():
            self.fig.draw_artist(artist)
        lerret.blit(self.fig.bbox)

    def tegn(self) -> None:
        """Oppdaterer figuren med gjeldende temperatur."""
        self._oppdater_data()
        numerisk, fourier = self._avlesning()
        self.les_numerisk_tekst.set_text(f"{numerisk:.3f} °C")
        self.les_fourier_tekst.set_text(f"{fourier:.3f} °C")
        self._blit()

    def _nytt_avlesningspunkt(self) -> None:
        self._flytt_merke()
        self.fig.canvas.draw_idle()
        self.tegn()

    # --- Tastatur og musehjul ------------------------------------------------

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

    # --- Tekstfelt -----------------------------------------------------------

    def _sett_felt(self, felt: TextBox, tekst: str) -> None:
        self._endrer_kontroller = True
        try:
            felt.set_val(tekst)
        finally:
            self._endrer_kontroller = False
        if not felt.capturekeystrokes:
            felt.cursor.set_visible(False)  # set_val viser skrivemerket selv om ingen skriver

    def _tallfelt(self, felt: TextBox, attributt: str, skala: float = 1.0, positiv: bool = True,
                  grenser: tuple[float, float] | None = None, etter=None):
        """Kobler et tekstfelt til en egenskap. Ugyldige tall settes tilbake, og
        tall utenfor `grenser` flyttes til nærmeste grense."""

        def ved_innsending(tekst: str) -> None:
            if self._endrer_kontroller:
                return
            gammel = getattr(self, attributt)
            try:
                verdi = les_tall(tekst) * skala
                if positiv and verdi <= 0:
                    raise ValueError
                if grenser is not None:
                    verdi = min(max(verdi, grenser[0]), grenser[1])
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

    # --- Kontroller ----------------------------------------------------------

    # Høyden (y) til radene med kontroller: tid, avlesning, materiale, rand og varmekilder
    RADER = (0.258, 0.198, 0.138, 0.078, 0.018)

    def _akse(self, x: float, y: float, bredde: float, hoyde: float = 0.039):
        return self.fig.add_axes((x, y, bredde, hoyde))

    def _knapp(self, akse, tekst: str, funksjon) -> Button:
        knapp = Button(akse, tekst)
        knapp.on_clicked(funksjon)
        self._widgets.append(knapp)
        return knapp

    def _tekstboks(self, akse, etikett: str, tekst: str = "") -> TextBox:
        felt = TextBox(akse, etikett, tekst, textalignment="right")
        self._widgets.append(felt)
        return felt

    def _bryter(self, akse, etikett: str, lav, hoy, **valg) -> Slider:
        bryter = Slider(akse, etikett, lav, hoy, **valg)
        self._widgets.append(bryter)
        return bryter

    def _bygg_kontroller(self, fourier_ledd: int) -> None:
        rad1, rad_avlesning, rad2, rad3, rad4 = self.RADER

        # Rad 1: avspilling og tid
        self.spill_knapp = self._knapp(self._akse(0.07, rad1, 0.09), "Spill av", self._spill_pause)
        self.tid_bryter = self._bryter(self._akse(0.19, rad1, 0.40), "", 0.0, 1.0, valinit=0.0)
        self.tid_bryter.valtext.set_visible(False)
        self.tid_bryter.drawon = False  # bryterne tegnes med blitting, ikke med full tegning
        self.tid_bryter.on_changed(
            lambda p: None if self._endrer_kontroller else self.vis_tid(self.tid_fra_bryter(p))
        )
        self.tid_felt = self._tekstboks(self._akse(0.66, rad1, 0.10), "t [s] ")
        self.tid_felt.on_submit(self._tid_skrevet)
        self.slutt_felt = self._tekstboks(self._akse(0.88, rad1, 0.10), "Slutt [s] ")
        self._tallfelt(self.slutt_felt, "t_slutt", etter=lambda: self._ny_tidsakse(behold_tid=True))

        # Avlesning av temperaturen i et valgt punkt
        self._avlesningsfelt = self._bygg_avlesningsfelt(rad_avlesning)
        midt = rad_avlesning + 0.0195
        self.fig.text(0.43, midt, "Numerisk:", va="center", fontsize=9)
        self.les_numerisk_tekst = self.fig.text(0.505, midt, "", va="center", fontsize=9, family="monospace")
        self.fig.text(0.68, midt, "Fourier:", va="center", fontsize=9)
        self.les_fourier_tekst = self.fig.text(0.74, midt, "", va="center", fontsize=9, family="monospace")

        # Rad 2: materiale, og bytte mellom stav og plate
        self.k_felt = self._tekstboks(self._akse(0.10, rad2, 0.08), "k ", formater_tall(self.k))
        self.rho_felt = self._tekstboks(self._akse(0.23, rad2, 0.08), "ρ ", formater_tall(self.rho))
        self.cp_felt = self._tekstboks(self._akse(0.37, rad2, 0.08), "c$_p$ ", formater_tall(self.c_p))
        for felt, navn in ((self.k_felt, "k"), (self.rho_felt, "rho"), (self.cp_felt, "c_p")):
            self._tallfelt(felt, navn, etter=self._nytt_materiale)
        self.alfa_tekst = self.fig.text(0.47, rad2 + 0.0195, "", va="center", fontsize=9)
        self.modus_knapp = self._knapp(self._akse(0.88, rad2, 0.10), self.BYTT_TIL, lambda _: self._bytt_modus())

        # Rad 3: randbetingelser
        self.rand_knapp = self._knapp(self._akse(0.07, rad3, 0.15), "", self._bytt_rand)
        self._randfelt = self._bygg_randfelt(rad3)
        self._vis_randvalg()

        # Rad 4: varmekilder, Fourier og fart
        self.kilde_felt = self._tekstboks(self._akse(0.15, rad4, 0.07), "Kilde [°C] ", formater_tall(self.T_kilde))
        self.bredde_felt = self._tekstboks(
            self._akse(0.33, rad4, 0.06), "Bredde [cm] ", formater_tall(self.bredde * 100)
        )
        self._tallfelt(self.kilde_felt, "T_kilde", positiv=False, etter=self._juster_temperaturskala)
        self._tallfelt(self.bredde_felt, "bredde", skala=0.01)
        self.ledd_bryter = self._bryter(
            self._akse(0.50, rad4, 0.15), "Fourier-ledd", 1, self.MAKS_FOURIER_LEDD,
            valinit=min(fourier_ledd, self.MAKS_FOURIER_LEDD), valstep=1,
        )
        self.ledd_bryter.drawon = False
        self.ledd_bryter.on_changed(lambda _: self.tegn())
        self.fart_bryter = self._bryter(
            self._akse(0.76, rad4, 0.08), "Fart", 0, len(FARTER) - 1, valinit=FARTER.index(1.0), valstep=1
        )
        self.fart_bryter.drawon = False
        self.fart_bryter.valtext.set_text("×1")
        self.fart_bryter.on_changed(self._fart_endret)
        self.nullstill_knapp = self._knapp(self._akse(0.88, rad4, 0.10), "Nullstill", self._nullstill)

        self._tekstfelt = [self.tid_felt, self.slutt_felt, self.k_felt, self.rho_felt, self.cp_felt,
                           self.kilde_felt, self.bredde_felt, *self._randfelt, *self._avlesningsfelt]
        self.timer = self.fig.canvas.new_timer(interval=15)
        self.timer.add_callback(self._tikk)

    def _bygg_avlesningsfelt(self, y: float) -> list:
        """Feltet for x-posisjonen der temperaturen leses av. Platen legger til et for y."""
        return [self._posisjonsfelt(0.20, y, "Les av ved x [cm] ", "les_x")]

    def _posisjonsfelt(self, x: float, y: float, etikett: str, attributt: str) -> TextBox:
        felt = self._tekstboks(self._akse(x, y, 0.06), etikett, formater_tall(getattr(self, attributt) * 100))
        self._tallfelt(felt, attributt, skala=0.01, positiv=False, grenser=(0.0, self.L),
                       etter=self._nytt_avlesningspunkt)
        return felt

    def _randtemperaturfelt(self, x: float, y: float, etikett: str, attributt: str) -> TextBox:
        """Et tekstfelt for en randtemperatur. Ny verdi starter simuleringen på nytt."""
        felt = self._tekstboks(self._akse(x, y, 0.07), etikett, formater_tall(getattr(self, attributt)))
        self._tallfelt(felt, attributt, positiv=False, etter=self.simuler)
        return felt


class Stavvisning(_Visning):
    """Temperaturen langs en stav: numerisk løsning og Fourier-rekke som to kurver."""

    BYTT_TIL = "Bytt til 2D"
    ANTALL_CELLER = 200
    T_AKSE = (-50.0, 400.0)  # temperaturområdet aksen minst skal dekke [°C]

    def __init__(self, fig, lengde: float, innstillinger: dict, bytt_modus):
        # Plotting: numerisk løsning i cellesentrene, forlenget ut til endeflatene
        # (se `_med_ender`), og Fourier-rekka på et fint rutenett fra 0 til L.
        rutenett = Rutenett.stav(lengde, self.ANTALL_CELLER)
        self.x_numerisk = np.concatenate([[0.0], rutenett.sentre[0], [lengde]])
        self.x_fourier = np.linspace(0.0, lengde, 1001)
        super().__init__(fig, rutenett, innstillinger, bytt_modus)

    def _startbetingelse(self) -> np.ndarray:
        return varmt_segment(self.rutenett, self.L / 2, self.bredde, self.T_omgivelse, self.T_kilde)

    def _rand(self):
        return ((self.T_venstre, self.T_hoyre),) if self.fast_rand else None

    def _bygg_randfelt(self, y: float) -> list:
        return [
            self._randtemperaturfelt(0.34, y, "Venstre [°C] ", "T_venstre"),
            self._randtemperaturfelt(0.53, y, "Høyre [°C] ", "T_hoyre"),
        ]

    def _med_ender(self, u: np.ndarray) -> np.ndarray:
        """Celleverdiene forlenget ut til endeflatene: til den faste temperaturen
        når enden holdes på den, ellers flatt (ingen varmestrøm over en isolert ende)."""
        if self.fast_rand:
            return np.concatenate([[self.T_venstre], u, [self.T_hoyre]])
        return np.concatenate([[u[0]], u, [u[-1]]])

    def _bygg_figur(self) -> None:
        L_cm = self.L * 100
        self.ax_temp = self.fig.add_axes((0.07, 0.445, 0.91, 0.53))
        self.ax_stav = self.fig.add_axes((0.07, 0.39, 0.91, 0.037), sharex=self.ax_temp)

        (self.numerisk_linje,) = self.ax_temp.plot(
            self.x_numerisk * 100, self._med_ender(self.u0), color="tab:red", lw=2, label="Numerisk"
        )
        (self.fourier_linje,) = self.ax_temp.plot(
            self.x_fourier * 100, np.zeros_like(self.x_fourier), color="tab:blue", lw=1.5, ls="--", label="Fourier"
        )
        (self.likevekt_linje,) = self.ax_temp.plot([0, L_cm], [0, 0], color="gray", lw=1, ls=":", label="Likevekt")
        self.ax_temp.legend(loc="upper right", fontsize=8, framealpha=0.8)
        # Avlesningspunktet: en loddrett strek, og en prikk som følger den numeriske kurven
        self.les_strek = self.ax_temp.axvline(self.les_x * 100, color="black", lw=0.8, alpha=0.35)
        (self.les_prikk,) = self.ax_temp.plot([self.les_x * 100], [0.0], "o", color="black", ms=4.5)
        self.likevekt_tekst = self.ax_temp.text(
            0.008, 0.975, "", transform=self.ax_temp.transAxes, va="top", fontsize=9, color="dimgray"
        )
        self.ax_temp.set_xlim(0, L_cm)
        self.ax_temp.set_ylabel("Temperatur [°C]")
        self.ax_temp.grid(alpha=0.3)
        self.ax_temp.tick_params(labelbottom=False)

        self.stav_bilde = self.ax_stav.imshow(
            self.u0[np.newaxis, :], aspect="auto", extent=(0, L_cm, 0, 1), cmap="inferno"
        )
        self.ax_stav.set_yticks([])
        self.ax_stav.set_xlabel("x [cm]", labelpad=1)

    def _ny_start(self) -> None:
        if self.fast_rand:
            # Faste ender: sinusrekke, og likevekten er en rett linje mellom endetemperaturene
            self.fourier = FourierSinusrekke(
                self.u0, self.rutenett, self.alfa, self.T_venstre, self.T_hoyre, self.MAKS_FOURIER_LEDD
            )
            self.likevekt_linje.set_ydata([self.T_venstre, self.T_hoyre])
        else:
            # Isolerte ender: cosinusrekke, og likevekten er snittemperaturen
            self.fourier = FourierCosinusrekke(self.u0, self.rutenett, self.alfa, self.MAKS_FOURIER_LEDD)
            self.likevekt_linje.set_ydata([self.fourier.a[0]] * 2)

    def _likevekt(self) -> tuple:
        if not self.fast_rand:
            tekst = f"Likevekt: {self.fourier.a[0]:.1f} °C"
        elif self.T_venstre == self.T_hoyre:
            tekst = f"Likevekt: {self.T_venstre:g} °C"
        else:
            tekst = f"Likevekt: rett linje fra {self.T_venstre:g} til {self.T_hoyre:g} °C"
        return self.likevekt_tekst, tekst

    def _juster_temperaturskala(self) -> None:
        """Temperaturaksen dekker minst T_AKSE, og utvides hvis noe er varmere/kaldere."""
        verdier = [self.u0.min(), self.u0.max(), self.T_kilde]
        if self.fast_rand:
            verdier += [self.T_venstre, self.T_hoyre]
        T_min = min(self.T_AKSE[0], *verdier)
        T_maks = max(self.T_AKSE[1], *verdier)
        margin = 0.03 * (T_maks - T_min)
        grenser = (T_min - margin, T_maks + margin)
        if tuple(self.ax_temp.get_ylim()) != grenser:
            self.ax_temp.set_ylim(grenser)
            self.stav_bilde.set_clim(T_min, T_maks)
            self.fig.canvas.draw_idle()  # aksen er en del av bakgrunnen

    def _data_artister(self) -> list:
        return [self.stav_bilde, self.likevekt_linje, self.numerisk_linje, self.fourier_linje, self.les_prikk]

    def _oppdater_data(self) -> None:
        ledd = int(self.ledd_bryter.val)
        self.numerisk_linje.set_ydata(self._med_ender(self.u))
        self.fourier_linje.set_ydata(self.fourier.evaluer(self.x_fourier, self.t, ledd))
        self.stav_bilde.set_data(self.u[np.newaxis, :])
        self.les_prikk.set_data([self.les_x * 100], [self._avlesning()[0]])

    def _avlesning(self) -> tuple[float, float]:
        # Numerisk: rett linje mellom cellesentrene, slik kurven er tegnet
        numerisk = np.interp(self.les_x, self.x_numerisk, self._med_ender(self.u))
        fourier = self.fourier.evaluer(np.array([self.les_x]), self.t, int(self.ledd_bryter.val))[0]
        return float(numerisk), float(fourier)

    def _flytt_merke(self) -> None:
        self.les_strek.set_xdata([self.les_x * 100] * 2)

    # --- Redigering med mus --------------------------------------------------

    def _mus_til_data(self, hendelse) -> tuple[float, float]:
        """Musposisjon i (x [m], T [°C]), også når musa er utenfor grafen."""
        x_cm, T = self.ax_temp.transData.inverted().transform((hendelse.x, hendelse.y))
        T_min, T_maks = self.ax_temp.get_ylim()
        return float(np.clip(x_cm / 100, 0.0, self.L)), float(np.clip(T, T_min, T_maks))

    def _trykk(self, hendelse) -> None:
        if hendelse.inaxes not in (self.ax_temp, self.ax_stav) or hendelse.button not in (VENSTRE, HOYRE):
            return
        self._start_redigering()
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
            self.simuler()

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
        else:
            return  # på selve staven kan man bare klikke, ikke tegne
        r["forrige"] = (x, T)
        self.simuler()

    def _slipp(self, hendelse) -> None:
        r = self._redigering
        if r is None:
            return
        if r["knapp"] == VENSTRE and not r["dratt"]:
            x, _ = r["forrige"]
            self.u0[np.abs(self.x - x) <= self.bredde / 2] = self.T_kilde
            self.u0[np.argmin(np.abs(self.x - x))] = self.T_kilde  # minst én celle
            self.simuler()
        super()._slipp(hendelse)

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


class Platevisning(_Visning):
    """Temperaturen i en kvadratisk plate: numerisk løsning og Fourier-rekke
    som to fargekart ved siden av hverandre."""

    BYTT_TIL = "Bytt til 1D"
    ANTALL_CELLER = 100  # per retning
    ANTALL_BILDER = 200
    MAKS_FOURIER_LEDD = 60  # per retning

    def __init__(self, fig, lengde: float, innstillinger: dict, bytt_modus):
        rutenett = Rutenett((lengde, lengde), (self.ANTALL_CELLER, self.ANTALL_CELLER))
        self.y = rutenett.sentre[1]
        self.forsterket = False
        # Cellesentrene pluss kantene, til avlesning mellom cellene (se `_avlesning`)
        self.x_med_kanter = np.concatenate([[0.0], rutenett.sentre[0], [lengde]])
        self.y_med_kanter = np.concatenate([[0.0], rutenett.sentre[1], [lengde]])
        super().__init__(fig, rutenett, innstillinger, bytt_modus)

    def _rand(self):
        return ((self.T_kant, self.T_kant),) * 2 if self.fast_rand else None

    def _bygg_randfelt(self, y: float) -> list:
        return [self._randtemperaturfelt(0.34, y, "Kant [°C] ", "T_kant")]

    def _bygg_avlesningsfelt(self, y: float) -> list:
        return [*super()._bygg_avlesningsfelt(y), self._posisjonsfelt(0.33, y, "y [cm] ", "les_y")]

    def _startbetingelse(self) -> np.ndarray:
        midten = (self.L / 2, self.L / 2)
        return varmt_segment(self.rutenett, midten, self.bredde, self.T_omgivelse, self.T_kilde)

    def _bygg_figur(self) -> None:
        L_cm = self.L * 100
        self.ax_numerisk = self.fig.add_axes((0.06, 0.38, 0.415, 0.553))
        self.ax_fourier = self.fig.add_axes(
            (0.505, 0.38, 0.415, 0.553), sharex=self.ax_numerisk, sharey=self.ax_numerisk
        )
        # Begge fargekartene deler fargeskala. Temperaturen lagres med x langs
        # første akse, mens imshow vil ha y som rader; derfor .T
        self.fargeskala = Normalize(0.0, 1.0)
        valg = dict(origin="lower", extent=(0, L_cm, 0, L_cm), cmap="inferno", norm=self.fargeskala,
                    interpolation="nearest")
        self.numerisk_bilde = self.ax_numerisk.imshow(self.u0.T, **valg)
        self.fourier_bilde = self.ax_fourier.imshow(self.u0.T, **valg)
        self.ax_numerisk.set_xlabel("x [cm]", labelpad=1)
        self.ax_fourier.set_xlabel("x [cm]", labelpad=1)
        self.ax_numerisk.set_ylabel("y [cm]")
        self.ax_fourier.tick_params(labelleft=False)
        self.fargesoyle = self.fig.colorbar(
            self.numerisk_bilde, cax=self.fig.add_axes((0.935, 0.38, 0.015, 0.553))
        )
        self.fargesoyle.ax.set_title("°C", fontsize=9)
        self.ax_numerisk.set_title("Numerisk", fontsize=10)
        self.ax_fourier.set_title("Fourier-rekke", fontsize=10)
        # Avlesningspunktet: et kryss i begge fargekartene
        merke = dict(marker="+", color="white", ms=11, mew=1.6, ls="none")
        (self.les_merke_numerisk,) = self.ax_numerisk.plot([self.les_x * 100], [self.les_y * 100], **merke)
        (self.les_merke_fourier,) = self.ax_fourier.plot([self.les_x * 100], [self.les_y * 100], **merke)

        self.farge_knapp = self._knapp(
            self._akse(0.72, self.RADER[2], 0.14), "Farger: lineære", self._bytt_fargeskala
        )

    def _ny_start(self) -> None:
        if self.fast_rand:
            # Faste kanter: dobbel sinusrekke, og likevekten er kanttemperaturen
            self.fourier = FourierSinusrekke2D(
                self.u0, self.rutenett, self.alfa, self.T_kant, self.MAKS_FOURIER_LEDD
            )
        else:
            # Isolerte kanter: dobbel cosinusrekke, og likevekten er snittemperaturen
            self.fourier = FourierCosinusrekke2D(self.u0, self.rutenett, self.alfa, self.MAKS_FOURIER_LEDD)

    def _likevekt(self) -> tuple:
        likevekt = self.T_kant if self.fast_rand else self.fourier.a[0, 0]
        return self.ax_numerisk.title, f"Numerisk   (likevekt: {likevekt:.1f} °C)"

    def _bytt_fargeskala(self, _=None) -> None:
        """Lineære farger: fargen følger temperaturen jevnt. Forsterkede: de lave
        temperaturene er strukket ut, så også svak varme synes godt."""
        self.forsterket = not self.forsterket
        self.farge_knapp.label.set_text("Farger: forsterket" if self.forsterket else "Farger: lineære")
        lav, hoy = self.fargeskala.vmin, self.fargeskala.vmax
        self.fargeskala = PowerNorm(0.4, lav, hoy) if self.forsterket else Normalize(lav, hoy)
        self.numerisk_bilde.set_norm(self.fargeskala)
        self.fourier_bilde.set_norm(self.fargeskala)
        self.fargesoyle.update_normal(self.numerisk_bilde)
        self.fig.canvas.draw_idle()  # fargesøylen er en del av bakgrunnen

    def _juster_temperaturskala(self) -> None:
        """Fargeskalaen går fra grunntemperaturen til den varmeste kilden (eller kanten)."""
        verdier = [self.T_omgivelse, self.u0.min(), self.u0.max(), self.T_kilde]
        if self.fast_rand:
            verdier.append(self.T_kant)
        lav, hoy = float(min(verdier)), float(max(verdier))
        if hoy - lav < 1e-9:
            hoy = lav + 1.0
        if (self.fargeskala.vmin, self.fargeskala.vmax) != (lav, hoy):
            self.fargeskala.vmin, self.fargeskala.vmax = lav, hoy
            self.fig.canvas.draw_idle()  # fargesøylen er en del av bakgrunnen

    def _data_artister(self) -> list:
        return [self.numerisk_bilde, self.fourier_bilde, self.les_merke_numerisk, self.les_merke_fourier]

    def _oppdater_data(self) -> None:
        ledd = int(self.ledd_bryter.val)
        self.numerisk_bilde.set_data(self.u.T)
        self.fourier_bilde.set_data(self.fourier.evaluer(self.x, self.y, self.t, ledd).T)

    def _avlesning(self) -> tuple[float, float]:
        # Numerisk: jevn overgang (bilineær interpolasjon) mellom cellesentrene. Ytterst
        # brukes kanttemperaturen når den er fast, ellers verdien i cellen innenfor.
        kant = dict(mode="constant", constant_values=self.T_kant) if self.fast_rand else dict(mode="edge")
        celler = RegularGridInterpolator((self.x_med_kanter, self.y_med_kanter), np.pad(self.u, 1, **kant))
        numerisk = celler((self.les_x, self.les_y))
        fourier = self.fourier.evaluer(
            np.array([self.les_x]), np.array([self.les_y]), self.t, int(self.ledd_bryter.val)
        )[0, 0]
        return float(numerisk), float(fourier)

    def _flytt_merke(self) -> None:
        for merke in (self.les_merke_numerisk, self.les_merke_fourier):
            merke.set_data([self.les_x * 100], [self.les_y * 100])

    # --- Redigering med mus --------------------------------------------------

    def _mus_til_data(self, hendelse, akse) -> tuple[float, float]:
        """Musposisjon i (x, y) [m] på platen, også når musa er utenfor fargekartet."""
        x_cm, y_cm = akse.transData.inverted().transform((hendelse.x, hendelse.y))
        return float(np.clip(x_cm / 100, 0.0, self.L)), float(np.clip(y_cm / 100, 0.0, self.L))

    def _trykk(self, hendelse) -> None:
        if hendelse.inaxes not in (self.ax_numerisk, self.ax_fourier) or hendelse.button not in (VENSTRE, HOYRE):
            return
        self._start_redigering()
        posisjon = self._mus_til_data(hendelse, hendelse.inaxes)
        self._redigering = {"knapp": hendelse.button, "akse": hendelse.inaxes, "forrige": posisjon}
        self._mal(posisjon, posisjon, hendelse.button)
        self.simuler()

    def _bevegelse(self, hendelse) -> None:
        r = self._redigering
        if r is None:
            return
        posisjon = self._mus_til_data(hendelse, r["akse"])
        self._mal(r["forrige"], posisjon, r["knapp"])
        r["forrige"] = posisjon
        self.simuler()

    def _mal(self, fra: tuple[float, float], til: tuple[float, float], knapp: int) -> None:
        """Maler et kvadrat med side `bredde` langs strekningen fra–til: varmekilde
        med venstre knapp, grunntemperatur med høyre."""
        verdi = self.T_kilde if knapp == VENSTRE else self.T_omgivelse
        halv_celle = self.rutenett.dx[0] / 2
        antall = int(np.hypot(til[0] - fra[0], til[1] - fra[1]) / halv_celle) + 1
        for s in np.linspace(0.0, 1.0, antall):
            x = fra[0] + s * (til[0] - fra[0])
            y = fra[1] + s * (til[1] - fra[1])
            i = np.abs(self.x - x) <= self.bredde / 2
            j = np.abs(self.y - y) <= self.bredde / 2
            i[np.argmin(np.abs(self.x - x))] = True  # minst én celle
            j[np.argmin(np.abs(self.y - y))] = True
            self.u0[np.ix_(i, j)] = verdi


class Varmesimulering:
    """Simuleringsvinduet. Starter som stav (1D) eller plate (2D), og kan
    byttes mellom dem med en knapp i vinduet."""

    VISNINGER = {"1D": Stavvisning, "2D": Platevisning}

    def __init__(
        self,
        materiale: Materiale,
        modus: str = "1D",
        lengde: float = 1.0,  # stavens lengde og platens sidekant [m]
        T_omgivelse: float = 0.0,  # grunntemperatur [°C]
        T_kilde: float = 300.0,  # temperatur for varmekilder [°C]
        bredde: float = 0.10,  # bredde for varmekilder [m]
        fast_rand: bool = False,  # False: isolerte kanter, True: fast temperatur på kantene
        T_rand: float = 0.0,  # temperaturen kantene holdes på når fast_rand er valgt [°C]
        fourier_ledd: int = 20,
    ):
        self.lengde = lengde
        self.fig = plt.figure(figsize=(10, 7.5))
        lerret = self.fig.canvas
        if lerret.manager is not None:
            lerret.manager.set_window_title("Varmesimulering")
            # Slå av matplotlibs hurtigtaster (f.eks. «l» = logaritmisk akse)
            lerret.mpl_disconnect(lerret.manager.key_press_handler_id)
        for egenskap in ("header_visible", "footer_visible", "toolbar_visible", "resizable"):
            if hasattr(lerret, egenskap):  # ipympl i Jupyter
                setattr(lerret, egenskap, False)

        self.visning = None
        self.sett_modus(modus, {
            "k": materiale.varmeledningsevne, "rho": materiale.tetthet, "c_p": materiale.spesifikk_varme,
            "T_omgivelse": T_omgivelse, "T_kilde": T_kilde, "bredde": bredde,
            "fast_rand": fast_rand, "T_venstre": T_rand, "T_hoyre": T_rand, "T_kant": T_rand,
            "les_x": lengde / 2, "les_y": lengde / 2,
            "fourier_ledd": fourier_ledd,
        })

    def sett_modus(self, modus: str, innstillinger: dict | None = None) -> None:
        """Bygger vinduet på nytt som stav ("1D") eller plate ("2D")."""
        if self.visning is not None:
            innstillinger = innstillinger or self.visning.innstillinger()
            self.visning.lukk()
        self.fig.clf()
        self.modus = modus
        self.visning = self.VISNINGER[modus](self.fig, self.lengde, innstillinger, self.bytt_modus)
        self.fig.canvas.draw_idle()

    def bytt_modus(self) -> None:
        """Bytter mellom stav og plate. Selve byttet skjer rett etterpå, slik at
        knappen som ble trykket rekker å bli ferdig før den fjernes."""
        ny = "2D" if self.modus == "1D" else "1D"
        self._bytte_timer = self.fig.canvas.new_timer(interval=1)
        self._bytte_timer.single_shot = True
        self._bytte_timer.add_callback(lambda: self.sett_modus(ny))
        self._bytte_timer.start()
