"""Kopierer et bilde til utklippstavlen, slik at det kan limes rett inn i et dokument.

På Windows brukes systemets egne funksjoner via ctypes, så det trengs ingen
ekstra pakker. På andre systemer lagres bildet i stedet som en PNG-fil.
"""

import struct
import sys
from pathlib import Path

import numpy as np


def kopier_bilde(rgba: np.ndarray) -> str:
    """Kopierer bildet (høyde × bredde × 4, øverste rad først) til utklippstavlen.
    Returnerer en kort melding om hva som ble gjort."""
    if sys.platform == "win32":
        _til_windows_utklippstavle(rgba)
        return "Kopiert!"
    import matplotlib.pyplot as plt

    fil = Path("varmesimulering_graf.png").resolve()
    plt.imsave(fil, rgba)
    print(f"Grafen er lagret som {fil}")
    return "Lagret som fil"


def _til_windows_utklippstavle(rgba: np.ndarray) -> None:
    """Legger bildet på utklippstavlen i formatet CF_DIB, som Word, PowerPoint,
    nettlesere og de fleste andre programmer kan lime inn."""
    import ctypes
    from ctypes import wintypes

    # Bildet som punktgrafikk (DIB): 24 bit per piksel i rekkefølgen blå-grønn-rød,
    # nederste rad først, og hver rad fylt opp til et helt antall 4 byte.
    hoyde, bredde = rgba.shape[:2]
    radlengde = (bredde * 3 + 3) // 4 * 4
    rader = np.zeros((hoyde, radlengde), dtype=np.uint8)
    rader[:, : bredde * 3] = rgba[::-1, :, 2::-1].reshape(hoyde, bredde * 3)
    hode = struct.pack("<IiiHHIIiiII", 40, bredde, hoyde, 1, 24, 0, rader.size, 3780, 3780, 0, 0)
    data = hode + rader.tobytes()

    kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
    user32 = ctypes.WinDLL("user32", use_last_error=True)
    kernel32.GlobalAlloc.argtypes = [wintypes.UINT, ctypes.c_size_t]
    kernel32.GlobalAlloc.restype = wintypes.HGLOBAL
    kernel32.GlobalLock.argtypes = [wintypes.HGLOBAL]
    kernel32.GlobalLock.restype = wintypes.LPVOID
    kernel32.GlobalUnlock.argtypes = [wintypes.HGLOBAL]
    kernel32.GlobalFree.argtypes = [wintypes.HGLOBAL]
    user32.OpenClipboard.argtypes = [wintypes.HWND]
    user32.SetClipboardData.argtypes = [wintypes.UINT, wintypes.HANDLE]
    user32.SetClipboardData.restype = wintypes.HANDLE
    GMEM_MOVEABLE, CF_DIB = 0x0002, 8

    minne = kernel32.GlobalAlloc(GMEM_MOVEABLE, len(data))
    if not minne:
        raise OSError("Fikk ikke satt av minne til bildet")
    ctypes.memmove(kernel32.GlobalLock(minne), data, len(data))
    kernel32.GlobalUnlock(minne)

    if not user32.OpenClipboard(None):
        kernel32.GlobalFree(minne)
        raise OSError("Fikk ikke åpnet utklippstavlen")
    try:
        user32.EmptyClipboard()
        if not user32.SetClipboardData(CF_DIB, minne):  # ved suksess eier systemet minnet
            kernel32.GlobalFree(minne)
            raise OSError("Fikk ikke lagt bildet på utklippstavlen")
    finally:
        user32.CloseClipboard()
