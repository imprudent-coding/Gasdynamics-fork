"""Prova del 01/03/1970, esercizio 1: diffrazione di un urto piano (M = 4) attorno a uno spigolo di 105°
(stessa traccia di B70-E01, 13/02/2023).

Riproduce tutti i numeri e le figure dello svolgimento.
Uso:  python calcoli.py   (dalla cartella dello svolgimento)

Convenzioni (gasdyn_guida.md §1): angoli in gradi, θ antiorario dall'asse x (verso dell'urto incidente e della
corrente 2, parallela alla parete superiore dello spigolo). La parete inclinata scende da C a −105°.
Stati: 1 gas in quiete davanti all'urto; 2 dietro l'urto piano; 3 dopo il ventaglio centrato in C (lungo il primo
tratto di SL, θ = −44°); 4 dopo il ventaglio centrato in F (θ = −61°); 5 dietro il piede B dell'urto diffratto;
6 dietro RC; ∞ campo asintotico allo spigolo.
Gas: aria, γ = 1,4, R = 287 J/(kg K) (la traccia dice "air").

Campo autosimile: ogni punto delle strutture si muove con velocità costante, proporzionale alla sua distanza da C
(le misure sull'immagine danno velocità, scala W_A / x_A).

Fuori modulo: misure sulla figura (coordinate in pixel lette su ingrandimenti con griglia, angolo della testa del
ventaglio dal profilo di luminosità), velocità normale di un punto di un'onda autosimile (W_n = h/t, h distanza di C
dalla tangente), cerchio acustico dello spigolo e punto D.
"""
import os
import sys
import math

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image
from scipy import ndimage

QUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(QUI, "..", ".."))
import gasdyn as gd  # noqa: E402

FIG = os.path.join(QUI, "figure")
os.makedirs(FIG, exist_ok=True)
G = gd.Gas(1.4, 287.0)
MS = 4.0
P1, T1 = 1e5, 293.0
RAPP_AB = 2.5
TH1, TH2 = 44.0, 17.0
SPIGOLO = 105.0


def riga(titolo):
    print(f"\n=== {titolo} ===")


# =============================================================================
# 1. STATO 1 E URTO PIANO (F §7: urto mobile u+a nel gas in quiete)
# =============================================================================
s1 = gd.State.da_pT(G, P1, T1, u=0.0, nome="1")
inc = gd.urto_mobile(s1, "u+a", Ms=MS, nome_valle="2")
s2 = inc.valle
WA = inc.w
riga("Stato 1 e urto piano")
print(f"  a1 = {s1.a:.1f} m/s, rho1 = {s1.rho:.3f} kg/m3")
print(f"  W = Ms a1 = {WA:.1f} m/s")
print(f"  stato 2: p2 = {s2.p:.4e} Pa (p2/p1 = {s2.p / P1:.3f}), rho2/rho1 = {s2.rho / s1.rho:.3f},"
      f" T2 = {s2.T:.1f} K (T2/T1 = {s2.T / T1:.4f})")
print(f"           u2 = {s2.u:.1f} m/s, a2 = {s2.a:.1f} m/s, M2 = u2/a2 = {s2.u / s2.a:.4f}")
print(f"           p02 = {s2.p0:.4e} Pa, T02 = {s2.T0:.1f} K, rho2 = {s2.rho:.3f} kg/m3")
print(f"  verifica R-H: residui = {inc.residui_rh()}")
print(f"  Mr a valle (u2-W)/a2 = {(s2.u - WA) / s2.a:.4f}  (|.|<1: urto subsonico rispetto al gas dietro)")
ds12 = gd.salto_entropia(s1, s2)["ds_R"]
print(f"  Δs/R attraverso l'urto piano = {ds12:.4f}")
# soglia di M_s oltre la quale la corrente dietro un urto in gas in quiete è supersonica
print(f"  limite |u2/a2| per Ms→∞ = {gd.M_valle_limite(G):.4f}")
ms_son = gd.urto_mobile(s1, "u+a", M_valle=1.0).Ms
print(f"  M2 = 1 per Ms = {ms_son:.3f}  ⇒ con Ms = 4 la corrente 2 è supersonica")

# =============================================================================
# 2. PIEDE B (autosimilarità + urto normale alla parete nel gas in quiete, F §7)
# =============================================================================
WB = WA / RAPP_AB
foot = gd.urto_mobile(s1, "u+a", w=WB, nome_valle="5")
s5 = foot.valle
riga("Piede B dell'urto diffratto")
print(f"  W_B = W_A / 2,5 = {WB:.1f} m/s, M_B = W_B/a1 = {foot.Ms:.4f}")
print(f"  stato 5: p5 = {s5.p:.4e} Pa (p5/p1 = {s5.p / P1:.3f}), T5 = {s5.T:.1f} K, u5 = {s5.u:.1f} m/s,"
      f" a5 = {s5.a:.1f} m/s, M5 = u5/a5 = {s5.u / s5.a:.3f}, rho5/rho1 = {s5.rho / s1.rho:.3f}")
ds15 = gd.salto_entropia(s1, s5)["ds_R"]
print(f"  Δs/R al piede = {ds15:.4f}  (contro {ds12:.4f} dell'urto piano)")
print(f"  u5/W_B = {s5.u / WB:.3f}  (rapporto atteso CE/CB se E si muove con il gas a parete)")

# =============================================================================
# 3. VENTAGLI DI PRANDTL–MEYER IN C E IN F (F §11, onde C+: la corrente ruota in senso orario)
# =============================================================================
s2d = gd.State.da_M(G, s2.u / s2.a, 0.0, p=s2.p, T=s2.T, nome="2")
eC = gd.espansione_pm(s2d, TH1, "C+", nome="3")
s3 = eC.valle
eF = gd.espansione_pm(s3, TH2, "C+", nome="4")
s4 = eF.valle
nu2, nu3, nu4 = (gd.prandtl_meyer(s.M, G) for s in (s2d, s3, s4))
riga("Ventagli di Prandtl–Meyer")
print(f"  ν(M2) = {nu2:.3f}°, μ2 = {gd.angolo_mach(s2d.M):.3f}°  (testa del ventaglio in C a θ+μ = {gd.angolo_mach(s2d.M):.2f}°)")
print(f"  stato 3: ν3 = ν2 + θ1 = {nu3:.3f}°, M3 = {s3.M:.4f}, p3 = {s3.p:.4e} Pa (p3/p2 = {s3.p / s2.p:.4f},"
      f" p3/p1 = {s3.p / P1:.3f}), T3 = {s3.T:.1f} K, a3 = {s3.a:.1f} m/s, V3 = {s3.V:.1f} m/s, θ3 = {s3.theta:.1f}°")
print(f"           μ3 = {s3.mu:.3f}°, coda del ventaglio in C a θ3+μ3 = {s3.theta + s3.mu:.2f}°")
print(f"  stato 4: ν4 = ν3 + θ2 = {nu4:.3f}°, M4 = {s4.M:.4f}, p4 = {s4.p:.4e} Pa (p4/p3 = {s4.p / s3.p:.4f},"
      f" p4/p2 = {s4.p / s2.p:.5f}), T4 = {s4.T:.1f} K, a4 = {s4.a:.1f} m/s, V4 = {s4.V:.1f} m/s, θ4 = {s4.theta:.1f}°")
print(f"           μ4 = {s4.mu:.3f}°, testa del ventaglio in F a {s3.theta + s3.mu:.2f}°, coda a θ4+μ4 = {s4.theta + s4.mu:.2f}°")
print(f"  p0 costante: p02 = {s2d.p0:.4e}, p03 = {s3.p0:.4e}, p04 = {s4.p0:.4e}; T0 = {s2d.T0:.1f}, {s3.T0:.1f}, {s4.T0:.1f}")
print(f"  V costante? no: V2 = {s2d.V:.1f}, V3 = {s3.V:.1f}, V4 = {s4.V:.1f} (V_max = {s2d.V_max:.1f})")
print(f"  ν_max = {gd.nu_max(G):.2f}°: deviazione massima possibile da 2 = {gd.nu_max(G) - nu2:.2f}° (> 105°)")
print(f"  deviazione totale di SL = {TH1 + TH2:.0f}° < {SPIGOLO:.0f}° (parete): SL resta staccata dalla parete")
att = gd.espansione_pm(s2d, SPIGOLO, "C+").valle
print(f"  ipotetica espansione aderente di 105°: M = {att.M:.2f}, p = {att.p:.3g} Pa, T = {att.T:.1f} K")

# =============================================================================
# 4. CAMPO ASINTOTICO: espansione da 2 fino a p1 (F §11 + isentropiche F §2)
# =============================================================================
Minf = gd.mach_da_rapporto(s2d.p0 / P1, G, "p0/p")
nuinf = gd.prandtl_meyer(Minf, G)
th_inf = nuinf - nu2
riga("Campo asintotico allo spigolo")
print(f"  p02/p1 = {s2d.p0 / P1:.3f} ⇒ M∞ = {Minf:.4f}, ν∞ = {nuinf:.3f}° ⇒ θ∞ = ν∞ − ν2 = {th_inf:.2f}°")

# =============================================================================
# 5. MISURE SULLA FIGURA (coordinate in pixel nel ritaglio della figura di E01.png, origine (200, 1300))
# =============================================================================
PX = dict(C=(584, 839), A=(1912, 51), B=(450, 1357), E=(523, 1080), F=(785, 1045), SLend=(905, 1248),
          Dvis=(1910, 420))
RC_PTS = np.array([(915, 1222), (960, 1178), (1050, 1127), (1150, 1071), (1215, 1010), (1245, 950)], float)
DS_PTS = np.array([(1912, 51), (1913, 360), (1885, 700), (1830, 800), (1760, 950), (1700, 1050), (1600, 1150),
                   (1450, 1260), (1300, 1385), (1180, 1455), (880, 1510), (600, 1450), (450, 1357)], float)
C = np.array(PX["C"], float)


def vec(p):
    """Vettore da C in un sistema con y verso l'alto (pixel)."""
    return np.array([p[0] - C[0], C[1] - p[1]], float)


xA = vec(PX["A"])[0]
CB = np.hypot(*vec(PX["B"]))
scala = WA / xA                       # (m/s) per pixel: velocità di un punto = distanza da C × scala
riga("Misure sulla figura")
print(f"  x_A − x_C = {xA:.0f} px, CB = {CB:.0f} px, rapporto = {xA / CB:.3f} (dato 2,5)")
ang_parete = math.degrees(math.atan2(*vec(PX["B"])[::-1]))
print(f"  direzione CB = {ang_parete:.1f}° (dato −105°)")
CF = vec(PX["F"]); FS = vec(PX["SLend"]) - CF
print(f"  SL: tratto CF a {math.degrees(math.atan2(CF[1], CF[0])):.1f}°, tratto dopo F a"
      f" {math.degrees(math.atan2(FS[1], FS[0])):.1f}° (dati −44° e −61°)")
CE = np.hypot(*vec(PX["E"]))
print(f"  CE/CB = {CE / CB:.3f}  ⇒ W_E ≈ {CE / CB * WB:.0f} m/s (u5 = {s5.u:.0f} m/s)")

# angolo della testa del ventaglio in C dal profilo di luminosità lungo i raggi da C
img = np.asarray(Image.open(os.path.join(FIG, "originale.png")).convert("L"), float)[1300:2951, 200:2400]
liscia = ndimage.gaussian_filter(img, 4)
angs = np.arange(30, 50.01, 0.5)
lum = []
for ang in angs:
    t = math.radians(ang)
    rs = np.arange(150, 650, 2.0)
    xs = C[0] + rs * math.cos(t); ys = C[1] - rs * math.sin(t)
    lum.append(ndimage.map_coordinates(liscia, [ys, xs]).mean())
lum = np.array(lum)
meta = 0.5 * (lum.min() + lum.max())
ang_testa = float(np.interp(meta, lum, angs))
print(f"  testa del ventaglio (bordo della zona scura da C): {ang_testa:.1f}° (μ2 = {s2d.mu:.1f}°)")

# punto D e cerchio acustico dello spigolo
u2, a2 = s2.u, s2.a
yD = math.sqrt(a2**2 - (WA - u2)**2)          # per unità di tempo
rT = math.sqrt(u2**2 - a2**2)                  # distanza di C dal punto di tangenza
print(f"  cerchio dello spigolo: centro u2 t = {u2:.0f} t, raggio a2 t = {a2:.0f} t")
print(f"  tangenza della testa rettilinea: a distanza {rT:.0f} t da C, ({rT * math.cos(math.radians(s2d.mu)):.0f},"
      f" {rT * math.sin(math.radians(s2d.mu)):.0f}) t")
print(f"  D previsto: y_D/x_A = sqrt(a2² − (W − u2)²)/W = {yD / WA:.3f};"
      f" D visibile sulla figura: y/x_A = {vec(PX['Dvis'])[1] / xA:.3f}")
print(f"  massima altezza del cerchio = a2 t ⇒ y/x_A = {a2 / WA:.3f}; centro a x/x_A = {u2 / WA:.3f}")


def normale_tangente(pts, i):
    """Tangente (differenze centrate) e normale verso valle (lontano da C) nel punto i, y verso l'alto."""
    p = np.array([vec(q) for q in pts])
    j0, j1 = max(i - 1, 0), min(i + 1, len(p) - 1)
    t = p[j1] - p[j0]; t /= np.hypot(*t)
    n = np.array([t[1], -t[0]])
    if n @ p[i] < 0:
        n = -n
    return p[i], n


riga("Urto di ricompressione RC")
rc = []
for i, (stato, nome) in [(0, (s4, "4")), (2, (s3, "3"))]:
    P, n = normale_tangente(RC_PTS, i)
    h = P @ n
    Wn = h * scala
    th_n = math.degrees(math.atan2(n[1], n[0]))
    phi = stato.theta - th_n
    Vn = stato.V * math.cos(math.radians(phi))
    Mrel = (Vn - Wn) / stato.a
    un = gd.urto_normale(Mrel, G)
    p6 = stato.p * un["p2_p1"]
    rc.append(dict(nome=nome, h=h, Wn=Wn, th_n=th_n, phi=phi, Mrel=Mrel, p6=p6, un=un, stato=stato))
    print(f"  tratto davanti allo stato {nome}: punto {RC_PTS[i]}, normale a {th_n:.1f}°,"
          f" h = {h:.0f} px (h/x_A = {h / xA:.3f}) ⇒ W_n = {Wn:.0f} m/s; angolo flusso–normale = {phi:.1f}°")
    print(f"      V·n = {Vn:.0f} m/s, M_rel = (V·n − W_n)/a = {Mrel:.3f}, p6/p{nome} = {un['p2_p1']:.3f},"
          f" p6 = {p6:.3e} Pa, M_rel valle = {un['M2']:.3f}, p0 valle/monte = {un['p02_p01']:.3f}")
    # stima semplificata (RC normale al flusso, W_n ≈ W_B)
    Ms_s = (stato.V - WB) / stato.a
    print(f"      con RC normale al flusso e W_n ≈ W_B: M_rel = {Ms_s:.3f}, p6 = {stato.p * gd.urto_normale(Ms_s, G)['p2_p1']:.3e} Pa")

riga("Urto diffratto DS: Mach locale (W_n = h/t)")
ds_loc = []
for i in range(len(DS_PTS)):
    P, n = normale_tangente(DS_PTS, i)
    h = P @ n
    Mloc = h * scala / s1.a
    th_n = math.degrees(math.atan2(n[1], n[0]))
    pr = gd.urto_mobile(s1, "u+a", Ms=Mloc).p2_p1 if Mloc > 1 else float("nan")
    ds_loc.append((th_n, Mloc, pr))
    print(f"  punto {DS_PTS[i]}: normale a {th_n:6.1f}°, M_locale ≈ {Mloc:.2f}, p/p1 dietro ≈ {pr:.2f}")



# =============================================================================
# 6. FIGURE
# =============================================================================
def px(v):
    """Da coordinate relative a C (y verso l'alto, pixel) a coordinate immagine."""
    v = np.asarray(v, float)
    return C[0] + v[..., 0], C[1] - v[..., 1]


def ray(ang, r0, r1):
    t = math.radians(ang)
    return px(np.array([[r0 * math.cos(t), r0 * math.sin(t)], [r1 * math.cos(t), r1 * math.sin(t)]]))


fig_img = Image.open(os.path.join(FIG, "originale.png")).convert("L").crop((200, 1300, 2400, 2951))
k = 1.0 / scala                    # pixel per (m/s · t): posizione = velocità × k
xc_c, rc_c = u2 * k, a2 * k        # cerchio acustico dello spigolo (pixel)
tang = rT * k * np.array([math.cos(math.radians(s2d.mu)), math.sin(math.radians(s2d.mu))])
yDp = yD * k


def base(ax):
    ax.imshow(fig_img, cmap="gray", vmin=0, vmax=255)
    ax.set_xlim(0, 2200); ax.set_ylim(1651, 0); ax.axis("off")


# --- 6a. figura annotata ------------------------------------------------------
fig, ax = plt.subplots(figsize=(9.6, 7.3))
base(ax)
ax.plot(DS_PTS[:, 0], DS_PTS[:, 1], "-", color="red", lw=1.6, alpha=.75)
ax.plot(RC_PTS[:, 0], RC_PTS[:, 1], "-", color="orange", lw=2.2)
ax.plot(*px(np.array([[0, 0], vec(PX["F"]), vec(PX["SLend"])])), "-", color="lime", lw=2)
# testa del ventaglio in C: tratto rettilineo fino alla tangenza, poi cerchio acustico fino a D
xs, ys = px(np.array([[0, 0], tang]))
ax.plot(xs, ys, "-", color="cyan", lw=1.8)
phT = math.atan2(tang[1], tang[0] - xc_c)
phD = math.atan2(yDp, xA - xc_c)
ph = np.linspace(phT, phD, 60)
ax.plot(*px(np.c_[xc_c + rc_c * np.cos(ph), rc_c * np.sin(ph)]), "-", color="cyan", lw=1.8)
ph2 = np.linspace(phD, -1.25, 80)
ax.plot(*px(np.c_[xc_c + rc_c * np.cos(ph2), rc_c * np.sin(ph2)]), ":", color="cyan", lw=1.2)
# code dei ventagli (C: −27,1°; F: testa −27,1°, coda −48,9°)
ax.plot(*ray(s3.theta + s3.mu, 0, 620), "--", color="cyan", lw=1)
Fv = vec(PX["F"])
for ang in (s3.theta + s3.mu, s4.theta + s4.mu):
    t = math.radians(ang)
    ax.plot(*px(np.array([Fv, Fv + 330 * np.array([math.cos(t), math.sin(t)])])), "--", color="deepskyblue", lw=1)
# CD (tracciato sulla figura, poi fino a E sotto la bolla)
CD_PTS = np.array([(1700, 650), (1640, 860), (1560, 1000), (1450, 1110), (1300, 1235), (1150, 1285),
                   (1000, 1335), (880, 1330), (760, 1230), (640, 1140), (523, 1080)], float)
ax.plot(CD_PTS[:, 0], CD_PTS[:, 1], "--", color="magenta", lw=1.6)
# punti
for nome, (x, y), dx, dy in [("C", PX["C"], -60, -25), ("A", PX["A"], 25, 45), ("B", PX["B"], -70, 10),
                             ("E", PX["E"], -70, 10), ("F", PX["F"], 25, -10)]:
    ax.plot(x, y, "o", ms=4, color="yellow")
    ax.text(x + dx, y + dy, nome, color="yellow", fontsize=11, weight="bold")
xD, yDpx = px(np.array([xA, yDp]))
ax.plot(xD, yDpx, "s", ms=6, mfc="none", mec="cyan", mew=1.5)
ax.text(xD - 120, yDpx + 10, "D", color="cyan", fontsize=11, weight="bold")
ax.plot(*PX["Dvis"], "s", ms=6, mfc="none", mec="yellow", mew=1.2)
# stati
for txt, (x, y) in [("1", (2080, 1200)), ("1", (1000, 1600)), ("2", (800, 250)), ("2", (300, 650)),
                    ("3", (1010, 1070)), ("4", (880, 1170)), ("5", (560, 1330)), ("6", (1240, 1190)),
                    ("bolla", (590, 1010))]:
    ax.text(x, y, txt, color="white", fontsize=11 if len(txt) == 1 else 9, weight="bold", ha="center",
            va="center", bbox=dict(boxstyle="circle" if len(txt) == 1 else "round", fc="black", ec="white",
                                   alpha=.75))
ax.text(1060, 520, f"{ang_testa:.0f}°", color="cyan", fontsize=10, weight="bold")
ax.text(1560, 560, "cerchio acustico\ndello spigolo", color="cyan", fontsize=8, ha="center")
fig.savefig(os.path.join(FIG, "annotata.png"), dpi=250, bbox_inches="tight", pad_inches=0.02)
plt.close(fig)


# --- 6b. linee di corrente -------------------------------------------------------
def bezier(P0, d0, P3, d3, n=60):
    L = np.hypot(*(P3 - P0))
    P1 = P0 + 0.4 * L * d0
    P2 = P3 - 0.4 * L * d3
    t = np.linspace(0, 1, n)[:, None]
    return (1 - t)**3 * P0 + 3 * (1 - t)**2 * t * P1 + 3 * (1 - t) * t**2 * P2 + t**3 * P3


DSv = np.array([vec(q) for q in DS_PTS])
ds_par = np.r_[0, np.cumsum(np.hypot(*np.diff(DSv, axis=0).T))]
s_fine = np.linspace(0, ds_par[-1], 400)
DS_f = np.c_[np.interp(s_fine, ds_par, DSv[:, 0]), np.interp(s_fine, ds_par, DSv[:, 1])]
tg = np.gradient(DS_f, axis=0); tg /= np.hypot(*tg.T)[:, None]
nDS = np.c_[tg[:, 1], -tg[:, 0]]
nDS[(nDS * DS_f).sum(1) < 0] *= -1
ang_nDS = np.degrees(np.arctan2(nDS[:, 1], nDS[:, 0]))
RCv = np.array([vec(q) for q in RC_PTS] + [vec((1255, 820)), vec((1240, 690))])


def verso_DS(P, th, piega=0.0):
    """Raccordo dal punto P (direzione th) fino a DS, dove la linea arriva normale all'urto.

    Il punto d'arrivo è dove la retta uscente da P con direzione th + piega incontra DS."""
    d = np.array([math.cos(math.radians(th + piega)), math.sin(math.radians(th + piega))])
    hit = interseca(np.array([P, P + 4000 * d]), DS_f)
    j = int(np.argmin(np.hypot(*(DS_f - hit[1]).T)))
    d0 = np.array([math.cos(math.radians(th)), math.sin(math.radians(th))])
    return bezier(P, d0, DS_f[j], nDS[j])


def interseca(path, poly):
    """Primo punto di path che attraversa la spezzata poly (None se nessuno)."""
    for i in range(len(path) - 1):
        p, r = path[i], path[i + 1] - path[i]
        for j in range(len(poly) - 1):
            q, s = poly[j], poly[j + 1] - poly[j]
            den = r[0] * s[1] - r[1] * s[0]
            if abs(den) < 1e-12:
                continue
            t = ((q - p)[0] * s[1] - (q - p)[1] * s[0]) / den
            u = ((q - p)[0] * r[1] - (q - p)[1] * r[0]) / den
            if 0 <= t <= 1 and 0 <= u <= 1:
                return i, p + t * r
    return None


def nel_cerchio(P):
    return (P[0] - xc_c)**2 + P[1]**2 <= rc_c**2


lines = []
for y0 in (12, 40, 110, 230, 380, 540, 650, 740):
    rh = y0 / math.sin(math.radians(s2d.mu))
    xh = rh * math.cos(math.radians(s2d.mu))
    seg = [np.array([-560.0, y0]), np.array([xh, y0])]
    if xh <= tang[0]:                    # entra nel ventaglio stazionario di C (centrato, ρ a r = cost.)
        path = []
        for th in np.linspace(0, -TH1, 200):
            st = gd.espansione_pm(s2d, -th, "C+") if th < 0 else None
            Mx = st.valle.M if st else s2d.M
            rho_a = (st.valle.rho * st.valle.a) if st else (s2d.rho * s2d.a)
            r = rh * s2d.rho * s2d.a / rho_a
            phi = math.radians(th + gd.angolo_mach(Mx))
            path.append([r * math.cos(phi), r * math.sin(phi)])
        path = np.array(path)
        thv = np.linspace(0, -TH1, 200)
        hit = interseca(path, RCv)
        dentro = [i for i, P in enumerate(path) if nel_cerchio(P)]
        if hit is not None and (not dentro or hit[0] <= dentro[0]):
            i_end, Pend = hit
        elif dentro:
            i_end = dentro[0]; Pend = path[i_end]
        else:                            # ventaglio completo: tratto rettilineo nello stato 3 fino a RC
            d3 = np.array([math.cos(math.radians(-TH1)), math.sin(math.radians(-TH1))])
            h3 = interseca(np.array([path[-1], path[-1] + 3000 * d3]), RCv)
            path = np.vstack([path, [h3[1]]])
            thv = np.r_[thv, -TH1]
            i_end, Pend = len(path) - 1, h3[1]
        th_end = float(thv[i_end])
        curva = np.vstack([np.array(seg), path[:i_end + 1], [Pend], verso_DS(Pend, th_end, -18.0 if y0 < 20 else -6.0)])
    else:                               # fuori dal ventaglio: orizzontale fino al cerchio o a DS
        if y0 < yDp:
            xcirc = xc_c - math.sqrt(rc_c**2 - y0**2)
            P = np.array([xcirc, y0])
            curva = np.vstack([seg[0], P, verso_DS(P, 0.0, -8.0)])
        else:
            curva = np.c_[np.linspace(-560.0, xA, 20), np.full(20, float(y0))]
    lines.append(curva)

fig, ax = plt.subplots(figsize=(9.6, 7.3))
base(ax)
for c in lines:
    xs, ys = px(c)
    ax.plot(xs, ys, "-", color="royalblue", lw=1.6)
    m = min(len(c) // 2, len(c) - 2)
    ax.annotate("", xy=(xs[m + 1], ys[m + 1]), xytext=(xs[m], ys[m]),
                arrowprops=dict(arrowstyle="-|>", color="royalblue", lw=1.4))
# linea di separazione e vortice
SLv = np.array([[0, 0], vec(PX["F"]), vec(PX["SLend"])])
vc = vec((850, 1300))
spir = np.array([vc + (60 - 9 * t) * np.array([math.cos(1.2 - 1.1 * t), math.sin(1.2 - 1.1 * t)])
                 for t in np.linspace(0, 5.5, 120)])
sl = np.vstack([np.array([[-560, 0]]), SLv, spir])
xs, ys = px(sl)
ax.plot(xs, ys, "-", color="lime", lw=1.8)
# ricircolo nella bolla (verso: giù lungo SL, su lungo la parete verso C)
loop = np.array([vec(q) for q in [(600, 880), (700, 980), (790, 1090), (840, 1230), (790, 1270), (700, 1210),
                                   (600, 1120), (545, 1000), (560, 900), (600, 880)]])
t_ = np.linspace(0, 1, len(loop)); tt = np.linspace(0, 1, 200)
lo = np.c_[np.interp(tt, t_, loop[:, 0]), np.interp(tt, t_, loop[:, 1])]
lo = ndimage.gaussian_filter1d(lo, 6, axis=0, mode="wrap")
xs, ys = px(lo)
ax.plot(xs, ys, "-", color="gold", lw=1.5)
for m in (40, 140):
    ax.annotate("", xy=(xs[m + 2], ys[m + 2]), xytext=(xs[m], ys[m]),
                arrowprops=dict(arrowstyle="-|>", color="gold", lw=1.3))
# lungo la parete inclinata da E verso B, e gas messo in moto da DS vicino al piede
pw = np.array([vec(PX["E"]) + 6 * np.array([math.cos(math.radians(-15)), math.sin(math.radians(-15))]),
               vec(PX["B"]) + 6 * np.array([math.cos(math.radians(-15)), math.sin(math.radians(-15))])])
xs, ys = px(pw)
ax.plot(xs, ys, "-", color="royalblue", lw=1.6)
ax.annotate("", xy=(xs[1], ys[1]), xytext=(xs[0] + .6 * (xs[1] - xs[0]), ys[0] + .6 * (ys[1] - ys[0])),
            arrowprops=dict(arrowstyle="-|>", color="royalblue", lw=1.4))
for xq in (700, 880):              # dietro DS la velocità è normale all'urto (gas davanti in quiete)
    j = int(np.argmin(np.abs(DS_f[:, 0] + C[0] - xq) + 1e4 * (DS_f[:, 1] > -200)))
    seg = np.array([DS_f[j] - 170 * nDS[j], DS_f[j]])
    xs, ys = px(seg)
    ax.plot(xs, ys, "-", color="royalblue", lw=1.3)
    ax.annotate("", xy=(xs[1], ys[1]), xytext=(xs[0] + .5 * (xs[1] - xs[0]), ys[0] + .5 * (ys[1] - ys[0])),
                arrowprops=dict(arrowstyle="-|>", color="royalblue", lw=1.2))
ax.plot(RC_PTS[:, 0], RC_PTS[:, 1], "-", color="orange", lw=1.6, alpha=.8)
ax.plot(CD_PTS[:, 0], CD_PTS[:, 1], "--", color="magenta", lw=1.2, alpha=.8)
fig.savefig(os.path.join(FIG, "linee_corrente.png"), dpi=250, bbox_inches="tight", pad_inches=0.02)
plt.close(fig)


# --- 6c. diagramma p–θ (log e lineare) -----------------------------------------
th_c = np.linspace(0, -105, 300)
p_c = [s2d.p if t == 0 else gd.espansione_pm(s2d, -t, "C+").valle.p for t in th_c]
punti = [("2", 0.0, s2d.p), ("3", s3.theta, s3.p), ("4", s4.theta, s4.p), ("∞", -th_inf, P1),
         ("5", -SPIGOLO, s5.p)]
for scala_y, nome in (("log", "p_theta_log.pdf"), ("linear", "p_theta_lin.pdf")):
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    ax.plot(th_c, np.array(p_c) / 1e5, "-", color="C0", lw=1.6, label="espansione $C^+$ da 2 (p$_0$ = p$_{02}$)")
    for r in rc:
        st = r["stato"]
        p6s = st.p * gd.urto_normale((st.V - WB) / st.a, G)["p2_p1"]
        ax.annotate("", xy=(st.theta, p6s / 1e5), xytext=(st.theta, st.p / 1e5),
                    arrowprops=dict(arrowstyle="-|>", color="C3", lw=1.4))
        ax.plot(st.theta, p6s / 1e5, "o", color="C3")
        ax.text(st.theta + 1.5, p6s / 1e5, f"6 (dietro RC, da {r['nome']})", color="C3", fontsize=8, va="bottom")
    for lab, t, p in punti:
        ax.plot(t, p / 1e5, "ko" if lab != "5" else "ks", ms=5)
        if lab == "∞":
            ax.text(t - 6, p / 1e5 * (0.45 if scala_y == "log" else 1) + (0 if scala_y == "log" else -0.9), lab,
                    fontsize=10)
        else:
            ax.text(t + 1.5, p / 1e5 * (1.12 if scala_y == "log" else 1) + (0 if scala_y == "log" else 0.25), lab,
                    fontsize=10)
    ax.axhline(P1 / 1e5, color="gray", ls=":", lw=1)
    ax.text(-104, P1 / 1e5 * (1.1 if scala_y == "log" else 1) + (0 if scala_y == "log" else 0.2),
            "$p_1$ (ambiente)", color="gray", fontsize=8)
    ax.axvline(-SPIGOLO, color="k", lw=0.8, ls="--")
    ax.text(-104, (s2d.p / 1e5) * (0.5 if scala_y == "log" else 0.9), "parete\ninclinata", fontsize=8)
    ax.set_yscale(scala_y)
    ax.set_xlabel(r"$\theta$ [°]"); ax.set_ylabel(r"$p$ [bar]")
    ax.set_xlim(-108, 5)
    ax.set_ylim(*((0, 20) if scala_y == "linear" else (5e-3, 40)))
    ax.grid(alpha=.3, which="both")
    ax.legend(loc="lower left" if scala_y == "log" else "upper left", fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, nome))
    plt.close(fig)
print("\nFigure salvate in figure/: annotata.png, linee_corrente.png, p_theta_log.pdf, p_theta_lin.pdf")
