"""Prova del 01/02/1970, esercizio 3: shadowgraph del getto di un ugello a spina troncata piano.

Riproduce tutte le misure sull'immagine, tutti i numeri e le figure dello svolgimento.
Uso:  python calcoli.py   (dalla cartella dello svolgimento)

L'immagine è una foto dello schermo: prima di misurare si raddrizza la prospettiva con
l'omografia che porta i quattro spigoli della cornice nera dello shadowgraph in un rettangolo.
Coordinate immagine (raddrizzata): x verso destra, y verso il basso, in pixel.
Coordinate fisiche: X = (x - x_base)/h_b, Y = (y_asse - y)/h_b (h_b = semialtezza della base),
angoli antiorari dall'asse (theta > 0 verso l'alto), come nel formulario (F §2D).
"""
import os
import sys
import math

import numpy as np
from PIL import Image, ImageOps, ImageDraw, ImageFont
from scipy.ndimage import gaussian_filter, map_coordinates
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

QUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(QUI, "..", ".."))
import gasdyn as gd  # noqa: E402

FIG = os.path.join(QUI, "figure")
aria = gd.Gas(1.4)
deg, rad = math.degrees, math.radians


def mu(M):
    return gd.angolo_mach(M)


def nu(M):
    return gd.prandtl_meyer(M, aria)


def M_da_nu(v):
    return gd.mach_da_pm(v, aria)


def p_p0(M):
    return 1.0 / gd.isentropiche(M, aria)["p0_p"]


# =============================================================================
# 1. RADDRIZZAMENTO DELLA FOTO
# =============================================================================
foto = ImageOps.exif_transpose(Image.open(os.path.join(FIG, "originale.jpg")))
g0 = np.asarray(foto.convert("L")).astype(float)
X0, Y0, X1, Y1 = 1250, 450, 3750, 3010          # riquadro che contiene la cornice nera
from scipy import ndimage  # noqa: E402
scuro = gaussian_filter(g0[Y0:Y1, X0:X1], 2) < 70
lab, n = ndimage.label(scuro)
k = 1 + int(np.argmax(ndimage.sum(scuro, lab, range(1, n + 1))))
ys, xs = np.nonzero(lab == k)
xs, ys = xs + X0, ys + Y0
spigoli = []
for f in (-(xs + ys), xs - ys, xs + ys, ys - xs):          # TL, TR, BR, BL
    i = int(np.argmax(f))
    spigoli.append((int(xs[i]), int(ys[i])))
src = np.array(spigoli, float)
Wr = 0.5 * (np.hypot(*(src[1] - src[0])) + np.hypot(*(src[2] - src[3])))
Hr = 0.5 * (np.hypot(*(src[3] - src[0])) + np.hypot(*(src[2] - src[1])))
dst = np.array([[0, 0], [Wr, 0], [Wr, Hr], [0, Hr]])


def omografia(s, d):
    A = []
    for (x, y), (u, v) in zip(s, d):
        A.append([x, y, 1, 0, 0, 0, -u * x, -u * y, -u])
        A.append([0, 0, 0, x, y, 1, -v * x, -v * y, -v])
    h = np.linalg.svd(np.array(A))[2][-1]
    return (h / h[-1]).reshape(3, 3)


Hinv = np.linalg.inv(omografia(src, dst))
rett = foto.transform((int(Wr), int(Hr)), Image.PERSPECTIVE, tuple(Hinv.flatten()[:8]), Image.BICUBIC)
rett.save(os.path.join(FIG, "raddrizzata.png"))
im = np.asarray(rett.convert("L")).astype(float)
# la foto dello schermo ha un reticolo di moiré di passo ~6 px: si filtra prima di misurare
g3 = gaussian_filter(im, 3.0)
g4 = gaussian_filter(im, 4.0)
print(f"Spigoli cornice (foto): {spigoli}; immagine raddrizzata {Wr:.0f} x {Hr:.0f} px "
      f"(rapporto altezza/larghezza {Hr / Wr:.3f})")

# =============================================================================
# 2. GEOMETRIA DELLA SPINA: BASE, RAMPE, ASSE
# =============================================================================
# base (troncamento): colonna di puntini chiari fra il corpo nero e la regione di base
prof_base = g3[1020:1330, 340:390].mean(0)
x_base = 340 + int(np.argmax(prof_base))
SOGLIA = 27.0        # corpo nero ~17, corrente >= 30 (immagine filtrata)


def bordo_sup(x):    # ultimo pixel chiaro sopra il corpo nero (rampa superiore)
    c = g4[760:1180, x]
    r = len(c) - 1
    while r > 0 and c[r - 1] <= SOGLIA:
        r -= 1
    return 760 + r


def bordo_inf(x):    # primo pixel chiaro sotto il corpo nero (rampa inferiore)
    c = g4[1180:1650, x]
    r = 0
    while r < len(c) - 1 and c[r + 1] <= SOGLIA:
        r += 1
    return 1180 + r


xr = np.arange(110, x_base - 18, 3)
ys_sup = np.array([bordo_sup(x) for x in xr], float)
ys_inf = np.array([bordo_inf(x) for x in xr], float)
c_sup = np.polyfit(xr, ys_sup, 2)
c_inf = np.polyfit(xr, ys_inf, 2)
# tratto finale delle rampe (ultimi ~100 px): retta, per l'inclinazione allo spigolo
kf = xr >= x_base - 120
l_sup = np.polyfit(xr[kf], ys_sup[kf], 1)
l_inf = np.polyfit(xr[kf], ys_inf[kf], 1)
y_cs = float(np.polyval(l_sup, x_base))          # spigolo superiore della base
y_ci = float(np.polyval(l_inf, x_base))          # spigolo inferiore
y_asse = 0.5 * (y_cs + y_ci)
h_b = 0.5 * (y_ci - y_cs)                        # semialtezza della base


def th_rampa(c, x):  # inclinazione della rampa [gradi, antiorario]
    return -deg(math.atan(np.polyval(np.polyder(c), x)))


th_rs, th_ri = -deg(math.atan(l_sup[0])), -deg(math.atan(l_inf[0]))
print(f"Base a x = {x_base} px; spigoli a y = {y_cs:.0f} e {y_ci:.0f} px; asse y = {y_asse:.0f} px; "
      f"semialtezza della base h_b = {h_b:.0f} px")
print(f"Rampe: residuo del profilo quadratico sup {np.std(ys_sup - np.polyval(c_sup, xr)):.1f} px, "
      f"inf {np.std(ys_inf - np.polyval(c_inf, xr)):.1f} px")
print(f"Inclinazione delle rampe nel tratto finale: sup {th_rs:.1f} gradi, inf {th_ri:+.1f} gradi; "
      f"a monte (x = 130): sup {th_rampa(c_sup, 130):.1f}, inf {th_rampa(c_inf, 130):+.1f}")


def coppia(a_sup, a_inf):
    """Strutture simmetriche: semiapertura media rispetto all'asse e inclinazione residua."""
    return 0.5 * (a_inf - a_sup), 0.5 * (a_inf + a_sup)


beta_r, rot_r = coppia(th_rs, th_ri)
print(f"Rampe allo spigolo: inclinazione media verso l'asse {beta_r:.1f} gradi "
      f"(asimmetria sup/inf {rot_r:+.1f} gradi)")
H_b = 2 * h_b


# =============================================================================
# 3. MISURE DELLE ONDE: rette di massima luminosità media (creste chiare)
# =============================================================================
def retta(x0, y0, x1, y1, xa, xb, modo="max", dy=14.0, dang=6.0, g=None):
    """Cerca la retta y = a + b x che massimizza (o minimizza) la luminosità media sul tratto
    [xa, xb], vicino alla retta di primo tentativo per (x0, y0)-(x1, y1)."""
    g = g3 if g is None else g
    ang0 = deg(math.atan2(y1 - y0, x1 - x0))
    xs_ = np.arange(xa, xb, 1.0)
    best = None
    for ang in np.arange(ang0 - dang, ang0 + dang + 1e-9, 0.1):
        t = math.tan(rad(ang))
        for off in np.arange(-dy, dy + 1e-9, 0.5):
            yy = y0 + off + t * (xs_ - x0)
            v = map_coordinates(g, [yy, xs_], order=1).mean()
            s_ = v if modo == "max" else -v
            if best is None or s_ > best[0]:
                best = (s_, t, y0 + off - t * x0, v)
    a, b, v = best[2], best[1], best[3]
    return a, b, -deg(math.atan(b)), v      # inclinazione fisica (antioraria)


def stampa(nome, r, xa, xb):
    a, b, th, v = r
    print(f"{nome}: inclinazione {th:+.1f} gradi, da ({xa}, {a + b * xa:.0f}) a ({xb}, {a + b * xb:.0f}) px, "
          f"luminosità media {v:.0f}")


# (1) contorni esterni del getto, dai labbri della carenatura
r1s = retta(200, 800, 550, 810, 200, 520)
r1i = retta(200, 1530, 550, 1525, 220, 520)
# (2) lip shock: tratto chiaro corto che parte dallo spigolo della base
r2s = retta(390, 955, 500, 960, 395, 480)
r2i = retta(390, 1380, 500, 1375, 395, 480)
# (3) contorni interni del getto (strati di scorrimento): primo tratto e tratto prima del collo
r3s_a = retta(400, 1010, 480, 1050, 400, 480, dang=10)
r3i_a = retta(400, 1330, 480, 1280, 400, 480, dang=10)
r3s_c = retta(560, 1090, 660, 1140, 560, 640, dang=10)
r3i_c = retta(560, 1230, 660, 1200, 560, 640, dang=10)
r3s_b = retta(480, 1050, 650, 1130, 480, 560, dang=10)
r3i_b = retta(480, 1280, 650, 1205, 480, 560, dang=10)
# (4) urti obliqui che partono dal collo
r4s = retta(800, 1100, 1390, 860, 850, 1300)
r4i = retta(800, 1220, 1400, 1430, 850, 1250)
for nome, r, xa, xb in [("(1) contorno esterno sup", r1s, 200, 520), ("(1) contorno esterno inf", r1i, 220, 520),
                        ("(2) lip shock sup", r2s, 395, 480), ("(2) lip shock inf", r2i, 395, 480),
                        ("(3) strato sup, primo tratto", r3s_a, 400, 480),
                        ("(3) strato inf, primo tratto", r3i_a, 400, 480),
                        ("(3) strato sup, tratto medio", r3s_b, 480, 560),
                        ("(3) strato inf, tratto medio", r3i_b, 480, 560),
                        ("(3) strato sup, prima del collo", r3s_c, 560, 640),
                        ("(3) strato inf, prima del collo", r3i_c, 560, 640),
                        ("(4) urto obliquo sup", r4s, 850, 1300), ("(4) urto obliquo inf", r4i, 850, 1250)]:
    stampa(nome, r, xa, xb)

# (5) collo della scia: ascissa dove i due contorni interni sono più vicini
coll = []
for x in range(600, 800, 4):
    c = g3[1090:1230, x]
    picchi = [i for i in range(2, len(c) - 2) if c[i] == c[i - 2:i + 3].max() and c[i] > 95]
    if len(picchi) >= 2:
        coll.append((picchi[-1] - picchi[0], x, 1090 + 0.5 * (picchi[0] + picchi[-1])))
gap_min, x_collo, y_collo = min(coll)
print(f"(5) collo della scia: x = {x_collo} px, y = {y_collo:.0f} px, distanza fra i contorni {gap_min} px; "
      f"distanza dalla base L = {x_collo - x_base} px = {(x_collo - x_base) / h_b:.2f} h_b")

# piedi degli urti (4) sui contorni esterni (1): incrocio delle due rette
def incrocio(r, s):
    x = (s[0] - r[0]) / (r[1] - s[1])
    return x, r[0] + r[1] * x


Ps = incrocio(r4s, r1s)
Pi = incrocio(r4i, r1i)
print(f"Gli urti (4) incontrano i contorni esterni a x = {Ps[0]:.0f} e {Pi[0]:.0f} px "
      f"(X = {(Ps[0] - x_base) / h_b:.1f} e {(Pi[0] - x_base) / h_b:.1f} h_b)")

# valori medi delle strutture simmetriche (semiaperture verso l'asse, positive)
beta_1, rot_1 = coppia(r1s[2], r1i[2])
beta_2, rot_2 = coppia(r2s[2], r2i[2])
beta_c, rot_c = coppia(r3s_a[2], r3i_a[2])
beta_b, rot_b = coppia(r3s_b[2], r3i_b[2])
beta_j, rot_j = coppia(r3s_c[2], r3i_c[2])
phi_4, rot_4 = coppia(-r4s[2], -r4i[2])          # l'urto si allontana dall'asse
for nome, b_, r_ in [("contorno esterno (1)", beta_1, rot_1), ("lip shock (2)", beta_2, rot_2),
                     ("strato, primo tratto (3)", beta_c, rot_c), ("strato, tratto medio (3)", beta_b, rot_b),
                     ("strato, prima del collo (3)", beta_j, rot_j)]:
    print(f"  {nome}: inclinato verso l'asse di {b_:.1f} gradi (asimmetria {r_:+.1f})")
print(f"  urto obliquo (4): inclinato di {phi_4:.1f} gradi sull'asse, verso l'esterno (asimmetria {rot_4:+.1f})")


# =============================================================================
# 4. VALORI LETTI (come con il goniometro: arrotondati al grado)
# =============================================================================
TH_R = -round(beta_r)          # rampa allo spigolo (metà superiore, verso l'asse)
TH_C = -round(beta_c)          # contorno interno subito dopo lo spigolo
TH_J = -round(beta_j)          # contorno interno subito prima del collo
PHI_4 = round(phi_4)           # urto obliquo (4) sull'asse
TH_1 = -round(beta_1)          # contorno esterno
L_HB = (x_collo - x_base) / h_b
print(f"\nValori usati: rampa {TH_R} gradi, contorno interno {TH_C} -> {TH_J} gradi, urto {PHI_4} gradi, "
      f"contorno esterno {TH_1} gradi, collo a L = {L_HB:.1f} h_b")


# =============================================================================
# 5. CALCOLO (metà superiore; la metà inferiore è speculare)
# =============================================================================
def catena(th_r, th_c, th_j, phi4, stampa_=False):
    """Mach e pressioni lungo il contorno interno: r (fine rampa) -> c (dopo il ventaglio allo
    spigolo, pressione di base) -> j (prima dell'urto) -> 2 (dopo l'urto, corrente assiale)."""
    theta = -th_j                       # deviazione imposta dall'urto (da th_j a 0)
    sigma = phi4 + theta                # angolo d'urto rispetto alla corrente a monte (F §12)
    Mj = gd.mach_da_theta_sigma(theta, sigma, aria)
    nu_j = nu(Mj)
    # onde C+ (nascono dal contorno interno, sotto la corrente): theta + nu costante (F §14)
    nu_c = nu_j + th_j - th_c
    nu_r = nu_j + th_j - th_r
    Mc, Mr = M_da_nu(nu_c), M_da_nu(nu_r)
    St_j = gd.State.da_totali(aria, Mj, p0=1.0, a0=1.0, theta=th_j, nome="j")
    u = gd.urto_obliquo(St_j, sigma=sigma, famiglia="C+", nome="2")
    out = dict(theta=theta, sigma=sigma, Mj=Mj, Mc=Mc, Mr=Mr, nu_j=nu_j, nu_c=nu_c, nu_r=nu_r,
               p_j=p_p0(Mj), p_c=p_p0(Mc), p_r=p_p0(Mr), urto=u, M2=u.valle.M, p2=u.valle.p,
               th2=u.valle.theta, p02=u.valle.p0, tm=gd.theta_max(Mj, aria))
    return out


C = catena(TH_R, TH_C, TH_J, PHI_4)
u4 = C["urto"]
print(f"Urto (4): deviazione theta = {C['theta']:.0f} gradi, sigma = {C['sigma']:.0f} gradi "
      f"-> M_j = {C['Mj']:.2f} (theta-sigma-M, F §12), ramo {u4.ramo}; theta_max(M_j) = {C['tm'][0]:.1f} gradi "
      f"(sigma = {C['tm'][1]:.1f})")
print(f"  M_n = {C['Mj'] * math.sin(rad(C['sigma'])):.2f}, p2/p_j = {u4.p2_p1:.2f}, M_2 = {C['M2']:.2f}, "
      f"theta_2 = {C['th2']:+.2f} gradi, p02/p0j = {u4.p02_p01:.3f}")
print(f"nu: j {C['nu_j']:.1f}, c {C['nu_c']:.1f}, r {C['nu_r']:.1f} gradi -> M_r = {C['Mr']:.2f}, "
      f"M_c = {C['Mc']:.2f}")
print(f"p/p0: r {C['p_r']:.4f}, c (base) {C['p_c']:.4f}, j {C['p_j']:.4f}, 2 {C['p2']:.4f}")
print(f"Rapporti: p_b/p_r = {C['p_c'] / C['p_r']:.2f}, p_j/p_b = {C['p_j'] / C['p_c']:.2f}, "
      f"p2/p_b = {C['p2'] / C['p_c']:.2f}")
print(f"mu: r {mu(C['Mr']):.1f}, c {mu(C['Mc']):.1f}, j {mu(C['Mj']):.1f}, 2 {mu(C['M2']):.1f} gradi")
# controllo geometrico: inclinazione media del contorno interno dalla base al collo
print(f"Testa del ventaglio allo spigolo: th_r + mu_r = {TH_R + mu(C['Mr']):.1f} gradi, coda th_c + mu_c = {TH_C + mu(C['Mc']):.1f}; linea 2 misurata {-beta_2:.1f} gradi")
print(f"Inclinazione media spigolo-collo: atan(h_b/L) = {deg(math.atan(1 / L_HB)):.1f} gradi "
      f"(fra {-TH_C} e {-TH_J}: {'sì' if -TH_J < deg(math.atan(1 / L_HB)) < -TH_C else 'no'})")
# lip shock: angolo con la corrente dopo il ventaglio dello spigolo (solo per la sintesi)
sig_ls = -TH_C - beta_2
print(f"(fuori PDF) lip shock: angolo con la corrente c = {sig_ls:.0f} gradi, mu(M_c) = {mu(C['Mc']):.1f}: "
      f"M_n = {C['Mc'] * math.sin(rad(sig_ls)):.2f} -> inclinazione non usata (tratto corto, nello strato viscoso)")

# --- sensibilità alle misure (solo chat): +-3 gradi sugli strati, +-2 su rampa e urto
print("\nSensibilità (solo sintesi):")
ris = []
for dr in (-2, 0, 2):
    for dc in (-3, 0, 3):
        for dj in (-3, 0, 3):
            for d4 in (-2, 0, 2):
                try:
                    c_ = catena(TH_R + dr, TH_C + dc, TH_J + dj, PHI_4 + d4)
                    ris.append((c_["Mj"], c_["Mc"], c_["Mr"], c_["p_c"], c_["p2"] / c_["p_c"], c_["M2"]))
                except gd.GasdynError:
                    pass
ris = np.array(ris)
for i, nome in enumerate(["M_j", "M_c", "M_r", "p_b/p0", "p2/p_b", "M_2"]):
    print(f"  {nome}: {ris[:, i].min():.3g} - {ris[:, i].max():.3g}")


# =============================================================================
# 6. FIGURE
# =============================================================================
plt.rcParams.update({"font.size": 10, "axes.grid": True, "grid.alpha": 0.3})


def punto(r, x):
    return x, r[0] + r[1] * x


# --- 6.1 figura annotata (immagine raddrizzata)
XA, XB, YA, YB = 0, 1650, 640, 1720
fig, ax = plt.subplots(figsize=(8.6, 5.7))
ax.imshow(im[YA:YB, XA:XB], cmap="gray", extent=(XA, XB, YB, YA), vmin=0, vmax=200)
ax.grid(False)
kw = dict(lw=1.6, alpha=0.9)
# asse e base
ax.plot([x_base, XB], [y_asse, y_asse], "w--", lw=0.8, alpha=0.7)
ax.text(XB - 10, y_asse - 10, "asse", color="w", ha="right", fontsize=8)
# rampe (profilo misurato) e base
xx = np.linspace(110, x_base, 50)
ax.plot(xx, np.polyval(c_sup, xx), color="tab:orange", **kw)
ax.plot(xx, np.polyval(c_inf, xx), color="tab:orange", **kw)
ax.plot([x_base, x_base], [y_cs, y_ci], color="tab:orange", **kw)
# (1) contorni esterni
for r, (xa, xb) in [(r1s, (170, 1395)), (r1i, (170, 1511))]:
    ax.plot(*zip(punto(r, xa), punto(r, xb)), color="deepskyblue", **kw)
# (2) lip shock
for r in (r2s, r2i):
    ax.plot(*zip(punto(r, x_base + 10), punto(r, 520)), color="yellow", **kw)
# (3) contorni interni: spezzata dei tre tratti fino al collo
for (ra, rb, rc), sgn in [((r3s_a, r3s_b, r3s_c), -1), ((r3i_a, r3i_b, r3i_c), +1)]:
    pts = [punto(ra, x_base), punto(ra, 480), punto(rb, 560), punto(rc, 640),
           (x_collo, y_collo + sgn * gap_min / 2)]
    ax.plot(*zip(*pts), color="lime", **kw)
# (4) urti obliqui dal collo ai contorni esterni
for r, P in [(r4s, Ps), (r4i, Pi)]:
    ax.plot([x_collo, P[0]], [y_collo, P[1]], color="red", **kw)
# (5) punto di ristagno al collo, (6) regione di base
ax.plot(x_collo, y_collo, "o", ms=7, mfc="none", mec="magenta", mew=2)
poli = [(x_base, y_cs), punto(r3s_a, 480), punto(r3s_b, 560), punto(r3s_c, 640), (x_collo, y_collo),
        punto(r3i_c, 640), punto(r3i_b, 560), punto(r3i_a, 480), (x_base, y_ci)]
ax.fill(*zip(*poli), color="magenta", alpha=0.18, lw=0)
etich = [("1", (300, 760), "deepskyblue"), ("2", (500, 930), "yellow"), ("3", (540, 1040), "lime"),
         ("4", (1150, 950), "red"), ("5", (730, 1215), "magenta"), ("6", (450, 1180), "magenta"),
         ("1", (300, 1590), "deepskyblue"), ("2", (500, 1415), "yellow"), ("3", (540, 1305), "lime"),
         ("4", (1150, 1355), "red")]
for t, (x, y), c in etich:
    ax.text(x, y, t, color=c, fontsize=12, fontweight="bold", ha="center", va="center",
            bbox=dict(boxstyle="circle,pad=0.15", fc="black", ec=c, alpha=0.75))
# stati e angoli misurati
for t, (x, y) in [("r", (300, 1000)), ("c", (430, 975)), ("j", (600, 1075)), ("2", (930, 1120)),
                  ("scia", (1000, 1175))]:
    ax.text(x, y, t, color="w", fontsize=10, style="italic", ha="center", va="center")
ax.text(60, 690, f"rampa {abs(TH_R)}°, contorno interno {abs(TH_C)}° → {abs(TH_J)}°,\n"
        f"urto obliquo {PHI_4}° sull'asse, contorno esterno {abs(TH_1)}°", color="w", fontsize=8.5,
        va="top", bbox=dict(fc="black", alpha=0.6, ec="none"))
ax.set_xlim(XA, XB)
ax.set_ylim(YB, YA)
ax.set_xticks([])
ax.set_yticks([])
fig.tight_layout()
fig.savefig(os.path.join(FIG, "annotata.png"), dpi=200)
plt.close(fig)

# --- 6.2 schema del campo (metà superiore, coordinate in unità h_b)
fig, ax = plt.subplots(figsize=(8.6, 4.3))
ax.grid(False)
Xr = lambda x: (x - x_base) / h_b              # noqa: E731
Yr = lambda y: (y_asse - y) / h_b              # noqa: E731
xx = np.linspace(110, x_base, 50)
# carenatura stilizzata e rampa
lip = (Xr(150), Yr(785))
ax.plot([lip[0] - 0.9, lip[0] - 0.9, lip[0]], [lip[1] + 1.2, lip[1] + 0.25, lip[1]], "k", lw=2)
ax.fill_between(Xr(xx), 0, Yr(np.polyval(c_sup, xx)), color="0.25")
ax.fill_between([Xr(110) - 0.6, Xr(110)], 0, Yr(np.polyval(c_sup, 110)), color="0.25")
ax.plot([0, 0], [0, 1], "k", lw=2)
ax.axhline(0, color="k", ls="-.", lw=0.8)
# ventaglio dal labbro: raggi verso la rampa
for xw in np.linspace(130, 330, 6):
    ax.plot([lip[0], Xr(xw)], [lip[1], Yr(np.polyval(c_sup, xw))], color="tab:blue", lw=0.6, alpha=0.6)
# contorno esterno (1) e urto (4) con riflessione in espansione
XP, YP = Xr(Ps[0]), Yr(Ps[1])
ax.plot([lip[0], XP], [lip[1], YP], color="deepskyblue", lw=2)
ax.plot([XP, XP + 3], [YP, YP + 3 * math.tan(rad(5))], color="deepskyblue", lw=2, ls="--")
Xn, Yn = Xr(x_collo), 0.0
ax.plot([Xn, XP], [Yn, YP], color="red", lw=2)
for d in (mu(C["M2"]), 34, 30):     # espansione riflessa (onde C-, testa a -mu_2)
    ax.plot([XP, XP + 2.2 * math.cos(rad(d))], [YP, YP - 2.2 * math.sin(rad(d))], color="tab:blue", lw=0.7)
# ventaglio allo spigolo (onde C+: testa a th_r + mu_r, coda a th_c + mu_c), la cui testa è la linea 2
testa, coda = TH_R + mu(C["Mr"]), TH_C + mu(C["Mc"])
for d in np.linspace(testa, coda, 4):
    ax.plot([0, 0.9 * math.cos(rad(d))], [1, 1 + 0.9 * math.sin(rad(d))], color="tab:blue", lw=0.7)
ax.plot([0, 0.95 * math.cos(rad(-beta_2))], [1, 1 - 0.95 * math.sin(rad(beta_2))], color="goldenrod", lw=1.8)
# contorno interno curvo: da TH_C allo spigolo fino al collo (parabola con pendenza iniziale TH_C)
s_ = np.linspace(0, 1, 60)
xs3 = Xn * s_
ys3 = 1 - (s_ * Xn) * math.tan(rad(abs(TH_C))) + (math.tan(rad(abs(TH_C))) * Xn - 1) * s_ ** 2
ax.plot(xs3, ys3, color="green", lw=2)
# onde di compressione C+ dal contorno curvo (inclinate di theta + mu), che rinforzano l'urto 4
for f in (0.45, 0.6, 0.75):
    i_ = int(f * 59)
    x0, y0 = xs3[i_], ys3[i_]
    th_loc = -deg(math.atan((ys3[i_ - 1] - ys3[i_ + 1]) / (xs3[i_ + 1] - xs3[i_ - 1])))
    d = th_loc + mu(C["Mj"])
    # fino all'incrocio con l'urto 4 (retta dal collo a PHI_4)
    t4, td = math.tan(rad(PHI_4)), math.tan(rad(d))
    xi = (y0 - td * x0 + t4 * Xn) / (t4 - td)
    ax.plot([x0, xi], [y0, y0 + td * (xi - x0)], color="red", lw=0.6, alpha=0.7)
ax.fill_between(xs3, 0, ys3, color="magenta", alpha=0.15)
ax.plot(Xn, 0, "o", ms=7, mfc="w", mec="magenta", mew=2)
ax.plot([Xn, Xn + 4.5], [0.02, 0.05], color="0.4", lw=3, alpha=0.5)
# etichette
for t, (x, y), c in [("1", (2.5, 1.95), "deepskyblue"), ("2", (1.05, 1.08), "goldenrod"),
                     ("3", (0.75, 0.62), "green"), ("4", (3.6, 0.95), "red"), ("5", (Xn, -0.22), "magenta"),
                     ("6", (0.55, 0.3), "magenta")]:
    ax.text(x, y, t, color=c, fontsize=12, fontweight="bold", ha="center", va="center")
for t, (x, y) in [("$r$", (-0.3, 1.2)), ("$c$", (0.82, 0.66)), ("$j$", (1.45, 0.55)), ("$2$", (3.8, 0.35)),
                  ("ambiente $p_a$", (2.0, 2.45)), ("scia", (4.6, 0.15)), ("ventaglio dal labbro", (-1.5, 2.48)),
                  ("ventaglio allo spigolo", (1.45, 1.22))]:
    ax.text(x, y, t, fontsize=9, ha="center", va="center")
ax.set_aspect("equal")
ax.set_xlim(-2.6, 7.2)
ax.set_ylim(-0.45, 2.7)
ax.set_xlabel("$x/h_b$")
ax.set_ylabel("$y/h_b$")
fig.tight_layout()
fig.savefig(os.path.join(FIG, "schema.pdf"))
plt.close(fig)

# --- 6.3 piano p-theta: percorso r -> c -> j -> 2 lungo il contorno interno (e speculare)
St_r = gd.State.da_totali(aria, C["Mr"], p0=1.0, a0=1.0, theta=TH_R, nome="r")
St_j = u4.monte
Ms_ = np.linspace(1.6, 4.3, 300)
th_is = TH_R - (np.array([nu(m) for m in Ms_]) - C["nu_r"])       # curva isentropica C+ per r
p_is = np.array([p_p0(m) for m in Ms_])
pol = gd.polare_p_theta(St_j, "C+", espansione=False)["urto"]
deb = pol["ramo"] == "debole"
fig, axs = plt.subplots(1, 2, figsize=(9.6, 4.2))
for ax, scala in zip(axs, ("log", "linear")):
    for sgn, ls in ((1, "-"), (-1, ":")):
        ax.plot(sgn * th_is, p_is, "tab:blue", ls=ls, lw=1.3)
        ax.plot(sgn * pol["theta"][deb], pol["p"][deb], "tab:red", ls=ls, lw=1.3)
        ax.plot(sgn * pol["theta"][~deb], pol["p"][~deb], "tab:red", ls=ls, lw=0.7, alpha=0.5)
    pts = [("r", TH_R, C["p_r"]), ("c", TH_C, C["p_c"]), ("j", TH_J, C["p_j"]), ("2", 0.0, C["p2"])]
    for t, th, p in pts:
        ax.plot(th, p, "ko", ms=4)
        ax.annotate(t, (th, p), textcoords="offset points", xytext=(5, 4), fontsize=10)
        if t != "2":
            ax.plot(-th, p, "o", mfc="w", mec="k", ms=4)
            ax.annotate(t + "'", (-th, p), textcoords="offset points", xytext=(5, 4), fontsize=9)
    ax.annotate("", xy=(TH_C, C["p_c"]), xytext=(TH_R, C["p_r"]),
                arrowprops=dict(arrowstyle="->", color="k", lw=0.8))
    ax.annotate("", xy=(TH_J, C["p_j"]), xytext=(TH_C, C["p_c"]),
                arrowprops=dict(arrowstyle="->", color="k", lw=0.8))
    ax.axvline(0, color="k", lw=0.6)
    ax.axhline(C["p_c"], color="magenta", lw=0.8, ls="--")
    ax.text(32, C["p_c"] * (1.08 if scala == "log" else 1.0) + (0 if scala == "log" else 0.004),
            "$p_b$", color="magenta", ha="right")
    ax.set_yscale(scala)
    ax.set_xlim(-45, 45)
    ax.set_ylim((0.008, 0.4) if scala == "log" else (0, 0.32))
    ax.set_xlabel(r"$\theta$ [°]")
    ax.set_ylabel("$p/p_0$")
    ax.set_title("scala logaritmica" if scala == "log" else "scala lineare", fontsize=10)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "p_theta.pdf"))
plt.close(fig)
print("\nFigure scritte in", FIG)
