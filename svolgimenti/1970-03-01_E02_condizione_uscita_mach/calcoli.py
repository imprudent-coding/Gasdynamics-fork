"""Prova del 01/03/1970, esercizio 2 (traccia identica a B70-E02, 13/02/2023): condizione al
contorno all'uscita di un condotto, per un flusso d'aria instazionario omoentropico subsonico,

    (gamma+1)/2 * (u/a)^(4/(gamma+1)) + (a/u)^(2(gamma-1)/(gamma+1)) = K,   K = cost. > 1,

ed effetto su un'onda semplice acustica che arriva all'uscita dall'interno del condotto.

Riproduce tutti i numeri e le figure dello svolgimento.
Uso:  python calcoli.py   (dalla cartella dello svolgimento)

Il problema e' letterale (nessun dato numerico): i numeri servono per l'aria (gamma = 1,4) e per un
esempio illustrativo (M_e = 0,3, p0 = a0 = 1, onde incidenti con p_i/p0 = 1,1 e 0,9).

Parti fuori modulo (scritte qui riusando gasdyn):
  - la funzione f(M) della traccia, il suo minimo e le radici di f(M) = K (brentq);
  - lo stato all'uscita dopo la riflessione con M imposto (R2 dall'interno + u = M a): formula
    chiusa, controllata con gd.onda_semplice (onda u-a da i) e con il limite M -> 0 (parete);
  - il diagramma x-t con caratteristiche rettilinee (semirette, come nella soluzione ufficiale
    del 14/09/2026) e il piano p-u con le curve del contorno.
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import brentq

QUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(QUI, "..", ".."))
import gasdyn as gd  # noqa: E402

FIG = os.path.join(QUI, "figure")
os.makedirs(FIG, exist_ok=True)
plt.rcParams.update({"font.size": 10, "axes.grid": True, "grid.alpha": 0.3})

G = gd.Gas(1.4)
g, d = G.gamma, G.delta                    # delta = (gamma-1)/2 (F §1)
E1 = 4 / (g + 1)                           # esponente di M nel primo termine
E2 = 2 * (g - 1) / (g + 1)                 # esponente di 1/M nel secondo termine


def f(M):
    """Primo membro della condizione della traccia, funzione del solo M = u/a."""
    M = np.asarray(M, float)
    return (g + 1) / 2 * M**E1 + M**(-E2)


def fp(M):
    """f'(M) = 2 M^-(3g-1)/(g+1) [M^2 - (g-1)/(g+1)]: si annulla solo in M_m."""
    return 2 * M**(-(3 * g - 1) / (g + 1)) * (M**2 - (g - 1) / (g + 1))


def coeff(M):
    """Coefficiente di riflessione con M imposto: dR1 = C dR2, C = (1 - delta M)/(1 + delta M)."""
    return (1 - d * M) / (1 + d * M)


# =============================================================================
# 1. LA CONDIZIONE E' M = cost: minimo di f, radici, casi esclusi
# =============================================================================
print("=" * 78)
print("1. Condizione f(M) = K")
Mm = np.sqrt((g - 1) / (g + 1))
fmin = (g + 1) / 2 * ((g + 1) / (g - 1))**((g - 1) / (g + 1))
f1 = (g + 3) / 2
print(f"   M_m = sqrt((g-1)/(g+1)) = {Mm:.4f}")
print(f"   f(M_m) = (g+1)/2 ((g+1)/(g-1))^((g-1)/(g+1)) = {fmin:.4f}  (diretto: {float(f(Mm)):.4f})")
print(f"   f(1) = (g+3)/2 = {f1:.4f}  (diretto: {float(f(1.0)):.4f})")
assert abs(float(f(Mm)) - fmin) < 1e-12 and abs(float(f(1.0)) - f1) < 1e-12
# controllo numerico della derivata analitica
for M in (0.2, 0.5, 0.9):
    h = 1e-6
    assert abs((f(M + h) - f(M - h)) / (2 * h) - fp(M)) < 1e-6
print("   f'(M) analitica controllata con le differenze finite: ok")
# f(M) >= f(M_m) > 1 per ogni M > 0: K > 1 e' necessaria, non sufficiente
Mg = np.linspace(0.01, 3, 3000)
assert f(Mg).min() >= fmin - 1e-9 and fmin > 1


def radici(K):
    """Radici subsoniche (0 < M < 1) di f(M) = K."""
    out = []
    if K > fmin:
        out.append(brentq(lambda M: f(M) - K, 1e-8, Mm))
        if K < f1:
            out.append(brentq(lambda M: f(M) - K, Mm, 1.0))
    elif abs(K - fmin) < 1e-12:
        out.append(Mm)
    return out


for K in (1.5, 1.7, 2.0, 2.5):
    print(f"   K = {K}: radici subsoniche {[round(r, 4) for r in radici(K)]}")

# =============================================================================
# 2. RIFLESSIONE: R1 = C R2 all'uscita (esatta), coefficiente per l'aria
# =============================================================================
print("=" * 78)
print("2. Coefficiente di riflessione C = (1 - delta M)/(1 + delta M)")
print(f"   M -> 0: C = 1 (parete);  M -> 1: C = (3-g)/(g+1) = {(3 - g) / (g + 1):.4f}")
assert abs(coeff(1.0) - (3 - g) / (g + 1)) < 1e-12
print("   M_e     C      1+C    1-C   (1+M)/(1-M)*C")
for M in (0.1, 0.2, 0.3, Mm, 0.6, 0.8):
    C = coeff(M)
    print(f"   {M:.3f}  {C:.3f}  {1 + C:.3f}  {1 - C:.3f}  {(1 + M) / (1 - M) * C:.3f}")

# =============================================================================
# 3. ESEMPIO: M_e = 0,3, p0 = a0 = 1; onda incidente u+a di compressione ed espansione
# =============================================================================
print("=" * 78)
ME = 0.3
KE = float(f(ME))
r = radici(KE)
print(f"3. Esempio M_e = {ME}: K = f(M_e) = {KE:.4f}, radici subsoniche {[round(x, 4) for x in r]}")
C = coeff(ME)
s0 = gd.State.da_pa(G, 1.0, 1.0, u=ME, nome="0")        # stato iniziale, u0 = M_e a0


def stato_uscita(si):
    """Stato e all'uscita dopo la riflessione: R2_e = R2_i (onda u-a), u_e = M_e a_e (F §3, §4)."""
    a_e = si.R2 / (1 + d * ME)
    e = gd.onda_semplice(si, "u-a", a=a_e, nome="e")
    assert abs(e.u - ME * e.a) < 1e-12 and abs(e.R2 - si.R2) < 1e-12
    return e


ris = {}
for nome, pr in (("compressione", 1.10), ("espansione", 0.90)):
    si = gd.onda_semplice(s0, "u+a", p=pr * s0.p, nome="i")   # R1 uniforme nell'onda incidente
    e = stato_uscita(si)
    dai, dar = si.a - s0.a, e.a - si.a
    dui, dur = si.u - s0.u, e.u - si.u
    print(f"   {nome}: p_i = {si.p:.4f}, u_i = {si.u:.4f}, a_i = {si.a:.5f}")
    print(f"      e: p_e = {e.p:.4f}, u_e = {e.u:.4f}, a_e = {e.a:.5f}, M_e = {e.u / e.a:.4f}")
    print(f"      da_r/da_i = {dar / dai:.4f} (C = {C:.4f}),  du_r/du_i = {dur / dui:.4f} (-C)")
    print(f"      dp_r/dp_i = {(e.p - si.p) / (si.p - s0.p):.4f},  dp_e/dp_i = {(e.p - s0.p) / (si.p - s0.p):.4f}"
          f"  (acustico: C = {C:.3f}, 1+C = {1 + C:.3f})")
    print(f"      relative: (dp_r/p_i)/(dp_i/p0) = {((e.p - si.p) / si.p) / ((si.p - s0.p) / s0.p):.4f}")
    assert abs(dar / dai - C) < 1e-12 and abs(dur / dui + C) < 1e-12   # esatte in a e u
    # confronto con altri contorni: pressione imposta (p_e = p0) e velocita' imposta (u_e = u0)
    ep = gd.onda_semplice(si, "u-a", p=s0.p, nome="p")
    eu = gd.onda_semplice(si, "u-a", u=s0.u, nome="u")
    print(f"      p imposta: u = {ep.u:.4f};  u imposta: p = {eu.p:.4f}")
    # limite M -> 0: la condizione diventa u = 0 (parete), C = 1
    ris[nome] = (si, e, ep, eu)

# limite M -> 0 controllato con il contorno di parete del modulo
q = gd.State.da_pa(G, 1.0, 1.0)
qi = gd.onda_semplice(q, "u+a", p=1.1)
w = gd.contorno_parete(qi, "destra")
aw = qi.R2 / (1 + d * 0.0)
print(f"   limite M -> 0: a_e = R2_i = {aw:.5f}, contorno di parete del modulo a = {w.valle.a:.5f}"
      " (il modulo riflette la compressione come urto: differenza del terzo ordine)")
assert abs(aw - w.valle.a) < 1e-4

# rapporto dei gradienti spaziali (onda riflessa piu' corta: velocita' a - u invece di a + u)
print(f"   gradienti: (R1)_x / (R2)_x = -C (1+M)/(1-M) = {-C * (1 + ME) / (1 - ME):.4f}")
print(f"   portata: dm/m = (g+1)/(g-1) da_e/a = (g+1)/(2g) dp_e/p")

# =============================================================================
# 4. FIGURE
# =============================================================================
# --- 4.1 la funzione f(M)
fig, ax = plt.subplots(figsize=(6.4, 4.2))
M = np.linspace(0.03, 1.6, 800)
ax.plot(M, f(M), "k", lw=1.6, label=r"$f(M)$")
ax.axhline(fmin, color="0.5", ls=":", lw=1)
ax.axhline(f1, color="0.5", ls=":", lw=1)
ax.axvline(1.0, color="0.3", ls="--", lw=1)
ax.plot([Mm], [fmin], "ko")
ax.annotate(rf"$M_m=\sqrt{{(\gamma-1)/(\gamma+1)}}\approx{Mm:.2f}$" + "\n" + rf"$f_{{min}}\approx{fmin:.2f}$",
            (Mm, fmin), (0.5, 1.22), fontsize=9, arrowprops=dict(arrowstyle="->", lw=0.8))
ax.annotate(rf"$f(1)=(\gamma+3)/2={f1:.1f}$", (1.0, f1), (1.12, 1.9), fontsize=9,
            arrowprops=dict(arrowstyle="->", lw=0.8))
ax.axhline(KE, color="C0", lw=1)
for x in r:
    ax.plot([x], [KE], "o", color="C0")
ax.text(0.68, KE + 0.03, rf"$K\approx{KE:.2f}$: due radici subsoniche", color="C0", fontsize=9)
ax.axhline(2.5, color="C3", lw=1)
r25 = radici(2.5)[0]
ax.plot([r25], [2.5], "o", color="C3")
r25s = brentq(lambda x: f(x) - 2.5, 1.0, 3.0)
ax.plot([r25s], [2.5], "x", color="C3")
ax.text(0.1, 2.58, r"$K>(\gamma+3)/2$: una sola radice subsonica", color="C3", fontsize=9)
ax.axhspan(1.0, fmin, color="0.85", alpha=0.6)
ax.text(0.9, 1.08, r"$1<K<f_{min}$: nessuna soluzione", fontsize=9)
ax.set_xlim(0, 1.6)
ax.set_ylim(1.0, 2.9)
ax.set_xlabel(r"$M=u/a$")
ax.set_ylabel(r"$f(M)=\frac{\gamma+1}{2}M^{4/(\gamma+1)}+M^{-2(\gamma-1)/(\gamma+1)}$")
ax.set_title(r"Condizione all'uscita, aria ($\gamma=1{,}4$)")
fig.tight_layout()
fig.savefig(os.path.join(FIG, "funzione_K.pdf"))
plt.close(fig)

# --- 4.2 coefficiente di riflessione
fig, ax = plt.subplots(figsize=(6.4, 3.8))
M = np.linspace(0, 1, 200)
ax.plot(M, coeff(M), "k", lw=1.8, label=r"$M$ imposto: $\frac{1-\delta M_e}{1+\delta M_e}$")
ax.axhline(1, color="C0", ls="--", lw=1.2, label=r"$u$ imposta (parete): $+1$")
ax.axhline(0, color="C2", ls="-.", lw=1.2, label=r"non riflettente: $0$")
ax.axhline(-1, color="C3", ls="--", lw=1.2, label=r"$p$ imposta (estremo aperto): $-1$")
ax.plot([ME], [C], "ko")
ax.annotate(rf"$M_e={ME}$: ${C:.2f}$", (ME, C), (0.36, 0.45), fontsize=9,
            arrowprops=dict(arrowstyle="->", lw=0.8))
ax.plot([1], [coeff(1)], "ks", mfc="w")
ax.annotate(r"$M_e\to1$: $(3-\gamma)/(\gamma+1)=2/3$", (1, coeff(1)), (0.55, 0.2), fontsize=9,
            arrowprops=dict(arrowstyle="->", lw=0.8))
ax.set_xlim(0, 1.02)
ax.set_ylim(-1.15, 1.15)
ax.set_xlabel(r"$M_e$ all'uscita")
ax.set_ylabel(r"$\Delta R_1^{rifl}/\Delta R_2^{inc}\ \ (=\Delta p_r/\Delta p_i$ acustico$)$")
ax.legend(loc="lower left", fontsize=8.5)
ax.set_title(r"Coefficiente di riflessione all'uscita, aria")
fig.tight_layout()
fig.savefig(os.path.join(FIG, "coeff_riflessione.pdf"))
plt.close(fig)

# --- 4.3 piano x-t (semirette), esempio di compressione
si, e, _, _ = ris["compressione"]
th, tt = 1.0, 1.5                         # arrivo all'uscita di testa e coda dell'onda incidente
lh, lt = s0.u + s0.a, si.u + si.a         # pendenze (dx/dt) dell'onda incidente
rh, rt = s0.u - s0.a, e.u - e.a           # pendenze dell'onda riflessa
xmin, tmax = -2.0, 3.0
fig, ax = plt.subplots(figsize=(6.0, 4.6))
ax.axvline(0, color="k", lw=2.5)
ax.text(0.04, 2.75, "uscita\n" + r"$M=M_e$", fontsize=9, va="top")
for k, tk in enumerate(np.linspace(th, tt, 6)):
    lam = lh + (lt - lh) * k / 5
    lr = rh + (rt - rh) * k / 5
    lw = 1.6 if k in (0, 5) else 0.7
    t0 = tk + xmin / lam                  # dove la caratteristica u+a esce dal riquadro
    ax.plot([xmin, 0], [t0, tk], color="C3", lw=lw)
    ax.plot([0, xmin], [tk, tk + xmin / lr], color="C0", lw=lw)
tP = (rh * th - lt * tt) / (rh - lt)
xP = rh * (tP - th)
ax.plot([xP], [tP], "k.", ms=6)
ax.text(xP - 0.12, tP + 0.02, "P", fontsize=9, ha="right")
ax.set_xlim(xmin, 0.45)
ax.set_ylim(0, tmax)


def ruota(pend):
    """Angolo sullo schermo di una retta con dx/dt = pend (per allineare le scritte)."""
    p1, p2 = ax.transData.transform([(0, 0), (pend, 1)])
    return np.degrees(np.arctan2(p2[1] - p1[1], p2[0] - p1[0])) % 180 - (180 if pend < 0 else 0)


ax.text(-0.3, 0.3, "0", fontsize=12, weight="bold")
xl = -1.25
ax.text(xl, 0.5 * (th + xl / lh + tt + xl / lt), r"onda incidente $u+a$", color="C3", fontsize=9,
        rotation=ruota(lh), ha="center", va="center", bbox=dict(fc="w", ec="none", pad=0.5))
ax.text(-1.7, 2.75, "1", fontsize=12, weight="bold")
ax.text(-0.16, 1.27, "Z", fontsize=10)
xl = -0.8
ax.text(xl, 0.5 * (th + xl / rh + tt + xl / rt), r"onda riflessa $u-a$", color="C0", fontsize=9,
        rotation=ruota(rh), ha="center", va="center", bbox=dict(fc="w", ec="none", pad=0.5))
ax.text(-0.35, 2.7, "2", fontsize=12, weight="bold")
ax.set_xlabel(r"$x$ (uscita in $x=0$)")
ax.set_ylabel(r"$t$")
ax.set_title(r"Piano $x$–$t$ (esempio $M_e=0{,}3$, compressione)")
fig.tight_layout()
fig.savefig(os.path.join(FIG, "xt.pdf"))
plt.close(fig)

# --- 4.4 piano p-u: compressione ed espansione incidenti, tre contorni a confronto
fig, ax = plt.subplots(figsize=(6.4, 4.8))
pp = np.linspace(0.8, 1.25, 300)
ax.plot(ME * s0.a * (pp / s0.p)**((g - 1) / (2 * g)), pp, "k", lw=1.8, label=r"contorno $u=M_e\,a$")
ax.axhline(s0.p, color="C3", ls="--", lw=1, label=r"$p$ imposta ($p=p_0$)")
ax.axvline(s0.u, color="C0", ls="--", lw=1, label=r"$u$ imposta ($u=u_0$)")
# onde isentropiche: u+a da 0 (R1 = R1_0), u-a da i (R2 = R2_i); a = a_ref (p/p_ref)^((g-1)/(2g))
a_p = lambda st, p: st.a * (p / st.p)**((g - 1) / (2 * g))
ax.plot(s0.u + (a_p(s0, pp) - s0.a) / d, pp, color="0.4", lw=1, label=r"onda $u+a$ da 0")
for nome, (si, e, ep, eu) in ris.items():
    pr = np.linspace(min(si.p, e.p, ep.p, eu.p) * 0.98, max(si.p, e.p, ep.p, eu.p) * 1.02, 200)
    ax.plot(si.u - (a_p(si, pr) - si.a) / d, pr, color="0.6", lw=1, ls="-.")
    for st, lab, col in ((si, "1", "k"), (e, "2", "k"), (ep, "", "C3"), (eu, "", "C0")):
        ax.plot([st.u], [st.p], "o", color=col, ms=5)
        if lab:
            ax.annotate(lab, (st.u, st.p), (6, -4 if lab == "1" else 4), textcoords="offset points", fontsize=11)
ax.plot([s0.u], [s0.p], "ko", ms=6)
ax.annotate("0", (s0.u, s0.p), (6, -12), textcoords="offset points", fontsize=11)
ax.text(0.37, 1.215, "compressione incidente", fontsize=9)
ax.text(0.22, 0.79, "espansione incidente", fontsize=9)
ax.set_xlim(0.15, 0.48)
ax.set_ylim(0.76, 1.25)
ax.set_xlabel(r"$u/a_0$")
ax.set_ylabel(r"$p/p_0$")
ax.legend(loc="lower right", fontsize=8.5)
ax.set_title(r"Piano $p$–$u$ (esempio $M_e=0{,}3$, $p_i/p_0=1{,}1$ e $0{,}9$)")
fig.tight_layout()
fig.savefig(os.path.join(FIG, "pu.pdf"))
plt.close(fig)
print("=" * 78)
print("Figure salvate in", FIG)
