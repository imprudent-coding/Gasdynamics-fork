"""Prova del 01/01/1975, esercizio 2: relazioni di salto per un flusso 1D non stazionario.

Dalle relazioni di salto di massa, quantita' di moto ed energia totale si ricavano quelle per
l'energia cinetica e per l'energia interna; quella per l'energia interna si riscrive con la
pressione come variabile dipendente; si mostra che il salto dell'energia cinetica e' il salto della
quantita' di moto moltiplicato per la velocita' media <u>.

Riproduce tutti i numeri e le figure dello svolgimento.
Uso:  python calcoli.py   (dalla cartella dello svolgimento)

Il problema e' letterale. I numeri servono solo come controllo e per i diagrammi: urto u+a con
Ms = 2 in aria (gamma = 1,4) in quiete, p1 = 1, a1 = 1 (rho1 = gamma), cioe' l'urto generato da un
pistone impulsivo che si muove a u2 = 1,25.

Parti fuori modulo (scritte qui riusando gasdyn):
  - i residui delle relazioni di salto derivate (cinetica, interna, pressione, Hugoniot, entalpia
    totale relativa), calcolati sugli stati monte/valle di gd.urto_mobile;
  - la curva di Hugoniot nel piano p-tau, l'isentropica e la ripartizione del lavoro del pistone
    (formule chiuse, controllate con gd.urto_mobile e gd.salto_entropia);
  - il diagramma x-t con il volume di controllo solidale con la discontinuita'.
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

QUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(QUI, "..", ".."))
import gasdyn as gd  # noqa: E402

FIG = os.path.join(QUI, "figure")
os.makedirs(FIG, exist_ok=True)
plt.rcParams.update({"font.size": 10, "axes.grid": True, "grid.alpha": 0.3})

G = gd.Gas(1.4)
g = G.gamma


def salto(a1, a2):
    """[q] = q_valle - q_monte (F §2, Rankine-Hugoniot generale)."""
    return a2 - a1


def media(a1, a2):
    """<q> = (q1 + q2)/2, come in S p.18 ([pu] = <p>[u] + <u>[p])."""
    return 0.5 * (a1 + a2)


def ei(st):
    """Energia interna specifica e_i = p/((gamma-1) rho) (F §1)."""
    return st.p / ((st.gas.gamma - 1) * st.rho)


def residui(mon, val, w):
    """Residui di tutte le relazioni di salto, nel riferimento solidale (v = u - w)."""
    gm = mon.gas.gamma
    v1, v2 = mon.u - w, val.u - w
    m = mon.rho * v1                                  # flusso di massa relativo (F §2)
    k1, k2 = 0.5 * mon.u**2, 0.5 * val.u**2
    e1, e2 = ei(mon) + k1, ei(val) + k2
    du, dp = salto(mon.u, val.u), salto(mon.p, val.p)
    um, pm = media(mon.u, val.u), media(mon.p, val.p)
    tau1, tau2 = 1 / mon.rho, 1 / val.rho
    r = {
        "massa      [m]=0": val.rho * v2 - m,
        "qdm        m[u]+[p]=0": m * du + dp,
        "energia    m[e]+[pu]=0": m * (e2 - e1) + salto(mon.p * mon.u, val.p * val.u),
        "cinetica   m[u^2/2]+<u>[p]=0": m * (k2 - k1) + um * dp,
        "interna    m[e_i]+<p>[u]=0": m * (ei(val) - ei(mon)) + pm * du,
        "pressione  -w[p]+<u>[p]+g<p>[u]=0": -w * dp + um * dp + gm * pm * du,
        "cin. rel.  m[v^2/2]+<v>[p]=0": m * (0.5 * v2**2 - 0.5 * v1**2) + media(v1, v2) * dp,
    }
    if abs(m) > 1e-12:
        r["Hugoniot   [e_i]+<p>[tau]=0"] = ei(val) - ei(mon) + pm * (tau2 - tau1)
        r["Hug. in p  <tau>[p]+g<p>[tau]=0"] = media(tau1, tau2) * dp + gm * pm * (tau2 - tau1)
        r["Rayleigh   [p]+m^2[tau]=0"] = dp + m**2 * (tau2 - tau1)
        r["entalpia   [h+v^2/2]=0"] = salto(mon.h + 0.5 * v1**2, val.h + 0.5 * v2**2)
    return r


def hugoniot(x, gm=g):
    """p2/p1 sulla curva di Hugoniot in funzione di x = tau2/tau1 (da <tau>[p] + g<p>[tau] = 0)."""
    return ((gm + 1) - (gm - 1) * x) / ((gm + 1) * x - (gm - 1))


# =============================================================================
# 1. Esempio: urto u+a con Ms = 2 in aria in quiete
# =============================================================================
q = gd.State.da_pa(G, 1.0, 1.0, nome="1")
urto = gd.urto_mobile(q, "u+a", Ms=2.0)
s1, s2, w = urto.monte, urto.valle, urto.w
print("=" * 72)
print("ESEMPIO (urto u+a, Ms = 2, gas in quiete, p1 = a1 = 1)")
print(urto)
print("residui R-H del modulo:", urto.residui_rh())
m = s1.rho * (s1.u - w)
du, dp = s2.u - s1.u, s2.p - s1.p
um, pm = media(s1.u, s2.u), media(s1.p, s2.p)
print(f"m = rho1 (u1 - w) = {m:.4f}   [u] = {du:.4f}   [p] = {dp:.4f}   <u> = {um:.4f}   <p> = {pm:.4f}")
print(f"tau1 = {1/s1.rho:.4f}   tau2 = {1/s2.rho:.4f}   rho2/rho1 = {s2.rho/s1.rho:.4f}")
dk = 0.5 * s2.u**2 - 0.5 * s1.u**2
dei = ei(s2) - ei(s1)
print(f"[u^2/2] = {dk:.4f}   m[u^2/2] = {m*dk:.4f}   <u>[p] = {um*dp:.4f}")
print(f"[e_i]   = {dei:.4f}   m[e_i]   = {m*dei:.4f}   <p>[u] = {pm*du:.4f}")
print(f"[pu] = {s2.p*s2.u - s1.p*s1.u:.4f} = <u>[p] + <p>[u] = {um*dp + pm*du:.4f}")
print(f"[pu] - m[u^2/2]... sorgente <p>[u] della cinetica = {pm*du:.4f}")
print(f"[e_i]/[u^2/2] = {dei/dk:.4f}   (p1+p2)/(p2-p1) = {(s1.p+s2.p)/(s2.p-s1.p):.4f}")
for k, v in residui(s1, s2, w).items():
    print(f"  {k:38s} residuo = {v:+.2e}")

# Lavoro compiuto dal pistone e sua ripartizione (per unita' di area e di tempo)
L = s2.p * s2.u
print(f"lavoro del pistone p2 u2 = {L:.4f}: cinetica |m|[u^2/2] = {abs(m)*dk:.4f} "
      f"({abs(m)*dk/L:.1%}), interna |m|[e_i] = {abs(m)*dei:.4f} ({abs(m)*dei/L:.1%})")

# Confronto con l'isentropica alla stessa densita'
tau1, tau2 = 1 / s1.rho, 1 / s2.rho
p_is = s1.p * (tau1 / tau2)**g
ei_is = p_is * tau2 / (g - 1)
lav_rev = ei_is - ei(s1)                   # = -int p dtau lungo l'isentropica
print(f"isentropica alla stessa densita': p = {p_is:.4f}, [e_i]_is = {lav_rev:.4f}; "
      f"Hugoniot [e_i] = {dei:.4f}; in piu' (riscaldamento a rho costante) = {dei - lav_rev:.4f}")
print(f"  quota irreversibile del lavoro del pistone = {abs(m)*(dei-lav_rev)/L:.1%}")
print("salto di entropia:", gd.salto_entropia(s1, s2))

# =============================================================================
# 2. Invarianza galileiana: lo stesso urto visto da un riferimento che si muove a -U0
# =============================================================================
U0 = 0.7
s1g, s2g, wg = s1.con(u=s1.u + U0), s2.con(u=s2.u + U0), w + U0
print("=" * 72)
print(f"STESSO URTO CON TUTTE LE VELOCITA' AUMENTATE DI {U0}")
dkg = 0.5 * s2g.u**2 - 0.5 * s1g.u**2
print(f"[u^2/2]: {dk:.4f} -> {dkg:.4f} (cambia);  [e_i]: {dei:.4f} -> {ei(s2g)-ei(s1g):.4f} (invariato)")
for k, v in residui(s1g, s2g, wg).items():
    print(f"  {k:38s} residuo = {v:+.2e}")

# =============================================================================
# 3. Contatto (m = 0): [p] = 0 e [u] = 0, [rho] e [e_i] liberi
# =============================================================================
print("=" * 72)
c1 = gd.State(G, 1.4, 0.5, 1.0, nome="1")
c2 = gd.State(G, 0.35, 0.5, 1.0, nome="2")
print("CONTATTO (rho1 = 1,4, rho2 = 0,35, u = 0,5, p = 1): [e_i] =", ei(c2) - ei(c1))
for k, v in residui(c1, c2, 0.5).items():
    print(f"  {k:38s} residuo = {v:+.2e}")

# =============================================================================
# 4. Hugoniot contro isentropica: entropia del terzo ordine e limite di urto forte
# =============================================================================
print("=" * 72)
print("HUGONIOT: controllo con gd.urto_mobile e ordine del salto di entropia")
for Ms in [1.01, 1.05, 1.1, 1.2, 1.5, 2.0, 5.0, 20.0]:
    u = gd.urto_mobile(q, "u+a", Ms=Ms)
    x = u.monte.rho / u.valle.rho
    ds = gd.salto_entropia(u.monte, u.valle)["ds_R"]
    dpp = u.p2_p1 - 1
    print(f"Ms={Ms:5.2f}  rho2/rho1={1/x:6.3f}  p2/p1={u.p2_p1:8.3f}  Hugoniot={hugoniot(x):8.3f}  "
          f"isentropica={x**-g:8.3f}  ds/R={ds:.3e}  ds/R/(dp/p)^3={ds/dpp**3:.4f}")
print(f"(gamma+1)/(12 gamma^2) = {(g+1)/(12*g*g):.4f} (limite per urti deboli)")
print(f"limite di urto forte rho2/rho1 -> (g+1)/(g-1) = {(g+1)/(g-1):.3f}")

# =============================================================================
# 5. Figure
# =============================================================================
# 5a. Piano x-t con il volume di controllo solidale con l'urto
fig, ax = plt.subplots(figsize=(6.0, 4.2))
T = 1.0
ax.plot([0, w * T], [0, T], "r-", lw=2.2, label="discontinuità $x=wt$ (urto)")
ax.plot([0, s2.u * T], [0, T], "k--", lw=1.4, label="pistone / traiettoria $dx/dt=u_2$")
for x0 in [0.4, 0.8, 1.2, 1.6]:
    tc = x0 / w                                   # istante d'arrivo dell'urto (u1 = 0)
    ax.plot([x0, x0], [0, tc], color="0.45", lw=0.9)
    ax.plot([x0, x0 + s2.u * (T - tc)], [tc, T], color="0.45", lw=0.9)
ax.plot([], [], color="0.45", lw=0.9, label="traiettorie delle particelle")
# volume di controllo: xi in [-0.25, 0.25] intorno all'urto, tau in [0.35, 0.75]
ta, tb, dxi = 0.35, 0.75, 0.25
P = np.array([[w * ta - dxi, ta], [w * ta + dxi, ta], [w * tb + dxi, tb], [w * tb - dxi, tb],
              [w * ta - dxi, ta]])
ax.fill(P[:, 0], P[:, 1], color="tab:blue", alpha=0.12)
ax.plot(P[:, 0], P[:, 1], color="tab:blue", lw=1.2, label="volume di controllo solidale")
ax.text(1.75, 0.35, "1 (monte)", fontsize=12, ha="center")
ax.text(1.0, 0.85, "2 (valle)", fontsize=12, ha="center")
ax.text(w * 0.55 + 0.27, 0.55, r"$\xi_1$", fontsize=10, color="tab:blue")
ax.text(w * 0.55 - 0.47, 0.55, r"$\xi_2$", fontsize=10, color="tab:blue")
ax.set_xlim(0, 2.2)
ax.set_ylim(0, T)
ax.set_xlabel("$x$")
ax.set_ylabel("$t$")
ax.legend(loc="upper left", fontsize=8)
ax.set_title(r"Piano $x$–$t$ (esempio: urto $u+a$, $M_s=2$, gas 1 in quiete)", fontsize=10)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "xt.pdf"))
plt.close(fig)

# 5b. Piano p-tau: Hugoniot, isentropica, retta di Rayleigh, trapezio [e_i]
x2 = tau2 / tau1
xl = (g - 1) / (g + 1)
fig, axs = plt.subplots(1, 2, figsize=(10.0, 4.2))
for ax, scala in zip(axs, ["lin", "log"]):
    xh = np.linspace(xl + 1e-3, 1.0, 600)
    xe = np.linspace(1.0, 1.6, 200)
    ax.plot(xh, hugoniot(xh), "r-", lw=1.8, label="Hugoniot (compressione)")
    ax.plot(xe, hugoniot(xe), "r:", lw=1.2, label="Hugoniot, ramo di espansione (escluso)")
    xi = np.linspace(xl + 1e-3 if scala == "log" else 0.22, 1.6, 600)
    ax.plot(xi, xi**-g, "b--", lw=1.4, label="isentropica per 1")
    ax.plot([x2, 1], [s2.p, 1], "k-", lw=1.0, label=r"retta di Rayleigh, pendenza $-m^2$")
    if scala == "lin":                            # l'area ha senso solo in scala lineare
        ax.fill([x2, 1, 1, x2], [0, 0, 1, s2.p], color="tab:orange", alpha=0.25,
                label=r"$[e_i]=-\langle p\rangle[\tau]$ (trapezio)")
    ax.plot([1, x2], [1, s2.p], "ko", ms=5)
    ax.annotate("1", (1, 1), xytext=(1.02, 1.15), fontsize=11)
    ax.annotate("2", (x2, s2.p), xytext=(x2 + 0.02, s2.p + 0.15 if scala == "lin" else s2.p * 1.1),
                fontsize=11)
    ax.plot([x2], [p_is], "bs", ms=4)
    ax.annotate(r"$2_s$", (x2, p_is), xytext=(x2 - 0.09, p_is * (0.75 if scala == "log" else 0.93)),
                fontsize=10, color="b")
    ax.axvline(xl, color="0.4", lw=0.8, ls="-.")
    ax.set_xlabel(r"$\tau/\tau_1=\rho_1/\rho$")
    ax.set_ylabel(r"$p/p_1$")
    if scala == "lin":
        ax.set_xlim(0.1, 1.6)
        ax.set_ylim(0, 7)
        ax.set_title("scala lineare", fontsize=10)
        ax.legend(fontsize=7.5, loc="upper right")
    else:
        ax.set_yscale("log")
        ax.set_xlim(0.1, 1.6)
        ax.set_ylim(0.5, 300)
        ax.text(xl + 0.01, 120, r"$\frac{\gamma-1}{\gamma+1}$", fontsize=11, color="0.3")
        ax.set_title("scala logaritmica", fontsize=10)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "p_tau.pdf"))
plt.close(fig)

# 5c. Ripartizione del lavoro del pistone (gas 1 in quiete) in funzione di Ms
Msv = np.linspace(1.0005, 6.0, 300)
f_cin, f_rev, f_irr = [], [], []
for Ms in Msv:
    u = gd.urto_mobile(q, "u+a", Ms=Ms)
    a_, b_ = u.monte, u.valle
    Lp = b_.p * b_.u
    mm = abs(a_.rho * (a_.u - u.w))
    k = 0.5 * b_.u**2
    e = ei(b_) - ei(a_)
    t1, t2 = 1 / a_.rho, 1 / b_.rho
    rev = a_.p * (t1 / t2)**g * t2 / (g - 1) - ei(a_)
    f_cin.append(mm * k / Lp)
    f_rev.append(mm * rev / Lp)
    f_irr.append(mm * (e - rev) / Lp)
f_cin, f_rev, f_irr = map(np.array, (f_cin, f_rev, f_irr))
assert np.allclose(f_cin + f_rev + f_irr, 1.0)
assert np.allclose(f_rev + f_irr, [(1 + gd.urto_mobile(q, "u+a", Ms=M).p2_p1)
                                   / (2 * gd.urto_mobile(q, "u+a", Ms=M).p2_p1) for M in Msv])
fig, ax = plt.subplots(figsize=(6.0, 3.8))
ax.fill_between(Msv, 0, f_cin, color="tab:green", alpha=0.45, label="energia cinetica")
ax.fill_between(Msv, f_cin, f_cin + f_rev, color="tab:blue", alpha=0.35,
                label="energia interna: compressione isentropica")
ax.fill_between(Msv, f_cin + f_rev, 1, color="tab:red", alpha=0.35,
                label="energia interna: riscaldamento a $\\rho$ costante")
ax.axhline(0.5, color="k", lw=0.8, ls=":")
ax.axvline(2.0, color="0.3", lw=0.8, ls="--")
ax.set_xlim(1, 6)
ax.set_ylim(0, 1)
ax.set_xlabel("$M_s$")
ax.set_ylabel(r"frazione del lavoro $p_2u_2$")
ax.legend(fontsize=8, loc="lower right")
ax.set_title("Gas 1 in quiete: dove va il lavoro del pistone", fontsize=10)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "ripartizione.pdf"))
plt.close(fig)
i2 = np.argmin(abs(Msv - 2.0))
print("=" * 72)
print(f"ripartizione a Ms=2: cin {f_cin[i2]:.3f}, rev {f_rev[i2]:.3f}, irr {f_irr[i2]:.3f}; "
      f"a Ms=6: cin {f_cin[-1]:.3f}, rev {f_rev[-1]:.3f}, irr {f_irr[-1]:.3f}")
print("figure scritte in", FIG)
