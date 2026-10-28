"""Prova del 01/01/1975, esercizio 3: derivate nel tempo di pressione e velocita' al contorno
d'ingresso di un flusso 1D instazionario, subsonico e omoentropico, con temperatura totale
costante assegnata.

Riproduce tutti i numeri e le figure dello svolgimento.
Uso:  python calcoli.py   (dalla cartella dello svolgimento)

Il problema e' letterale (nessun dato numerico). I numeri servono per l'aria (gamma = 1,4) e per un
esempio illustrativo: ingresso a sinistra (x = 0), M1 = 0,3, grandezze totali p_0 = a_0 = 1, onda
incidente u-a con p2/p1 = 1,1 (compressione) e 0,9 (espansione).

Parti fuori modulo (scritte qui riusando gasdyn):
  - le formule chiuse delle derivate al contorno (u_t, p_t, variabile d'onda entrante), controllate
    sostituendole nelle tre equazioni di Eulero in forma caratteristica (S p.46) con stati casuali;
  - lo stato all'ingresso dopo la riflessione (R1 dall'interno + a^2 + delta u^2 = a_0^2): formula
    chiusa, controllata con gd.contorno_serbatoio(..., modello='isentropico');
  - il diagramma x-t con caratteristiche rettilinee (semirette) e il piano p-u con le curve del contorno.
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
g, d = G.gamma, G.delta                    # delta = (gamma-1)/2 (F §1)


def derivate_ingresso(rho, u, a, Y):
    """Derivate al contorno sinistro con s e T0 costanti (S p.46, F §3).

    Y = p_x - rho a u_x e' la variabile d'onda della u-a, uscente, nota dall'interno.
    Ritorna u_t, p_t, rho_t e la variabile d'onda entrante X = p_x + rho a u_x (u+a).
    """
    u_t = -(a - u) / (a + u) * Y / rho
    p_t = -rho * u * u_t                   # T0 e s costanti: dp = -rho u du
    rho_t = p_t / a**2                     # omoentropico
    X = ((a - u) / (a + u))**2 * Y
    return u_t, p_t, rho_t, X


def coeff(M):
    """Coefficiente di riflessione all'ingresso: dR2 = C dR1, C = -(1-M)/(1+M)."""
    return -(1 - M) / (1 + M)


# =============================================================================
# 1. CONTROLLO DELLE FORMULE NELLE EQUAZIONI DI EULERO IN FORMA CARATTERISTICA
# =============================================================================
print("=" * 78)
print("1. Derivate al contorno d'ingresso: controllo nelle equazioni (S p.46)")
rng = np.random.default_rng(1)
for _ in range(5):
    rho, a = rng.uniform(0.5, 2.0), rng.uniform(0.5, 2.0)
    u = rng.uniform(0.01, 0.99) * a
    Y = rng.uniform(-1, 1)
    Z = 0.0                                # omoentropico: rho_x - p_x/a^2 = 0
    u_t, p_t, rho_t, X = derivate_ingresso(rho, u, a, Y)
    e_rho = rho_t + u * Z + (u + a) / (2 * a**2) * X + (u - a) / (2 * a**2) * Y
    e_u = u_t + (u + a) / (2 * rho * a) * X - (u - a) / (2 * rho * a) * Y
    e_p = p_t + (u + a) / 2 * X + (u - a) / 2 * Y
    # condizione T0 costante: a^2 + delta u^2 = a_0^2; omoentropico: a_t = delta a p_t/(gamma p)
    a_t = d * a * p_t / (rho * a**2)
    e_T0 = a * a_t + d * u * u_t           # d(a^2 + delta u^2)/dt / 2
    res = max(abs(e_rho), abs(e_u), abs(e_p), abs(e_T0))
    assert res < 1e-12, res
print("   tre equazioni + condizione T0 = cost soddisfatte (residuo < 1e-12) per 5 stati casuali")

# stessa cosa con gli invarianti di Riemann: (R1)_t = (a-u)(R1)_x, a a_t + delta u u_t = 0
rho, a, u = 1.3, 1.1, 0.4
R1x = 0.37
Y = rho * a / d * R1x                      # p_x - rho a u_x = (rho a/delta)(a - delta u)_x
u_t, p_t, _, X = derivate_ingresso(rho, u, a, Y)
M = u / a
u_t_R = -(a - u) * R1x / (d * (1 + M))
a_t_R = M * (a - u) * R1x / (1 + M)
assert abs(u_t - u_t_R) < 1e-12 and abs(p_t - rho * a / d * a_t_R) < 1e-12
print("   stesse derivate con gli invarianti di Riemann: ok")
print("   (R2)_x/(R1)_x = ((1-M)/(1+M))^2, (R2)_t/(R1)_t = -(1-M)/(1+M)")

# =============================================================================
# 2. COEFFICIENTI PER L'ARIA
# =============================================================================
print("=" * 78)
print("2. Coefficienti all'ingresso in funzione di M (aria)")
print("   M      C=dR2/dR1  dp/dp_i=2M/(1+M)  du/du_i=2/(1+M)  gradienti ((1-M)/(1+M))^2  "
      "(dm/m)/(dp_i/p)")
TAB = (0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9)
for M in TAB:
    C = coeff(M)
    print(f"   {M:.2f}   {C:+.3f}      {2 * M / (1 + M):.3f}             {2 / (1 + M):.3f}"
          f"            {C**2:.3f}                    {-2 * (1 - M) / (g * M):+.2f}")

# =============================================================================
# 3. ESEMPIO: M1 = 0,3, p_0 = a_0 = 1; onda incidente u-a di compressione ed espansione
# =============================================================================
print("=" * 78)
M1 = 0.3
P0, A0 = 1.0, 1.0                          # grandezze totali (serbatoio a monte dell'ingresso)
s1 = gd.State.da_totali(G, M1, p0=P0, a0=A0, nome="1")
print(f"3. Esempio M1 = {M1}: p1 = {s1.p:.4f}, a1 = {s1.a:.4f}, u1 = {s1.u:.4f}, "
      f"rho1 = {s1.rho:.4f}  (T0, p0 = 1)")
C1 = coeff(M1)


def stato_ingresso(s2):
    """Stato 3 all'ingresso: a3 - delta u3 = R1_2 (dall'interno), a3^2 + delta u3^2 = a_0^2.

    delta(1+delta) u^2 + 2 delta R1 u + R1^2 - a_0^2 = 0, radice con u > 0.
    """
    R1 = s2.R1
    u3 = (-d * R1 + np.sqrt(d * (1 + d) * A0**2 - d * R1**2)) / (d * (1 + d))
    s3 = gd.onda_semplice(s2, "u+a", u=u3, nome="3")
    assert abs(s3.a**2 + d * s3.u**2 - A0**2) < 1e-12 and abs(s3.R1 - R1) < 1e-12
    # controllo con il contorno di serbatoio del modulo (onda interna isentropica)
    r = gd.contorno_serbatoio(s2, "sinistra", P0, a0=A0, modello="isentropico")
    assert r["regime"] == "ingresso" and abs(r["stato"].u - s3.u) < 1e-10 and abs(r["stato"].p - s3.p) < 1e-10
    return s3


ris = {}
for nome, pr in (("compressione", 1.10), ("espansione", 0.90)):
    s2 = gd.onda_semplice(s1, "u-a", p=pr * s1.p, nome="2")   # R2 uniforme nell'incidente
    s3 = stato_ingresso(s2)
    dR1, dR2 = s2.R1 - s1.R1, s3.R2 - s2.R2
    dpi, dpr = s2.p - s1.p, s3.p - s2.p
    print(f"   {nome}: p2 = {s2.p:.4f}, u2 = {s2.u:.4f}, a2 = {s2.a:.4f}, M2 = {s2.M:.3f}")
    print(f"      T0_2/T0 = {s2.a0**2 / A0**2:.4f}, p0_2/p0 = {s2.p0 / P0:.4f}  (nel condotto T0 e p0 non sono costanti)")
    print(f"      3: p3 = {s3.p:.4f}, u3 = {s3.u:.4f}, a3 = {s3.a:.4f}, M3 = {s3.Mu:.3f}, "
          f"T0_3/T0 = {s3.a0**2 / A0**2:.6f}, p0_3 = {s3.p0:.6f}")
    print(f"      dR2/dR1 = {dR2 / dR1:+.4f} (locale C(M1) = {C1:+.4f}),  dp_r/dp_i = {dpr / dpi:+.4f}")
    print(f"      (p3 - p1)/dp_i = {(s3.p - s1.p) / dpi:.4f} (acustico 2M/(1+M) = {2 * M1 / (1 + M1):.4f})")
    print(f"      (u3 - u1)/du_i = {(s3.u - s1.u) / (s2.u - s1.u):.4f} (acustico 2/(1+M) = {2 / (1 + M1):.4f})")
    dm = s3.rho * s3.u / (s1.rho * s1.u) - 1
    print(f"      portata: m3/m1 - 1 = {dm:+.4f}  (acustico -2(1-M)/(gamma M) dp_i/p1 = "
          f"{-2 * (1 - M1) / (g * M1) * dpi / s1.p:+.4f})")
    assert s3.u > 0 and s3.Mu < 1                      # resta ingresso subsonico
    assert abs(s3.p0 - P0) < 1e-12                     # p0 costante: omoentropico + T0
    # confronto: pressione imposta (p3 = p1) e velocita' imposta (u3 = u1)
    ep = gd.onda_semplice(s2, "u+a", p=s1.p, nome="p")
    eu = gd.onda_semplice(s2, "u+a", u=s1.u, nome="u")
    print(f"      p imposta: u = {ep.u:.4f};  u imposta: p = {eu.p:.4f}")
    ris[nome] = (s2, s3, ep, eu)

# limite acustico: onda incidente molto debole -> coefficienti locali
s2 = gd.onda_semplice(s1, "u-a", p=1.0001 * s1.p)
s3 = stato_ingresso(s2)
rC = (s3.R2 - s2.R2) / (s2.R1 - s1.R1)
print(f"   onda debole (p2/p1 = 1,0001): dR2/dR1 = {rC:+.5f}, C = {C1:+.5f}")
assert abs(rC - C1) < 1e-3

# compressione che annulla la velocita' all'ingresso (fine della validita' della condizione)
# u3 = 0 => a3 = a_0, R1_2 = a_0 => a2 = a_0 - delta u2 ... onda u-a da 1: a2 + delta u2 = R2_1
u2z = (s1.R2 - A0) / (2 * d)               # da a2 - delta u2 = a_0 e a2 + delta u2 = R2_1
s2z = gd.onda_semplice(s1, "u-a", u=u2z)
print(f"   compressione che porta u3 = 0: p2/p1 = {s2z.p / s1.p:.4f}  (oltre, il gas esce: resta solo p = p_0)")
# espansione che porta l'ingresso a M = 1: a3 = u3 = a* = a_0 sqrt(2/(g+1))
ast = A0 * np.sqrt(2 / (g + 1))
R1s = ast * (1 - d)
s2s = gd.onda_semplice(s1, "u-a", a=(s1.R2 + R1s) / 2)
print(f"   espansione che porta l'ingresso a M = 1: p2/p1 = {s2s.p / s1.p:.4f}, a* = {ast:.4f}")

# =============================================================================
# 4. FIGURE
# =============================================================================
# --- 4.1 coefficienti in funzione di M
fig, ax = plt.subplots(figsize=(6.4, 4.0))
M = np.linspace(0, 1, 200)
ax.plot(M, coeff(M), "k", lw=1.8, label=r"$T_0$, $s$ imposti: $-\frac{1-M}{1+M}$")
ax.plot(M, 2 * M / (1 + M), color="C1", lw=1.4, label=r"$\Delta p_{ingr}/\Delta p_i=\frac{2M}{1+M}$")
ax.axhline(1, color="C0", ls="--", lw=1.1, label=r"$u$ imposta (parete): $+1$")
ax.axhline(0, color="C2", ls="-.", lw=1.1, label=r"non riflettente: $0$")
ax.axhline(-1, color="C3", ls="--", lw=1.1, label=r"$p$ imposta: $-1$")
ax.plot([M1], [C1], "ko")
ax.annotate(rf"$M_1={M1}$: ${C1:.2f}$".replace(".", "{,}"), (M1, C1), (0.38, -0.85), fontsize=9,
            arrowprops=dict(arrowstyle="->", lw=0.8))
ax.set_xlim(0, 1.02)
ax.set_ylim(-1.15, 1.15)
ax.set_xlabel(r"$M$ all'ingresso")
ax.set_ylabel(r"$\Delta R_2^{rifl}/\Delta R_1^{inc}\ \ (=\Delta p_r/\Delta p_i$ acustico$)$")
ax.legend(loc="lower right", fontsize=8.3, framealpha=1)
ax.set_title(r"Riflessione all'ingresso con $T_0$ e $s$ costanti")
fig.tight_layout()
fig.savefig(os.path.join(FIG, "coeff_riflessione.pdf"))
plt.close(fig)

# --- 4.2 piano x-t (semirette), esempio di compressione incidente
s2, s3, _, _ = ris["compressione"]
th, tt = 1.0, 1.5                          # arrivo all'ingresso di testa e coda dell'incidente
lh, lt = s1.u - s1.a, s2.u - s2.a          # pendenze dx/dt dell'incidente u-a (< 0)
rh, rt = s2.u + s2.a, s3.u + s3.a          # pendenze della riflessa u+a (> 0)
xmax, tmax = 2.0, 3.0
fig, ax = plt.subplots(figsize=(6.0, 4.6))
ax.axvline(0, color="k", lw=2.5)
ax.text(-0.04, 2.75, "ingresso\n" + r"$T_0,\ s$", fontsize=9, va="top", ha="right")
for k, tk in enumerate(np.linspace(th, tt, 6)):
    lam = lh + (lt - lh) * k / 5
    lr = rh + (rt - rh) * k / 5
    lw = 1.6 if k in (0, 5) else 0.7
    ax.plot([xmax, 0], [tk + xmax / lam, tk], color="C0", lw=lw)
    ax.plot([0, xmax], [tk, tk + xmax / lr], color="C3", lw=lw)
tP = (rh * th - lt * tt) / (rh - lt)
xP = rh * (tP - th)
ax.plot([xP], [tP], "k.", ms=6)
ax.text(xP + 0.08, tP - 0.02, "P", fontsize=9)
ax.set_xlim(-0.45, xmax)
ax.set_ylim(0, tmax)


def ruota(pend):
    """Angolo sullo schermo di una retta con dx/dt = pend (per allineare le scritte)."""
    p1, p2 = ax.transData.transform([(0, 0), (pend, 1)])
    ang = np.degrees(np.arctan2(p2[1] - p1[1], p2[0] - p1[0]))
    return (ang + 90) % 180 - 90


ax.text(0.25, 0.3, "1", fontsize=12, weight="bold")
xl = 0.6
ax.text(xl, 0.5 * (th + xl / lh + tt + xl / lt), r"onda incidente $u-a$", color="C0", fontsize=9,
        rotation=ruota(lh), ha="center", va="center", bbox=dict(fc="w", ec="none", pad=0.5))
ax.text(1.55, 1.2, "2", fontsize=12, weight="bold")
ax.text(0.1, 1.2, "Z", fontsize=10)
xl = 0.95
ax.text(xl, 0.5 * (th + xl / rh + tt + xl / rt), r"onda riflessa $u+a$", color="C3", fontsize=9,
        rotation=ruota(rh), ha="center", va="center", bbox=dict(fc="w", ec="none", pad=0.5))
ax.text(0.2, 2.7, "3", fontsize=12, weight="bold")
ax.set_xlabel(r"$x$ (ingresso in $x=0$, condotto per $x>0$)")
ax.set_ylabel(r"$t$")
ax.set_title(r"Piano $x$–$t$ (esempio $M_1=0{,}3$, compressione incidente)")
fig.tight_layout()
fig.savefig(os.path.join(FIG, "xt.pdf"))
plt.close(fig)

# --- 4.3 piano p-u: compressione ed espansione incidenti, tre contorni a confronto
fig, ax = plt.subplots(figsize=(6.4, 4.8))
uu = np.linspace(0.12, 0.5, 300)
ax.plot(uu, P0 * (1 - d * uu**2 / A0**2)**(g / (g - 1)), "k", lw=1.8,
        label=r"contorno: $a^2+\delta u^2=a_0^2$, $s$ cost.")
ax.axhline(s1.p, color="C3", ls="--", lw=1, label=r"$p$ imposta ($p=p_1$)")
ax.axvline(s1.u, color="C0", ls="--", lw=1, label=r"$u$ imposta ($u=u_1$)")
a_p = lambda st, p: st.a * (p / st.p)**((g - 1) / (2 * g))
pp = np.linspace(0.8, 1.05, 300)
ax.plot(s1.u - (a_p(s1, pp) - s1.a) / d, pp, color="0.4", lw=1, label=r"onda $u-a$ da 1")
for nome, (s2, s3, ep, eu) in ris.items():
    pr = np.linspace(min(s2.p, s3.p, ep.p, eu.p) * 0.98, max(s2.p, s3.p, ep.p, eu.p) * 1.02, 200)
    ax.plot(s2.u + (a_p(s2, pr) - s2.a) / d, pr, color="0.6", lw=1, ls="-.")
    for st, lab, col in ((s2, "2", "k"), (s3, "3", "k"), (ep, "", "C3"), (eu, "", "C0")):
        ax.plot([st.u], [st.p], "o", color=col, ms=5)
        if lab:
            ax.annotate(lab, (st.u, st.p), (8, -3), textcoords="offset points", fontsize=11)
ax.plot([s1.u], [s1.p], "ko", ms=6)
ax.annotate("1", (s1.u, s1.p), (6, 4), textcoords="offset points", fontsize=11)
ax.text(0.245, 1.045, "compressione incidente", fontsize=9)
ax.text(0.385, 0.81, "espansione incidente", fontsize=9)
ax.set_xlim(0.12, 0.5)
ax.set_ylim(0.8, 1.06)
ax.set_xlabel(r"$u/a_0$")
ax.set_ylabel(r"$p/p_0$")
ax.legend(loc="lower left", fontsize=8.3)
ax.set_title(r"Piano $p$–$u$ (esempio $M_1=0{,}3$, $p_2/p_1=1{,}1$ e $0{,}9$)")
fig.tight_layout()
fig.savefig(os.path.join(FIG, "pu.pdf"))
plt.close(fig)
print("=" * 78)
print("Figure salvate in", FIG)
