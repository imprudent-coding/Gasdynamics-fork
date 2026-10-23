"""Prova del 01/02/1970, esercizio 1: riflessione di un urto obliquo su una lastra piana con corrente a M = 1,3
(stessa traccia di B69-E01, 25/01/2023).

Riproduce tutti i numeri e le figure dello svolgimento.
Uso:  python calcoli.py   (dalla cartella dello svolgimento)

Convenzioni (gasdyn_guida.md §1): angoli in gradi, θ antiorario dall'asse x (verso della corrente 1, parallela
alla lastra). L'urto incidente viene dall'alto: famiglia C−, devia la corrente di −θ_w verso la lastra; l'urto
riflesso nasce dalla lastra (parete sotto la corrente): famiglia C+, riporta θ = 0. Stati: 1 monte, 2 dopo
l'incidente, 3 dopo il riflesso, 4 dietro il gambetto (metodo_docente §b). Grandezze adimensionali: p1 = 1, a1 = 1.
La traccia non dà la deviazione θ_w dell'urto incidente: è il parametro della discussione (0 < θ_w < θmax(M1)).
Gas ideale: aria, γ = 1,4 (non indicato nella traccia).
"""
import os
import sys
import math

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
G = gd.Gas(1.4)
M1 = 1.3
s1 = gd.State.da_M(G, M1, p=1.0, a=1.0, nome="1")
TW_ES = 2.0          # deviazione dell'esempio di riflessione regolare
TW_IRR = 5.0         # deviazione dell'esempio oltre il distacco


def pn(M):
    """p2/p1 dell'urto normale (F §12 con σ = 90°)."""
    return gd.urto_normale(M, G)["p2_p1"]


def p_forte_1(theta):
    """Pressione sul ramo forte della polare di 1 alla deviazione |θ| (polare simmetrica)."""
    t = abs(theta)
    if t < 1e-9:
        return pn(M1)
    return gd.urto_obliquo(s1, theta=t, famiglia="C+", ramo="forte").valle.p


# =============================================================================
# 1. POLARE DELLO STATO 1 (F §12–13)
# =============================================================================
mu1 = gd.angolo_mach(M1)
tmax1, smax1 = gd.theta_max(M1, G)
ts1, ss1 = gd.sonico(M1, G)
un1 = gd.urto_normale(M1, G)
pmax1 = 1 + 2 * G.gamma / (G.gamma + 1) * (M1**2 * math.sin(math.radians(smax1))**2 - 1)
ps1 = 1 + 2 * G.gamma / (G.gamma + 1) * (M1**2 * math.sin(math.radians(ss1))**2 - 1)
print("=== Polare dello stato 1 (M1 = 1.3, γ = 1.4) ===")
print(f"  μ1 = {mu1:.2f}°  (onda nulla, p/p1 = 1)")
print(f"  θmax = {tmax1:.3f}° a σ = {smax1:.2f}°, p/p1 = {pmax1:.4f}")
print(f"  punto sonico: θs = {ts1:.3f}° a σ = {ss1:.2f}°, p/p1 = {ps1:.4f}")
print(f"  urto normale: p/p1 = {un1['p2_p1']:.4f}, M = {un1['M2']:.4f}, p02/p01 = {un1['p02_p01']:.4f}")

# =============================================================================
# 2. RIFLESSIONE REGOLARE DI ESEMPIO, θ_w = 2° (F §12, §14)
# =============================================================================
print(f"\n=== Riflessione regolare, θ_w = {TW_ES}° ===")
rr = gd.riflessione_regolare(s1, TW_ES)
inc, rif = rr["incidente"], rr["riflesso"]
st2, st3 = rr["2"], rr["3"]
tmax2 = gd.theta_max(st2.M, G)[0]
print(f"  incidente: σi = {inc.sigma:.2f}°, Mn1 = {inc.Mn1:.4f}, p2/p1 = {st2.p:.4f}, M2 = {st2.M:.4f},"
      f" θmax(M2) = {tmax2:.3f}°, Δs/R = {gd.salto_entropia(s1, st2)['ds_R']:.2e}")
print(f"  riflesso:  σr = {rif.sigma:.2f}° (rispetto alla corrente 2), angolo con la lastra φr = "
      f"{rr['angolo_riflesso_parete']:.2f}°, Mn2 = {rif.Mn1:.4f}")
print(f"             p3/p2 = {st3.p / st2.p:.4f}, p3/p1 = {st3.p:.4f}, M3 = {st3.M:.4f}, θ3 = {st3.theta:.2e}°,"
      f" Δs/R = {gd.salto_entropia(st2, st3)['ds_R']:.2e}")
print(f"  confronto: p2/p1 = {st2.p:.4f} < p3/p2 = {st3.p / st2.p:.4f} (riflesso più forte),"
      f" σi = {inc.sigma:.2f}° ≠ φr = {rr['angolo_riflesso_parete']:.2f}°")
rr_f = gd.riflessione_regolare(s1, TW_ES, ramo_riflesso="forte")
print(f"  RR con riflesso forte: p3'/p1 = {rr_f['3'].p:.4f}, M3' = {rr_f['3'].M:.4f}, "
      f"φr' = {rr_f['angolo_riflesso_parete']:.2f}°  (solo con contropressione imposta)")

# Verifica di limite acustico (S: onde deboli): Δp/p ≈ γM²θ/√(M²−1) per onda; RR ⇒ p3 − 1 ≈ 2(p2 − 1).
tw_small = 0.5
rr_s = gd.riflessione_regolare(s1, tw_small)
lin = G.gamma * M1**2 / math.sqrt(M1**2 - 1) * math.radians(tw_small)
print(f"  limite debole θ_w = {tw_small}°: p2 − 1 = {rr_s['2'].p - 1:.4f} (lineare {lin:.4f}), "
      f"p3 − 1 = {rr_s['3'].p - 1:.4f} ≈ 2(p2 − 1) = {2 * (rr_s['2'].p - 1):.4f}")

# =============================================================================
# 3. CRITERI: DISTACCO, SONICO, VON NEUMANN (F §14)
# =============================================================================
print("\n=== Criteri di transizione ===")
tD = gd.theta_distacco(M1, G)
inc_D = gd.urto_obliquo(s1, theta=tD, famiglia="C-")
rr_D = gd.riflessione_regolare(s1, tD * (1 - 1e-7))
print(f"  distacco θD = {tD:.3f}°: M2 = {inc_D.valle.M:.4f}, θmax(M2) = {gd.theta_max(inc_D.valle.M, G)[0]:.3f}°,"
      f" σi = {inc_D.sigma:.2f}°, p2 = {inc_D.valle.p:.4f}, p3 = {rr_D['3'].p:.4f}, M3 = {rr_D['3'].M:.4f},"
      f" φr = {rr_D['angolo_riflesso_parete']:.2f}°")
# Criterio sonico (fuori modulo: radice di M3(θ_w) = 1 con la RR del modulo)
tS = brentq(lambda tw: gd.riflessione_regolare(s1, tw)["3"].M - 1.0, 1.0, tD * (1 - 1e-9), xtol=1e-12)
print(f"  sonico (M3 = 1) θS = {tS:.3f}°, p3 = {gd.riflessione_regolare(s1, tS)['3'].p:.4f}")
print(f"  von Neumann: p3(θ_w) massimo (a θD) = {rr_D['3'].p:.4f} < p normale = {un1['p2_p1']:.4f} ⇒ θN non esiste")
try:
    gd.theta_von_neumann(M1, G)
except gd.NessunaSoluzione as e:
    print("  theta_von_neumann:", e)

print("\n  Tabella RR (riflesso debole):")
print("   θw     σi     p2/p1   M2     θmax(M2)  φr     p3/p1   M3")
tab = [1.0, 2.0, 3.0, tS, tD * (1 - 1e-7)]
for tw in tab:
    r = gd.riflessione_regolare(s1, tw)
    print(f"  {tw:5.2f}  {r['incidente'].sigma:6.2f}  {r['2'].p:.4f}  {r['2'].M:.4f}  {gd.theta_max(r['2'].M, G)[0]:6.3f}"
          f"   {r['angolo_riflesso_parete']:6.2f}  {r['3'].p:.4f}  {r['3'].M:.4f}")

# =============================================================================
# 4. OLTRE IL DISTACCO: TEORIA DEI TRE URTI (F §14)
# =============================================================================
print(f"\n=== Teoria dei tre urti, θD < θ_w < θs(M1) ===")
n_mr = 0
grid = np.linspace(0.1, ts1 - 0.01, 80)
for tw in grid:
    try:
        gd.riflessione_mach(s1, tw)
        n_mr += 1
    except gd.NessunaSoluzione:
        pass
print(f"  θ_w su {len(grid)} valori in (0, θs): soluzioni a tre urti trovate = {n_mr}")
# Polare riflessa tutta interna alla polare di 1: confronto punto per punto
marg = []
for tw in grid:
    st2_ = gd.urto_obliquo(s1, theta=tw, famiglia="C-").valle
    pol = gd.polare_p_theta(st2_, "C+", n=200, espansione=False)["urto"]
    marg.append(min(p_forte_1(t) - p for t, p in zip(pol["theta"], pol["p"])))
print(f"  min [p_forte,1(θ) − p_riflessa(θ)] su tutti i θ_w = {min(marg):.4f} > 0 ⇒ polare riflessa interna")
inc_I = gd.urto_obliquo(s1, theta=TW_IRR, famiglia="C-")
st2_I = inc_I.valle
tmax2_I = gd.theta_max(st2_I.M, G)[0]
print(f"  esempio θ_w = {TW_IRR}°: σi = {inc_I.sigma:.2f}°, M2 = {st2_I.M:.4f}, p2 = {st2_I.p:.4f},"
      f" θmax(M2) = {tmax2_I:.3f}° ⇒ polare riflessa in θ ∈ [{-TW_IRR - tmax2_I:.2f}, {-TW_IRR + tmax2_I:.2f}]°,"
      f" p max = {st2_I.p * pn(st2_I.M):.4f}; polare 1 (ramo forte) a θ = −{TW_IRR}°: p = {p_forte_1(TW_IRR):.4f}")
print(f"  θs(M1) = {ts1:.3f}° < θ_w < θmax(M1) = {tmax1:.3f}°: M2 < 1, nessun urto riflesso possibile")

# =============================================================================
# 5. FIGURE
# =============================================================================
C_INC, C_RIF = "#1f4e9c", "#c0392b"


def disegna_polare(ax, st, fam, colore, lab, lw=1.6, entrambe=False, ls="-"):
    fams = ["C+", "C-"] if entrambe else [fam]
    for i, f in enumerate(fams):
        u = gd.polare_p_theta(st, f, n=500, espansione=False)["urto"]
        ax.plot(u["theta"], u["p"], color=colore, lw=lw, ls=ls, label=lab if i == 0 else None)


def punto(ax, th, p, lab, dx=0.12, dy=0.015, c="k"):
    ax.plot(th, p, "o", color=c, ms=4.5, zorder=5)
    ax.annotate(lab, (th, p), (th + dx, p + dy), fontsize=10, color=c)


# --- Figura 1: piano p–θ, tre casi ---------------------------------------------------------
fig, axs = plt.subplots(1, 3, figsize=(13.5, 4.6), sharey=True)
casi = [(TW_ES, f"(a) $\\theta_I={TW_ES:.0f}^\\circ<\\theta_D$: riflessione regolare"),
        (tD, f"(b) $\\theta_I=\\theta_D={tD:.2f}^\\circ$: distacco".replace(".", "{,}")),
        (TW_IRR, f"(c) $\\theta_I={TW_IRR:.0f}^\\circ>\\theta_D$: nessuna RR né MR")]
for ax, (tw, tit) in zip(axs, casi):
    disegna_polare(ax, s1, None, C_INC, "polare di 1", entrambe=True)
    st2_ = gd.urto_obliquo(s1, theta=tw * (1 - 1e-7), famiglia="C-").valle
    disegna_polare(ax, st2_, None, C_RIF, "polare di 2", lw=1.2, entrambe=True, ls=":")
    disegna_polare(ax, st2_, "C+", C_RIF, "polare di 2, ramo $C^+$")
    ax.axvline(0, color="0.4", lw=0.8)
    punto(ax, 0, 1, "1", dx=0.15, dy=-0.05)
    punto(ax, -tw, st2_.p, "2", dx=-0.55, dy=0.0)
    ax.plot(0, un1["p2_p1"], "s", color=C_INC, ms=4)
    ax.annotate("urto normale\ndi 1", (0, un1["p2_p1"]), (0.4, un1["p2_p1"] + 0.02), fontsize=8, color=C_INC)
    if tw < tD:
        r = gd.riflessione_regolare(s1, tw)
        rf = gd.riflessione_regolare(s1, tw, ramo_riflesso="forte")
        punto(ax, 0, r["3"].p, "3", dx=0.15, dy=-0.04, c=C_RIF)
        punto(ax, 0, rf["3"].p, "3'", dx=0.15, dy=-0.04, c=C_RIF)
    elif abs(tw - tD) < 1e-6:
        punto(ax, 0, rr_D["3"].p, "3 (tangenza)", dx=0.15, dy=-0.04, c=C_RIF)
    ax.set_title(tit, fontsize=10)
    ax.set_xlabel(r"$\theta$ [°]")
    ax.set_xlim(-8.5, 8.5)
    ax.grid(alpha=0.3)
axs[0].set_ylabel(r"$p/p_1$")
axs[0].set_ylim(0.95, 1.9)
axs[0].legend(fontsize=8, loc="lower right")
fig.tight_layout()
fig.savefig(os.path.join(FIG, "p_theta.pdf"))
plt.close(fig)

# --- Figura 2: pressioni a parete in funzione di θ_w ---------------------------------------
tws = np.linspace(0.01, tD * (1 - 1e-6), 120)
p2s, p3s, p3f = [], [], []
for tw in tws:
    r = gd.riflessione_regolare(s1, tw)
    p2s.append(r["2"].p)
    p3s.append(r["3"].p)
    p3f.append(gd.riflessione_regolare(s1, tw, ramo_riflesso="forte")["3"].p)
fig, ax = plt.subplots(figsize=(6.4, 4.2))
ax.plot(tws, p2s, color=C_INC, label=r"$p_2/p_1$ (dietro l'incidente)")
ax.plot(tws, p3s, color=C_RIF, label=r"$p_3/p_1$ (RR, riflesso debole)")
ax.plot(tws, p3f, color=C_RIF, ls="--", label=r"$p_{3'}/p_1$ (RR, riflesso forte)")
ax.axhline(un1["p2_p1"], color="k", lw=0.9, ls=":", label=r"urto normale di 1 ($p/p_1=1{,}805$)")
ax.axvline(tS, color="0.5", lw=0.8, ls="-.")
ax.axvline(tD, color="0.3", lw=1.0)
ax.text(tS - 0.05, 1.02, r"$\theta_S$", ha="right", fontsize=9)
ax.text(tD + 0.05, 1.02, r"$\theta_D$", fontsize=9)
ax.set_xlim(0, 3.6)
ax.set_xlabel(r"$\theta_I$ [°]")
ax.set_ylabel(r"$p/p_1$")
ax.grid(alpha=0.3)
ax.legend(fontsize=8, loc="upper left")
fig.tight_layout()
fig.savefig(os.path.join(FIG, "p_parete.pdf"))
plt.close(fig)

# --- Figura 3: diagramma θ_w–M1 ------------------------------------------------------------
Ms = np.linspace(1.05, 4.0, 60)
tmx = [gd.theta_max(M, G)[0] for M in Ms]
tsn = [gd.sonico(M, G)[0] for M in Ms]
d = gd.diagramma_rr_mr(Ms, G)
fig, ax = plt.subplots(figsize=(6.4, 4.4))
ax.plot(Ms, tmx, color="0.2", lw=1.2, label=r"$\theta_{\max}(M_1)$: urto incidente staccato oltre")
ax.plot(Ms, tsn, color="0.2", lw=1.0, ls=":", label=r"$\theta_s(M_1)$: $M_2<1$ oltre")
ax.plot(Ms, d["theta_D"], color=C_RIF, lw=1.6, label=r"$\theta_D$ (distacco)")
ax.plot(Ms, d["theta_N"], color=C_INC, lw=1.6, label=r"$\theta_N$ (von Neumann)")
ax.axvline(M1, color="k", lw=0.8, ls="--")
ax.plot([M1, M1], [tD, tD], "o", color=C_RIF, ms=4)
ax.annotate(rf"$M_1=1{{,}}3$: $\theta_D={tD:.2f}^\circ$".replace(".", "{,}"), (M1, tD), (1.45, 1.0), fontsize=8,
            arrowprops=dict(arrowstyle="-", lw=0.6))
ax.text(2.3, 7.0, "RR", fontsize=10)
ax.text(1.75, 17.0, "riflessione\nirregolare", fontsize=9, ha="center")
ax.text(3.45, 22.4, "RR o MR", fontsize=7.5, rotation=6, ha="center")
ax.text(2.0, 33.0, "urto incidente staccato", fontsize=8, ha="center")
ax.set_xlabel(r"$M_1$")
ax.set_ylabel(r"$\theta_I$ [°]")
ax.set_xlim(1.0, 4.0)
ax.set_ylim(0, 40)
ax.grid(alpha=0.3)
ax.legend(fontsize=7.5, loc="lower right")
fig.tight_layout()
fig.savefig(os.path.join(FIG, "diagramma_M_theta.pdf"))
plt.close(fig)
print(f"\n  θN(M1) esiste da M1 ≈ {Ms[np.where(~np.isnan(d['theta_N']))[0][0]]:.2f} (griglia)")

# --- Figura 4: piano fisico ----------------------------------------------------------------
fig, axs = plt.subplots(1, 2, figsize=(13, 5.0))
H = 1.0
# (a) RR esatta a θ_w = 2°: punto di riflessione nell'origine
ax = axs[0]
si, fr = inc.sigma, rr["angolo_riflesso_parete"]
xa = -H / math.tan(math.radians(si))
xb = H / math.tan(math.radians(fr))
ax.fill_between([-1.6, 1.6], -0.08, 0, color="0.75")
ax.plot([-1.6, 1.6], [0, 0], "k", lw=1.5)
ax.plot([xa, 0], [H, 0], color=C_INC, lw=2.2)
ax.plot([0, xb], [0, H], color=C_RIF, lw=2.2)
# linee di corrente
for y0 in (0.25, 0.55, 0.85):
    xi = -y0 / math.tan(math.radians(si))
    yr = None
    # tratto 2: inclinato di −θ_w fino al riflesso
    k = -math.tan(math.radians(TW_ES))
    # intersezione con il riflesso y = x·tan(fr)
    xr = (y0 - k * xi) / (math.tan(math.radians(fr)) - k)
    yr = xr * math.tan(math.radians(fr))
    ax.plot([-1.6, xi], [y0, y0], color="0.45", lw=0.8)
    ax.plot([xi, xr], [y0, yr], color="0.45", lw=0.8)
    ax.plot([xr, 1.6], [yr, yr], color="0.45", lw=0.8)
ax.annotate("", (-1.15, 0.4), (-1.5, 0.4), arrowprops=dict(arrowstyle="->", lw=1.2))
ax.text(-1.45, 0.45, r"$M_1=1{,}3$", fontsize=10)
ax.text(-0.95, 0.12, "1", fontsize=13)
ax.text(-0.05, 0.62, "2", fontsize=13)
ax.text(1.05, 0.12, "3", fontsize=13)
ax.text(xa - 0.05, H + 0.05, r"urto incidente ($C^-$)", fontsize=9, color=C_INC, ha="center")
ax.text(xb + 0.05, H + 0.05, r"urto riflesso ($C^+$)", fontsize=9, color=C_RIF, ha="center")
ang = np.radians(np.linspace(180 - si, 180, 30))
ax.plot(0.22 * np.cos(ang), 0.22 * np.sin(ang), "k", lw=0.7)
ax.text(-0.42, 0.06, rf"$\sigma_i={si:.1f}^\circ$".replace(".", "{,}"), fontsize=9)
ang = np.radians(np.linspace(0, fr, 30))
ax.plot(0.22 * np.cos(ang), 0.22 * np.sin(ang), "k", lw=0.7)
ax.text(0.25, 0.06, rf"$\varphi_r={fr:.1f}^\circ$".replace(".", "{,}"), fontsize=9)
ax.text(0.0, -0.2, "lastra, $\\theta=0$", ha="center", fontsize=9)
ax.set_title(rf"(a) Riflessione regolare, $\theta_I={TW_ES:.0f}^\circ$ (angoli in scala)", fontsize=10)
ax.set_xlim(-1.6, 1.6)
ax.set_ylim(-0.3, 1.15)
ax.set_aspect("equal")
ax.axis("off")

# (b) riflessione irregolare (schema): urto incidente che si incurva fino a incidere normalmente la lastra
ax = axs[1]
ax.fill_between([-1.6, 1.6], -0.08, 0, color="0.75")
ax.plot([-1.6, 1.6], [0, 0], "k", lw=1.5)
si5 = inc_I.sigma
xa = -0.55 - H / math.tan(math.radians(si5))
# tratto rettilineo fino a y = 0.45, poi raccordo che arriva normale alla parete in x = 0
y_c = 0.62
x_c = xa + (H - y_c) / math.tan(math.radians(si5))
ax.plot([xa, x_c], [H, y_c], color=C_INC, lw=2.2)
tt = np.linspace(0, 1, 50)
# Bezier quadratica tangente al tratto rettilineo in (x_c, y_c) e verticale in (x_f, 0)
x_f = x_c + 0.16
P0 = np.array([x_c, y_c])
P2 = np.array([x_f, 0.0])
P1 = np.array([x_f, y_c - (x_f - x_c) * math.tan(math.radians(si5))])
B = ((1 - tt)[:, None] ** 2) * P0 + (2 * (1 - tt) * tt)[:, None] * P1 + (tt[:, None] ** 2) * P2
ax.plot(B[:, 0], B[:, 1], color=C_INC, lw=2.2)
# linea sonica (tratteggiata) che chiude la sacca subsonica dietro il tratto quasi normale
th = np.linspace(0, 1, 40)
xs = x_c + 0.02 + (x_f + 0.75 - x_c - 0.02) * th
ys = (y_c - 0.04) * (1 - th) ** 1.6
ax.plot(xs, ys, "k--", lw=0.9)
ax.annotate("", (-1.15, 0.25), (-1.5, 0.25), arrowprops=dict(arrowstyle="->", lw=1.2))
ax.text(-1.45, 0.3, r"$M_1=1{,}3$", fontsize=10)
ax.text(-1.1, 0.08, "1", fontsize=13)
ax.text(0.45, 0.75, "2", fontsize=13)
ax.text(x_f + 0.1, 0.1, "4", fontsize=13)
ax.text(x_f + 0.45, 0.3, "linea sonica", fontsize=8)
ax.text(x_c - 0.05, 0.98, "urto incidente", fontsize=9, color=C_INC, ha="right")
ax.text(x_f - 0.06, -0.2, "urto quasi normale\n a parete", fontsize=8, color=C_INC, ha="center")
ax.set_title(r"(b) $\theta_D<\theta_I<\theta_s(M_1)$: nessuna riflessione, urto curvo (schema)", fontsize=10)
ax.set_xlim(-1.6, 1.6)
ax.set_ylim(-0.3, 1.15)
ax.set_aspect("equal")
ax.axis("off")
fig.tight_layout()
fig.savefig(os.path.join(FIG, "campo.pdf"))
plt.close(fig)
print("\nFigure salvate in", FIG)
