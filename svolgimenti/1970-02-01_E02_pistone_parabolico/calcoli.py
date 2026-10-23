"""Prova del 01/02/1970, esercizio 2 (traccia identica a B69-E02, 25/01/2023): pistone a sinistra
di un condotto molto lungo con legge x = t^2 per 0 <= t <= 1 e x = 2t - 1 per t > 1, gas in quiete
con p = 1, rho = gamma (quindi a0 = 1 per ogni gas). Evoluzione del campo e dipendenza da gamma e R.

Riproduce tutti i numeri e le figure dello svolgimento.
Uso:  python calcoli.py   (dalla cartella dello svolgimento)

Parti fuori modulo (scritte qui riusando gasdyn):
  - inviluppo delle caratteristiche (forma analitica, controllata con pistone_legge);
  - stima della traiettoria dell'urto fra la formazione e l'assorbimento dell'ultima
    caratteristica di compressione: in ogni istante la velocita' dell'urto e' quella del problema
    di Riemann fra lo stato dell'onda semplice che lo raggiunge e il gas in quiete (pseudo-Riemann
    locale; trascura le onde riflesse). Controllo: nel limite di urto debole si ritrova la regola
    della media delle velocita' caratteristiche, che ha soluzione analitica;
  - istanti di arrivo delle onde stimati con caratteristiche rettilinee (come nella soluzione
    ufficiale del 14/09/2026).
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

QUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(QUI, "..", ".."))
import gasdyn as gd  # noqa: E402

FIG = os.path.join(QUI, "figure")
os.makedirs(FIG, exist_ok=True)
plt.rcParams.update({"font.size": 10, "axes.grid": True, "grid.alpha": 0.3})

VP = 2.0           # velocita' finale del pistone
T1 = 1.0           # fine della fase parabolica


def xp(t):
    t = np.asarray(t, float)
    return np.where(t <= T1, t**2, 2 * t - 1)


# =============================================================================
# 1. ANALISI PER UN GAS GENERICO (funzione di gamma), F §4, §6, §7, §8, §9
# =============================================================================
def analisi(gamma, stampa=False, nome_gas=""):
    G = gd.Gas(gamma)
    d = G.delta
    s0 = gd.State.da_pa(G, 1.0, 1.0, nome="0")          # p = 1, a = 1  =>  rho = gamma
    assert abs(s0.rho - gamma) < 1e-12
    r = {"gamma": gamma, "G": G, "s0": s0}

    # --- onda semplice u+a: R1 = a - delta u = a0 uniforme (F §4)
    # caratteristica emessa in tp: u = 2tp, a = 1 + 2 delta tp, pendenza 1 + (gamma+1) tp
    lam = lambda tp: 1 + (gamma + 1) * tp
    r["lam"] = lam

    # --- blow-up (F §6): t_BU = a0/((1+delta) xpp(0)); controllo con pistone_legge
    t_bu = 1.0 / ((1 + d) * 2.0)
    P = gd.pistone_legge(lambda t: t**2, s0, T1)
    assert abs(P.t_bu - t_bu) < 1e-8 and abs(P.x_bu - t_bu) < 1e-8
    r["t_bu"], r["x_bu"], r["P"] = t_bu, t_bu, P

    # --- stato 1: plateau dietro la compressione, pistone a velocita' costante VP
    s1 = gd.onda_semplice(s0, "u+a", u=VP, nome="1")
    assert abs(s1.a - gamma) < 1e-12                     # a1 = 1 + 2 delta = gamma
    r["s1"] = s1

    # --- pseudo-Riemann 1 | 0 (F §9) quando l'urto ha inglobato tutta la compressione
    sol = gd.interazione(s1.con(nome="1"), s0.con(nome="4"))
    assert sol.configurazione == "RCS"
    r["rp"] = sol
    s2, s3 = sol.stati["2"].con(nome="2"), sol.stati["3"].con(nome="3")
    r["s2"], r["s3"] = s2, s3
    r["w_rp"] = sol.onda_dx.w

    # --- stato 4: riflessione dell'espansione u-a sul pistone (u = VP, R1 portato dalla u-a)
    a4 = s2.a - d * (s2.u - VP)
    s4 = gd.isentropica(s2, a=a4).con(u=VP, nome="4")
    assert abs(s4.R1 - s2.R1) < 1e-12
    r["s4"] = s4

    # --- regime asintotico: urto da pistone impulsivo VP (F §7) + gas isentropico lato pistone
    imp = gd.pistone_impulsivo(s0, VP)
    w_inf = ((1 + d) * VP + np.sqrt((1 + d) ** 2 * VP**2 + 4)) / 2
    assert abs(imp.w - w_inf) < 1e-12
    s6 = imp.valle.con(nome="6")
    p_inf = 1 + s0.rho * w_inf * VP                     # quantita' di moto: p6 - p0 = rho0 w u6
    assert abs(s6.p - p_inf) < 1e-10
    s5 = gd.isentropica(s0, p=s6.p).con(u=VP, nome="5")   # stessa entropia del gas iniziale
    r["imp"], r["w_inf"], r["s5"], r["s6"] = imp, w_inf, s5, s6
    r["ds_rp"] = gd.salto_entropia(s0, s3)
    r["ds_inf"] = gd.salto_entropia(s0, s6)
    r["ds_5"] = gd.salto_entropia(s0, s5)

    if stampa:
        print(f"\n=== Gas {nome_gas}: gamma = {gamma:.4f}, delta = {d:.4f} ===")
        print(f"stato 0: p = 1, rho = {s0.rho:.4f}, a = {s0.a:.4f}, R1 = R2 = {s0.R1:.4f}")
        print(f"t_BU = x_BU = 1/(gamma+1) = {t_bu:.4f} (pistone_legge: {P.t_bu:.6f}, {P.x_bu:.6f}),"
              f" x_p(t_BU) = {t_bu**2:.4f}")
        print(f"stato 1: u = {s1.u:.3f}, a = {s1.a:.4f}, p = gamma^(2g/(g-1)) = {s1.p:.4f},"
              f" rho = {s1.rho:.4f}, T1/T0 = a1^2 = {s1.a**2:.4f},"
              f" M1 = {s1.M:.4f}, u-a = {s1.u - s1.a:.4f}, u+a = {s1.u + s1.a:.4f}")
        it = sol.iterazioni
        print("pseudo-Riemann 1|0, iterazioni (k, u, p2, p3):")
        for k in it:
            print(f"   k={k[0]}: u = {k[1]:.5f}, p2 = {k[2]:.4f}, p3 = {k[3]:.4f}, |p2-p3| = {abs(k[2]-k[3]):.2e}")
        print(sol)
        v = sol.velocita()
        print(f"   velocita': {v}")
        print(f"   Ms = {sol.onda_dx.Ms:.4f}, p3/p0 = {s3.p:.4f}, rho3/rho0 = {s3.rho/s0.rho:.4f},"
              f" T3/T0 = {s3.a**2:.4f}, s3/s0 = {r['ds_rp']['s2_s1']:.4f}")
        print(f"   T2/T0 = {s2.a**2:.4f}, rho2/rho0 = {s2.rho/s0.rho:.4f}  (impedenze rho a: "
              f"2 -> {s2.rho*s2.a:.3f}, 3 -> {s3.rho*s3.a:.3f})")
        print(f"stato 4 (pistone, dopo la riflessione): u = {s4.u:.3f}, a = {s4.a:.4f}, p = {s4.p:.4f},"
              f"  confronto p_inf = {s6.p:.4f}")
        print(f"asintotico: w_inf = {w_inf:.4f}, Ms = {imp.Ms:.4f}, p6 = 1 + rho0 w u = {s6.p:.4f},"
              f" rho6/rho0 = w/(w-u) = {s6.rho/s0.rho:.4f}, T6/T0 = {s6.a**2:.4f}, s6/s0 = {r['ds_inf']['s2_s1']:.4f}")
        print(f"            stato 5 (lato pistone): u = {s5.u}, p = {s5.p:.4f}, a = {s5.a:.4f},"
              f" rho5/rho0 = {s5.rho/s0.rho:.4f}, T5/T0 = {s5.a**2:.4f}")
        print(f"   Delta s/R: urto pseudo-RP {r['ds_rp']['ds_R']:.4f}, urto asintotico {r['ds_inf']['ds_R']:.4f}")
        print(f"   impedenze: 5 -> {s5.rho*s5.a:.3f}, 6 -> {s6.rho*s6.a:.3f}")
    return r


aria = analisi(1.4, stampa=True, nome_gas="aria")
G = aria["G"]
g, d = G.gamma, G.delta
lam = aria["lam"]
s0, s1, s2, s3, s4, s5, s6 = (aria[k] for k in ("s0", "s1", "s2", "s3", "s4", "s5", "s6"))

# =============================================================================
# 2. GEOMETRIA DELLE CARATTERISTICHE E INVILUPPO (aria)
# =============================================================================
# x = tp^2 + (1 + (g+1) tp)(t - tp);  dx/dtp = 0  =>  t_env = (1 + 2 g tp)/(g+1)
X = lambda tp, t: tp**2 + lam(tp) * (t - tp)
t_env = lambda tp: (1 + 2 * g * tp) / (g + 1)
print("\n=== Inviluppo (aria) ===")
for tp in (0.0, 0.5, 1.0):
    print(f"tp = {tp}: t_env = {t_env(tp):.4f}, x_env = {X(tp, t_env(tp)):.4f}")
# controllo numerico: due caratteristiche vicine si incontrano in t_env
for tp in (0.2, 0.7):
    h = 1e-6
    tc = (X(tp + h, 0) - X(tp, 0) - 0) / (lam(tp) - lam(tp + h))   # X(tp,t) = X(tp,0) + lam t
    assert abs(tc - t_env(tp)) < 1e-4
print("testa x = t raggiunta dal pistone x = t^2 in t = 1 (mai: l'urto si forma prima, t_BU < 1)")
print(f"coda: caratteristica da (1, 1) con pendenza u1 + a1 = 2 + gamma = {lam(1.0):.3f}")

# =============================================================================
# 3. STIMA DELLA TRAIETTORIA DELL'URTO (fuori modulo)
# =============================================================================
def w_locale(tp, G=G, s0=s0):
    """Velocita' dell'urto fra lo stato dell'onda semplice emesso in tp e il gas in quiete."""
    if tp < 1e-9:
        return 1.0
    st = gd.onda_semplice(s0, "u+a", u=2 * tp)
    return gd.riemann(st, s0).onda_dx.w


def traiettoria_urto(G, s0, n=200):
    gm = G.gamma
    lamg = lambda tp: 1 + (gm + 1) * tp
    Xg = lambda tp, t: tp**2 + lamg(tp) * (t - tp)
    # t(tp) con dt/dtp = (dX/dtp)/(w - lam); dX/dtp = (g+1) t - 1 - 2 g tp
    def rhs(tp, y):
        t = y[0]
        return [((gm + 1) * t - 1 - 2 * gm * tp) / (w_locale(tp, G, s0) - lamg(tp))]
    eps = 1e-3
    t_start = (1 + 4 * gm * eps / 3) / (gm + 1)          # avvio con la soluzione d'urto debole
    sol = solve_ivp(rhs, (eps, 1.0), [t_start], rtol=1e-9, atol=1e-11,
                    dense_output=True, max_step=0.01)
    tp = np.concatenate([[0.0], np.linspace(eps, 1.0, n)])
    t = np.concatenate([[1 / (gm + 1)], sol.sol(tp[1:])[0]])
    x = Xg(tp, t)
    return tp, t, x


tp_s, t_s, x_s = traiettoria_urto(G, s0)
tA, xA = t_s[-1], x_s[-1]
# controllo: Lax (1 < w < u+a dietro) e posizione prima dell'inviluppo
for tp, tt in zip(tp_s[1:], t_s[1:]):
    assert 1 < w_locale(tp) < lam(tp) and tt < t_env(tp)
# confronto con la regola dell'urto debole (media delle velocita' caratteristiche, analitica):
tA_deb = (1 + 4 * g / 3) / (g + 1)
xA_deb = tA_deb + g / 3
print("\n=== Traiettoria dell'urto (stima, fuori modulo) ===")
print(f"punto A (urto raggiunge la coda): tA = {tA:.4f}, xA = {xA:.4f};  w(A-) = {w_locale(1.0):.4f}")
print(f"regola d'urto debole (analitica): tA = {tA_deb:.4f}, xA = {xA_deb:.4f}, w = (g+3)/2 = {(g+3)/2:.3f}"
      " (sottostima: Ms ~ 2,9 non e' debole)")
print(f"pistone in tA: x = {float(xp(tA)):.4f}; coda: x = 1 + (2+g)(tA-1) = {1 + lam(1)*(tA-1):.4f}")

# =============================================================================
# 4. ISTANTI DI ARRIVO (caratteristiche rettilinee) DOPO IL PUNTO A
# =============================================================================
sol = aria["rp"]
v = sol.velocita()
w_rp, uC = sol.onda_dx.w, sol.u_star
# testa e coda dell'espansione u-a raggiungono il pistone x = 2t - 1
tP_h = (xA + 1 - v["testa_sx"] * tA) / (VP - v["testa_sx"])
tP_c = (xA + 1 - v["coda_sx"] * tA) / (VP - v["coda_sx"])
# testa dell'espansione riflessa u+a: nel gas 2 (u2 + a2) fino al contatto, poi nel gas 3 fino all'urto
lam_r2 = s2.u + s2.a
xPh = VP * tP_h - 1
tC = (xPh - lam_r2 * tP_h - (xA - uC * tA)) / (uC - lam_r2)
xC = xA + uC * (tC - tA)
lam_r3 = s3.u + s3.a
tS = (xC - lam_r3 * tC - (xA - w_rp * tA)) / (w_rp - lam_r3)
xS = xA + w_rp * (tS - tA)
print("\n=== Istanti (stima con semirette) ===")
print(f"fronti dell'espansione u-a da A: testa {v['testa_sx']:.4f}, coda {v['coda_sx']:.4f}; contatto {uC:.4f}; urto {w_rp:.4f}")
print(f"testa al pistone: t = {tP_h:.3f}, x = {xPh:.3f};  coda al pistone: t = {tP_c:.3f}, x = {VP*tP_c-1:.3f}")
print(f"riflessa u+a: pendenze {lam_r2:.4f} (gas 2), {lam_r3:.4f} (gas 3); al contatto t = {tC:.3f}, x = {xC:.3f};"
      f" all'urto t = {tS:.3f}, x = {xS:.3f}")

# =============================================================================
# 5. DIPENDENZA DA GAMMA E DA R
# =============================================================================
gas_tab = [("elio (He)", 5 / 3, 2077.0), ("aria", 1.4, 287.0), ("SF$_6$", 1.1, 56.9)]
print("\n=== Dipendenza dal gas ===")
print("gas, gamma, t_BU, a1, p1, M1, pRP, uRP, wRP, MsRP, w_inf, p_inf, rho6/rho0, T6/T0, T5/T0, p4, tA, xA")
righe = []
for nome, gm, R in gas_tab:
    rr = analisi(gm, stampa=(nome != "aria"), nome_gas=nome)
    tps, ts, xs = traiettoria_urto(rr["G"], rr["s0"], n=60)
    rr["tA"], rr["xA"] = ts[-1], xs[-1]
    righe.append((nome, rr))
    print(f"{nome:10s} {gm:.3f} {rr['t_bu']:.3f} {rr['s1'].a:.3f} {rr['s1'].p:.2f} {rr['s1'].M:.3f}"
          f" {rr['rp'].p_star:.2f} {rr['rp'].u_star:.3f} {rr['w_rp']:.3f} {rr['rp'].onda_dx.Ms:.3f}"
          f" {rr['w_inf']:.3f} {rr['s6'].p:.2f} {rr['s6'].rho/rr['s0'].rho:.3f} {rr['s6'].a**2:.3f}"
          f" {rr['s5'].a**2:.3f} {rr['s4'].p:.2f} {rr['tA']:.3f} {rr['xA']:.3f}")
    # T0 = p0/(rho0 R) = 1/(gamma R) nelle unita' della traccia: R cambia solo la scala di T
    print(f"           T0 = 1/(gamma R) (unita' della traccia); a0 a 293 K: {np.sqrt(gm*R*293.15):.0f} m/s")
print("limite gamma -> 1: p1 -> e^2 =", f"{np.e**2:.3f}", "; t_BU -> 1/2; w_inf -> 1+sqrt2 =", f"{1+np.sqrt(2):.3f}")

# =============================================================================
# 6. FIGURE
# =============================================================================
C_URTO, C_COMP, C_ESP, C_CONT, C_PIST = "crimson", "tab:blue", "tab:green", "tab:purple", "k"

# --- traiettoria della particella che l'urto investe in BU (limite del gas mai urtato)
def u_onda(x, t):
    """u nell'onda semplice fra pistone e urto (prima della coda), 2 nel plateau dopo la coda."""
    tp_sh = np.interp(t, t_s, tp_s) if t < tA else 1.0
    if t > 1 and x <= 1 + lam(1.0) * (t - 1):
        return VP
    # dietro l'urto restano le caratteristiche emesse fra tp_sh (appena raggiunta dall'urto) e min(t, 1)
    f = lambda tp: X(tp, t) - x
    lo, hi = tp_sh, min(t, 1.0)
    if f(lo) <= 0:
        return 2 * lo
    if f(hi) >= 0:
        return 2 * hi
    return 2 * brentq(f, lo, hi)


def u_part(t, y):
    x = y[0]
    if t < tA:
        return [u_onda(x, t)]
    # dopo A: gas 2 (u*) a destra della testa dell'espansione u-a, plateau 1 (u = 2) a sinistra;
    # dopo la riflessione sul pistone il gas 4 torna a u = 2
    if x > xA + v["testa_sx"] * (t - tA) and t < tP_h + (x - xPh) / lam_r2:
        return [uC]
    return [VP]


t_bu = 1 / (g + 1)
sp = solve_ivp(u_part, (t_bu + 1e-6, 6.0), [t_bu], max_step=0.005, rtol=1e-8)
tpart, xpart = sp.t, sp.y[0]
iA = np.searchsorted(tpart, tA)
print(f"\nparticella urtata in BU: x(tA) = {xpart[iA]:.3f} (xA = {xA:.3f}); x(t=3) = {np.interp(3, tpart, xpart):.3f}")


def strato_entropico(ax, tmax):
    tt = tpart[tpart <= tmax]
    xx = xpart[tpart <= tmax]
    dx_dx = np.where(tt < tA, np.interp(tt, t_s, x_s), xA + uC * (tt - tA))
    ax.fill_betweenx(tt, xx, dx_dx, color=C_CONT, alpha=0.15, lw=0)
    ax.plot(xx, tt, color=C_CONT, lw=1.0, ls=":", label="particella urtata in BU")


# --- 6.1 x-t: formazione dell'urto
fig, ax = plt.subplots(figsize=(6.6, 5.4))
TF = 1.8
tt = np.linspace(0, TF, 400)
ax.plot(xp(tt), tt, color=C_PIST, lw=2.2, label="pistone")
for tp in np.linspace(0, 1, 11):
    t_hit = np.interp(tp, tp_s, t_s) if tp > 0 else t_bu
    tv = np.linspace(tp, t_hit, 30)
    ax.plot(X(tp, tv), tv, color=C_COMP, lw=0.9)
    # prolungamento a piu' valori (senza l'urto) fino all'inviluppo
    tv = np.linspace(t_hit, min(t_env(tp), TF), 30)
    ax.plot(X(tp, tv), tv, color=C_COMP, lw=0.6, ls=":", alpha=0.8)
tpe = np.linspace(0, 1, 100)
ax.plot(X(tpe, t_env(tpe)), t_env(tpe), color="gray", ls="--", lw=1.2, label="inviluppo (soluzione a più valori)")
ax.plot(x_s, t_s, color=C_URTO, lw=2.2, label="urto (stima)")
tv = np.linspace(tA, TF, 20)
ax.plot(xA + w_rp * (tv - tA), tv, color=C_URTO, lw=2.2)
ax.plot(xA + uC * (tv - tA), tv, color=C_CONT, lw=1.4, ls="--", label="contatto (pseudo-Riemann)")
strato_entropico(ax, TF)
for lamv, tP in ((v["testa_sx"], tP_h), (v["coda_sx"], tP_c)):
    tv = np.array([tA, min(tP, TF)])
    ax.plot(xA + lamv * (tv - tA), tv, color=C_ESP, lw=1.0)
ax.plot([], [], color=C_ESP, lw=1.0, label="espansione $u-a$")
ax.plot(t_bu, t_bu, "o", color=C_URTO, ms=6)
ax.annotate("BU", (t_bu, t_bu), xytext=(0.62, 0.30), color=C_URTO,
            arrowprops=dict(arrowstyle="-", color=C_URTO, lw=0.7))
ax.plot(xA, tA, "s", color="k", ms=5)
ax.annotate("A", (xA, tA), xytext=(xA + 0.18, tA - 0.12))
ax.plot(1, 1, "o", color="k", ms=4)
ax.annotate("testa $x=t$", (0.3, 0.3), xytext=(0.02, 0.62), color=C_COMP, fontsize=8,
            arrowprops=dict(arrowstyle="-", color=C_COMP, lw=0.6))
ax.annotate("coda $x=1+(2+\\gamma)(t-1)$", (1.6, 1 + 0.6 / lam(1.0)), xytext=(1.75, 0.98), color=C_COMP,
            fontsize=8, arrowprops=dict(arrowstyle="-", color=C_COMP, lw=0.6))
ax.text(2.7, 0.55, "0", fontsize=13, weight="bold")
ax.text(1.55, 1.22, "1", fontsize=12, weight="bold")
ax.set_xlim(0, 3.5)
ax.set_ylim(0, TF)
ax.set_xlabel("$x$")
ax.set_ylabel("$t$")
ax.set_title("Formazione dell'urto (aria, $\\gamma=1{,}4$)")
ax.legend(loc="upper left", fontsize=7.2)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "xt_formazione.pdf"))
plt.close(fig)

# --- 6.2 x-t: campo dopo A
fig, ax = plt.subplots(figsize=(7.0, 5.6))
TM = 3.5
tt = np.linspace(0, TM, 400)
ax.plot(xp(tt), tt, color=C_PIST, lw=2.2, label="pistone")
for tp in np.linspace(0, 1, 6):
    t_hit = np.interp(tp, tp_s, t_s) if tp > 0 else t_bu
    tv = np.linspace(tp, t_hit, 20)
    ax.plot(X(tp, tv), tv, color=C_COMP, lw=0.9)
ax.plot(x_s, t_s, color=C_URTO, lw=2.2, label="urto")
tv = np.linspace(tA, tS, 10)
ax.plot(xA + w_rp * (tv - tA), tv, color=C_URTO, lw=2.2)
# dopo S l'urto rallenta verso w_inf (raccordo qualitativo)
tv = np.linspace(tS, TM, 50)
tau = 1.0
wv = aria["w_inf"] + (w_rp - aria["w_inf"]) * np.exp(-(tv - tS) / tau)
xv = xS + np.concatenate([[0], np.cumsum(0.5 * (wv[1:] + wv[:-1]) * np.diff(tv))])
ax.plot(xv, tv, color=C_URTO, lw=2.2)
tv = np.linspace(tA, TM, 20)
ax.plot(xA + uC * (tv - tA), tv, color=C_CONT, lw=1.4, ls="--", label="contatto (pseudo-Riemann)")
strato_entropico(ax, TM)
for lamv, tP in ((v["testa_sx"], tP_h), (v["coda_sx"], tP_c)):
    ax.plot([xA, xA + lamv * (tP - tA)], [tA, tP], color=C_ESP, lw=1.0)
ax.plot([], [], color=C_ESP, lw=1.0, label="espansione $u-a$ (da A)")
ax.plot([xPh, xC, xS], [tP_h, tC, tS], color=C_ESP, lw=1.0, ls="-.", label="testa dell'espansione riflessa $u+a$")
for p_, lab, off in (((xA, tA), "A", (0.15, -0.1)), ((xC, tC), "C", (0.15, -0.12)),
                     ((xS, tS), "S", (0.15, -0.12)), ((xPh, tP_h), "P", (-0.3, 0.03))):
    ax.plot(*p_, "o", color="k", ms=4)
    ax.annotate(lab, p_, xytext=(p_[0] + off[0], p_[1] + off[1]), fontsize=9)
ax.text(5.5, 1.3, "0", fontsize=13, weight="bold")
ax.text(1.55, 1.27, "1", fontsize=11, weight="bold")
ax.annotate("2", (3.12, 1.85), xytext=(2.2, 2.55), fontsize=12, weight="bold",
            arrowprops=dict(arrowstyle="->", lw=0.7))
ax.annotate("3", (4.1, 2.1), xytext=(5.0, 1.85), fontsize=12, weight="bold",
            arrowprops=dict(arrowstyle="->", lw=0.7))
ax.annotate("4", (5.6, 3.1), xytext=(4.2, 3.3), fontsize=12, weight="bold",
            arrowprops=dict(arrowstyle="->", lw=0.7))
ax.annotate("strato entropico", (2.9, 2.2), xytext=(1.3, 3.0), fontsize=8, color=C_CONT,
            arrowprops=dict(arrowstyle="->", lw=0.6, color=C_CONT))
ax.set_xlim(1.2, 8.2)
ax.set_ylim(0.9, TM)
ax.set_xlabel("$x$")
ax.set_ylabel("$t$")
ax.set_title("Piano $x$–$t$ (aria): dopo A caratteristiche rettilinee, raccordo dell'urto qualitativo")
ax.title.set_fontsize(9.5)
ax.legend(loc="lower right", fontsize=7)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "xt_campo.pdf"))
plt.close(fig)

# --- 6.3 p-u (lineare e logaritmica)
def pu(ax, log):
    pp = np.linspace(1, 13, 600)
    # compressione isentropica u+a da 0 (onda semplice): u = (a0/delta)[(p)^q - 1]
    u_is = (1 / d) * (pp ** G.q - 1)
    ax.plot(u_is, pp, color=C_COMP, lw=1.6, label="isentropica $u+a$ da 0 (compressione 0→1)")
    ax.plot(gd.curva_onda(s0, "u+a", pp), pp, color=C_URTO, lw=1.6, label="Hugoniot $u+a$ da 0 (urto)")
    pl = np.linspace(6.0, s1.p, 200)
    ax.plot(gd.curva_onda(s1, "u-a", pl), pl, color=C_ESP, lw=1.6, label="espansione $u-a$ da 1")
    p4l = np.linspace(s4.p - 0.6, s2.p, 100)
    ax.plot(gd.curva_onda(s2, "u+a", p4l), p4l, color=C_ESP, lw=1.3, ls="-.",
            label="espansione $u+a$ da 2 (riflessione sul pistone)")
    ax.axvline(VP, color="k", lw=0.8, ls=":")
    for st, lab, off in ((s0, "0", (0.05, 0.15)), (s1, "1", (-0.18, 0.15)), (s2, "2=3", (0.05, 0.1)),
                         (s4, "4", (-0.2, -0.45)), (s6, "5=6", (0.06, -0.5))):
        ax.plot(st.u, st.p, "o", color="k", ms=4)
        ax.text(st.u + off[0], st.p * (1 + 0.04 * np.sign(off[1])) if log else st.p + off[1], lab, fontsize=9)
    ax.set_xlabel("$u$")
    ax.set_ylabel("$p$")
    if log:
        ax.set_yscale("log")
        ax.set_ylim(0.9, 14)
    else:
        ax.set_ylim(0, 13)
    ax.set_xlim(-0.05, 2.5)


fig, axs = plt.subplots(1, 2, figsize=(10.5, 4.6))
pu(axs[0], False)
pu(axs[1], True)
axs[0].set_title("scala lineare")
axs[1].set_title("scala logaritmica")
axs[0].legend(loc="upper left", fontsize=7.5)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "pu.pdf"))
plt.close(fig)

# zoom attorno all'intersezione
fig, ax = plt.subplots(figsize=(5.2, 4.2))
pp = np.linspace(7.5, 11, 400)
ax.plot((1 / d) * (pp ** G.q - 1), pp, color=C_COMP, lw=1.6, label="isentropica $u+a$ da 0")
ax.plot(gd.curva_onda(s0, "u+a", pp), pp, color=C_URTO, lw=1.6, label="Hugoniot $u+a$ da 0")
ax.plot(gd.curva_onda(s1, "u-a", pp[pp <= s1.p]), pp[pp <= s1.p], color=C_ESP, lw=1.6, label="espansione $u-a$ da 1")
p4l = np.linspace(7.5, s2.p, 100)
ax.plot(gd.curva_onda(s2, "u+a", p4l), p4l, color=C_ESP, lw=1.3, ls="-.", label="espansione $u+a$ da 2")
ax.axvline(VP, color="k", lw=0.8, ls=":")
it0 = sol.iterazioni[0]
ax.plot([it0[1], it0[1]], [it0[3], it0[2]], color="gray", lw=3, alpha=0.4)
ax.text(it0[1] + 0.01, 0.5 * (it0[2] + it0[3]), "$u^{(0)}=2$:\n$p_2^{(0)}-p_3^{(0)}$", fontsize=8, color="gray")
for st, lab, off in ((s1, "1", (-0.06, 0.1)), (s2, "2=3", (0.01, 0.1)), (s4, "4", (0.012, -0.12)),
                     (s6, "5=6", (0.012, 0.05))):
    ax.plot(st.u, st.p, "o", color="k", ms=4)
    ax.text(st.u + off[0], st.p + off[1], lab, fontsize=9)
ax.set_xlim(1.8, 2.3)
ax.set_ylim(7.5, 11)
ax.set_xlabel("$u$")
ax.set_ylabel("$p$")
ax.set_title("Piano $p$–$u$ vicino a $u=2$ (aria)")
ax.legend(loc="lower left", fontsize=7)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "pu_zoom.pdf"))
plt.close(fig)

# --- 6.4 dipendenza da gamma
gv = np.linspace(1.05, 5 / 3, 40)
tb, wr, wi, p1v, prp, pinf, r6 = [], [], [], [], [], [], []
for gm in gv:
    rr = analisi(gm)
    tb.append(rr["t_bu"]); wr.append(rr["w_rp"]); wi.append(rr["w_inf"])
    p1v.append(rr["s1"].p); prp.append(rr["rp"].p_star); pinf.append(rr["s6"].p)
    r6.append(rr["s6"].rho / rr["s0"].rho)
fig, axs = plt.subplots(1, 3, figsize=(11, 3.6))
axs[0].plot(gv, tb, color="k")
axs[0].set_ylabel("$t_{BU}=x_{BU}$")
axs[1].plot(gv, wr, color=C_URTO, ls="--", label="$w$ pseudo-Riemann")
axs[1].plot(gv, wi, color=C_URTO, label="$w_\\infty$ asintotica")
axs[1].set_ylabel("velocità dell'urto")
axs[1].legend(fontsize=8)
axs[2].plot(gv, p1v, color=C_COMP, label="$p_1$ (isentropica)")
axs[2].plot(gv, prp, color=C_ESP, ls="--", label="$p_2=p_3$ (pseudo-Riemann)")
axs[2].plot(gv, pinf, color=C_URTO, label="$p_5=p_6$ (asintotica)")
axs[2].set_ylabel("$p$")
axs[2].legend(fontsize=8)
for ax in axs:
    ax.set_xlabel("$\\gamma$")
    for nome, gm, R in gas_tab:
        ax.axvline(gm, color="gray", lw=0.7, ls=":")
for nome, gm, R in gas_tab:
    axs[0].text(gm, max(tb) * 1.005, nome.split(" ")[0], fontsize=7, ha="center", va="bottom")
fig.tight_layout()
fig.savefig(os.path.join(FIG, "gas.pdf"))
plt.close(fig)
print("\nfigure salvate in", FIG)
