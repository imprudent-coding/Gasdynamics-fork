"""Prova del 01/01/1975, esercizio 1: serbatoio di grande capacita' collegato all'ambiente da un condotto a
sezione costante, chiuso a meta' da un setto rimosso istantaneamente a t0. Evoluzione del campo nel condotto
fino al completamento del transitorio.

Dati: aria, gamma = 1,4, R = 287 J/(kg K) (il testo scrive "KJ/kg K": refuso).
Serbatoio (a sinistra): P0 = 0,04 atm, T0 = 98 K.  Ambiente (a destra): p = 1 atm, T = 287,7 K.

Riproduce tutti i numeri e le figure dello svolgimento.
Uso:  python calcoli.py   (dalla cartella dello svolgimento)

Convenzioni: x verso destra, serbatoio in x = 0, ambiente in x = L, setto in x = L/2. La lunghezza L non e'
data: posizioni in x/L, tempi in tau = L/a4 (a4 = velocita' del suono dell'ambiente).

Parti fuori modulo (scritte qui riusando gasdyn):
  - testa dell'onda riflessa dall'ingresso dentro il ventaglio centrato (integrazione analitica di
    dx/dt = u - a nel ventaglio), controllata con un'integrazione numerica che usa gd.Ventaglio.stato_xi;
  - condizione d'ingresso dall'ambiente (a^2 + delta u^2 = a4^2, s = s4, R2 dall'interno) e suo limite
    sonico; gd.contorno_serbatoio non limita l'ingresso a M = 1 (vedi nota in fondo all'output);
  - stima asintotica della durata del transitorio (u + a = 3 R2 - 2 R1 vicino allo stato finale);
  - controllo numerico indipendente dell'intero transitorio (schema ai volumi finiti HLLC con le
    condizioni al contorno di ingresso e uscita): serve solo a confermare il quadro qualitativo e lo
    stato finale, non compare nello svolgimento.
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon

QUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(QUI, "..", ".."))
import gasdyn as gd  # noqa: E402

FIG = os.path.join(QUI, "figure")
os.makedirs(FIG, exist_ok=True)
plt.rcParams.update({"font.size": 10, "axes.grid": True, "grid.alpha": 0.3,
                     "axes.formatter.use_locale": False})
from matplotlib.ticker import FuncFormatter  # noqa: E402
VIRGOLA = FuncFormatter(lambda v, pos: f"{v:g}".replace(".", ","))

ATM = gd.ATM
G = gd.Gas(1.4, 287.0)
g, d = G.gamma, G.delta                      # delta = (gamma-1)/2 (F §1)
ESP = 2 * g / (g - 1)                        # p ~ a^(2 gamma/(gamma-1)) (F §4)


def riga(titolo):
    print("\n" + "=" * 78 + "\n" + titolo + "\n" + "=" * 78)


# =============================================================================
riga("1. Dati e stati iniziali")
# =============================================================================
s1 = gd.State.da_pT(G, 0.04 * ATM, 98.0, nome="1")        # gas del serbatoio, meta' sinistra
s4 = gd.State.da_pT(G, 1.00 * ATM, 287.7, nome="4")       # aria ambiente, meta' destra
a1, a4 = s1.a, s4.a
print(f"stato 1: p = {s1.p/ATM:.3f} atm, T = {s1.T:.1f} K, rho = {s1.rho:.4f} kg/m3, a = {a1:.1f} m/s")
print(f"stato 4: p = {s4.p/ATM:.3f} atm, T = {s4.T:.1f} K, rho = {s4.rho:.4f} kg/m3, a = {a4:.1f} m/s")
print(f"p4/p1 = {s4.p/s1.p:.1f}, a4/a1 = {a4/a1:.4f}, rho1 a1 = {s1.rho*a1:.1f}, rho4 a4 = {s4.rho*a4:.1f}")
ds14 = gd.salto_entropia(s4, s1)["ds_R"]
print(f"entropia: (s1 - s4)/R = {ds14:+.4f}  (gas del serbatoio piu' freddo: entropia minore)")

# =============================================================================
riga("2. Problema di Riemann alla rimozione del setto (F §8)")
# =============================================================================
q = (g - 1) / (2 * g)
z = (a4 / a1) * (s1.p / s4.p) ** q                        # innesco isentropico (F §8)
u0 = (z * (a1 + d * s1.u) - (a4 - d * s4.u)) / (d * (1 + z))
sol = gd.riemann(s1, s4)
k0, u0m, p20, p30 = sol.iterazioni[0]
print(f"z = {z:.4f}, u0 = {u0:.1f} m/s (modulo: {u0m:.1f})")
print(f"primo tentativo: p2 = {p20/ATM:.4f} atm (urto), p3 = {p30/ATM:.4f} atm (espansione), "
      f"|1 - p2/p3| = {abs(1-p20/p30):.3f}")
for k, uk, p2k, p3k in sol.iterazioni:
    print(f"  k = {k}: u = {uk:9.3f} m/s, p2 = {p2k/ATM:.5f}, p3 = {p3k/ATM:.5f} atm")
print(sol)
s2, s3 = sol.stati["2"], sol.stati["3"]
urto, vent = sol.onda_sx, sol.onda_dx
ustar, pstar = sol.u_star, sol.p_star
for s in (s2, s3):
    print(f"stato {s.nome}: p = {s.p/ATM:.4f} atm, u = {s.u:.1f} m/s, rho = {s.rho:.4f}, T = {s.T:.1f} K, "
          f"a = {s.a:.1f} m/s, M = {s.Mu:+.3f}, R1 = {s.R1:.1f}, R2 = {s.R2:.1f}, "
          f"T0 = {s.T0:.1f} K, p0 = {s.p0/ATM:.3f} atm")
print(f"urto u-a: w = {urto.w:.1f} m/s, Ms = {urto.Ms:.3f}, p2/p1 = {urto.p2_p1:.3f}, "
      f"rho2/rho1 = {urto.rho2_rho1:.3f}, Mr valle = {urto.Mr_valle:.3f}")
ds_urto = gd.salto_entropia(s1, s2)["ds_R"]
print(f"  (s2 - s1)/R = {ds_urto:+.4f}   s2/s1 = {s2.s/s1.s:.4f}")
print(f"contatto: w = u* = {ustar:.1f} m/s; rho2/rho3 = {s2.rho/s3.rho:.5f}, T2/T3 = {s2.T/s3.T:.5f}, "
      f"(s2 - s3)/R = {gd.salto_entropia(s3, s2)['ds_R']:+.5f}")
print(f"ventaglio u+a: testa {vent.xi_testa:.1f} m/s, coda u3 + a3 = {vent.xi_coda:.2f} m/s, "
      f"transonico = {vent.transonico}")

# Temperatura del serbatoio che annulla esattamente il contatto (solo per il commento in chat).
from scipy.optimize import brentq  # noqa: E402
T_null = brentq(lambda T: (lambda r: r.stati["2"].rho - r.stati["3"].rho)(
    gd.riemann(gd.State.da_pT(G, 0.04 * ATM, T), s4)), 60.0, 150.0)
print(f"T del serbatoio che annulla il contatto: {T_null:.2f} K (dato: 98 K)")

# Punto sonico del ventaglio (u + a = 0 sul raggio x = L/2), F §5 con u4 = 0
a_S = 2 * a4 / (g + 1)
sS = vent.stato_xi(0.0, nome="S")
print(f"punto sonico: a_S = 2 a4/(g+1) = {a_S:.1f} m/s (modulo {sS.a:.1f}), u_S = {sS.u:.1f}, "
      f"p_S = {sS.p/ATM:.4f} atm = p4 (2/(g+1))^(2g/(g-1)) = {(2/(g+1))**ESP:.4f} atm, T_S = {sS.T:.1f} K")

# =============================================================================
riga("3. Tempi caratteristici (tau = L/a4)")
# =============================================================================
t_urto = 0.5 * a4 / abs(urto.w)              # urto in x = 0
t_cont = 0.5 * a4 / abs(ustar)               # contatto in x = 0
t_testa = 0.5                                # testa del ventaglio in x = L
xiT = vent.xi_coda / a4                      # pendenza della coda (adimensionale)
print(f"urto all'uscita (x=0):     t = {t_urto:.3f} tau")
print(f"contatto all'uscita:       t = {t_cont:.3f} tau")
print(f"testa ventaglio all'ingresso (x=L): t = {t_testa:.3f} tau")

# Testa della riflessa u-a nel ventaglio (fuori modulo): con X = x - L/2, xi = X/t,
#   u - a = [(1-delta) xi - 2 a4]/(1+delta)  =>  dX/dt = k X/t - c,  k = (1-d)/(1+d), c = 2a4/(1+d)
#   =>  X = C t^k - (a4/delta) t   (adimensionale: a4 = L = 1)
k = (1 - d) / (1 + d)
C = (0.5 + 0.5 / d) / 0.5**k
t_esce = (C / (1 / d + xiT)) ** (1 / (1 - k))          # incontro con la coda X = xiT t
x_esce = 0.5 + xiT * t_esce
t_son = (C * d) ** (1 / (1 - k))                       # attraversa il raggio sonico X = 0
lam3 = (s3.u - s3.a) / a4
t_rif0 = t_esce + x_esce / abs(lam3)
print(f"testa riflessa nel ventaglio: X = C t^k - t/delta, k = {k:.4f}, C = {C:.4f}")
print(f"  attraversa il raggio sonico x = L/2 a t = {t_son:.3f} tau")
print(f"  esce dal ventaglio (coda) a t = {t_esce:.3f} tau, x = {x_esce:.3f} L")
print(f"  poi in 3 con u3 - a3 = {s3.u - s3.a:.1f} m/s: arriva in x = 0 a t = {t_rif0:.3f} tau")
# controllo numerico con lo stato del ventaglio del modulo
tt, xx = 0.5, 1.0                                     # parte dall'ingresso a t = L/(2 a4)
dt = 1e-5
while True:
    st = vent.stato_xi((xx - 0.5) / tt * a4)
    if (xx - 0.5) / tt * a4 <= vent.xi_coda:
        break
    xx += (st.u - st.a) / a4 * dt
    tt += dt
print(f"  controllo (integrazione con gd.Ventaglio.stato_xi): esce a t = {tt:.3f}, x = {xx:.3f}")
t_coda0 = 0.5 / abs(xiT)
print(f"coda del ventaglio come semiretta: arriverebbe in x = 0 a t = {t_coda0:.2f} tau (stima per eccesso)")

# =============================================================================
riga("4. Contorni")
# =============================================================================
print(f"uscita (x=0) dopo l'urto: u2 + a2 = {s2.u + s2.a:.1f} m/s < 0 -> uscita supersonica, "
      f"nessuna condizione; getto sottoespanso p2/p_serb = {s2.p/s1.p:.2f}")


def ingresso(R2, a0=a4):
    """Stato all'ingresso dall'ambiente (fuori modulo): a = R2 - delta u, a^2 + delta u^2 = a0^2."""
    A, B, Cc = d * (1 + d), -2 * d * R2, R2**2 - a0**2
    u = (-B - np.sqrt(B * B - 4 * A * Cc)) / (2 * A)
    a = R2 - d * u
    return u, a


a_st = a4 * np.sqrt(2 / (g + 1))
R2_st, R1_st = (1 - d) * a_st, (1 + d) * a_st
u_chk, a_chk = ingresso(R2_st)
print(f"ingresso sonico: a* = a4 sqrt(2/(g+1)) = {a_st:.1f} m/s, R2* = (1-d) a* = {R2_st:.1f}, "
      f"R1* = (1+d) a* = {R1_st:.1f}  (controllo: u = {u_chk:.1f}, a = {a_chk:.1f})")
# raggio del ventaglio con R2 = R2*: R2 = a4 + 2 delta (xi - a4)/(1+delta)
xi_st = a4 + (R2_st - a4) * (1 + d) / (2 * d)
print(f"raggio del ventaglio con R2 = R2*: xi = {xi_st:.1f} m/s; raggio sonico R2 = {sS.R2:.1f}; coda R2 = {s3.R2:.1f}")
for R2 in (a4, 320.0, 300.0, 280.0, 260.0, R2_st):
    u, a = ingresso(R2)
    print(f"  R2 = {R2:6.1f}: u = {u:7.1f}, a = {a:6.1f}, M = {u/a:+.3f}, R1 = {a - d*u:6.1f}, "
          f"p = {(a/a4)**ESP:.4f} atm")

# Stato finale F: corrente sonica uniforme
cr = gd.critiche(G)
sF = gd.State.da_pa(G, cr["p*_p0"] * s4.p, a_st, u=-a_st).con(nome="F")
mdot = sF.rho * a_st
print(f"stato F: M = {sF.Mu:+.3f}, a = {sF.a:.1f} m/s, u = {sF.u:.1f} m/s, T = {sF.T:.2f} K, "
      f"p = {sF.p/ATM:.4f} atm, rho = {sF.rho:.4f} kg/m3, p0 = {sF.p0/ATM:.3f} atm, T0 = {sF.T0:.1f} K")
print(f"portata per unita' di area: rho* a* = {mdot:.1f} kg/(s m2); "
      f"p_uscita / p_serbatoio = {sF.p/s1.p:.1f}")
print(f"salto di R1 dovuto alla riflessa: da {a4:.1f} a {R1_st:.1f} m/s (+{R1_st - a4:.1f})")

# Intensita' massima della compressione riflessa nella zona 3 (R2 = R2,3 uniforme, R1 da a4 a R1*)
a_c = (R1_st + s3.R2) / 2
u_c = (s3.R2 - R1_st) / (2 * d)
pr = (a_c / s3.a) ** ESP
Ms_c = gd.mach_da_urto(pr, G, "p2/p1")
print(f"compressione u-a in 3 fino a R1 = R1*: a = {a_c:.1f}, u = {u_c:.1f}, M = {u_c/a_c:+.3f}, "
      f"p/p3 = {pr:.3f} -> urto equivalente Ms = {Ms_c:.3f}")
print(f"u - a = 2 R2 - 3 R1, u + a = 3 R2 - 2 R1 (gamma = 1,4): controllo su F: u+a = "
      f"{3*R2_st - 2*R1_st:.2e}")

# Stima della durata: vicino a F, u + a = 3 (R2 - R2*): una u+a con deficit dR2 attraversa il condotto
# in L/(3 dR2) -> al tempo t il deficit residuo e' ~ L/(3 t) (in tau: a4/(3 t)).
for t in (5, 10, 30):
    dR2 = a4 / (3 * t)
    a_ = (R1_st + R2_st - dR2) / 2
    u_ = (R2_st - dR2 - R1_st) / (2 * d)
    print(f"  t = {t:3d} tau: dR2 ~ {dR2:5.1f} m/s -> M ~ {u_/a_:+.3f}, p ~ {(a_/a4)**ESP:.3f} atm")

# =============================================================================
riga("5. Figure")
# =============================================================================
COL = {"urto": "C3", "cont": "C2", "vent": "C0", "rif": "C1", "son": "0.4"}

# --- schema ---
fig, ax = plt.subplots(figsize=(7.0, 2.6))
ax.axis("off")
ax.add_patch(Rectangle((-2.2, -1.6), 2.2, 3.2, fc="#dfe9f5", ec="k", lw=1.5))
ax.add_patch(Rectangle((0, -0.35), 4, 0.7, fc="white", ec="none"))
ax.add_patch(Rectangle((0, -0.35), 2, 0.7, fc="#dfe9f5", ec="none"))
ax.add_patch(Rectangle((2, -0.35), 2, 0.7, fc="#fde9d9", ec="none"))
ax.add_patch(Rectangle((4, -1.6), 1.7, 3.2, fc="#fde9d9", ec="none"))
ax.plot([0, 4], [0.35, 0.35], "k", lw=1.5)
ax.plot([0, 4], [-0.35, -0.35], "k", lw=1.5)
ax.plot([0, 0], [0.35, 1.6], "w", lw=2.5)
ax.plot([0, 0], [-0.35, -1.6], "k", lw=1.5)
ax.plot([0, 0], [0.35, 1.6], "k", lw=1.5)
ax.plot([0, 0], [-0.35, 0.35], color="#dfe9f5", lw=3)
ax.plot([2, 2], [-0.35, 0.35], "k", lw=3)
ax.text(-1.1, 0.9, "serbatoio", ha="center")
ax.text(-1.1, -0.9, r"$p_1=0{,}04$ atm" + "\n" + r"$T_1=98$ K", ha="center", fontsize=9)
ax.text(1.0, 0.0, "1", ha="center", va="center", fontsize=12)
ax.text(3.0, 0.0, "4", ha="center", va="center", fontsize=12)
ax.text(4.85, 0.9, "ambiente", ha="center")
ax.text(4.85, -0.9, r"$p_4=1$ atm" + "\n" + r"$T_4=287{,}7$ K", ha="center", fontsize=9)
ax.text(2, -0.55, "setto", ha="center", va="top", fontsize=9)
ax.annotate("", xy=(4.0, -1.35), xytext=(0, -1.35), arrowprops=dict(arrowstyle="<->"))
ax.text(2, -1.3, "$L$", ha="center", va="bottom")
ax.text(0, 0.45, "$x=0$", ha="left", fontsize=8)
ax.text(4, 0.45, "$x=L$", ha="right", fontsize=8)
ax.set_xlim(-2.3, 5.8)
ax.set_ylim(-1.7, 1.7)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "schema.pdf"))
plt.close(fig)

# --- piano x-t ---
TMAX = 1.6
fig, ax = plt.subplots(figsize=(6.6, 6.0))
ax.plot([0.5, 0], [0, t_urto], color=COL["urto"], lw=2.2)
ax.plot([0.5, 0], [0, t_cont], color=COL["cont"], lw=1.2, ls="--")
# testa riflessa (curva esatta nel ventaglio + semiretta in 3)
tr = np.linspace(0.5, t_esce, 200)
xr = 0.5 + C * tr**k - tr / d
ax.plot(xr, tr, color=COL["rif"], lw=2.0)
ax.plot([x_esce, 0], [t_esce, t_rif0], color=COL["rif"], lw=2.0)
# raggi del ventaglio fino alla testa riflessa
for xi in np.linspace(xiT, 1.0, 9):
    if abs(xi) < 1e-9:
        continue
    # incontro raggio X = xi t con la testa X = C t^k - t/d
    tc = (C / (1 / d + xi)) ** (1 / (1 - k))
    tc = min(tc, TMAX)
    lw = 1.8 if xi in (xiT, 1.0) else 0.7
    ax.plot([0.5, 0.5 + xi * tc], [0, tc], color=COL["vent"], lw=lw)
# raggio sonico verticale
ax.plot([0.5, 0.5], [0, t_son], color=COL["son"], lw=1.4, ls=":")
ax.text(0.505, 0.62, "$M=-1$", fontsize=8, color=COL["son"], rotation=90, va="center")
# coda dopo l'interazione (semiretta, stima)
ax.plot([x_esce, x_esce + xiT * (TMAX - t_esce)], [t_esce, TMAX], color=COL["vent"], lw=1.2, ls="--")
# riflesse u-a successive (schematiche)
for te, sl in ((0.62, -1.25), (0.78, -1.45), (0.95, -1.6), (1.12, -1.7), (1.3, -1.78), (1.45, -1.8)):
    tt_ = np.array([te, TMAX])
    ax.plot(1 + sl * (tt_ - te), tt_, color=COL["rif"], lw=0.7, ls="--")
ax.fill_between([0.5 + xiT * TMAX, 1.0], [TMAX, TMAX], [TMAX, TMAX], color="none")
ax.text(0.22, 0.08, "1", fontsize=13, ha="center")
ax.text(0.80, 0.12, "4", fontsize=13, ha="center")
ax.text(0.12, 0.40, "2", fontsize=13, ha="center")
ax.text(0.30, 0.60, "3", fontsize=13, ha="center")
ax.text(0.62, 0.30, "E", fontsize=12, ha="center", color=COL["vent"])
ax.text(0.80, 0.95, "I", fontsize=12, ha="center", color="k")
ax.text(0.15, 1.30, "C", fontsize=12, ha="center", color=COL["rif"])
ax.text(0.93, 1.50, r"$\to$F", fontsize=12, ha="center")
ax.text(0.17, 0.24, "urto", color=COL["urto"], rotation=-50, fontsize=8)
ax.text(0.30, 0.42, "contatto", color=COL["cont"], rotation=-55, fontsize=8)
ax.text(0.16, 1.07, "testa riflessa", color=COL["rif"], rotation=-27, fontsize=8)
ax.annotate("coda (semiretta)", xy=(0.375, 1.42), xytext=(0.17, 1.52), fontsize=7, color=COL["vent"],
            arrowprops=dict(arrowstyle="->", color=COL["vent"], lw=0.7))
ax.annotate("coda", xy=(0.47, 0.30), xytext=(0.36, 0.20), fontsize=8, color=COL["vent"],
            arrowprops=dict(arrowstyle="->", color=COL["vent"], lw=0.7))
ax.annotate("testa", xy=(0.93, 0.43), xytext=(0.95, 0.25), fontsize=8, color=COL["vent"],
            arrowprops=dict(arrowstyle="->", color=COL["vent"], lw=0.7))
for tm, lab in ((t_urto, f"{t_urto:.2f}"), (t_cont, f"{t_cont:.2f}"), (t_rif0, f"{t_rif0:.2f}")):
    ax.plot([0], [tm], "k.", ms=5)
    ax.text(0.01, tm, lab.replace(".", ","), fontsize=8, va="bottom")
ax.plot([1], [0.5], "k.", ms=5)
ax.text(0.99, 0.5, "0,5", fontsize=8, ha="right", va="bottom")
ax.plot([x_esce], [t_esce], "k.", ms=5)
ax.set_xlim(0, 1)
ax.set_ylim(0, TMAX)
ax.set_xlabel("$x/L$  (serbatoio a sinistra, ambiente a destra)")
ax.set_ylabel(r"$t\,a_4/L$")
fig.tight_layout()
ax.xaxis.set_major_formatter(VIRGOLA)
ax.yaxis.set_major_formatter(VIRGOLA)
fig.savefig(os.path.join(FIG, "xt.pdf"))
plt.close(fig)

# --- piano p-u (lineare e logaritmico) ---
pp_s = np.linspace(s1.p, 1.3 * ATM, 300)                    # urto u-a da 1
u_s = [gd.curva_onda(s1, "u-a", p) for p in pp_s]
pp_r = np.linspace(0.02 * ATM, s4.p, 300)                    # espansione u+a da 4
u_r = [gd.curva_onda(s4, "u+a", p) for p in pp_r]
uu_in = np.linspace(-a_st, 0, 200)                           # stati d'ingresso ammissibili
p_in = s4.p * ((a4**2 - d * uu_in**2) / a4**2) ** (ESP / 2)
uu_in2 = np.linspace(-1.25 * a_st, -a_st, 50)
p_in2 = s4.p * ((a4**2 - d * uu_in2**2) / a4**2) ** (ESP / 2)
uu_son = np.linspace(-380, -150, 200)                        # luogo sonico M = -1 a entropia s4
p_son = s4.p * (np.abs(uu_son) / a4) ** ESP
punti = [(s1, "1"), (s4, "4"), (s2, "2=3"), (sS, "S"), (sF, "F")]


def disegna_pu(ax, log):
    ax.plot(u_s, np.array(pp_s) / ATM, color=COL["urto"], lw=1.8, label="urto $u-a$ da 1")
    ax.plot(u_r, pp_r / ATM, color=COL["vent"], lw=1.8, label="espansione $u+a$ da 4")
    ax.plot(uu_in, p_in / ATM, color="k", lw=1.8, label=r"ingresso: $a^2+\delta u^2=a_4^2$, $s=s_4$")
    ax.plot(uu_in2, p_in2 / ATM, color="k", lw=1.0, ls=":")
    ax.plot(uu_son, p_son / ATM, color=COL["son"], lw=1.0, ls="--", label=r"stati sonici $M=-1$, $s=s_4$")
    for s, n in punti:
        ax.plot(s.u, s.p / ATM, "o", color="k", ms=4)
        dx_, dy_ = {"F": (-16, 4), "S": (4, 6), "2=3": (-26, -12)}.get(n, (8, 4))
        ax.annotate(n, (s.u, s.p / ATM), textcoords="offset points", xytext=(dx_, dy_), fontsize=10)
    ax.set_xlabel("$u$ [m/s]")
    ax.set_ylabel("$p$ [atm]")
    if log:
        ax.set_yscale("log")
        ax.set_ylim(0.02, 1.5)
    else:
        ax.set_ylim(0, 1.3)
    ax.set_xlim(-400, 100)
    ax.xaxis.set_major_formatter(VIRGOLA)
    ax.yaxis.set_major_formatter(VIRGOLA)


fig, axs = plt.subplots(1, 2, figsize=(10.5, 4.4))
disegna_pu(axs[0], False)
disegna_pu(axs[1], True)
axs[0].set_title("scala lineare")
axs[1].set_title("scala logaritmica")
axs[0].legend(fontsize=7.5, loc="upper left")
fig.tight_layout()
fig.savefig(os.path.join(FIG, "pu.pdf"))
plt.close(fig)

# --- profili lungo il condotto prima delle riflessioni e stato finale ---
x = np.linspace(0, 1, 801)
fig, axs = plt.subplots(2, 1, figsize=(6.6, 5.4), sharex=True)
for tau, st_ in ((0.25, "-"), (0.5, "-.")):
    # fino a t = 0,5 L/a4 nessuna onda riflessa e' rientrata: il campo e' quello del problema di Riemann
    # (l'urto uscito nel serbatoio non riflette nulla perche' l'uscita e' supersonica)
    c = sol.campiona((x - 0.5) * 1.0, tau / a4)           # L = 1 m, t in s = tau L/a4
    axs[0].plot(x, c["p"] / ATM, st_, color="k", lw=1.4, label=fr"$t={tau:.2f}\,L/a_4$".replace(".", "{,}"))
    axs[1].plot(x, c["M"], st_, color="k", lw=1.4)
axs[0].axhline(sF.p / ATM, color=COL["rif"], ls="--", lw=1.4, label=r"stato finale F ($t\to\infty$)")
axs[1].axhline(-1, color=COL["rif"], ls="--", lw=1.4)
axs[0].set_ylabel("$p$ [atm]")
axs[1].set_ylabel("$M=u/a$")
axs[1].set_xlabel("$x/L$")
axs[0].legend(fontsize=8, loc="center left")
fig.tight_layout()
for a_ in axs:
    a_.xaxis.set_major_formatter(VIRGOLA)
    a_.yaxis.set_major_formatter(VIRGOLA)
fig.savefig(os.path.join(FIG, "profili.pdf"))
plt.close(fig)
print("figure salvate in", FIG)


# =============================================================================
riga("6. Controllo numerico indipendente del transitorio (solo verifica, fuori modulo)")
# =============================================================================
def simulazione(N=400, tout=(0.5, 1.0, 1.5, 2.0, 5.0, 10.0, 20.0)):
    """Volumi finiti HLLC, primo ordine. Contorni: x = 0 uscita nel serbatoio (p = p1 se subsonica,
    estrapolazione se supersonica); x = L ingresso dall'ambiente (a^2 + delta u^2 = a4^2, s = s4,
    R2 dall'interno, limitato a M = -1)."""
    L = 1.0
    dx = L / N
    xc = (np.arange(N) + 0.5) * dx
    rho = np.where(xc < 0.5, s1.rho, s4.rho)
    u = np.zeros(N)
    p = np.where(xc < 0.5, s1.p, s4.p)

    def cons(r, u, p):
        return np.array([r, r * u, p / (g - 1) + 0.5 * r * u * u])

    def flux(r, u, p):
        E = p / (g - 1) + 0.5 * r * u * u
        return np.array([r * u, r * u * u + p, u * (E + p)])

    def hllc(rl, ul, pl, rr, ur, pr):
        al, ar = np.sqrt(g * pl / rl), np.sqrt(g * pr / rr)
        SL, SR = np.minimum(ul - al, ur - ar), np.maximum(ul + al, ur + ar)
        Ss = (pr - pl + rl * ul * (SL - ul) - rr * ur * (SR - ur)) / (rl * (SL - ul) - rr * (SR - ur))
        UL, UR, FL, FR = cons(rl, ul, pl), cons(rr, ur, pr), flux(rl, ul, pl), flux(rr, ur, pr)

        def star(r, u, p, S, U):
            cc = r * (S - u) / (S - Ss)
            return np.array([cc, cc * Ss, cc * (U[2] / r + (Ss - u) * (Ss + p / (r * (S - u))))])
        ULs, URs = star(rl, ul, pl, SL, UL), star(rr, ur, pr, SR, UR)
        return np.where(SL >= 0, FL, np.where(Ss >= 0, FL + SL * (ULs - UL),
                                              np.where(SR > 0, FR + SR * (URs - UR), FR)))

    U = cons(rho, u, p)
    t, kout, out = 0.0, 0, []
    while kout < len(tout):
        rho, u = U[0], U[1] / U[0]
        p = (g - 1) * (U[2] - 0.5 * rho * u * u)
        a = np.sqrt(g * p / rho)
        dt = 0.8 * dx / np.max(np.abs(u) + a)
        if t + dt * a4 / L > tout[kout]:
            dt = (tout[kout] - t) * L / a4
        gl = (rho[0], u[0], p[0]) if u[0] + a[0] < 0 else (rho[0], u[0], s1.p)
        if u[-1] < 0:
            ub, ab = ingresso(a[-1] + d * u[-1])
            if ub + ab < 0:
                ab, ub = a_st, -a_st
            gr = (s4.rho * (ab / a4) ** (2 / (g - 1)), ub, s4.p * (ab / a4) ** ESP)
        else:
            gr = (rho[-1], u[-1], s4.p)
        R_ = np.concatenate([[gl[0]], rho, [gr[0]]])
        Uu = np.concatenate([[gl[1]], u, [gr[1]]])
        P = np.concatenate([[gl[2]], p, [gr[2]]])
        F = hllc(R_[:-1], Uu[:-1], P[:-1], R_[1:], Uu[1:], P[1:])
        U = U - dt / dx * (F[:, 1:] - F[:, :-1])
        t += dt * a4 / L
        if abs(t - tout[kout]) < 1e-12:
            rho, u = U[0], U[1] / U[0]
            p = (g - 1) * (U[2] - 0.5 * rho * u * u)
            M = u / np.sqrt(g * p / rho)
            out.append((t, p[0] / ATM, M[0], p[N // 2] / ATM, M[N // 2], p[-1] / ATM, M[-1]))
            kout += 1
    return out


print("   t/tau | uscita x=0: p [atm], M | mezzeria: p, M | ingresso x=L: p, M")
for t, p0_, M0_, pm, Mm, pL, ML in simulazione():
    print(f"  {t:6.2f} | {p0_:.3f} {M0_:+.3f} | {pm:.3f} {Mm:+.3f} | {pL:.3f} {ML:+.3f}")
print(f"  (stato finale atteso: p = {sF.p/ATM:.3f} atm, M = -1)")
print("\nNota sul modulo: gd.contorno_serbatoio(stato3, 'destra', ...) restituisce un ingresso con |M| > 1\n"
      "(es. M = 1,32 da uno stato supersonico): un ingresso da un serbatoio in quiete in un condotto a sezione\n"
      "costante non puo' superare M = 1. Qui l'ingresso e' scritto a parte con il limite sonico.")
