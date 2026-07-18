"""Build the talk figures from real Cosmic.jl output and a matched CLASS run.

Inputs (all produced by the runs described in the talk):
  distances.tsv        z, d_L, d_A, d_C            [Mpc]      -- Cosmic.jl
  cls.tsv              l, D_l^TT, D_l^EE, D_l^TE   [muK^2]    -- Cosmic.jl
  pk.tsv               k, P_lin, P_halofit         [Mpc^3]    -- Cosmic.jl
  class_cl.dat         CLASS unlensed cl output
  class_pk.dat         CLASS P(k) output

Nothing here is synthetic: if a file is missing the corresponding figure is
skipped rather than filled in with made-up numbers.
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

GOLD, DEEP, INK, RED, BLUE = "#f3a914", "#b06e0a", "#28282d", "#be2828", "#1e50a0"

plt.rcParams.update({
    "font.size": 11, "axes.labelsize": 12, "axes.titlesize": 12,
    "axes.edgecolor": INK, "axes.labelcolor": INK,
    "xtick.color": INK, "ytick.color": INK, "text.color": INK,
    "axes.grid": True, "grid.alpha": 0.25, "grid.linewidth": 0.6,
    "figure.dpi": 200, "savefig.bbox": "tight", "savefig.facecolor": "white",
})

HERE = os.path.dirname(os.path.abspath(__file__))
p = lambda f: os.path.join(HERE, f)


def have(*files):
    missing = [f for f in files if not os.path.exists(p(f))]
    if missing:
        print("  skip (missing %s)" % ", ".join(missing))
        return False
    return True


# --- 1. distances -----------------------------------------------------------
if have("distances.tsv"):
    d = np.loadtxt(p("distances.tsv"))
    z, dL, dA, dC = d[:, 0], d[:, 1], d[:, 2], d[:, 3]
    fig, ax = plt.subplots(figsize=(6.2, 3.6))
    ax.loglog(z, dL, color=BLUE, lw=2, label=r"$d_L$  luminosity")
    ax.loglog(z, dC, color=DEEP, lw=2, label=r"$d_C$  comoving")
    ax.loglog(z, dA, color=RED, lw=2, label=r"$d_A$  angular diameter")
    imax = int(np.argmax(dA))
    ax.plot(z[imax], dA[imax], "o", color=RED, ms=6)
    ax.annotate("turnover\n$z\\simeq%.1f$" % z[imax],
                xy=(z[imax], dA[imax]), xytext=(z[imax] * 1.7, dA[imax] * 0.42),
                color=RED, fontsize=10,
                arrowprops=dict(arrowstyle="->", color=RED, lw=1.2))
    ax.set_xlabel("redshift $z$")
    ax.set_ylabel("distance  [Mpc]")
    ax.legend(frameon=False, loc="upper left", fontsize=10)
    fig.savefig(p("fig_distances.png"))
    plt.close(fig)
    print("  fig_distances.png")


# --- 2. CMB TT --------------------------------------------------------------
if have("cls.tsv"):
    c = np.loadtxt(p("cls.tsv"))
    l, tt = c[:, 0], c[:, 1]
    m = l >= 2
    fig, ax = plt.subplots(figsize=(6.0, 3.7))
    ax.plot(l[m], tt[m], color=BLUE, lw=1.8)
    ax.set_xlabel(r"multipole $\ell$")
    ax.set_ylabel(r"$\ell(\ell+1)C_\ell^{TT}/2\pi\ \ [\mu\mathrm{K}^2]$")
    ax.set_xlim(2, l.max())
    # mark the first acoustic peak actually present in the data
    band = m & (l > 150) & (l < 300)
    if band.any():
        ipk = np.argmax(np.where(band, tt, -np.inf))
        ax.annotate("first acoustic peak", xy=(l[ipk], tt[ipk]),
                    xytext=(l[ipk] * 1.9, tt[ipk] * 0.80), color=DEEP, fontsize=10,
                    arrowprops=dict(arrowstyle="->", color=DEEP, lw=1.2))
    ax.annotate("Silk damping", xy=(1500, np.interp(1500, l, tt)),
                xytext=(560, np.interp(1500, l, tt) + 0.42 * tt[m].max()),
                color=RED, fontsize=10,
                arrowprops=dict(arrowstyle="->", color=RED, lw=1.2))
    fig.savefig(p("fig_cls.png"))
    plt.close(fig)
    print("  fig_cls.png")


# --- 3. Cosmic vs CLASS, TT -------------------------------------------------
if have("cls.tsv", "class_cl.dat"):
    c = np.loadtxt(p("cls.tsv"))
    l, tt = c[:, 0], c[:, 1]
    cl = np.loadtxt(p("class_cl.dat"))
    lc, ttc = cl[:, 0], cl[:, 1]
    # CLASS writes l(l+1)C_l/2pi in dimensionless units of T_cmb^2
    Tcmb_muK = 2.7255e6
    ttc = ttc * Tcmb_muK ** 2
    lo, hi = max(l.min(), lc.min()), min(l.max(), lc.max())
    g = np.linspace(max(lo, 2), hi, 800)
    a = np.interp(g, l, tt)
    b = np.interp(g, lc, ttc)

    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(6.4, 4.4), sharex=True,
        gridspec_kw=dict(height_ratios=[2.4, 1], hspace=0.08))
    ax1.plot(g, b, color=INK, lw=3.0, alpha=0.35, label="CLASS 3.3.4")
    ax1.plot(g, a, color=BLUE, lw=1.5, label="Cosmic.jl")
    ax1.set_ylabel(r"$\ell(\ell+1)C_\ell^{TT}/2\pi\ [\mu\mathrm{K}^2]$")
    ax1.legend(frameon=False, fontsize=10)
    ax2.axhline(0.0, color=INK, lw=0.8)
    ax2.plot(g, 100 * (a / b - 1.0), color=DEEP, lw=1.4)
    ax2.set_ylabel("diff [%]", fontsize=10)
    ax2.set_xlabel(r"multipole $\ell$")
    ax2.set_xlim(2, hi)
    ax2.set_ylim(-3, 3)
    fig.savefig(p("fig_cls_compare.png"))
    plt.close(fig)
    print("  fig_cls_compare.png  (max |diff| = %.2f%% over 30<l<2000)"
          % np.nanmax(np.abs(100 * (a / b - 1))[(g > 30)]))


# --- 4. matter power --------------------------------------------------------
if have("pk.tsv"):
    d = np.loadtxt(p("pk.tsv"))
    k, plin, pnl = d[:, 0], d[:, 1], d[:, 2]
    fig, ax = plt.subplots(figsize=(6.0, 3.7))
    if os.path.exists(p("class_pk.dat")):
        # CLASS writes k in h/Mpc and P in (Mpc/h)^3; Cosmic uses 1/Mpc and Mpc^3
        H = 0.6766
        cp = np.loadtxt(p("class_pk.dat"))
        ax.loglog(cp[:, 0] * H, cp[:, 1] / H ** 3, color=INK, lw=3.0, alpha=0.30,
                  label="CLASS linear")
    ax.loglog(k, plin, color=BLUE, lw=1.6, label="Cosmic.jl linear")
    ax.loglog(k, pnl, color=RED, lw=1.6, label="Cosmic.jl halofit")
    ipk = int(np.argmax(plin))
    ax.annotate("turnover: matter--radiation equality",
                xy=(k[ipk], plin[ipk]),
                xytext=(k[ipk] * 2.2, plin[ipk] * 3.0),
                color=DEEP, fontsize=9.5,
                arrowprops=dict(arrowstyle="->", color=DEEP, lw=1.2))
    ax.annotate("non-linear growth", xy=(1.0, np.interp(1.0, k, pnl)),
                xytext=(0.12, np.interp(1.0, k, pnl) * 6.0),
                color=RED, fontsize=9.5,
                arrowprops=dict(arrowstyle="->", color=RED, lw=1.2))
    ax.set_xlabel(r"$k$  [Mpc$^{-1}$]")
    ax.set_ylabel(r"$P(k)$  [Mpc$^3$]")
    ax.legend(frameon=False, fontsize=10, loc="lower left")
    fig.savefig(p("fig_pk.png"))
    plt.close(fig)
    print("  fig_pk.png")

print("done")
