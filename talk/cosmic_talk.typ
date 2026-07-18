#import "@preview/pinit:0.2.2": *

// ---------------------------------------------------------------- style ----
#let gold  = rgb("#f3a914")
#let deep  = rgb("#b06e0a")
#let ink   = rgb("#28282d")
#let red   = rgb("#be2828")
#let blue  = rgb("#1e50a0")

#set page(
  width: 25.4cm, height: 19.05cm, margin: (x: 1.6cm, y: 1.5cm),
  footer: context [
    #set text(8pt, fill: ink.lighten(35%))
    #grid(columns: (1fr, 1fr, 1fr),
      align: (left, center, right),
      [Cosmic.jl], [Kazi Abu Rousan], [#counter(page).display()]
    )
  ],
)
#set text(font: "New Computer Modern", size: 15pt, fill: ink)
#set par(justify: true, leading: 0.68em)
#show heading.where(level: 1): it => block(
  width: 100%, inset: (x: 10pt, y: 7pt), radius: 4pt, fill: gold,
  text(19pt, weight: "bold", fill: ink, it.body),
)
#show heading.where(level: 2): it => text(15pt, weight: "bold", fill: deep, it.body)
#let hl(body) = text(fill: deep, weight: "bold", body)
#let note(body) = block(
  width: 100%, inset: 8pt, radius: 3pt,
  fill: gold.lighten(80%), stroke: (left: 2pt + gold), text(13pt, body),
)
#let slide(title, body) = { pagebreak(weak: true); heading(level: 1, title); v(4pt); body }

// ---------------------------------------------------------------- title ----
#align(center + horizon)[
  #block(width: 100%, inset: 16pt, radius: 6pt, fill: gold)[
    #text(30pt, weight: "bold")[Cosmic.jl]
    #v(2pt)
    #text(15pt)[A cosmology and Einstein--Boltzmann engine in pure Julia]
  ]
  #v(18pt)
  #text(17pt, weight: "bold")[Kazi Abu Rousan]
  #v(4pt)
  #text(12pt, fill: deep)[`github.com/aburousan/Cosmic.jl`]
  #v(10pt)
  #text(11pt, fill: ink.lighten(30%))[Equations follow D. Baumann, _Cosmology_ (CUP, 2021)]
]

// =========================================================================
#slide[What this talk is about][
  I made a code called #hl[Cosmic.jl]. It is written fully in #hl[Julia].
  It does the same job as CLASS and CAMB, but in one language.

  In this talk I will show:

  - What a Boltzmann code actually solves (little bit of physics and maths).
  - Why I wrote one more code, when CLASS and CAMB are already very good.
  - Some features, with small code examples.
  - How close my results are to CLASS and CAMB.
  - What is still missing. I will say this part honestly.

  #note[
    One thing I want to say early. CLASS and CAMB are the *references* here.
    My code tries to reproduce them, not the other way. Where my numbers are
    different, I assume my code is wrong until I can prove otherwise.
  ]
]

// =========================================================================
#slide[The physics: two problems, one after another][
  #grid(columns: (1fr, 1fr), gutter: 16pt,
    [
      #text(14pt)[*1. Smooth universe*]

      How fast is the universe expanding? This comes from the Friedmann
      equation. If we write every component as a density, then

      $ H^2(a) = H_0^2 sum_s (rho_s (a)) / rho_(c,0) $

      Each component $s$ (photons, baryons, dark matter, neutrinos, dark
      energy) only needs to know its own $rho_s (a)$.
    ],
    [
      #text(14pt)[*2. Small ripples*]

      On top of the smooth part there are small perturbations. Photons are
      described by a distribution function, and it follows the Boltzmann
      equation

      $ (dif f) / (dif eta) = C[f] $

      where $C[f]$ is Thomson scattering off free electrons.
      #text(11pt, fill: deep)[(Baumann, App. B.1)]
    ]
  )
  #v(6pt)
  #note[
    Everything a Boltzmann code prints --- the CMB spectrum, the matter power
    spectrum --- comes from these two problems, plus one projection integral.
  ]
]

// =========================================================================
#slide[From ripples to what we observe][
  We expand the photon temperature field in multipoles $Theta_ell (k, eta)$.
  After decoupling the photons just free-stream. Then the observed angular
  power spectrum is one integral over $k$:

  #v(4pt)
  $ C_ell = 4 pi integral (dif k) / k #pin(1)P_cal(R)(k)#pin(2)
        underbrace(|Theta_ell (k, eta_0)|^2, "transfer function") $
  #v(4pt)

  // transparentize so the highlighted symbol still reads through the box
  #pinit-highlight(1, 2, fill: gold.transparentize(65%))
  #pinit-point-from(2, offset-dx: 6pt, offset-dy: 46pt, fill: deep)[
    #text(12pt, fill: deep)[primordial spectrum from inflation (Baumann §8.3)]
  ]

  #v(46pt)
  So there are exactly three pieces:

  + *Initial conditions* --- what inflation gives us.
  + *Evolution* --- the Boltzmann + Einstein system, until today.
  + *Projection* --- turning a 3D $k$ into a 2D angle $ell$
    #text(11pt, fill: deep)[(Baumann §7.3.1, App. B.2)].

  A Boltzmann code is basically a careful machine for these three steps.
]

// =========================================================================
#slide[Why write one more code?][
  CLASS (C) and CAMB (Fortran) are excellent and very well tested. I did not
  write Cosmic.jl to replace them. I wrote it for three reasons.

  #v(6pt)
  #grid(columns: (1fr,), gutter: 10pt,
    note[
      *1. One language from top to bottom.* In CLASS and CAMB the fast part is
      C or Fortran, and the part you actually use is a Python wrapper. Julia is
      fast enough to be both. So the code you read is the code that runs.
    ],
    note[
      *2. Composition should be generic.* A cosmology here is just a *list of
      species*. To add a new component you write one function, its $rho(a)$.
      You do not touch ten other files.
    ],
    note[
      *3. Compute, do not fit.* Many codes use fitting formulae for hard
      integrals. Where I can, I compute the actual integral, and keep the fit
      only as a cross-check.
    ],
  )
]

// =========================================================================
#slide[The design, in one screen][
  ```julia
  using Cosmic

  c = cosmology()                    # Planck 2018 LCDM
  age(c)                             # 13.787 Gyr
  luminosity_distance(c, 1.0)        # Mpc

  rec = recombination(c)             # thermal history
  sp  = cmb_spectra(c, rec; lmax=1500, nk=20000)
  mp  = matter_power_spectrum(c, rec; kmin=1e-3, kmax=10.0)
  ```
  #v(6pt)
  That is the whole interface for a standard run. `cosmology()` builds the
  species list, `recombination` solves the thermal history, and the last two
  do the perturbations and the projection.
]

// =========================================================================
#slide[Distances: the easy check first][
  #align(center)[#image("fig_distances.png", width: 68%)]
  #v(-2pt)
  Luminosity distance $d_L$, angular diameter distance $d_A$, comoving
  distance $d_C$. Note $d_A$ turns over near $z approx 1.6$: far away objects
  start looking *bigger* again, because the universe was smaller when the
  light left.

  #text(12pt, fill: deep)[Against CLASS these agree to about $10^(-9)$.]
]

// =========================================================================
#slide[The CMB spectrum, and what the bumps mean][
  #grid(columns: (1.35fr, 1fr), gutter: 12pt,
    align(center + horizon)[#image("fig_cls.png", width: 100%)],
    [
      #text(13pt)[
        Before decoupling, photons and baryons are one fluid. Gravity pulls it
        in, radiation pressure pushes it out. This gives *sound waves*
        #text(11pt, fill: deep)[(Baumann §7.4)].

        - The *first peak* is the mode that had just enough time for one
          compression.
        - Peak *positions* tell us the geometry.
        - Peak *heights* tell us how much baryon and dark matter there is
          #text(11pt, fill: deep)[(§7.5.1--7.5.2)].
        - The fall at high $ell$ is *Silk damping*: photons diffuse out of
          small overdensities #text(11pt, fill: deep)[(§7.4.4)].
      ]
    ]
  )
]

// =========================================================================
#slide[Cosmic.jl compared with CLASS][
  #align(center)[#image("fig_cls_compare.png", width: 56%)]
  #v(-2pt)
  Same cosmological parameters on both sides, same $Y_p$, and no reionization
  on either side. Top: the two curves sit on top of each other. Bottom: the
  difference in percent. It stays under #hl[0.5%] below $ell approx 1000$, and
  under 1.1% up to $ell = 1500$.

  #note[
    #text(13pt)[
      *First make it converge, then compare.* My first try at this plot had 8%
      ringing near the first peak. It was not physics. The code itself had
      warned me that the $k$ grid was too coarse ($Delta k = 7.9 times 10^(-5)$,
      but $1.2 times 10^(-5)$ was needed). After raising `nk` to 20000 the
      ringing was gone. So always converge both codes before you call a
      difference a bug.
    ]
  ]
]

// =========================================================================
#slide[The matter power spectrum][
  #grid(columns: (1.35fr, 1fr), gutter: 12pt,
    align(center + horizon)[#image("fig_pk.png", width: 100%)],
    [
      #text(13pt)[
        $P(k)$ tells us how clumpy matter is at scale $k$
        #text(11pt, fill: deep)[(Baumann §5.3.4)].

        - It rises, then turns over. The turnover is the scale that entered
          the horizon at matter--radiation equality.
        - Modes entering earlier grow more slowly, because radiation stops
          them. So small scales are suppressed.
        - The small wiggles are the same sound waves as in the CMB.
        - At large $k$ gravity is no longer linear, so we need a non-linear
          correction (halofit here).
      ]
    ]
  )
]

// =========================================================================
#slide[Some features I like][
  #grid(columns: (1fr, 1fr), gutter: 14pt,
    note[
      *BBN is solved, not looked up.*
      CLASS reads $Y_p$ from a table, CAMB uses a fitting formula. Here the
      12-nuclide, 63-reaction network is integrated in the same background,
      so the helium that goes into recombination is the helium your own
      $omega_b$ produced.
    ],
    note[
      *Curvature done exactly.*
      Hyperspherical Bessel functions are solved from their own recurrences at
      every $nu$. CLASS switches to rescaled flat Bessels above $nu = 4000$;
      here that switch never happens.
    ],
    note[
      *Spectral distortions from first principles.*
      Exact Karzas--Latter Gaunt factor, exact double-Compton emissivity, and
      the exact Klein--Nishina $times.o$ Maxwell--Jüttner kernel instead
      of the Kompaneets approximation.
    ],
    note[
      *Early-universe gravitational waves.*
      Scalar-induced GWs and primordial black hole abundances, including
      non-Gaussianity ($f_"NL"$, $g_"NL"$).
    ],
  )
]

// =========================================================================
#slide[CLASS, CAMB and Cosmic.jl side by side][
  #set text(12.5pt)
  #table(
    columns: (auto, 1fr, 1fr, 1fr),
    stroke: none,
    inset: 6pt,
    fill: (_, y) => if y == 0 { gold } else if calc.odd(y) { gold.lighten(88%) },
    table.header([], [*CAMB*], [*CLASS*], [*Cosmic.jl*]),
    [First release], [2000], [2011], [in development],
    [Language], [Fortran 90], [C], [Julia],
    [Interface], [Python wrapper], [`classy` wrapper], [native],
    [Recombination], [RECFAST, HyRec, CosmoRec], [RECFAST, HyRec], [RECFAST-class, HyRec-2020],
    [BBN], [fitting formula], [interpolated table], [solved network],
    [Curvature], [exact], [flat approx. above $nu = 4000$], [exact at all $nu$],
    [Non-linear], [halofit, HMcode], [halofit, HMcode], [halofit, HMcode, EFT],
    [Distortions], [---], [fit based], [first principles],
    [Maturity], [very high], [very high], [#text(fill: red)[young]],
  )
]

// =========================================================================
#slide[How close are the numbers?][
  #set text(13pt)
  Cross-checks at *matched settings*, against live CLASS 3.3.4 and CAMB runs:
  #v(4pt)
  #table(
    columns: (1.5fr, 1fr, 1fr),
    stroke: none, inset: 5.5pt,
    fill: (_, y) => if y == 0 { gold } else if calc.odd(y) { gold.lighten(88%) },
    table.header([*Quantity*], [*Reference*], [*Agreement*]),
    [Distances], [CLASS], [$10^(-9)$],
    [$P(k)$, $sigma_8$], [CAMB], [0.04%, 0.007%],
    [$C_ell^(T T)$ unlensed], [CLASS], [$approx 0.4%$],
    [lensed $T T$ / $E E$], [CLASS], [0.3% / 0.4%],
    [$C_ell^(phi phi)$, $ell <= 1000$], [CLASS], [$<= 0.9%$],
    [halofit $P_"nl" \/ P_"lin"$], [CLASS], [0.25%],
    [EFT 1-loop], [FAST-PT], [$< 0.5%$],
    [$Y_p$], [CLASS (PRIMAT)], [$< 0.1%$],
    [$Y_p$, D/H, $""^7$Li], [PRyMordial], [$-0.004%$, $+0.03%$, $-1.1%$],
    [$mu$, $y$ distortions], [CLASS], [0.02%],
    [open / closed $C_ell^(T T)$], [CLASS / CAMB], [$<= 0.6%$ / $<= 0.4%$],
  )
]

// =========================================================================
#slide[What is missing, and what is approximate][
  I think a code is easier to trust if it says its own weak points.

  #v(6pt)
  #note[
    *Not implemented yet:* vector modes; number-count and shear spectra;
    late-time $y$-distortion; scalar-field and interacting dark energy
    perturbations; a full second-order (bispectrum) solver.
  ]
  #v(6pt)
  #note[
    *A known approximation.* For tensor and vector modes I use the "optimal"
    Boltzmann hierarchy, which re-expands them on scalar normal modes. Pitrou,
    Pereira & Lesgourgues (PRD *102*, 023511) proved that the separability this
    depends on #hl[fails when the space is curved]. For scalars the error is
    tiny, but for *tensor polarisation* it is about $50 |Omega_K| %$.

    So curved-space tensor $E E$ and $B B$ in my code are not exact. I keep
    this written down instead of hiding it.
  ]
]

// =========================================================================
#slide[Summary][
  - Cosmic.jl solves the full chain in one language: background $arrow.r$ BBN
    $arrow.r$ recombination $arrow.r$ perturbations $arrow.r$ CMB and $P(k)$.

  - It agrees with CLASS and CAMB at the few $times 10^(-3)$ level on
    distances, $C_ell$, $P(k)$, BBN abundances and spectral distortions.

  - What makes it different: species-based design, exact curved geometry, a
    real BBN network, and distortion physics computed instead of fitted.

  - It is still young. CLASS and CAMB stay the references, and the open
    approximations are written down, not hidden.

  #v(14pt)
  #align(center)[
    #block(inset: 12pt, radius: 5pt, fill: gold)[
      #text(15pt, weight: "bold")[github.com/aburousan/Cosmic.jl]
    ]
    #v(6pt)
    #text(13pt)[Kazi Abu Rousan]
  ]
]
