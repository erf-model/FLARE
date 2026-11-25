# Single span after Ma, Bocchini and Christou 2020

## Purpose

One conductor span in a crosswind, set up after the single-span study of Ma, Bocchini and
Christou 2020 (Section 3), to compare the conductor model with a published case:

1. **Mean tension.** In a steady uniform wind, does the run's end tension match a closed-form
   catenary and an independent lumped chain at rest? (`run_check.py`, CTest)
2. **Peak tension.** Under a turbulent wind, how do the run's peak and minimum tensions over 10 min
   compare with the FEM values in Ma et al.'s Table 2? Ma et al.'s intensity measure V_ms is the
   highest 1-min mean wind at 10 m. Here the stripe is the inflow's wind at 10 m, a long mean,
   which is lower than V_ms in the same wind by a gust factor that is not applied.
3. **Fragility.** How does a curve fitted to the runs (multiple stripes, Baker 2015) compare with
   Ma et al.'s curves?

The setup, from Ma et al. Section 3 and Table 3:

| Item | Value | Source |
|---|---|---|
| Conductor | ACSR Drake 26/7, 15.95 N/m, rated strength 140.18 kN | Ma 2020 Table 3 |
| Supports | hinged, at the same level | Ma 2020 Section 3 |
| Sag at rest | span / 50 | Ma 2020 Section 3 |
| Drag coefficient, air density | 1.0, 1.226 kg/m³ | Ma 2020 Section 3 |
| Wind | normal to the span, log law with z0 = 0.05 m | Ma 2020 Section 3 |
| Conductor height | 40 m, read here as the span's reference height, two thirds of the sag below the supports (supports at 44 m for 300 m) | Ma 2020 Section 3; the reading is ours |
| Span | 300 m | Ma 2020 Fig. 9 uses 300 m; the span behind Table 2 is not stated next to the table, and 300 m is assumed |
| Diameter | 0.0281 m | Ma 2020 Table 3: 26 aluminium strands of 4.4 mm and 7 steel strands of 3.5 mm give 3 × 3.5 + 4 × 4.4 = 28.1 mm |
| Axial stiffness EA | 3.0e7 N | the conductor model's tests; **not a supplier figure** (Ma 2020 Table 3 gives none) |
| Segments | 80 | Ma 2020 use at least 80 elements |

## Files

- `static_reference.py`: the unstretched length that gives a sag of span / 50 under gravity, and
  the end tension in a steady uniform wind three ways: the elastic catenary
  (`Structures/Tools/broken_wire`, standard library) under the resultant of weight and full drag;
  a lumped chain at rest in 3-D whose drag takes only the wind normal to each link (numpy); and
  the same chain with MoorDyn's conventions (drag on the unstretched length, weight less the
  displaced air's; MoorDyn-C 2.7.1 `source/Line.cpp`).
- `make_case.py`: writes a deck. `--mode prescribed` is still air with a uniform wind handed to
  every node (`conductors.prescribed_velocity`). `--mode flow` is a neutral boundary layer over
  flat ground (k-equation RANS, surface layer with z0, log-law inflow at x-low) with random gusts
  (`conductors.gust_type = random`) and the wind at 10 m recorded beside the peak tension
  (`conductors.intensity_measure = wind`, 600 s average).
- `run_check.py`: runs the prescribed decks at 20, 60 and 100 m/s and checks the mean tension
  (below). Used by the CTest.
- `make_plan.py`: writes the plan of the fragility ensemble driver (one stripe per wind, one
  member per gust seed, a shared spin-up per stripe) and of the fit.

## Data and acknowledgements

No data set is downloaded or redistributed. The setup and the published values compared with
here are from

- Ma, L., Bocchini, P. and Christou, V. (2020). Fragility models of electrical conductors in power
  transmission networks subjected to hurricanes. *Structural Safety* 82, 101890.
  doi:10.1016/j.strusafe.2019.101890

We thank the authors for the detail of their Section 3 and Table 2, which make this comparison
possible. Values read from their figures are labelled as read by eye.

The fit follows Baker, J. W. (2015). Efficient analytical fragility function fitting using dynamic
structural analysis. *Earthquake Spectra* 31(1), 579-599. doi:10.1193/021113EQS025M

## How to run

Build FLARE with the real MoorDyn-C (the stub has no dynamics). The mean-tension check:

```bash
python3 run_check.py --exe path/to/flare --out check
```

Add `--without-fft` to `run_check.py` and `make_plan.py` for an ERF built without FFT.

The stripes (each member is 3 flow-throughs of spin-up, shared per stripe, then 660 s):

```bash
python3 make_plan.py --out stripes --winds 40 50 60 70 80 --seeds 1 2 3 4 5
python3 ../../Tools/fragility_ensemble/fragility_ensemble.py run --plan stripes/plan.json --exe path/to/flare --ranks 2
python3 ../../Tools/fragility_ensemble/fragility_ensemble.py collect --plan stripes/plan.json
python3 ../../Tools/fragility_fit/fragility_fit.py --plan stripes/fit_plan.json
```

## Expected results

### Mean tension in a steady wind

Drake on 300 m, unstretched length 300.0204 m, EA 3.0e7 N, 80 segments, MoorDyn-C 2.7.1 on one
rank (November 24, 2026). End tension averaged over the last 20 s of 120 s:

| Wind (m/s) | Run (N) | Chain, MoorDyn conventions (N) | Difference | Catenary (N) | Difference | Swing, run / atan(f/w) |
|---|---|---|---|---|---|---|
| 20 | 31766.5 | 31766.5 | 0.0001 % | 31787.7 | -0.07 % | 23.37° / 23.36° |
| 60 | 76506.3 | 76506.1 | 0.0002 % | 76991.5 | -0.63 % | 75.57° / 75.58° |
| 100 | 148786.7 | 148786.0 | 0.0005 % | 150713.9 | -1.28 % | 84.70° / 84.71° |

The run reproduces the chain with MoorDyn's conventions. The catenary takes the whole wind as
normal to the cable. The chain takes off the component along each sloping link, which on its own
lowers the tension by 1.28 % at 100 m/s, and puts its drag on the stretched length, which restores
about 0.33 %; net, the chain is 0.03 % (20 m/s) to 1.0 % (100 m/s) below the catenary. MoorDyn's
drag on the unstretched length and the air's buoyancy take another 0.03 % to 0.33 %.

**Sensitivity to EA** (chain with MoorDyn conventions, `static_reference.py --EA`):

| Wind (m/s) | EA 1.5e7 N | EA 3.0e7 N | EA 6.0e7 N |
|---|---|---|---|
| 40 | 43.7 kN | 47.8 kN | 51.6 kN |
| 60 | 66.2 kN | 76.5 kN | 87.2 kN |
| 80 | 93.3 kN | 110.9 kN | 130.2 kN |

Halving or doubling EA moves the tension by 8 % to 17 %. Without a supplier's EA for Drake this
is the largest open input of the comparison.

### Peak tension and fragility under a turbulent wind

Drake on 300 m, 5 stripes (inflow wind at 10 m of 40 to 80 m/s) by 5 gust seeds, each member
600 s of statistics after the shared spin-up and 60 s, MoorDyn-C 2.7.1 on 2 ranks (November 25,
2026; 46 min in all). Peak tension is the driver's `peak_demand` (every step); minimum and mean
tension and the wind at mid-span are from `conductors/S.dat` (every 1 s). The mean wind at
mid-span matches the log law at 40 m (50.46 m/s against 50.43 m/s at 40 m/s).

**Peak and minimum tension against Ma et al.'s Table 2**, with the stripe set equal to Table 2's
V_ms:

| Wind (m/s) | Peak, mean of 5 (range) (kN) | Table 2 peak (kN) | Difference | Minimum, mean of 5 (kN) | Table 2 minimum (kN) | Mean tension (kN) |
|---|---|---|---|---|---|---|
| 40 | 79.3 (77.2-81.4) | 81.6 | -2.8 % | 48.2 | 53.0 | 62.1 |
| 50 | 106.6 (105.7-108.2) | 109.1 | -2.3 % | 63.1 | 65.9 | 82.1 |
| 60 | 135.7 (134.3-137.4) | 138.8 | -2.2 % | 78.4 | 79.7 | 104.2 |
| 70 | 166.5 (164.2-168.2) | 170.9 | -2.6 % | 96.2 | 95.2 | 127.6 |
| 80 | 198.6 (195.9-200.7) | 205.3 | -3.3 % | 114.2 | 111.3 | 152.4 |

The peaks lie 2 % to 3 % below Table 2 at every wind. Table 2's caption names V_ms, the highest
1-min mean wind at 10 m, while the runs' stripe is a long mean. The table matches the runs when
the two are taken as equal. Our reading, not stated in the paper, is that Table 2's FEM runs took
V_ms as the mean wind without the conversion to a 10-min mean that the paper describes for its
fragility analysis.

**Fragility.** Against Drake's rated strength (140.18 kN), no member exceeds it at 60 m/s and all
five do at 70 m/s, so the Baker 2015 fit is separated (`fragility_fit.py` reports a median between
60 and 70 m/s). A straight line through log peak tension against log stripe (25 members, peak
proportional to the wind to the power 1.33, residual standard deviation 0.011) reaches 140.18 kN at
61.5 m/s, with a log standard deviation of 0.009 from the gusts alone. Table 2's peaks, interpolated
the same way, reach 140.18 kN at 60.4 m/s.

Ma et al.'s Fig. 11 (Drake, 0° yaw; read by eye) puts the median near V_ms = 72 m/s, about
1e-8 near 57.5 m/s and about 1 near 75 m/s. The ratio of the medians, 72 / 61.5 = 1.17, has the size
of a gust factor between a 1-min and a long mean, which the paper converts with the ESDU method
but does not report; the comparison of medians therefore needs that factor. The runs' spread is
narrower than Fig. 11's for three reasons:

- the capacity here is the rated strength, a single value, where Ma et al. take a distribution
  (their Fig. 10 peaks near 143 to 144 kN, read by eye);
- each member has one fully coherent gust along the span, where Ma et al. use a Kaimal spectrum
  with Davenport coherence; a coherent gust raises the peaks and lowers the minima of the span's
  load;
- 5 seeds per stripe.

**What is not settled.** The gust factor between V_ms and the runs' wind; EA (see above; halving or
doubling it moves the mean tension by 8 % to 17 %); the span behind Table 2 (300 m assumed); and
the height of the reference wind (40 m read as the span's reference height; reading it as the
supports' height lowers the drag by about 3 %).

## CTest

`SingleSpanMa2020_MeanTension` (labels `regression`, `conductors`; MoorDyn-C only) runs
`run_check.py` on one rank (with `--without-fft` when `ERF_ENABLE_FFT` is off) and passes when it
prints `SingleSpanMa2020: all checks hold`; the runs are removed after a pass. The end
tensions must lie within 0.01 % of the chain with MoorDyn's conventions, within 1.5 % of the
catenary, and the swing within 0.05° of atan(f/w). The chain's values are pinned in
`run_check.py`; with numpy present they are recomputed and the pins checked to 0.001 %. It takes
about 14 s.

**The check fails when it should.** Two deliberate changes to `make_case.py`:

| Change | What the check reports |
|---|---|
| Drag coefficient in the deck 0.5 % high | 0.05 % off the chain at 20 m/s (and the swing 0.11° off), 0.33 % at 100 m/s |
| Mass per length from g = 9.80665 instead of 9.81 | 0.02 % off the chain at 20 m/s |
