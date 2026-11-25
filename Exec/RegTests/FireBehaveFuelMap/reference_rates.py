#!/usr/bin/env python3
"""FireBehaveFuelMap: the closed-form reference rates of check_behave_fuel_map.py.

    python3 reference_rates.py

Multi-class Rothermel (1972) surface spread with Albini's (1976) size-class
weighting, as Andrews (2018, RMRS-GTR-371) gives it, written from those
papers and not from the FLARE code:

  * characteristic SAV sigma = sum_i f_i sigma_i over the dead and the live
    category, f_i the category area shares, sigma_i = sum_j f_ij sigma_ij;
  * net load per category from Albini's size-class weights g_ij;
  * live extinction moisture 2.9 W' (1 - M'_f/Mx_dead) - 0.226, with
    W' = sum w e^(-138/sigma) (dead) / sum w e^(-500/sigma) (live) and M'_f
    the dead moisture weighted by w e^(-138/sigma), per size class;
  * heat sink rho_b sum_i f_i sum_j f_ij e^(-138/sigma_ij) (250 + 1116 M_ij);
  * the dynamic herbaceous transfer (Scott and Burgan 2005): a share
    (1.20 - M_lh)/(1.20 - 0.30) of the live herbaceous load moves to the
    dead class, keeping its SAV and taking the 1-h moisture;
  * Rothermel's wind cap U <= 0.9 I_R (erf.fire.use_wind_limit = true).

Fuel parameters are Scott and Burgan (2005, RMRS-GTR-153) SH5 (145) and GR2
(102). Prints R0 and the rate at 3 m/s and 6.5 m/s midflame wind for the
deck's moistures, 3/4/5 % dead and 60 % live.
"""
import math

TPA = 2000.0 / 43560.0          # t/ac -> lb/ft^2
FT_MIN_TO_M_S = 0.00508
M_S_TO_FT_MIN = 196.85

# code: (w1, w10, w100, w_live_herb, w_live_woody [t/ac],
#        sigma_1h, sigma_live_herb, sigma_live_woody [1/ft], depth [ft], Mx_dead, h [Btu/lb])
FUELS = {
    145: ("SH5", (3.60, 2.10, 0.0, 0.0, 2.90, 750.0, 1500.0, 1600.0, 6.0, 0.15, 8000.0)),
    102: ("GR2", (0.10, 0.0, 0.0, 1.00, 0.0, 2000.0, 1800.0, 9999.0, 1.0, 0.15, 8000.0)),
}
SIGMA_10H, SIGMA_100H = 109.0, 30.0
RHO_P, S_T, S_E = 32.0, 0.0555, 0.010


def size_class(s):
    """Albini's size classes by SAV [1/ft]."""
    for k, edge in enumerate((1200.0, 192.0, 96.0, 48.0, 16.0)):
        if s >= edge:
            return k
    return 5


def eta_m(m, mx):
    r = min(1.0, m / mx)
    return 1.0 - 2.59 * r + 5.11 * r * r - 3.52 * r ** 3


def rate(p, m1, m10, m100, m_lh, m_lw, u_ms, cap=True):
    w1, w10, w100, wlh, wlw, s1, slh, slw, depth, mx_dead, h = p
    t = min(1.0, max(0.0, (1.20 - m_lh) / (1.20 - 0.30))) if wlh > 0 else 0.0
    dead = [(w1 * TPA, s1, m1), (w10 * TPA, SIGMA_10H, m10), (w100 * TPA, SIGMA_100H, m100),
            (t * wlh * TPA, slh, m1)]
    live = [((1.0 - t) * wlh * TPA, slh, m_lh), (wlw * TPA, slw, m_lw)]
    cats = [[c for c in dead if c[0] > 0], [c for c in live if c[0] > 0]]

    area = [[w * s / RHO_P for w, s, _ in cat] for cat in cats]
    a_cat = [sum(a) for a in area]
    a_tot = sum(a_cat)
    f_cat = [a / a_tot for a in a_cat]
    f_cls = [[a / a_cat[k] for a in area[k]] if a_cat[k] > 0 else [] for k in range(2)]
    sigma = sum(f_cat[k] * sum(f * c[1] for f, c in zip(f_cls[k], cats[k])) for k in range(2))

    w0 = sum(c[0] for cat in cats for c in cat)
    rho_b = w0 / depth
    beta = rho_b / RHO_P

    w_net, m_cat = [], []
    for k in range(2):
        g = {}
        for f, c in zip(f_cls[k], cats[k]):
            g[size_class(c[1])] = g.get(size_class(c[1]), 0.0) + f
        w_net.append(sum(g[size_class(c[1])] * c[0] * (1.0 - S_T) for c in cats[k]))
        m_cat.append(sum(f * c[2] for f, c in zip(f_cls[k], cats[k])))

    hn_dead = sum(w * math.exp(-138.0 / s) for w, s, _ in cats[0])
    hn_live = sum(w * math.exp(-500.0 / s) for w, s, _ in cats[1])
    mf_dead = sum(w * math.exp(-138.0 / s) * m for w, s, m in cats[0]) / hn_dead
    w_prime = hn_dead / hn_live if hn_live > 0 else 0.0
    mx = [mx_dead, max(mx_dead, 2.9 * w_prime * (1.0 - mf_dead / mx_dead) - 0.226)]

    eta_s = min(1.0, 0.174 * S_E ** -0.19)
    beta_op = 3.348 * sigma ** -0.8189
    gamma_max = sigma ** 1.5 / (495.0 + 0.0594 * sigma ** 1.5)
    a = 133.0 * sigma ** -0.7913
    ratio = beta / beta_op
    gamma = gamma_max * ratio ** a * math.exp(a * (1.0 - ratio))
    i_r = gamma * sum(w_net[k] * h * eta_m(m_cat[k], mx[k]) * eta_s for k in range(2) if cats[k])

    xi = math.exp((0.792 + 0.681 * sigma ** 0.5) * (beta + 0.1)) / (192.0 + 0.2595 * sigma)
    q_sink = sum(f_cat[k] * sum(f * math.exp(-138.0 / c[1]) * (250.0 + 1116.0 * c[2])
                                for f, c in zip(f_cls[k], cats[k])) for k in range(2))
    r0 = i_r * xi / (rho_b * q_sink)

    u = u_ms * M_S_TO_FT_MIN
    if cap:
        u = min(u, 0.9 * i_r)
    c_w = 7.47 * math.exp(-0.133 * sigma ** 0.55)
    b_w = 0.02526 * sigma ** 0.54
    e_w = 0.715 * math.exp(-3.59e-4 * sigma)
    phi_w = c_w * u ** b_w * ratio ** -e_w
    return r0 * FT_MIN_TO_M_S, r0 * (1.0 + phi_w) * FT_MIN_TO_M_S, sigma


if __name__ == "__main__":
    print("fuel  sigma[1/ft]  R0 [m/s]              R(3 m/s)              R(6.5 m/s)")
    for code, (name, p) in FUELS.items():
        r0, _, sigma = rate(p, 0.03, 0.04, 0.05, 0.60, 0.60, 0.0)
        r3 = rate(p, 0.03, 0.04, 0.05, 0.60, 0.60, 3.0)[1]
        r65 = rate(p, 0.03, 0.04, 0.05, 0.60, 0.60, 6.5)[1]
        print(f"{name} ({code})  {sigma:7.1f}  {r0!r:22s} {r3!r:22s} {r65!r}")
