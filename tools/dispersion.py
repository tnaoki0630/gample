import numpy as np
import matplotlib.pyplot as plt
from scipy.special import iv  # Modified Bessel function I_n

EPS0 = 8.8541878128e-12


def sine_over_x(x, t):
    """
    Return sin(x t) / x with the x -> 0 limit handled safely.
    """
    x = np.asarray(x, dtype=float)
    out = np.empty_like(x)

    mask = np.abs(x) < 1e-14
    out[mask] = t
    out[~mask] = np.sin(x[~mask] * t) / x[~mask]
    return out


def one_minus_cos_over_x2(x, t):
    """
    Return (1 - cos(x t)) / x^2 with the x -> 0 limit handled safely.

    Limit x -> 0:
        (1 - cos(x t)) / x^2 -> t^2 / 2
    """
    x = np.asarray(x, dtype=float)
    out = np.empty_like(x)

    mask = np.abs(x) < 1e-14
    out[mask] = 0.5 * t**2
    out[~mask] = (1.0 - np.cos(x[~mask] * t)) / (x[~mask] ** 2)
    return out


def compute_b_s(ky, rho_Ts):
    return (ky * rho_Ts) ** 2


def compute_n_res(omega, ky, uy, omega_cs):
    if np.abs(omega_cs) < 1e-30:
        raise ValueError("omega_cs is too close to zero.")
    return (omega - ky * uy) / omega_cs


def nmax_guess_from_physics(omega, ky, uy, omega_cs, rho_Ts, margin=20):
    """
    Physics-informed initial guess:
        n_max ~ max(b_s, |n_res|) + margin
    """
    b_s = compute_b_s(ky, rho_Ts)
    n_res = compute_n_res(omega, ky, uy, omega_cs)
    guess = int(np.ceil(max(b_s, abs(n_res)) + margin))
    return max(1, guess)


def friction_series_chunk(
    n_start,
    n_end,
    t,
    omega,
    ky,
    uy,
    omega_cs,
    rho_Ts,
    average=False,
):
    """
    Compute terms for n = n_start, ..., n_end inclusive, for either:

    instantaneous:
        sum_{n=-inf}^{inf} (n^2 / b_s) e^{-b_s} I_n(b_s) sin(Delta_n t)/Delta_n

    cumulative average:
        (1/t) * integral_0^t instantaneous(s) ds
      = sum_{n=-inf}^{inf} (n^2 / b_s) e^{-b_s} I_n(b_s) * (1-cos(Delta_n t)) / (t Delta_n^2)

    by explicitly combining +n and -n contributions through Delta_{+n}, Delta_{-n},
    while using I_{-n}(b)=I_n(b).
    """
    if n_end < n_start:
        return np.array([], dtype=int), np.array([], dtype=float)

    b_s = compute_b_s(ky, rho_Ts)
    if b_s == 0.0:
        raise ValueError("b_s = (k_y rho_Ts)^2 is zero; choose nonzero ky and rho_Ts.")

    n_vals = np.arange(n_start, n_end + 1, dtype=int)
    I_n = iv(n_vals, b_s)
    pref = (n_vals**2 / b_s) * np.exp(-b_s) * I_n

    delta_p = omega - ky * uy - n_vals * omega_cs   # +n
    delta_m = omega - ky * uy + n_vals * omega_cs   # -n

    if average:
        if np.abs(t) < 1e-30:
            # cumulative average at t=0 equals instantaneous value at t=0, which is 0
            terms = np.zeros_like(pref)
        else:
            terms = pref * (
                one_minus_cos_over_x2(delta_p, t) / t
                + one_minus_cos_over_x2(delta_m, t) / t
            )
    else:
        terms = pref * (sine_over_x(delta_p, t) + sine_over_x(delta_m, t))

    return n_vals, terms


def friction_value(
    t,
    Ehat_y,
    q_s,
    ky,
    lambda_Ds,
    omega,
    uy,
    omega_cs,
    rho_Ts,
    margin=20,
    atol=1e-12,
    rtol=1e-10,
    patience=20,
    chunk_size=50,
    hard_max_n=20000,
    average=False,
    return_diagnostics=False,
):
    """
    Compute either

    instantaneous:
        <δE_y δn_s>

    or its cumulative time average:
        (1/t) * integral_0^t <δE_y δn_s>(s) ds

    using:
      1) Initial upper bound guess from physics:
         n_max_guess ~ max(b_s, |n_res|) + margin
      2) Adaptive continuation:
         continue summing in chunks until the full term is sufficiently small
         for 'patience' consecutive n.
    """
    if np.abs(ky) < 1e-30:
        raise ValueError("ky is too close to zero.")
    if lambda_Ds <= 0.0:
        raise ValueError("lambda_Ds must be positive.")
    if chunk_size < 1:
        raise ValueError("chunk_size must be >= 1.")
    if patience < 1:
        raise ValueError("patience must be >= 1.")

    prefactor = -(EPS0 / (2.0 * q_s * ky * lambda_Ds**2)) * (np.abs(Ehat_y) ** 2)

    b_s = compute_b_s(ky, rho_Ts)
    n_res = compute_n_res(omega, ky, uy, omega_cs)
    n_guess = nmax_guess_from_physics(
        omega=omega, ky=ky, uy=uy, omega_cs=omega_cs, rho_Ts=rho_Ts, margin=margin
    )

    partial_sum = 0.0
    all_n = []
    all_terms = []

    small_count = 0
    used_n = 0

    n_start = 1
    n_end = min(n_guess, hard_max_n)

    while True:
        n_vals, terms = friction_series_chunk(
            n_start=n_start,
            n_end=n_end,
            t=t,
            omega=omega,
            ky=ky,
            uy=uy,
            omega_cs=omega_cs,
            rho_Ts=rho_Ts,
            average=average,
        )

        for n, term in zip(n_vals, terms):
            partial_sum += term
            all_n.append(n)
            all_terms.append(term)
            used_n = int(n)

            thresh = atol + rtol * abs(partial_sum)

            if abs(term) < thresh:
                small_count += 1
            else:
                small_count = 0

            if small_count >= patience:
                value = prefactor * partial_sum
                if return_diagnostics:
                    return value, {
                        "used_n": used_n,
                        "partial_sum": partial_sum,
                        "last_term": term,
                        "b_s": b_s,
                        "n_res": n_res,
                        "n_guess": n_guess,
                        "threshold": thresh,
                        "n_history": np.array(all_n, dtype=int),
                        "term_history": np.array(all_terms, dtype=float),
                    }
                return value

        if n_end >= hard_max_n:
            value = prefactor * partial_sum
            if return_diagnostics:
                return value, {
                    "used_n": used_n,
                    "partial_sum": partial_sum,
                    "last_term": all_terms[-1] if all_terms else np.nan,
                    "b_s": b_s,
                    "n_res": n_res,
                    "n_guess": n_guess,
                    "warning": "Reached hard_max_n before convergence.",
                    "n_history": np.array(all_n, dtype=int),
                    "term_history": np.array(all_terms, dtype=float),
                }
            return value

        n_start = n_end + 1
        n_end = min(n_end + chunk_size, hard_max_n)


def friction_vs_time(
    t_array,
    Ehat_y,
    q_s,
    ky,
    lambda_Ds,
    omega,
    uy,
    omega_cs,
    rho_Ts,
    margin=20,
    atol=1e-12,
    rtol=1e-10,
    patience=20,
    chunk_size=50,
    hard_max_n=20000,
    average=False,
):
    vals = np.empty_like(t_array, dtype=float)
    used_ns = np.empty_like(t_array, dtype=int)

    for i, t in enumerate(t_array):
        val, diag = friction_value(
            t=t,
            Ehat_y=Ehat_y,
            q_s=q_s,
            ky=ky,
            lambda_Ds=lambda_Ds,
            omega=omega,
            uy=uy,
            omega_cs=omega_cs,
            rho_Ts=rho_Ts,
            margin=margin,
            atol=atol,
            rtol=rtol,
            patience=patience,
            chunk_size=chunk_size,
            hard_max_n=hard_max_n,
            average=average,
            return_diagnostics=True,
        )
        vals[i] = val
        used_ns[i] = diag["used_n"]

    return vals, used_ns


if __name__ == "__main__":
    # -----------------------------
    # physical parameters from PIC
    # -----------------------------
    pi = np.atan(1)*4
    params = {
        "Ehat_y": 1.5e8,          # V/m
        "q_s": -1.602176634e-19,  # C
        "ky": 0.8e4/(2*pi),       # 1/m
        "lambda_Ds": 2.5e-2/512,  # m (using grid size in PIC alternate of debye length)
        "omega": 4e7,             # rad/s
        "uy": 3.0e6,              # m/s
        "omega_cs": 1.0e8,        # rad/s
        "rho_Ts": 2.5e-3,         # m
    }

    b_s = compute_b_s(params["ky"], params["rho_Ts"])
    n_res = compute_n_res(params["omega"], params["ky"], params["uy"], params["omega_cs"])
    n_guess = nmax_guess_from_physics(
        omega=params["omega"],
        ky=params["ky"],
        uy=params["uy"],
        omega_cs=params["omega_cs"],
        rho_Ts=params["rho_Ts"],
        margin=20,
    )

    print(f"b_s   = {b_s:.6g}")
    print(f"n_res = {n_res:.6g}")
    print(f"initial n_max guess = {n_guess}")

    # -----------------------------
    # Time sweep
    # -----------------------------
    ncyc = 4 # number of free cyclotron motion in chamber
    tmax = ncyc*(2*pi/params["omega_cs"])
    t_array = np.linspace(0.0, tmax, 1000)

    R_vals, used_ns_inst = friction_vs_time(
        t_array=t_array,
        margin=20,
        atol=1e-13,
        rtol=1e-9,
        patience=25,
        chunk_size=50,
        hard_max_n=20000,
        average=False,
        **params,
    )

    R_avg_vals, used_ns_avg = friction_vs_time(
        t_array=t_array,
        margin=20,
        atol=1e-13,
        rtol=1e-9,
        patience=25,
        chunk_size=50,
        hard_max_n=20000,
        average=True,
        **params,
    )

    plt.figure(figsize=(8, 5))
    plt.plot(t_array*params["omega_cs"]/(2*pi), R_vals, label=r"$\langle \delta E_y \delta n_s \rangle$")
    plt.plot(t_array*params["omega_cs"]/(2*pi), R_avg_vals, label="cumulative average")
    plt.xlabel(r"$t\ \frac{\omega_{ce}}{2\pi} [-]$")
    plt.ylabel("value")
    plt.title("Instantaneous vs cumulative-average response")
    plt.legend()
    plt.tight_layout()
    plt.show()

    plt.figure(figsize=(8, 4))
    plt.plot(t_array, used_ns_inst, label="instantaneous")
    plt.plot(t_array, used_ns_avg, label="cumulative average")
    plt.xlabel("t [s]")
    plt.ylabel("used n cutoff")
    plt.title("Adaptive harmonic cutoff")
    plt.legend()
    plt.tight_layout()
    plt.show()