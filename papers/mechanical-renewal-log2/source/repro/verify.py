#!/usr/bin/env python3
"""Independent finite audit of the mechanical and Qin--Sun counting ledgers.

Python standard library only. Default execution reads but never writes files.
No computation here certifies period rank, persistence, or limit cycles.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import math
from fractions import Fraction as F
from pathlib import Path

if not __debug__:
    raise RuntimeError("Run without -O or -OO; optimized execution is forbidden.")


def check(ok: bool, label: str) -> None:
    if not ok:
        raise RuntimeError(label)


def prefixes(n: int) -> list[int]:
    """Build by reading the binary word; cross-check the floor formula."""
    check(n >= 2, "half-degree must be at least 2")
    word = bin(n)[2:]
    values = [1]
    for bit in word[1:]:
        values.append(2 * values[-1] + int(bit))
    r = len(word) - 1
    check(values == [n // 2 ** (r - i) for i in range(r + 1)], "prefix floor identity")
    return values


def mechanical(n: int, enumerate_modes: bool = False) -> dict:
    ns = prefixes(n)
    r, total, steps = len(ns) - 1, 0, []
    weighted_square_sum = weighted_main = weighted_B = 0
    for i, (s, t) in enumerate(zip(ns, ns[1:]), 1):
        eta, birth, squared, linear = t - 2 * s, 0, 0, 0
        check(eta in (0, 1), "binary digit")
        check(max(4 * s, 2 * t) == 2 * t, "Hamiltonian composition/padding degree")
        check(2 * (2 * s - 1) + 1 <= 2 * t - 1, "old form degree")
        layers = []
        for k in range(i):
            A = 2 ** k
            d, e = divmod(s, A)
            check(d == ns[i - 1 - k], "birth half-degree is a prefix")
            modes = d * (t - A * d)
            check(d >= 1 and modes >= 1, "positive canonical dimension")
            check(A * modes == s * (s + eta) - e * e - eta * e,
                  "copied dimension identity")
            if enumerate_modes:
                counted, exponents = 0, set()
                for q in range(d):
                    coefficient_degree = 2 * q
                    for _ in range(k):
                        coefficient_degree = 2 * coefficient_degree + 1
                    check(coefficient_degree == A * (2 * q + 1) - 1,
                          "iterated auxiliary coefficient degree")
                    p = 0
                    while p + coefficient_degree <= t - 2:
                        source_degree = (2 * p + 1) + 2 * coefficient_degree + 1
                        check(source_degree <= 2 * t - 2 < 2 * t - 1,
                              "final source total degree")
                        exponent = 2 * d * p + d + 2 * q + 1
                        check(exponent not in exponents, "distinct scaled period exponents")
                        exponents.add(exponent)
                        counted += 1
                        p += 1
                    check(p == t - A * (2 * q + 1), "exact p-range")
                check(counted == modes, "enumerated canonical modes equal closed count")
            birth += 4 * A * (modes - 1)
            squared += e * e
            linear += e
            layers.append({"k": k, "A": A, "d": d, "e": e, "modes": modes})
        closed_birth = i * (t * t - eta) - 4 * (squared + eta * linear) - 4 * (2 ** i - 1)
        check(birth == closed_birth and birth >= 0, "closed birth identity/nonnegativity")
        check(3 * squared <= s * s, "binary squared-remainder bound")
        check(linear <= 2 ** i - 1 - i, "binary linear-remainder bound")
        B_i = i * eta + 4 * eta * linear + 4 * (2 ** i - 1)
        check(3 * birth >= (3 * i - 1) * t * t - 3 * B_i,
              "strong per-step lower estimate")
        weighted_square_sum += 4 ** (r - i) * squared
        weighted_main += i * 4 ** (r - i) * t * t
        weighted_B += 4 ** (r - i) * B_i
        total = 4 * total + birth
        steps.append({"i": i, "s": s, "t": t, "eta": eta,
                      "field_degree": 2 * t - 1, "birth": birth,
                      "total": total, "layers": layers})
    intermediate = F(n * n * (3 * r * r + r) - 20 * n * 2 ** r - 24 * 4 ** r, 6)
    uniform = F(n * n * (3 * r * r + r - 44), 6)
    check(total >= intermediate >= uniform, "uniform finite estimates")
    if n == 2 ** r:
        closed = F(n * n * r * (r + 1), 2) - F(8 * 4 ** r, 3) + 4 * 2 ** r - F(4, 3)
        check(total == closed, "dyadic closed form")
    alpha = F(n, 2 ** r)
    profile_F = sum((F((n * 2 ** ell) % 2 ** r, 2 ** r) ** 2 / 4 ** (ell + 1)
                     for ell in range(r)), F(0))
    beta = F(1, 2) - 4 * profile_F / (alpha * alpha)
    Q = r * profile_F - F(weighted_square_sum, 4 ** r)
    E0 = F(weighted_main, 4 ** r) - alpha * alpha * r * (r + 1) / 2
    B = F(weighted_B, 4 ** r)
    check(-4 * alpha <= E0 <= 0 and 0 <= B <= 4, "profile auxiliary errors")
    reconstructed = alpha * alpha * r * (r + 1) / 2 - 4 * r * profile_F + E0 + 4 * Q - B
    check(F(total, 4 ** r) == reconstructed, "exact finite profile decomposition")
    check(0 <= Q <= F(19, 9), "finite binary-profile remainder")
    check(F(1, 6) <= beta <= F(1, 2), "binary-profile beta interval")
    profile = n * n * (F(r * r, 2) + beta * r)
    check(-8 * n * n <= total - profile <= F(76, 9) * 4 ** r,
          "uniform binary-profile error")
    return {"n": n, "binary": bin(n)[2:], "prefixes": ns,
            "depth": r, "degree": 2 * n - 1, "certificate": total,
            "uniform_lower": str(uniform), "intermediate_lower": str(intermediate),
            "profile_F": str(profile_F), "profile_beta": str(beta), "steps": steps}


def qs_seeds(D: int):
    """The seed theorem is an input; only its displayed integers are checked."""
    r, q = 2, 9
    while 4 * q - 5 <= D:
        yield r, 4 * q - 5, 4 * r * r * (q // 3) ** 2 - (q - 1) ** 2
        r, q = r + 1, 3 * q


def qs_envelope(D: int) -> dict:
    direct, best, chosen = 0, 0, None
    for r, d, cycles in qs_seeds(D):
        direct = max(direct, cycles)
        m = (D + 1) // (d + 1)
        if m * m * cycles > best:
            best, chosen = m * m * cycles, {"r": r, "seed_degree": d,
                "seed_cycles": cycles, "m": m, "realized_degree": m * (d + 1) - 1}
    check(best >= direct, "replication dominates direct seed selection")
    if chosen:
        check(chosen["realized_degree"] <= D, "replication degree threshold")
    return {"direct_seed_bound": direct, "replication_bound": best, "chosen": chosen}


def compare(D: int) -> dict:
    n = (D + 1) // 2
    row = mechanical(n)
    check(row["degree"] <= D < row["degree"] + 2, "target-degree half-degree threshold")
    out = {"D": D, "mechanical_n": n, "mechanical_degree": row["degree"],
           "mechanical_bound": row["certificate"], **qs_envelope(D)}
    scale = D * D * math.log(D) ** 2
    out["ratios_descriptive_only"] = {key: format(out[key] / scale, ".12g")
        for key in ("mechanical_bound", "direct_seed_bound", "replication_bound")}
    return out


def audit_replication(max_degree: int) -> None:
    check(list(qs_seeds(319)) == [(2, 31, 80), (3, 103, 2240), (4, 319, 40256)],
          "seed fixtures")
    # Independent order: enumerate realizations, then form the monotone closure.
    realized = [0] * (max_degree + 1)
    for _, d, cycles in qs_seeds(max_degree):
        for m in range(1, (max_degree + 1) // (d + 1) + 1):
            degree = m * (d + 1) - 1
            realized[degree] = max(realized[degree], m * m * cycles)
        for a in range(1, 13):
            for b in range(1, 13):
                once = a * (d + 1) - 1
                check(b * (once + 1) - 1 == a * b * (d + 1) - 1,
                      "composed replication degree")
                check(b * b * (a * a * cycles) == (a * b) ** 2 * cycles,
                      "composed replication count")
    previous = 0
    for D, value in enumerate(realized):
        previous = max(previous, value)
        check(previous == qs_envelope(D)["replication_bound"], "independent degree envelope")
        if D >= 3:
            n = (D + 1) // 2
            check(2 * n - 1 <= D < 2 * n + 1, "all target-degree thresholds")


def dataset(max_n: int, direct_through: int, max_degree: int) -> dict:
    keep = {2, 3, 4, 5, 7, 8, 15, 16, 31, 32, 64, 100, 256, 1000, 10000, 20000}
    examples = []
    for n in range(2, max_n + 1):
        row = mechanical(n, enumerate_modes=n <= direct_through)
        if n in keep:
            examples.append(row)
    audit_replication(max_degree)
    constants = {"mechanical_all_degree": 1 / (8 * math.log(2) ** 2),
        "QinSun_subsequence": 1 / (36 * math.log(3) ** 2),
        "QinSun_direct_all_degree": 1 / (324 * math.log(3) ** 2),
        "QinSun_replication_all_degree": 1 / (36 * math.log(3) ** 2)}
    constants["mechanical_divided_by_QinSun_subsequence"] = (
        constants["mechanical_all_degree"] / constants["QinSun_subsequence"])
    return {"schema": 1, "status": "PASS", "scope": "finite discrete identities only; no analytic certification",
        "range": {"n_min": 2, "n_max": max_n, "explicit_pq_enumeration_through": min(max_n, direct_through),
                  "replication_degree_min": 0, "replication_degree_max": max_degree},
        "constants_approximate_descriptive_only": {k: format(v, ".15g") for k, v in constants.items()},
        "mechanical_examples": examples,
        "degree_comparison": [compare(D) for D in (3, 5, 7, 9, 13, 15, 31, 63, 100, 103,
            319, 1000, 1023, 10000, 10 ** 6, 10 ** 9, 10 ** 12, 10 ** 18)]}


def table(data: dict) -> str:
    output = io.StringIO()
    fields = ["D", "mechanical_n", "mechanical_degree", "mechanical_bound",
              "direct_seed_bound", "replication_bound", "seed_r", "copies_m", "replicated_degree"]
    writer = csv.DictWriter(output, fieldnames=fields, lineterminator="\n")
    writer.writeheader()
    for row in data["degree_comparison"]:
        chosen = row["chosen"] or {}
        writer.writerow({**{key: row[key] for key in fields[:6]},
            "seed_r": chosen.get("r", ""), "copies_m": chosen.get("m", ""),
            "replicated_degree": chosen.get("realized_degree", "")})
    return output.getvalue()


def exact_payload(value):
    """Approximate display values never participate in a validation decision."""
    if isinstance(value, dict):
        return {key: exact_payload(item) for key, item in value.items()
                if not key.endswith("descriptive_only")}
    if isinstance(value, list):
        return [exact_payload(item) for item in value]
    return value


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-n", type=int, default=20000)
    parser.add_argument("--direct-through", type=int, default=256)
    parser.add_argument("--max-degree", type=int, default=20000)
    parser.add_argument("--write-data", action="store_true", help="explicitly regenerate adjacent data.json")
    parser.add_argument("--table", action="store_true", help="print CSV to stdout after verification")
    args = parser.parse_args()
    if args.max_n < 2 or args.direct_through < 2 or args.max_degree < 0:
        parser.error("require max-n>=2, direct-through>=2, max-degree>=0")
    data = dataset(args.max_n, args.direct_through, args.max_degree)
    target = Path(__file__).with_name("data.json")
    if args.write_data:
        target.write_text(json.dumps(data, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    elif target.is_file():
        saved = json.loads(target.read_text(encoding="utf-8"))
        if saved["range"] == data["range"]:
            check(exact_payload(saved) == exact_payload(data), "saved exact data equal recomputation")
    if args.table:
        print(table(data), end="")
    else:
        print(f"PASS: mechanical n=2..{args.max_n}; explicit p,q enumeration through "
              f"{min(args.max_n, args.direct_through)}; replication D=0..{args.max_degree}.")
        print("PASS: degree thresholds, exact births, dyadic/profile identities, uniform finite bounds.")
        print("Scope: finite discrete checks only; analytic and dynamical claims are not certified.")


if __name__ == "__main__":
    main()
