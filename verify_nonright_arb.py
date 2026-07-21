"""Rigorous Arb certificate for the non-right-angled numerical data.

The weighted nerve is (2,3,3,2,3), in the order
    (A_1 A_2, A_1 C_1, A_1 C_2, C_1 C_2, C_1 O_c).

For each requested p, this program first uses a 256-bit Krawczyk
inclusion to certify an exact solution of the five edge relations near
the printed decimal data (a, c, and w_c are held at their printed decimal
values).  The fifth coordinate o of O_c is an additional unknown.  The
program then checks every non-edge by Arb ball arithmetic and certifies a
nonzero 6-by-6 pole minor.

With no command-line arguments, the program certifies p = 2, 3, and 4.
No network access is used while the program runs.
"""

import sys

from flint import arb, arb_mat, ctx

ctx.prec = 256
N = 5
M = 5


DATA = {
    2: {
        "fixed": {
            "a": "-0.6268817510226523",
            "c": "0.5487408449899542",
            "w": "0",
        },
        # x_a, x_c, y_c, z_c, o
        "center": [
            "-1.401554716067853470030676095811538950854765709159313015550226158491128256496580418191794542811477826",
            "-0.6476210307741850815508051702909336278248988667881674448748027442863469636697556417635114604121346237",
            "-0.4705242208881140487392931205531704259453020382949411402482043595330349440535111453765058127394534191",
            "0.7868299486133950741177542241560532017915468140143337204055661269737953831550700042113036088650644047",
            "1.2785906461197737040079695439130898140424808182648477310558262814212761851615962968768737605591127742890494422",
        ],
    },
    3: {
        "fixed": {
            "a": "0.48414887858672406",
            "c": "-0.60870282888116",
            "w": "0",
        },
        "center": [
            "-1.574018170609196835110080873645721885249266041239479257143383578723507525725412862587621725334034779",
            "-0.5654007505721593803086442232470231588189012716737036274208823665183934497391908399726154421972104636",
            "-0.4107876906568251714973844016419541370136743205056705790445308955699881099373183992452078082159215927",
            "0.7237662032437768771294339977016883706182465192791178069382628277154701067131275011378254528212999306",
            "-1.2563383350869236388545380472908446520318897567578880092089734358459042418554865260380081917415479992854580558",
        ],
    },
    4: {
        "fixed": {
            "a": "0.574393980597",
            "c": "-0.354053027688",
            "w": "0.320143790304",
        },
        "center": [
            "1.472547826708362146810734554709705594465915647848468832252742223617789948590729264614229394949123023",
            "0.5635165586028417620214297481440625490719200761461305606159987135528915001935757617552391120282962491",
            "0.4094187450601897497338999002691085524133728703422732619769196064362360340986375850203161847951369032",
            "0.8121101333467932146400457844938095007994974244557607382986864441659548757663847353385597665871821876",
            "-1.670750761584002775292220840606283507754307212678575288338452578989812",
        ],
    },
}

SEPARATION_THRESHOLDS = {
    2: arb("1.45"),
    3: arb("1.10"),
    4: arb("1.09"),
}


class AD:
    def __init__(self, value, derivative):
        self.value = value
        self.derivative = derivative

    def __add__(self, other):
        if not isinstance(other, AD):
            other = AD(arb(other), [arb(0)] * N)
        return AD(
            self.value + other.value,
            [self.derivative[j] + other.derivative[j] for j in range(N)],
        )

    __radd__ = __add__

    def __neg__(self):
        return AD(-self.value, [-x for x in self.derivative])

    def __sub__(self, other):
        return self + (-other)

    def __rsub__(self, other):
        return other + (-self)

    def __mul__(self, other):
        if not isinstance(other, AD):
            other = AD(arb(other), [arb(0)] * N)
        return AD(
            self.value * other.value,
            [
                self.derivative[j] * other.value
                + self.value * other.derivative[j]
                for j in range(N)
            ],
        )

    __rmul__ = __mul__


def constant(value):
    return AD(arb(value), [arb(0)] * N)


def variables(values):
    result = []
    for i, value in enumerate(values):
        derivative = [arb(0)] * N
        derivative[i] = arb(1)
        result.append(AD(value, derivative))
    return result


def lorentz(u, v):
    return sum(u[j] * v[j] for j in range(5)) - u[5] * v[5]


def rotate(v, k, p):
    alpha = 2 * arb.pi() * k / M
    beta = 2 * arb.pi() * k / (M * p)
    ca, sa = alpha.cos(), alpha.sin()
    cb, sb = beta.cos(), beta.sin()
    return [
        ca * v[0] + sa * v[1],
        -sa * v[0] + ca * v[1],
        cb * v[2] + sb * v[3],
        -sb * v[2] + cb * v[3],
        v[4],
        v[5],
    ]


def equations_and_jacobian(p, fixed, values):
    xa, xc, yc, zc, o = variables(values)
    zero, one = constant(0), constant(1)
    a, c, wc = (constant(fixed[key]) for key in ("a", "c", "w"))
    A = [xa, zero, zero, zero, a, one]
    C = [xc, yc, zc, wc, c, one]
    O = [zero, zero, zero, zero, o, one]
    SA = rotate(A, 1, p)
    SC = rotate(C, 1, p)
    qa = lorentz(A, A)
    qc = lorentz(C, C)
    ac1 = lorentz(A, C)
    ac2 = lorentz(A, SC)
    eqs = [
        lorentz(A, SA),              # H(A_1,A_2)=0
        4 * ac1 * ac1 - qa * qc,     # H(A_1,C_1)=1/4
        4 * ac2 * ac2 - qa * qc,     # H(A_1,C_2)=1/4
        lorentz(C, SC),              # H(C_1,C_2)=0
        4 * lorentz(C, O) * lorentz(C, O)
        - qc * lorentz(O, O),        # H(C_1,O_c)=1/4
    ]
    return (
        [eq.value for eq in eqs],
        [[eq.derivative[j] for j in range(N)] for eq in eqs],
    )


def column(entries):
    return arb_mat([[x] for x in entries])


def identity(n):
    return arb_mat(
        [[arb(1 if i == j else 0) for j in range(n)] for i in range(n)]
    )


def is_edge(left, right, p):
    left_label, i, _ = left
    right_label, j, _ = right
    if left_label > right_label:
        return is_edge(right, left, p)
    if left_label == right_label == "A":
        return (i - j) % M in (1, M - 1)
    if left_label == right_label == "C":
        return (i - j) % (M * p) in (1, M * p - 1)
    if (left_label, right_label) == ("A", "C"):
        return (j - i) % M in (0, 1)
    return (left_label, right_label) == ("C", "O")


def certify(p, datum):
    fixed = datum["fixed"]
    center = [arb(value) for value in datum["center"]]
    f0, j0 = equations_and_jacobian(p, fixed, center)
    inverse = arb_mat(j0).inv()
    radius = "1e-40"
    box = [arb(value, radius) for value in datum["center"]]
    _, jbox = equations_and_jacobian(p, fixed, box)
    krawczyk = (
        column(center)
        - inverse * column(f0)
        + (identity(N) - inverse * arb_mat(jbox)) * (column(box) - column(center))
    )
    if not all(box[i].contains_interior(krawczyk[i, 0]) for i in range(N)):
        raise RuntimeError("p=%d: Krawczyk inclusion failed" % p)

    xa, xc, yc, zc, o = box
    a, c, wc = (arb(fixed[key]) for key in ("a", "c", "w"))
    A = [xa, arb(0), arb(0), arb(0), a, arb(1)]
    C = [xc, yc, zc, wc, c, arb(1)]
    O = [arb(0), arb(0), arb(0), arb(0), o, arb(1)]
    qa, qc, qo = (lorentz(v, v) for v in (A, C, O))
    if not (qa > arb(0) and qc > arb(0) and qo > arb(0)):
        raise RuntimeError("p=%d: a pole is not certified spacelike" % p)
    if not (
        lorentz(A, C) < arb(0)
        and lorentz(A, rotate(C, 1, p)) < arb(0)
    ):
        raise RuntimeError("p=%d: wrong branch for an A-C angle" % p)
    if not (lorentz(C, O) < arb(0)):
        raise RuntimeError("p=%d: wrong branch for the C_1 O_c angle" % p)

    vertices = [("A", i, rotate(A, i, p)) for i in range(M)]
    vertices.extend(("C", i, rotate(C, i, p)) for i in range(M * p))
    vertices.append(("O", 0, O))

    minor_columns = [A, rotate(A, 1, p), C, rotate(C, 1, p), rotate(C, 2, p), O]
    pole_minor = arb_mat(
        [[vector[row] for vector in minor_columns] for row in range(6)]
    ).det()
    if pole_minor.contains(arb(0)):
        raise RuntimeError("p=%d: the selected pole minor contains zero" % p)

    failures = []
    minimum = None
    family_minima = {}
    checked = 0
    for i, left in enumerate(vertices):
        for right in vertices[i + 1 :]:
            if is_edge(left, right, p):
                continue
            inner = lorentz(left[2], right[2])
            if not (inner < arb(0)):
                raise RuntimeError(
                    "p=%d: non-edge sign failed for %s and %s: %s"
                    % (p, left[0:2], right[0:2], inner)
                )
            h_value = inner**2 / (
                lorentz(left[2], left[2]) * lorentz(right[2], right[2])
            )
            checked += 1
            if minimum is None or h_value.mid() < minimum[2].mid():
                minimum = (left[0:2], right[0:2], h_value)
            family = left[0] + right[0]
            if (
                family not in family_minima
                or h_value.mid() < family_minima[family][2].mid()
            ):
                family_minima[family] = (left[0:2], right[0:2], h_value)
            if not (h_value > arb(1)):
                failures.append((left[0:2], right[0:2], h_value))

    a_o = lorentz(A, O) ** 2 / (qa * qo)
    c_o = lorentz(C, O) ** 2 / (qc * qo)
    print("m=%d, p=%d" % (M, p))
    print("  Krawczyk inclusion: PASS (unique zero in radius-1e-40 box)")
    print("  spacelike poles and intended angle branches: PASS")
    print("  nonzero pole minor det(A_1,A_2,C_1,C_2,C_3,O_c):", pole_minor)
    print("  non-edges checked:", checked)
    print("  all non-edge inner products are strictly negative: PASS")
    print("  H(C_1,O_c) edge:", c_o)
    print("  H(A_1,O_c):", a_o)
    print("  smallest non-edge interval:", minimum)
    print("  minima by unordered orbit family:")
    for family in sorted(family_minima):
        print("   ", family, family_minima[family])
    print("  non-edges not certified above 1:", len(failures))
    if p in SEPARATION_THRESHOLDS:
        threshold = SEPARATION_THRESHOLDS[p]
        if not (minimum[2] > threshold):
            raise RuntimeError(
                "p=%d: claimed separation threshold %s was not certified"
                % (p, threshold)
            )
        print("  all non-edges satisfy H > %s: PASS" % threshold)
    return failures


requested_p = tuple(int(value) for value in sys.argv[1:]) or (2, 3, 4)
if any(value not in DATA for value in requested_p):
    raise SystemExit("supported values of p are 2, 3, and 4")

all_failures = []
for p_value in requested_p:
    all_failures.extend((p_value, item) for item in certify(p_value, DATA[p_value]))

if all_failures:
    raise SystemExit("AUDIT FAILED: a non-edge inequality was not certified")

print("All non-edges satisfy H > 1: PASS")
