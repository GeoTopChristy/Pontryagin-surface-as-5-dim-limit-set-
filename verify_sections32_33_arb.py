"""Arb certificate for the finite candidates in Subsections 3.2 and 3.3.

The program treats

    p = 2,  m = 5,...,14,
    p = 3,  m = 6.

For each pair it certifies the intended root of the three remaining edge
equations, spacelikeness, strict negativity of every non-edge inner product,
the strict hyperparallelism inequality H > 1 for every non-edge, and a
nonzero 6-by-6 pole minor.  Consequently, the normalized Gram matrix has
non-positive off-diagonal entries and inertia

    (5, 1, (2*p + 1)*m - 5).

The computation uses 256-bit Arb balls through python-flint.  No network
access is used.
"""

from flint import arb, arb_mat, ctx


ctx.prec = 256
NVAR = 3


P2_DATA = {
    5: {
        "c": "0.799017656465542408620140381",
        "xc": "0.122273943504802098569920749",
        "center": (
            "0.256610919062989110001238515408447848433888719840894108224521177953209494316849889281225659",
            "-0.592493888387420168789735254386298458271675840304697732081178293523980543380217423818388142",
            "-0.254978775728845090738675627169753935318549313332070831415943065158322662062122698190804237",
        ),
    },
    6: {
        "c": "0.752207433154698587777682664",
        "xc": "0.295340841312926983316869000",
        "center": (
            "0.363985065970052014639277212315960396880382841158660137398440428827103709501842598635024746",
            "-0.526566055023335379943943003796005418162108163024400546588056277019403795723810374162039722",
            "-0.311817671928598973334391351303762001781513015972230633135622972606925666029392932688209969",
        ),
    },
    7: {
        "c": "0.701421169267411435561004691",
        "xc": "0.420692468577898393449324242",
        "center": (
            "0.416080364969525631082726734247665726881919482306499406776460691019000531925432753262309190",
            "-0.476935715148980222843728675820368087651149852275912561012693956276461176340652003121371325",
            "-0.306755363087187764457644280630695305526011864683559038054286747474715057632926777551681078",
        ),
    },
    8: {
        "c": "0.651811674204430000938609478",
        "xc": "0.520053035253210257751731310",
        "center": (
            "0.435291969867236946927939848121252647830081558747659740666998602913478220191883144896099151",
            "-0.433745563011742360357425874943701371999029513606080225018205097479667423398964080957805633",
            "-0.287010563561798073209905676920059012924352435328896480332440089318389790320075862425392439",
        ),
    },
    9: {
        "c": "0.605581569857028708543867298",
        "xc": "0.605581569857028708543867298",
        "center": (
            "0.439851669558000814855856579250312065372511937969897365139200916955766590073876477961467186",
            "-0.381430666678514989194727364314692647503096184083047567594870654075279862362499739795142339",
            "-0.267853511289241861029648800792160351402667671286296905265215795060758821231092430591125744",
        ),
    },
    10: {
        "c": "0.56347508022852538271539092",
        "xc": "0.66350870473171691330300603",
        "center": (
            "0.428846317454309366566990459305802788719551285901906496273325342452905062486597921970404503",
            "-0.357760149177176824297864486041081584124651188735441093346817518804139467271226912330448928",
            "-0.242260270880108919470808099932196096214856404614998338115956704245009945959035750993237202",
        ),
    },
    11: {
        "c": "0.52553978753537131089536249",
        "xc": "0.71507437010670445314418016",
        "center": (
            "0.416051264961175805451308109344748686769901573048322245292840640378189408716054022568934506",
            "-0.324318887928020485878207703136567525682169085609778232488381585430244252747337687065747940",
            "-0.221574760705922492101762758475563559092011599332455454172668268053265771723338601173713557",
        ),
    },
    12: {
        "c": "0.49151494147922154573416078",
        "xc": "0.75695338565955465860498795",
        "center": (
            "0.400983390010308034769607539721510201452691522615958522140401098731744595757810268711791338",
            "-0.293770653273517316888207580807645788840372793235373277111743714570749533207156941865396891",
            "-0.202457619347883871307250914755255604948880729321633090657723681965697716929897540430355508",
        ),
    },
    13: {
        "c": "0.46102461398697999394927073",
        "xc": "0.79121807855499177484077895",
        "center": (
            "0.385106737122972988774030850612146586177716817752413746891642121369058059555491230503996577",
            "-0.265944079524448015092487469135946680967687643467897883162605405165949472457540808791464885",
            "-0.184767996326696318969744743335361859088443906491765227031685361635033213571454001473286084",
        ),
    },
    14: {
        "c": "0.43367051369784541920121363",
        "xc": "0.81947258888977858851569865",
        "center": (
            "0.369248574155755040009093066151462235326856740049873488873042141688243478572236349076094125",
            "-0.240622385245369589735133970006646227994441547461654736699611875029370385305347351400412097",
            "-0.168278974156987446545499956821986144395305747584869123913156038906773703752282020972004085",
        ),
    },
}


P3_DATA = {
    6: {
        "a": "-0.4",
        "b": "0.05",
        "c": "0.6424527258674142103774534942",
        "xc": "0.3489363914190385393938022941",
        "center": (
            "0.500395202593635920076550913459585013110777966039778635774555432458625790443071827886437324",
            "-0.590703366810916717927348010052355323793567717534518819972381922891619040838686437579514422",
            "-0.279275471132337364875414758800705922594837655653259184894963195924916226964621955732622253",
        ),
    }
}


class AD:
    def __init__(self, value, derivative):
        self.value = value
        self.derivative = derivative

    def __add__(self, other):
        if not isinstance(other, AD):
            other = AD(arb(other), [arb(0)] * NVAR)
        return AD(
            self.value + other.value,
            [
                self.derivative[j] + other.derivative[j]
                for j in range(NVAR)
            ],
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
            other = AD(arb(other), [arb(0)] * NVAR)
        return AD(
            self.value * other.value,
            [
                self.derivative[j] * other.value
                + self.value * other.derivative[j]
                for j in range(NVAR)
            ],
        )

    __rmul__ = __mul__


def variables(values):
    result = []
    for i, value in enumerate(values):
        derivative = [arb(0)] * NVAR
        derivative[i] = arb(1)
        result.append(AD(value, derivative))
    return result


def lorentz(u, v):
    return sum(u[j] * v[j] for j in range(5)) - u[5] * v[5]


def rotate(v, k, p, m):
    alpha = 2 * arb.pi() * k / m
    beta = 2 * arb.pi() * k / (p * m)
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


def base_poles(p, m, a, b, c, xc, root):
    cosine = (2 * arb.pi() / m).cos()
    xa = ((1 - a * a) / cosine).sqrt()
    xb = (1 - a * b) * (cosine / (1 - a * a)).sqrt()
    yb = (
        (1 - a * b)
        * (arb.pi() / m).tan()
        * (cosine / (1 - a * a)).sqrt()
    )
    numerator = (
        2 * (a * b - 1) ** 2 * cosine**2
        - (1 - a * a) * (1 - b * b) * (1 + cosine)
    )
    denominator = (
        (2 * arb.pi() / (p * m)).cos()
        * (cosine + 1)
        * (a * a - 1)
    )
    zb = -(numerator / denominator).sqrt()
    y, z, w = root
    return {
        "A": [xa, arb(0), arb(0), arb(0), a, arb(1)],
        "B": [xb, yb, zb, arb(0), b, arb(1)],
        "C": [xc, y, z, w, c, arb(1)],
        "O": [arb(0), arb(0), arb(0), arb(0), 1 / c, arb(1)],
    }


def equations_and_jacobian(p, m, a, b, c, xc, values):
    poles = base_poles(p, m, a, b, c, xc, variables(values))
    B, C = poles["B"], poles["C"]
    equations = [
        lorentz(B, C),
        lorentz(B, rotate(C, 1, p, m)),
        lorentz(C, rotate(C, 1, p, m)),
    ]
    return (
        [equation.value for equation in equations],
        [
            [equation.derivative[j] for j in range(NVAR)]
            for equation in equations
        ],
    )


def column(entries):
    return arb_mat([[x] for x in entries])


def identity(n):
    return arb_mat(
        [[arb(1 if i == j else 0) for j in range(n)] for i in range(n)]
    )


def is_edge(left, right, p, m):
    left_label, i, _ = left
    right_label, j, _ = right
    order = {"A": 0, "B": 1, "C": 2, "O": 3}
    if order[left_label] > order[right_label]:
        return is_edge(right, left, p, m)
    if left_label == right_label == "A":
        return (i - j) % m in (1, m - 1)
    if left_label == right_label and left_label in ("B", "C"):
        return (i - j) % (p * m) in (1, p * m - 1)
    if (left_label, right_label) == ("A", "B"):
        return (j - i) % m in (0, 1)
    if (left_label, right_label) == ("B", "C"):
        return (j - i) % (p * m) in (0, 1)
    return (left_label, right_label) == ("C", "O")


def certify(p, m, datum):
    if p == 2:
        a = -(arb.pi() / (2 * m)).sin()
        b = -a
    else:
        a, b = arb(datum["a"]), arb(datum["b"])
    c, xc = arb(datum["c"]), arb(datum["xc"])
    center = [arb(value) for value in datum["center"]]
    f0, j0 = equations_and_jacobian(p, m, a, b, c, xc, center)
    inverse = arb_mat(j0).inv()
    radius = "1e-30"
    box = [arb(value, radius) for value in datum["center"]]
    _, jbox = equations_and_jacobian(p, m, a, b, c, xc, box)
    krawczyk = (
        column(center)
        - inverse * column(f0)
        + (identity(NVAR) - inverse * arb_mat(jbox))
        * (column(box) - column(center))
    )
    if not all(
        box[i].contains_interior(krawczyk[i, 0]) for i in range(NVAR)
    ):
        raise RuntimeError("p=%d, m=%d: Krawczyk inclusion failed" % (p, m))

    poles = base_poles(p, m, a, b, c, xc, box)
    vertices = [
        ("A", i, rotate(poles["A"], i, p, m)) for i in range(m)
    ]
    vertices.extend(
        ("B", i, rotate(poles["B"], i, p, m)) for i in range(p * m)
    )
    vertices.extend(
        ("C", i, rotate(poles["C"], i, p, m)) for i in range(p * m)
    )
    vertices.append(("O", 0, poles["O"]))

    for label in "ABCO":
        if not (lorentz(poles[label], poles[label]) > arb(0)):
            raise RuntimeError(
                "p=%d, m=%d: a pole is not certified spacelike" % (p, m)
            )

    checked_non_edges = 0
    closest_to_zero = None
    minimum_h = None
    for i, left in enumerate(vertices):
        for right in vertices[i + 1 :]:
            if is_edge(left, right, p, m):
                continue
            inner = lorentz(left[2], right[2])
            if not (inner < arb(0)):
                raise RuntimeError(
                    "p=%d, m=%d: non-edge sign failed for %s and %s: %s"
                    % (p, m, left[0:2], right[0:2], inner)
                )
            left_norm = lorentz(left[2], left[2])
            right_norm = lorentz(right[2], right[2])
            h_value = inner * inner / (left_norm * right_norm)
            if not (h_value > arb(1)):
                raise RuntimeError(
                    "p=%d, m=%d: H > 1 failed for %s and %s: %s"
                    % (p, m, left[0:2], right[0:2], h_value)
                )
            checked_non_edges += 1
            if (
                closest_to_zero is None
                or inner.mid() > closest_to_zero[2].mid()
            ):
                closest_to_zero = (left[0:2], right[0:2], inner)
            if minimum_h is None or h_value.mid() < minimum_h[2].mid():
                minimum_h = (left[0:2], right[0:2], h_value)

    minor_columns = [
        poles["A"],
        rotate(poles["A"], 1, p, m),
        poles["B"],
        rotate(poles["B"], 1, p, m),
        poles["C"],
        poles["O"],
    ]
    pole_minor = arb_mat(
        [[vector[row] for vector in minor_columns] for row in range(6)]
    ).det()
    if pole_minor.contains(arb(0)):
        raise RuntimeError(
            "p=%d, m=%d: selected pole minor contains zero" % (p, m)
        )

    vector_count = (2 * p + 1) * m + 1
    nullity = vector_count - 6
    print("p=%d, m=%d, N=%d" % (p, m, vector_count))
    print("  Krawczyk inclusion: PASS")
    print("  spacelike poles: PASS")
    print("  non-edges checked:", checked_non_edges)
    print("  all non-edge inner products are strictly negative: PASS")
    print("  closest non-edge inner-product interval:", closest_to_zero)
    print("  all non-edges satisfy H > 1: PASS")
    print("  smallest non-edge H interval:", minimum_h)
    print(
        "  nonzero pole minor det(A_1,A_2,B_1,B_2,C_1,O_c):",
        pole_minor,
    )
    print("  normalized Gram inertia: (5,1,%d): PASS" % nullity)


for m_value in sorted(P2_DATA):
    certify(2, m_value, P2_DATA[m_value])
for m_value in sorted(P3_DATA):
    certify(3, m_value, P3_DATA[m_value])

print("All Subsection 3.2 and 3.3 candidates: PASS")
