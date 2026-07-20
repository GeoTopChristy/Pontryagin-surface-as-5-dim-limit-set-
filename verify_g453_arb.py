"""Rigorous Arb certificate for the A-B-C-D-E-F-O_c realization.

This file requires python-flint 0.6 or later.  It performs a 256-bit
Krawczyk inclusion and then checks every non-edge using Arb ball arithmetic.
No network access is used while the certificate runs.
"""

from flint import arb, arb_mat, ctx

ctx.prec = 256
N = 16

FIXED = {
    "xC": "0.06620923121333717465",
    "c": "0.37871672186762506129",
    "xD": "-0.07722547088248804537",
    "d": "0.63500241411989255536",
    "xE": "-0.22245295031144918241",
    "e": "0.81317731860166608815",
    "xF": "-0.17391442856664743322",
    "f": "0.93028963143779908052",
}

# xa, xb,yb,zb, yC,zC,wC, yD,zD,wD, yE,zE,wE, yF,zF,wF
CENTERS = [
    "1.7989074399478672722612275836242556383275449836252138705954054322966074398005331",
    "0.5558929702514211719920480478975692506910578149970370615776507310375082433606461",
    "0.40387988390687641267543535018388981475887243710327903705249865349443092534868527",
    "0.94765814150418754521609983343350862488091396979942633297159815277248946649108994",
    "0.20377106091432812103428517583297702460618709740489630765945188658044221404926705",
    "0.92955010124102065291915784832944224636146058266090614970378585583014880446242229",
    "0.14722627253866626618623436755464636020770033929193657392964854660825226196473711",
    "0.150709639288625046328496260937497127463550593638301122387157139538965248935566",
    "0.75453785267906673624523483424184668461182395756145319278545519996212887577956293",
    "0.22099320594624878838852050640397138879044884065835423615009208921251092913202197",
    "0.15935754342818882821809735698822571225883617187383859819841351295493932113616403",
    "0.50540331391642384391544417836882614977514224392735680493651976110864945990076987",
    "0.27642821133104210787368076092603665468613417626462918913558414606612599270594354",
    "-0.024165230772630758431635624879506161863819109477136680428789232249511121153120805",
    "0.31022447960072457765107245724362730269159691769540886094086757470520514373823294",
    "0.1876951189690994718670860402918496806348824126559974344185272764725320228375304",
]


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


def rotate_once(v):
    alpha = 2 * arb.pi() / 5
    beta = 2 * arb.pi() / 20
    ca, sa = alpha.cos(), alpha.sin()
    cb, sb = beta.cos(), beta.sin()
    return [
        v[0] * ca + v[1] * sa,
        v[0] * (-sa) + v[1] * ca,
        v[2] * cb + v[3] * sb,
        v[2] * (-sb) + v[3] * cb,
        v[4],
        v[5],
    ]


def equations_and_jacobian(values):
    v = variables(values)
    xa, xb, yb, zb = v[0:4]
    yC, zC, wC, yD, zD, wD = v[4:10]
    yE, zE, wE, yF, zF, wF = v[10:16]
    z0, z1 = constant(0), constant(1)
    fixed = {key: constant(value) for key, value in FIXED.items()}

    A = [xa, z0, z0, z0, z0, z1]
    B = [xb, yb, zb, z0, z0, z1]
    C = [fixed["xC"], yC, zC, wC, fixed["c"], z1]
    D = [fixed["xD"], yD, zD, wD, fixed["d"], z1]
    E = [fixed["xE"], yE, zE, wE, fixed["e"], z1]
    F = [fixed["xF"], yF, zF, wF, fixed["f"], z1]

    eqs = [
        lorentz(A, rotate_once(A)),
        lorentz(A, B),
        lorentz(A, rotate_once(B)),
        lorentz(B, rotate_once(B)),
    ]
    for outer, inner in ((B, C), (C, D), (D, E), (E, F)):
        eqs.extend(
            [
                lorentz(outer, inner),
                lorentz(outer, rotate_once(inner)),
                lorentz(inner, rotate_once(inner)),
            ]
        )
    return (
        [eq.value for eq in eqs],
        [[eq.derivative[j] for j in range(N)] for eq in eqs],
    )


def column(entries):
    return arb_mat([[x] for x in entries])


def identity(n):
    return arb_mat([[arb(1 if i == j else 0) for j in range(n)] for i in range(n)])


# Krawczyk inclusion.
center = [arb(s) for s in CENTERS]
F0, J0 = equations_and_jacobian(center)
Cinv = arb_mat(J0).inv()
radius = "1e-40"
box = [arb(s, radius) for s in CENTERS]
_, Jbox = equations_and_jacobian(box)

center_vector = column(center)
box_vector = column(box)
K = (
    center_vector
    - Cinv * column(F0)
    + (identity(N) - Cinv * arb_mat(Jbox)) * (box_vector - center_vector)
)

if not all(box[i].contains_interior(K[i, 0]) for i in range(N)):
    raise RuntimeError("Krawczyk inclusion failed")
print("Krawczyk inclusion: PASS")
print("Unique zero certified in the radius-1e-40 box")


def pole_vectors(root_box):
    f = {key: arb(value) for key, value in FIXED.items()}
    xa, xb, yb, zb = root_box[0:4]
    yC, zC, wC, yD, zD, wD = root_box[4:10]
    yE, zE, wE, yF, zF, wF = root_box[10:16]
    return {
        "A": [xa, arb(0), arb(0), arb(0), arb(0), arb(1)],
        "B": [xb, yb, zb, arb(0), arb(0), arb(1)],
        "C": [f["xC"], yC, zC, wC, f["c"], arb(1)],
        "D": [f["xD"], yD, zD, wD, f["d"], arb(1)],
        "E": [f["xE"], yE, zE, wE, f["e"], arb(1)],
        "F": [f["xF"], yF, zF, wF, f["f"], arb(1)],
        "O": [arb(0), arb(0), arb(0), arb(0), arb(1) / f["f"], arb(1)],
    }


def rotate(v, k):
    alpha = 2 * arb.pi() * k / 5
    beta = 2 * arb.pi() * k / 20
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


def is_edge(left, right):
    L, i, _ = left
    M, j, _ = right
    if L > M:
        return is_edge(right, left)
    if L == M:
        period = 5 if L == "A" else 20
        return L != "O" and (i - j) % period in (1, period - 1)
    if (L, M) == ("A", "B"):
        return (j - i) % 5 in (0, 1)
    chain = "BCDEF"
    if L in chain and M in chain and chain.index(M) == chain.index(L) + 1:
        return (j - i) % 20 in (0, 1)
    return (L, M) == ("F", "O")


base = pole_vectors(box)
vertices = [("A", i, rotate(base["A"], i)) for i in range(5)]
vertices.extend(
    (label, i, rotate(base[label], i)) for label in "BCDEF" for i in range(20)
)
vertices.append(("O", 0, base["O"]))

norm_lower_bounds = {}
for label in "ABCDEFO":
    q = lorentz(base[label], base[label])
    if not (q > arb(0)):
        raise RuntimeError("A pole is not certified spacelike")
    norm_lower_bounds[label] = q.lower()

threshold = arb("1.0627")
minimum_midpoint = None
minimum_pair = None
checked_non_edges = 0
for i, left in enumerate(vertices):
    for right in vertices[i + 1 :]:
        if is_edge(left, right):
            continue
        inner = lorentz(left[2], right[2])
        qleft = lorentz(left[2], left[2])
        qright = lorentz(right[2], right[2])
        H = inner**2 / (qleft * qright)
        if not (H > threshold):
            raise RuntimeError(
                "Non-edge separation failed for %s and %s: %s"
                % (left[0:2], right[0:2], H)
            )
        checked_non_edges += 1
        midpoint = float(H.mid())
        if minimum_midpoint is None or midpoint < minimum_midpoint:
            minimum_midpoint = midpoint
            minimum_pair = (left[0:2], right[0:2], H)

print("Spacelike poles: PASS")
print("Norm lower bounds:", norm_lower_bounds)
print("Certified non-edges checked:", checked_non_edges)
print("All non-edges satisfy H > 1.0627: PASS")
print("Smallest interval encountered:", minimum_pair)
