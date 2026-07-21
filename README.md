# Ancillary Analytic Note and Verification Certificates

This directory contains the ancillary material accompanying the manuscript:

1. a self-contained analytic note proving the uniform separation estimate
   for the right-angled (p=2,3) construction; and
2. three machine-checkable certificate programs for the finite and
   computer-assisted constructions.

The computations use rigorous Arb ball arithmetic through
<code>python-flint</code>; they do not rely on ordinary floating-point
rounding. No manuscript source or other unpublished material is needed to
run them.

For the finite candidates, the newly added (H>1) branch was also
independently rerun with Sage <code>RealBallField(256)</code>, which is backed
by Arb; the resulting enclosures are included in the supplied log.

Each Python program has two stages:

1. a Krawczyk inclusion certifies the existence and local uniqueness of an
   exact solution near the displayed numerical approximation; and
2. interval arithmetic verifies spacelikeness, the prescribed Coxeter edge
   relations, every strict non-edge separation inequality, the negative
   sign of every non-edge inner product, and a nonzero \(6\times6\) pole
   minor used to certify the Gram-matrix signature.

The programs perform no network requests and require no input files other
than the files in this directory.

## Files

| File | Description |
|---|---|
| [uniform_separation_analytic_supplement.tex](uniform_separation_analytic_supplement.tex) | Source of the analytic verification of the uniform separation estimates. |
| [uniform_separation_analytic_supplement.pdf](uniform_separation_analytic_supplement.pdf) | Compiled analytic verification cited in the manuscript. |
| [requirements_certificate.txt](requirements_certificate.txt) | Pinned Python dependency used to reproduce the certificates. |
| [verify_g453_arb.py](verify_g453_arb.py) | Certificate for the $A$-$B$-$C$-$D$-$E$-$F$-$O_c$ realization of $G_{4,5,3}$. |
| [verify_g453_arb.log](verify_g453_arb.log) | Recorded successful output of the $G_{4,5,3}$ certificate. |
| [verify_nonright_arb.py](verify_nonright_arb.py) | Certificate for the non-right-angled realizations with $m=5$, $p=2,3,4$, and weights $(2,3,3,2,3)$. |
| [verify_nonright_arb.log](verify_nonright_arb.log) | Recorded successful output of the non-right-angled certificate. |
| [verify_sections32_33_arb.py](verify_sections32_33_arb.py) | Certificate for the finite right-angled candidates in Subsections 3.2 and 3.3. |
| [verify_sections32_33_arb.log](verify_sections32_33_arb.log) | Recorded certificate output, augmented by an independent 256-bit Sage/Arb verification of every newly added (H>1) line. |

The supplied <code>.log</code> files are included for reference only. The
Python programs recompute every enclosure and inequality independently.

## Requirements

- Python 3;
- <code>python-flint==0.6.0</code>, as pinned in
  [requirements_certificate.txt](requirements_certificate.txt).

For a clean installation, create a virtual environment and install the
pinned dependency:

~~~bash
python -m venv .venv
~~~

On Linux or macOS:

~~~bash
source .venv/bin/activate
python -m pip install -r requirements_certificate.txt
~~~

On Windows PowerShell:

~~~powershell
.\.venv\Scripts\Activate.ps1
python -m pip install -r requirements_certificate.txt
~~~

## Running the certificates

Run the $G_{4,5,3}$ certificate with

~~~bash
python verify_g453_arb.py
~~~

Run all three non-right-angled cases with

~~~bash
python verify_nonright_arb.py
~~~

Individual non-right-angled cases may also be checked separately:

~~~bash
python verify_nonright_arb.py 2
python verify_nonright_arb.py 3
python verify_nonright_arb.py 4
~~~

Run the finite candidates from Subsections 3.2 and 3.3 with

~~~bash
python verify_sections32_33_arb.py
~~~

To create the corresponding release log on Windows PowerShell, run

~~~powershell
python .\verify_sections32_33_arb.py 2>&1 |
    Tee-Object .\verify_sections32_33_arb.log
~~~

A successful run exits with status code <code>0</code> and ends with
<code>PASS</code>. If a Krawczyk inclusion, spacelikeness condition, edge
relation, sign condition, or non-edge bound cannot be certified, the program
raises an error or exits with a nonzero status.

## Certified results

### The finite candidates in Subsections 3.2 and 3.3

The program [verify_sections32_33_arb.py](verify_sections32_33_arb.py)
certifies the candidates \((p,m)=(2,5),\ldots,(2,14)\) and \((3,6)\).
For the complete pole collection

$$
\{A_i\}_{i=1}^{m}\cup
\{B_j\}_{j=1}^{pm}\cup
\{C_q\}_{q=1}^{pm}\cup\{O_c\},
$$

the number of vectors is \(N=(2p+1)m+1\).  In every case, the certificate
evaluates every non-edge on the full root enclosure and verifies both

$$
\langle n(X),n(Y)\rangle<0
\qquad\text{and}\qquad
H(n(X),n(Y))>1.
$$

It also verifies that all normalized off-diagonal Gram entries are
non-positive and that the Gram inertia is

$$
(5,1,N-6)=(5,1,(2p+1)m-5).
$$

### The $G_{4,5,3}$ realization

The program [verify_g453_arb.py](verify_g453_arb.py) uses 256-bit Arb balls
and a 16-variable Krawczyk operator. It certifies:

- a unique exact zero in the radius-$10^{-40}$ box around the stored
  numerical center;
- spacelikeness of all pole vectors;
- all 5,240 non-edge inequalities and their negative inner-product signs;
- the nonvanishing of the \(A_1,B_1,C_1,D_1,E_1,F_1\) pole minor; and
- the uniform bound
  $$
  H>1.0627.
  $$

The smallest interval encountered is associated with the pair
$(B_0,D_1)$:

$$
H(B_0,D_1)
=
[1.0627365391573188754714869416872504002
\pm4.21\times10^{-38}].
$$

### The non-right-angled realizations

The program [verify_nonright_arb.py](verify_nonright_arb.py) treats $m=5$
and the edge weights

$$
(A_1A_2,A_1C_1,A_1C_2,C_1C_2,C_1O_c)
=(2,3,3,2,3).
$$

For each $p\in\{2,3,4\}$, it uses 256-bit Arb balls and a five-variable
Krawczyk operator. The certified separation results are:

| $p$ | Non-edges checked | Certified bound | Smallest interval |
|---:|---:|---:|---|
| 2 | 75 | $H>1.45$ | $H(C_1,C_6)\in[1.4569476819993073,1.4569476819993074]$ |
| 3 | 145 | $H>1.10$ | $H(C_1,C_6)\in[1.1082946366272705,1.1082946366272706]$ |
| 4 | 240 | $H>1.09$ | $H(C_1,C_6)\in[1.0929469998988836,1.0929469998988838]$ |

The program also certifies that all relevant poles are spacelike, that the
intended signed branches of the weight-three relations are used, that every
non-edge inner product is negative, that the
\(A_1,A_2,C_1,C_2,C_3,O_c\) pole minor is nonzero, and that

$$
H(C_1,O_c)=\cos^2(\pi/3)=\frac14.
$$

## Certification method

Let $F:\mathbb R^n\to\mathbb R^n$ be the system of edge equations and
$x_0$ the stored high-precision center. For
<code>verify_g453_arb.py</code> and <code>verify_nonright_arb.py</code>, the
root box is

$$
X=x_0+[-10^{-40},10^{-40}]^n.
$$

For <code>verify_sections32_33_arb.py</code>, it is

$$
X=x_0+[-10^{-30},10^{-30}]^3.
$$

Using an approximate inverse $C$ of $DF(x_0)$, each program evaluates
the Krawczyk enclosure

$$
K(X)=x_0-CF(x_0)+(I-CDF(X))(X-x_0)
$$

with outward-rounded Arb balls and verifies

$$
K(X)\subset\mathrm{int}(X).
$$

This inclusion certifies an exact zero in $X$, locally unique in the
certified box. All subsequent norm, edge, sign, and non-edge calculations
are evaluated on the entire root enclosure $X$, rather than only at its
midpoint.

For spacelike poles $u,v$, the separation quantity used by the programs is

$$
H(u,v)=
\frac{\langle u,v\rangle^2}
{\langle u,u\rangle\langle v,v\rangle}.
$$

Thus $H=\cos^2(\pi/k)$ certifies an edge of weight $k$, while $H>1$
certifies that a non-edge pair is hyperparallel.

## Reproducibility notes

- Decimal constants are parsed directly as Arb values; they are not first
  converted to binary floating-point numbers.
- The precision is set explicitly to 256 bits in all three programs.
- The programs enumerate all non-edges determined by the corresponding
  combinatorics.
- The pinned dependency records the software version used for the supplied
  logs.
- No manuscript text or external data are transmitted when any certificate
  is run.

## Citation

The manuscript cites this repository as a single ancillary record rather
than citing each Python file separately:

> <https://github.com/GeoTopChristy/Pontryagin-surface-as-5-dim-limit-set->.

For a submitted or published version, replace the moving repository URL by
a tagged release URL or a Zenodo DOI whenever possible. The Krawczyk method,
Arb, and <code>python-flint</code> are cited separately in the manuscript's
bibliography; the individual <code>.py</code> files are not separate
bibliography entries.
