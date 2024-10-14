# A Python Interface \& Extension to [Singular](https://www.singular.uni-kl.de/)

[![CI Lint](https://github.com/GDeLaurentis/syngular/actions/workflows/ci_lint.yml/badge.svg)](https://github.com/GDeLaurentis/syngular/actions/workflows/ci_lint.yml)
[![CI Test](https://github.com/GDeLaurentis/syngular/actions/workflows/ci_test.yml/badge.svg)](https://github.com/GDeLaurentis/syngular/actions/workflows/ci_test.yml)
[![Coverage](https://img.shields.io/badge/Coverage-86%25-greenyellow?labelColor=2a2f35)](https://github.com/GDeLaurentis/syngular/actions)
[![Docs](https://github.com/GDeLaurentis/syngular/actions/workflows/cd_docs.yml/badge.svg?label=Docs)](https://gdelaurentis.github.io/syngular/)
[![PyPI](https://img.shields.io/pypi/v/syngular.svg?label=PyPI)](https://pypi.org/project/syngular/)
[![PyPI Downloads](https://img.shields.io/pypi/dm/syngular.svg?label=PyPI%20downloads)](https://pypistats.org/packages/syngular)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/GDeLaurentis/syngular/HEAD)
[![DOI](https://zenodo.org/badge/378157462.svg)](https://zenodo.org/doi/10.5281/zenodo.11113680)
[![Python](https://img.shields.io/pypi/pyversions/syngular?label=Python)](https://pypi.org/project/syngular/)

The `syngular` library is a Python 3 package for algebraic geometry computations. It provides an intuitive and object-oriented interface to [Singular](https://www.singular.uni-kl.de/). Furthermore, it extends the numerical capabilities of Singular, providing a numerical solver for arbitrary systems of polynomial equations in tandem with [pyadic](https://github.com/GDeLaurentis/pyadic), and its applicaibility to physics computations, where generic algorithms may be insufficient.

## Interface

Python classes for 'Ideal', 'Ring' and 'QuotientRing'. Several related functions accessible as attributes or methods. Intuitive operations through magic methods, e.g. Ideal addition '+' and intersection '&'.

## Extension

### Points on varieties over $\mathbb{F}_p$, $\mathbb{Q}_p$ and $\mathbb{C}$.

The function `ideal.point_on_variety` allows to obtain numerical solutions to arbirary systems of polynomial equations in arbitrary polynomial quotient rings, over any of the three above mentioned fields. The system of equations may be underconstrained, i.e. the ideal may have any dimension. The $p$-adic and complex solutions can be requested as not exact. I.e. the point may lie close to but not exactly on the associated variety. This is essential for numerical computations where otherwise division-by-zero erros may occur. The limitation is that Singular must be able to compute an indepednent set for the ideal, in order to reduce it to a zero dimensional sub-variety.

### Primality test (lighter than a primary decomposition).

The function `ideal.test_primality` allows to test whether an ideal is prime or not, without performing a full primary decomposition. The algorithm can run also with successively looser degree bounds. It returns True if the idea is prime, False if it is not, or raises an `Inconclusive` exception if it cannot decide. The latter case happens, for instance, if the ideal is not radical, since the algorithm may not be able to find a linear projection. Inconclusive cases include when the ideal is primary but not prime.


## Requirements
```
numpy, sympy, Singular
```


## Installation
```
pip install -e path/to/repo
```

## Testing

```
pytest --cov syngular/ --cov-report html tests/ --verbose
```

## Quick Start

Define an ideal over a ring in two variables
```
from syngular import Ideal, Ring
I = Ideal(Ring('0', ('x1', 'x2'), 'dp'), ['x1*x2'])
```
You can now inspect `I` to see what methods and attributes are available.

## Solving arbitrary systems of polynomial equations

Generate a $p$-adic solution to a system of 2 polynomial equations in 3 variables, controlling the precision to which they are solved.
```
field = Field("padic", 2 ** 31 - 1, 8)
ring = Ring('0', ('x', 'y', 'z', ), 'dp')
I = Ideal(ring, ['x*y^2+y^3-z^2', 'x^3+y^3-z^2', ])
```

The variety associated to `I` has 3 branches. In other words, the system of equations has 3 types of solutions.
```
(Q1, P1), (Q2, P2), (Q3, P3) = I.primary_decomposition
```

Generate a solution on the first branch
```
numerical_point = Q1.point_on_variety(field=field, directions=I.generators, valuations=(1, 1, ), ) 
```
is a dictionary of numerical values for each variable in the ring.

These are small with valuations (1, 1)
```
list(map(lambda string: Polynomial(string, field).subs(numerical_point), Q1.generators))
```

while these are O(1) with valuations (0, 0)
```
list(map(lambda string: Polynomial(string, field).subs(numerical_point), Q2.generators))
```

See [arXiv:2207.10125](https://arxiv.org/pdf/2207.10125) Fig. 1 for a graphical depiction.

## Citation

If you found this library useful, please consider citing it and [Singular](https://www.singular.uni-kl.de/)


```bibtex
@inproceedings{DeLaurentis:2023qhd,
    author = "De Laurentis, Giuseppe",
    title = "{Lips: $p$-adic and singular phase space}",
    booktitle = "{21th International Workshop on Advanced Computing and Analysis Techniques in Physics Research}: {AI meets Reality}",
    eprint = "2305.14075",
    archivePrefix = "arXiv",
    primaryClass = "hep-th",
    reportNumber = "PSI-PR-23-14",
    month = "5",
    year = "2023"
}
```

```bibtex
@misc {DGPS,
 title = {{\sc Singular} {4-3-0} --- {A} computer algebra system for polynomial computations},
 author = {Decker, Wolfram and Greuel, Gert-Martin and Pfister, Gerhard and Sch\"onemann, Hans},
 year = {2022},
 howpublished = {\url{http://www.singular.uni-kl.de}},
}
```
