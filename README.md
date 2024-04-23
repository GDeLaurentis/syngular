# A Python Interface \& Extension to [Singular](https://www.singular.uni-kl.de/)

[![Continuous Integration Status](https://github.com/GDeLaurentis/syngular/actions/workflows/continuous_integration.yml/badge.svg)](https://github.com/GDeLaurentis/syngular/actions)
[![Coverage](https://img.shields.io/badge/Coverage-86%25-greenyellow?labelColor=2a2f35)](https://github.com/GDeLaurentis/syngular/actions)
[![PyPI Downloads](https://img.shields.io/pypi/dm/syngular.svg?label=PyPI%20downloads)](https://pypi.org/project/syngular/)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/GDeLaurentis/syngular/HEAD)

## Interface

- Ideals, Rings, Quotient Rings functions

## Extension

### Points on varieties over $\mathbb{F}_p$, $\mathbb{Q}_p$ and $\mathbb{C}$.

The function `ideal.point_on_variety` allows generate solutions to arbirary systems of polynomial equations in arbitrary polynomial quotient rings, over any of the three above mentioned fields. The system of equations may be underconstrained, i.e. the ideal may have any dimension. The $p$-adic and complex solutions can be requested as not exact. I.e. the point may lie close to but not exactly on the associated variety. This is essential for numerical computations where otherwise division-by-zero erros may occur. The limitation is that Singular must be able to compute an indepednent set for the ideal, in order to reduce it to a zero dimensional sub-variety.

### Primality test (lighter than a primary decomposition).

The function `ideal.test_primality` allows to test whether an ideal is prime or not, without performing a full primary decomposition. The algorithm can run also with successively looser degree bounds. It returns True if the idea is prime, False if it is not, or raises an `Inconclusive` exception if it cannot decide. The latter case happens, for instance, if the ideal is not radical, since the algorithm may not be able to find a linear projection. Inconclusive cases include when the ideal is primary but not prime.


## Requirements
```
numpy, sympy, mutableint, Singular
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

```
from syngular import Ideal, Ring
I = Ideal(Ring('0', ('x1', 'x2'), 'dp'), ['x1*x2'])
```
