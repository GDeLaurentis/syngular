# A Python Interface \& Extension to [Singular](https://www.singular.uni-kl.de/)

[![Continuous Integration Status](https://github.com/GDeLaurentis/syngular/actions/workflows/continuous_integration.yml/badge.svg)](https://github.com/GDeLaurentis/syngular/actions)
[![Coverage](https://img.shields.io/badge/Coverage-88%25-greenyellow?labelColor=2a2f35)](https://github.com/GDeLaurentis/syngular/actions)
[![PyPI Downloads](https://img.shields.io/pypi/dm/syngular.svg?label=PyPI%20downloads)](https://pypi.org/project/syngular/)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/GDeLaurentis/syngular/HEAD)


## Extension

- Points on varieties over $\mathbb{F}_p$, $\mathbb{Q}_p$ and $\mathbb{C}$.
- Primality test (lighter than a primary decomposition).

## Interface

- Ideals, Rings, Quotient Rings functions


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
