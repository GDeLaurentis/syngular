# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- Added warning `to point_on_variety` if `directions` are specified when the field is a finite field (and thus the directions cannot be imposed).

### Changed

### Fixed

- Fixed issue where TIMEOUT was not being correctly set after removing mutableint.
- Fixed issue in README where the example could not be copy-pasted into code (missing quotation marks).

### Deprecated

## [0.2.2] - 2024-06-06

### Added

- Custom `Monomial` and `Polynomial` classes, to better handle polynomials with coefficients in any of $\mathbb{C}, \mathbb{F}_p, \mathbb{Q}_p$.
- DOI.
- `Ideal.point_on_variety` now checks whether the provided `directions` are consistent with the given variety.
- Added support for rational field, `Field('rational', 0, 0)`.
- Tentatively introducing shorter aliases for field names, such as `C`, `Qp` and `Fp`.
- Continuous depolyment of Documentation via GitHub Pages.
- CI checks several python versions.

### Changed

- `Ideal.point_on_variety` no longer relies on `sympy` for parsing the polynomials, resulting in a drastic performance boost for complicated systems of polynomials. Interface should be unchanged.
- Splitting CI Test and CI Lint.
- <span style="color:red"> Breaking: </span> `syngular.TIMEOUT` and `syngular.DEGBOUND` are now normal integers (instead of mutableints) (`.set()` will not work any longer). Removed `mutableint` because it breaks with python3.12 and might not be needed in the first place.

### Fixed

- Improved logic for call to `Singular` with long commands requiring the writing of a file. As a consequence, parallelised calls should now have much better stability.

### Deprecated

- `Field.random_element` is now deprecated. Use `Field.random` instead.

## [0.2.1] - 2024-05-04

### Added

- `QRing` is an alias for `QuotientRing`.
- `test_primality` accepts optional `kwarg`, `iterated_degbound_computation`, with defult value `False`. Behaviour is unchanged on default setting. If set to `True`, the codimensions of ideals involving f-poly factors are computed with subsequently looser degree bounds, until either a given degree (18) is exceeded or the codim is greater than that of the original ideal. Use is recommended only for ideals that are expected to be prime and for which the computation without degree bound fails.

### Changed

- Default `max_tries` for `Ideal.point_on_variety` increased to 100 from 10: needed when generating a large number of points.

### Fixed

- Occasinally, when generating a `point_on_variety` with multi-precision complex (`mpc`) the `lex_groebner_basis` computation can return the unit ideal, due to precision loss. Introduced a new `RootPrecisionError` exception and added this to the retry logic.

## [0.2.0] - 2023-12-27

### Added

- Ideal method to generate points of varieties, `Ideal.point_on_variety`.
- Field object (moved from [lips](https://github.com/GDeLaurentis/lips))

### Fixed

- Syntax consistency when printing polynomials over rings with single-character variables (forced `short=0` always).
- Import path in Singular for `poly.lib` is now `polylib.lib`.

## [0.1.3] - 2023-03-02

### Added

- Basic interface functions for `Ring`, `QuotientRing` and `Ideal`.
- Primality test, `Ideal.test_primality`.

[unreleased]: https://github.com/GDeLaurentis/syngular/compare/v0.2.2...HEAD
[0.2.2]: https://github.com/GDeLaurentis/syngular/compare/v0.2.1...v0.2.2
[0.2.1]: https://github.com/GDeLaurentis/syngular/compare/v0.2.0...v0.2.1
[0.2.0]: https://github.com/GDeLaurentis/syngular/compare/v0.1.3...v0.2.0
[0.1.3]: https://github.com/GDeLaurentis/syngular/releases/tag/v0.1.3
