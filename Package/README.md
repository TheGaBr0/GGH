# GGH-crypto
GGH-crypto is a Python package implementing the Goldreich-Goldwasser-Halevi (GGH) public key cryptosystem and its optimization, GGH-HNF by Micciancio. This package is designed for educational and research purposes, offering insights into lattice-based cryptography.
This project was developed as part of a 3-year degree program at the Università degli Studi di Milano (University of Milan). It explores the resilience of lattice-based cryptography against quantum threats and introduces an hybrid variant.

# Features

- Implementation of the original GGH cryptosystem (1997)
- Implementation of the GGH-HNF optimization (2002)
- Utility functions for lattice-based cryptography
- Algorithms for solving the Closest Vector Problem (CVP)
- Lattice reduction algorithms

# Usage and details
For detailed installation, usage, examples and documentation, please visit the [GitHub repository](https://github.com/TheGaBr0/GGH).


# Note
Both the original GGH cryptosystem and its GGH-HNF optimization have known security vulnerabilities. This implementation is not intended for production use.

# Changelog

## 1.1.0
### Performance
- `Utils.babai_rounding` now solves the linear system `x · basis = point` instead of computing `point * basis.inv()`, removing a full exact-rational matrix inversion (and an extra matrix multiply) from every `decrypt()` call in both `GGHCryptosystem` and `GGHHNFCryptosystem`. On large dimensions this is the dominant cost of decryption — roughly an order-of-magnitude speedup at n ≈ 200, with bit-identical output.
- `GGHCryptosystem.decrypt` applies the same change to the final `CVP * public_basis.inv()` step.
- `GGHHNFCryptosystem.generate_keys_from_R` and `GGHCryptosystem.generate_keys_from_R` no longer compute an unused `R.inv()`, avoiding an O(n³) inversion when a private basis is supplied.

### Bug fixes
- Constructing `GGHCryptosystem` with a supplied `private_basis` was broken: `__init__` called the non-existent `generate_keys_from_basis()` (renamed to `generate_keys_from_R`). Fixed the call site.
- `GGHCryptosystem.generate_sigma` was inverting the basis twice on the keys-from-R path, so `sigma` was derived from `R` instead of `R⁻¹` and disagreed with normal key generation. It now receives the basis, like every other call site.
- `GGHHNFCryptosystem.generate_keys_from_R` computed `R_rho` only when `debug=True`, so building from a private basis with `debug=False` left it `None` and broke `generate_error()`. The computation is now unconditional.

### Breaking changes
- `GGHHNFCryptosystem.private_key` is now the private basis matrix (`fmpz_mat`), not a `(R_inv, R)` tuple, matching the other key-generation path and `GGHCryptosystem.private_key`. Replace `private_key[1]` with `private_key`.

## 1.0.5
- Added `nguyen_fix` parameter to `GGHCryptosystem` (default `False`). When enabled, implements the Mandangan et al. (2020) countermeasure against Nguyen's attack: error entries are drawn from {σ-2, σ-1, σ, σ+1} instead of {-σ, +σ}, preserving ||e|| = σ√n while breaking the elimination stage of the attack.
- Key generation with `nguyen_fix=True` automatically retries until a basis yielding σ > 2 is found (required by the countermeasure). Raises `ValueError` after 100 failed attempts with a suggestion to increase the dimension.

## 1.0.4
- Initial stable release.    