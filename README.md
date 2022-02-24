# STMMRunner

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaRemoteSensing.github.io/STMMRunner.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaRemoteSensing.github.io/STMMRunner.jl/dev)
[![Build Status](https://github.com/JuliaRemoteSensing/STMMRunner.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/JuliaRemoteSensing/STMMRunner.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/JuliaRemoteSensing/STMMRunner.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaRemoteSensing/STMMRunner.jl)

General runner for multiple STMM models.

Currently supported STMM models:

- [MSTM v4](https://github.com/dmckwski/MSTM)
  - You **do not need** to compile MSTM v4 locally since there is already the `MSTM_jll` package.
  - Lattice mode not supported yet.
- [MSTM v3](https://www.eng.auburn.edu/~dmckwski/scatcodes/)
  - You **need** to compile a parallel version of MSTM v3.

Working in process:

- [FaSTMM](https://bitbucket.org/planetarysystemresearch/fastmm_v1.0)
- [CELES](https://github.com/disordered-photonics/celes.git)
- [TERMS](https://github.com/nano-optics/terms.git)
