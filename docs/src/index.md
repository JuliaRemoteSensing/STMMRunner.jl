# STMMRunner

Documentation for [STMMRunner](https://github.com/JuliaRemoteSensing/STMMRunner.jl).

Currently supported STMM models:

- [MSTM v4](https://github.com/dmckwski/MSTM)
  - You **do not need** to compile `MSTM v4` locally since there is already the `MSTM_jll` package.
  - Lattice mode not supported yet.
- [MSTM v3](https://www.eng.auburn.edu/~dmckwski/scatcodes/)
  - You **need** to compile a parallel version of `MSTM v3`.
- [FaSTMM](https://bitbucket.org/planetarysystemresearch/fastmm_v1.0)
  - You **need** a compiled version of `FaSTMM`.
- [SMUTHI](https://gitlab.com/AmosEgel/smuthi)
  - `SMUTHI` will be automatically installed via `CondaPkg`.

Working in progress:

- [CELES](https://github.com/disordered-photonics/celes.git)
- [TERMS](https://github.com/nano-optics/terms.git)
