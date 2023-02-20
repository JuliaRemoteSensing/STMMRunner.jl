# STMMRunner

[STMMRunner](https://github.com/JuliaRemoteSensing/STMMRunner.jl) provides a common configuration interface for several STMM (Superposition T-Matrix Method) codes. Since each runner package (they are listed below) reexports all the structs and methods defined in `STMMRunner`, you do not need to specify `STMMRunner` as a dependency. Instead, you should use the corresponding STMM code package directly.

Currently supported STMM models:

- [MSTM4Runner](@ref): wrapping [MSTM v4](https://github.com/dmckwski/MSTM)
  - You **do not need** to compile `MSTM v4` locally since there is already the `MSTM_jll` package.
  - Lattice mode not supported yet.
- [MSTM3Runner](@ref): wrapping [MSTM v3](https://www.eng.auburn.edu/~dmckwski/scatcodes/)
  - You **need** to compile a parallel version of `MSTM v3`.
- [FaSTMMRunner](@ref): wrapping [FaSTMM](https://bitbucket.org/planetarysystemresearch/fastmm_v1.0)
  - You **need** a compiled version of `FaSTMM`.
- [SMUTHIRunner](@ref): wrapping [SMUTHI](https://gitlab.com/AmosEgel/smuthi)
  - `SMUTHI` will be automatically installed via `CondaPkg`.

Working in progress:

- [CELES](https://github.com/disordered-photonics/celes.git)
- [TERMS](https://github.com/nano-optics/terms.git)
