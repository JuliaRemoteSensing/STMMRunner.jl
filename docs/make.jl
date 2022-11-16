using STMMRunner
using Documenter

DocMeta.setdocmeta!(STMMRunner, :DocTestSetup, :(using STMMRunner); recursive = true)

makedocs(;
         modules = [STMMRunner],
         authors = "Gabriel Wu <wuzihua@pku.edu.cn> and contributors",
         repo = "https://github.com/JuliaRemoteSensing/STMMRunner.jl/blob/{commit}{path}#{line}",
         sitename = "STMMRunner.jl",
         format = Documenter.HTML(;
                                  prettyurls = get(ENV, "CI", "false") == "true",
                                  canonical = "https://JuliaRemoteSensing.github.io/STMMRunner.jl",
                                  assets = String[]),
         pages = [
             "Home" => "index.md",
             "Configurations" => "config.md",
             "MSTM v4" => "mstmv4.md",
             "MSTM v3" => "mstmv3.md",
         ])

deploydocs(; repo = "github.com/JuliaRemoteSensing/STMMRunner.jl", devbranch = "main")
