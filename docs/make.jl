using STMMRunner
using Documenter

DocMeta.setdocmeta!(STMMRunner, :DocTestSetup, :(using STMMRunner); recursive = true)

makedocs(;
    modules = [STMMRunner],
    authors = "Gabriel Wu <wuzihua@pku.edu.cn> and contributors",
    repo = "https://github.com/lucifer1004/STMMRunner.jl/blob/{commit}{path}#{line}",
    sitename = "STMMRunner.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://lucifer1004.github.io/STMMRunner.jl",
        assets = String[],
    ),
    pages = ["Home" => "index.md"],
)

deploydocs(; repo = "github.com/lucifer1004/STMMRunner.jl", devbranch = "main")
