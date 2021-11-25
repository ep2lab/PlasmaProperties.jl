using Documenter
using PlasmaProperties

push!(LOAD_PATH,"../src/")
makedocs(sitename="PlasmaProperties.jl Documentation",
         pages = [
            "Index" => "index.md", 
         ],
         format = Documenter.HTML(prettyurls = false)
) 
deploydocs(
    repo = "github.com/ep2lab/PlasmaProperties.jl.git",
    devbranch = "main"
)
