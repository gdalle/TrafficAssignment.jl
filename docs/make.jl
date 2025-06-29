using Documenter
using TrafficAssignment
ENV["DATADEPS_ALWAYS_ACCEPT"] = true

TrafficAssignmentProblem(TransportationNetworks, "SiouxFalls")
TrafficAssignmentProblem(Unified, "San Francisco")

cp(joinpath(@__DIR__, "..", "README.md"), joinpath(@__DIR__, "src", "index.md"); force=true)

makedocs(;
    modules=[TrafficAssignment],
    authors="Changhyun Kwon, Guillaume Dalle",
    sitename="TrafficAssignment.jl",
    format=Documenter.HTML(),
    pages=["Home" => "index.md", "tutorial.md", "api.md"],
)

deploydocs(; repo="github.com/gdalle/TrafficAssignment.jl", devbranch="master")
