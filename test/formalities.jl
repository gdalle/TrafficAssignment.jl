using TestItems

@testitem "Aqua" begin
    using Aqua
    Aqua.test_all(TrafficAssignment; undocumented_names=true)
end

@testitem "Doctests" begin
    using Documenter
    DocMeta.setdocmeta!(
        TrafficAssignment, :DocTestSetup, :(using TrafficAssignment); recursive=true
    )
    Documenter.doctest(TrafficAssignment)
end
