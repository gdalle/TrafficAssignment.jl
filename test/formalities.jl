using TestItems

@testitem "Aqua" begin
    using Aqua
    Aqua.test_all(TrafficAssignment; undocumented_names=true)
end

@testitem "JET" begin
    using JET
    JET.test_package(TrafficAssignment; target_defined_modules=true)
end

@testitem "JuliaFormatter" begin
    using JuliaFormatter
    @test JuliaFormatter.format(TrafficAssignment; overwrite=false)
end

@testitem "Doctests" begin
    using Documenter
    DocMeta.setdocmeta!(
        TrafficAssignment, :DocTestSetup, :(using TrafficAssignment); recursive=true
    )
    Documenter.doctest(TrafficAssignment)
end
