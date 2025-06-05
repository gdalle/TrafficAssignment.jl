using TestItems

@testitem "Aqua" begin
    using Aqua
    Aqua.test_all(TrafficAssignment)
end

@testitem "JET" begin
    using JET
    JET.test_package(TrafficAssignment; target_defined_modules=true)
end

@testitem "JuliaFormatter" begin
    using JuliaFormatter
    @test JuliaFormatter.format(TrafficAssignment; overwrite=false)
end

@testitem "Undocumented names" begin
    if isdefined(Base.Docs, :undocumented_names)
        @test isempty(Base.Docs.undocumented_names(TrafficAssignment))
    end
end

@testitem "Doctests" begin
    using Documenter
    Documenter.doctest(TrafficAssignment)
end
