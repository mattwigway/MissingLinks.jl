@testitem "Break lines" begin
    using ArchGDAL
    import MissingLinks: break_long_line
    include("util.jl")

    # first, simple: a 475m-long straight line
    broken = break_long_line(ArchGDAL.createlinestring([[0.0, 0.0], [475.0, 0.0]]), 250)
    @test length(broken) == 2
    @test all(get_coords_xy(broken[1]) .≈ [[0.0, 0.0], [475 / 2, 0.0]])
    @test all(get_coords_xy(broken[2]) .≈ [[475 / 2, 0.0], [475.0, 0.0]])

    # Diagonal line from 0, 0 to 200, 200 - 280m long
    broken = break_long_line(ArchGDAL.createlinestring([[0.0, 0.0], [200.0, 200.0]]), 250)
    @test length(broken) == 2
    @test all(get_coords_xy(broken[1]) .≈ [[0.0, 0.0], [100.0, 100.0]])
    @test all(get_coords_xy(broken[2]) .≈ [[100.0, 100.0], [200.0, 200.0]])

    # Complex line
    broken = break_long_line(ArchGDAL.createlinestring([
        [0.0, 0.0],       # cumulative distance
        [100.0, 0.0],     # 100m
        [100.0, 50.0],    # 150m
        [100.0, 150.0],   # 250m
        [210.0, 150.0]    # 360m
        ]),
    100)
    @test length(broken) == 4
    @test all(get_coords_xy(broken[1]) .≈ [[0.0, 0.0], [90.0, 0.0]]) # 90m
    @test all(get_coords_xy(broken[2]) .≈ [[90.0, 0.0], [100.0, 0.0], [100.0, 50.0], [100.0, 80.0]]) # 90m
    @test all(get_coords_xy(broken[3]) .≈ [[100.0, 80.0], [100.0, 150.0], [120.0, 150.0]]) # 90m
    @test all(get_coords_xy(broken[4]) .≈ [[120.0, 150.0], [210.0, 150.0]]) # 90m

    # short line
    coords = [[10.0, 5.0], [15.0, 4.0], [17.2, 9.3]]
    broken = break_long_line(ArchGDAL.createlinestring(coords), 250)
    @test length(broken) == 1
    @test all(get_coords_xy(broken[1]) .≈ coords)
end