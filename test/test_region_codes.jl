@testset "region codes" begin
    @test region_code("1") == 1
    @test region_code("1.3") == 103
    @test region_code("1.3.1") == 10301
    @test region_code("2.12.45") == 21245
    @test region_level.((1, 99, 103, 9999, 10301)) == (1, 1, 2, 2, 3)
    @test region_path(10301) == "1.3.1"
    @test region_path(region_code("12.3.99")) == "12.3.99"
    @test_throws ArgumentError region_code("1.100")
    @test_throws ArgumentError region_path(1300)
    @test region_ancestor(10301, 1) == 1
    @test region_ancestor(10301, 2) == 103
    @test region_ancestor(103, 2) == 103
    codes = [0 1 103; 10301 104 10401]
    @test in_region(codes, 1) == [false true true; true true true]
    @test in_region(codes, 103) == [false false true; true false false]
    @test in_region(codes, 10301) == [false false false; true false false]
    a = region_flag_attrib([1, 103], ["Northern Hemisphere", "Greenland"])
    @test a[2] == ("flag_meanings" => "Northern_Hemisphere Greenland")
end
