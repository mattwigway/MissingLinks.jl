@testitem "add_unless_typemax" begin
    import MissingLinks: add_unless_typemax
    for T in [Int8, UInt8, Int16, UInt16, Int32, UInt32, Int64, UInt64]
        # No overflow
        @test add_unless_typemax(one(T), one(T)) == convert(T, 2)
        @test add_unless_typemax(typemax(T), one(T)) == typemax(T)
        @test_throws OverflowError println(T, add_unless_typemax(typemax(T) - one(T), convert(T, 42)))
    end
end