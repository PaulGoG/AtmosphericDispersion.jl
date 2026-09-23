@testset "PasquillClass" begin
    @test length(PASQUILL_CLASSES) == 6
    @test PASQUILL_CLASSES[1] === PASQUILL_A
    @test PASQUILL_CLASSES[6] === PASQUILL_F
    @test pasquill('D') === PASQUILL_D
    @test pasquill('d') === PASQUILL_D
    @test pasquill("F") === PASQUILL_F
    @test pasquill(PASQUILL_B) === PASQUILL_B
    @test_throws ArgumentError pasquill('G')
    @test_throws ArgumentError pasquill("DE")
    for (i, c) in enumerate(PASQUILL_CLASSES)
        @test classindex(c) == i
        @test pasquill(letter(c)) === c
    end
end
