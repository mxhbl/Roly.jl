@testset "encoder" begin
    n, s = 4, 3
    enc = Roly.PolyEncoder([[[1], [2], [3]], [[4], [5], [6]], [[7], [8], [9]], [[10], [11], [12]]])
    for i in 1:n, j in 1:s
        @test enc.bwd[enc.fwd[i][j][1]] == [i, j, 1]
    end

    perm = [12, 4, 6, 3, 2, 1, 5, 9, 10, 7, 11, 8]
    permute!(enc, perm)
    for i in 1:n, j in 1:s
        @test enc.bwd[enc.fwd[i][j][1]] == [i, j, 1]
    end
end