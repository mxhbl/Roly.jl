@testset "geometry" begin

    # Triangle Offsets
    geom = PolygonGeometry(3, 1.)

    @test all(Roly.attachment_offset(1, 1, geom, geom) .≈ ([0., -inv(√3)], Roly.Angle(1.)))
    @test all(Roly.attachment_offset(2, 2, geom, geom) .≈ ([-0.5, inv(√3)/2], Roly.Angle(1.)))
    @test all(Roly.attachment_offset(3, 3, geom, geom) .≈ ([0.5, inv(√3)/2], Roly.Angle(1.)))
    @test all(Roly.attachment_offset(1, 2, geom, geom) .≈ ([0., -inv(√3)], Roly.Angle(1. + 2/3)))
    @test all(Roly.attachment_offset(2, 3, geom, geom) .≈ ([-0.5, inv(√3)/2], Roly.Angle(1. + 2/3)))
    @test all(Roly.attachment_offset(3, 1, geom, geom) .≈ ([0.5, inv(√3)/2], Roly.Angle(1. + 2/3)))
    @test all(Roly.attachment_offset(1, 3, geom, geom) .≈ ([0., -inv(√3)], Roly.Angle(1. - 2/3)))
    @test all(Roly.attachment_offset(2, 1, geom, geom) .≈ ([-0.5, inv(√3)/2], Roly.Angle(1. - 2/3)))
    @test all(Roly.attachment_offset(3, 2, geom, geom) .≈ ([0.5, inv(√3)/2], Roly.Angle(1. - 2/3)))

    # Inverse Offsets
    @test only(Roly.face_pairs(Roly.Point(0., -inv(√3)), Roly.Angle(0.), Roly.Angle(1.), geom, geom)) == (1, 1)
    @test only(Roly.face_pairs(Roly.Point(0., -inv(√3)), Roly.Angle(0.), Roly.Angle(1. + 2/3), geom, geom)) == (1, 2)
    @test only(Roly.face_pairs(Roly.Point(0., -inv(√3)), Roly.Angle(0.), Roly.Angle(1. - 2/3), geom, geom)) == (1, 3)
    @test only(Roly.face_pairs(Roly.Point(-0.5, inv(√3)/2), Roly.Angle(0.), Roly.Angle(1.), geom, geom)) == (2, 2)
    @test only(Roly.face_pairs(Roly.Point(-0.5, inv(√3)/2), Roly.Angle(0.), Roly.Angle(1. + 2/3), geom, geom)) == (2, 3)
    @test only(Roly.face_pairs(Roly.Point(-0.5, inv(√3)/2), Roly.Angle(0.), Roly.Angle(1. - 2/3), geom, geom)) == (2, 1)
    @test only(Roly.face_pairs(Roly.Point(0.5, inv(√3)/2), Roly.Angle(0.), Roly.Angle(1.), geom, geom)) == (3, 3)
    @test only(Roly.face_pairs(Roly.Point(0.5, inv(√3)/2), Roly.Angle(0.), Roly.Angle(1. + 2/3), geom, geom)) == (3, 1)
    @test only(Roly.face_pairs(Roly.Point(0.5, inv(√3)/2), Roly.Angle(0.), Roly.Angle(1. - 2/3), geom, geom)) == (3, 2)
end