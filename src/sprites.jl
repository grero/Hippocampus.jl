using Makie
using GeometryBasics
using Rotations
using Colors
using CoordinateTransformations

struct Sprite{T<:RGB,T2<:Integer,T3<:Point3, T4<:Point2,T5<:Vec3}
    points::Vector{T3}
    faces::Vector{TriangleFace{T2}}
    normals::Vector{T5}
    uv::Vector{T4}
    img::Matrix{T}
end

"""
    (tt::CoordinateTransformations.LinearMap{T})(sp::Sprite) where T <: Rotation

Apply rotation to `sprite`.
"""
function (tt::CoordinateTransformations.LinearMap{T})(sp::Sprite) where T <: Rotation
    mp = mean(sp.points)
    new_points = tt.(sp.points .- mp)
    new_points .+= convert(eltype(new_points), mp)
    new_points = eltype(sp.points).(new_points)
    new_normals = eltype(sp.normals).(tt.(sp.normals))
    Sprite(new_points, sp.faces, new_normals, sp.uv, sp.img)
end

function (tt::CoordinateTransformations.LinearMap{T})(sp::Sprite) where T <: Translation 
    new_points = eltype(sp.points).([point + tt.linear.translation for point in sp.points])
    Sprite(new_points, sp.faces, sp.normals, sp.uv, sp.img)
end

"""
Map `img` onto a flat 3D surface
"""
function sprite(img,rect=Rect2(-0.5, -0.5, 1.0, 1.0))
    points = decompose(Point3f, rect)
    normals = decompose(GeometryBasics.Normal(Vec3f), rect)
    faces = decompose(GLTriangleFace, rect)
    uv = decompose(UV(Point2f), rect)

    Sprite(points, faces, normals, uv, img)
end

function Makie.convert_arguments(::Type{<:AbstractPlot}, sp::Sprite)
    gb_mesh = Hippocampus.GeometryBasics.Mesh(Hippocampus.GeometryBasics.meta(sp.points; uv=sp.uv, sp.normals), sp.faces)
    PlotSpec(Makie.Mesh, gb_mesh, color=sp.img)
end