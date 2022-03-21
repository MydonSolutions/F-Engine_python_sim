struct NArray{N} <: AbstractArray{Int, N}
	data::Array{Int, N}
end
struct Scalar
	data::Int
end

NArray(data::Integer) = Scalar(data)

Base.size(A::NArray) = size(A.data)
Base.size(A::Scalar) = size(A.data)

Base.@inline function Base.getindex(A::NArray{N}, i::Vararg{Int, N}) where {N}
	@boundscheck checkbounds(A.data, i...)
	@inbounds NArray(getindex(A.data, i...))
end

Base.@inline function Base.setindex!(A::NArray{N}, v, i::Vararg{Int, N}) where {N}
	@boundscheck checkbounds(A.data, i...)
	@inbounds setindex!(A.data, v, i...)
end

Base.convert(::Type{<:Integer}, arr::Scalar) = arr.data

function Base.show(io::IO, v::NArray)
	print(v.data)
end

using Random

a = NArray{3}(Random.rand(Int, 3, 3, 3))
println("a:\n", a, "\n")
println("a.data[3,3,3]:\n", a.data[3,3,3], "\n")
println("a[3,3,3]:\n", a[3,3,3], "\n")
println("a.data[3,3,:]:\n", a.data[3,3,:], "\n")
println("a[3,3,:]:\n", a[3,3,:], "\n")
println("a[3,1:2:3,:] .+ 1:\n", a[3,1:2:3,:] .+ 1, "\n")
println("eachindex(a):\n", eachindex(a), "\n")
println("map(println, a):\n")
map(println, a)
println("\n")

println("hcat(a[1:2, :, :], Random.rand(Int, 2, 3, 3)):\n", hcat(a[1:2, :, :], Random.rand(Int, 2, 3, 3)), "\n")