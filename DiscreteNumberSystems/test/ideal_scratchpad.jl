struct NArray{N} <: AbstractArray{Integer, N}
	data::Array{Integer, N}
end
struct Scalar
	data::Integer
end

struct CScalar
	real::Integer
	imag::Integer
end
struct CNArray{N} <:AbstractArray{Complex{<:Integer},N}
	real :: NArray{N}
	imag :: NArray{N}
	function CNArray{N}(real :: NArray{N}, imag :: NArray{N}) where {N}
		return new{N}(real,imag)
	end
end

Base.size(A::NArray) = size(A.data)
Base.size(A::Scalar) = size(A.data)
Base.length(A::Scalar) = length(A.data)
Base.size(CA::CNArray) = size(CA.real.data)
Base.size(CS::CScalar) = size(CS.real)
Base.length(CS::CScalar) = length(CS.real)

Base.@inline function Base.getindex(CA::CNArray{N}, i::Vararg{Integer,N}) where {N}
	#Assume that real.data has the same bounds as imag.data
	@boundscheck checkbounds(CA.real.data, i...)
	@inbounds rdata, idata = getindex(CA.real, i...), getindex(CA.imag, i...)
	length(rdata) == 1 ? CScalar(rdata, idata) : CNArray(rdata, idata)
end
Base.@inline function Base.setindex!(CA::CNArray{N},v, i::Vararg{Int,N}) where{N}
	#Assume that real.data has the same bounds as imag.data
	@boundscheck checkbounds(CA.real.data,i...)
	@inbounds CNArray{N}(setindex!(CA.real,v.real,i...),setindex!(CA.imag,v.imag,i...))
end

Base.@inline function Base.getindex(A::NArray{N}, i::Vararg{Int, N}) where {N}
	@boundscheck checkbounds(A.data, i...)
	@inbounds data = getindex(A.data, i...)
	length(data) == 1 ? Scalar(data) : NArray(data)
end

Base.@inline function Base.setindex!(A::NArray{N}, v, i::Vararg{Int, N}) where {N}
	@boundscheck checkbounds(A.data, i...)
	@inbounds setindex!(A.data, v, i...)
end

Base.convert(::Type{<:Integer}, arr::Scalar) = arr.data
Base.convert(::Type{<:Complex{<:Integer}}, arr::CScalar) = arr.real + 1im* arr.imag

function Base.show(io::IO, v::NArray)
	print(v.data)
end
function Base.show(io::IO, v::CNArray)
	print(v.real,v.imag)
end

using Random

a = NArray{3}(Random.rand(Int, 3, 3, 3))
s = NArray{1}(Random.rand(Int,2))

ca = CNArray{3}(a,a)
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
println("s: ",s[:])
a[3,1:2,3] = s[:]
println("a[3,1:2,3] = s[:]",a[3,1:2,3])
println("a[3,3,2:end]",a[3,3,2:end])
println("hcat(a[1:2, :, :], Random.rand(Int, 2, 3, 3)):\n", hcat(a[1:2, :, :], Random.rand(Int, 2, 3, 3)), "\n")
println("ca:\n",ca,"\n")
println("ca.real.data[3,3,3]:\n", ca.real.data[3,3,3], "\n")
println("ca.imag.data[3,3,3]:\n", ca.imag.data[3,3,3], "\n")
println("ca[3,3,3]:\n", ca[3,3,3], "\n")
println("ca.real[3,3,:]:\n", ca.real[3,3,:], "\n")
println("ca[3,3,:]:\n", ca[3,3,:], "\n")
println("ca[3,1:2:3,:] .+ 1:\n", ca[3,1:2:3,:] .+ 1, "\n")
println("eachindex(ca):\n", eachindex(ca), "\n")
println("map(println, ca):\n")
map(println, ca)
println("\n")
println("hcat(ca[1:2, :, :], Random.rand(Complex{Int}, 2, 3, 3)):\n", hcat(ca[1:2, :, :], Random.rand(Complex{Int}, 2, 3, 3)), "\n")
