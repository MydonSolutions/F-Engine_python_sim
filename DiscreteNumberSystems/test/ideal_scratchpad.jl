struct NArray{N} <: AbstractArray{Int, N}
	data::Array{Int, N}
end
struct Scalar <: NArray{1}
	data::Int
end
# NArray(data::Array{<:Integer, N}) where {N} = NArray{N}(data)
NArray(data::Integer) = Scalar(data)
# NArray(arr::NArray) = NArray(NArray.data)

Base.size(A::NArray) = size(A.data)
Base.size(A::Scalar) = size(A.data)

Base.@inline function Base.getindex(A::NArray{N}, i::Vararg{Int, N}) where {N}
	@boundscheck checkbounds(A.data, i...)
	@inbounds NArray(getindex(A.data, i...))
	# @inbounds data = getindex(A.data, i...)
	# if length(data) > 1
	# 	NArray(data)
	# else
	# 	data
	# end
end

Base.@inline function Base.setindex!(A::NArray{N}, v, i::Vararg{Int, N}) where {N}
	@boundscheck checkbounds(A.data, i...)
	@inbounds setindex!(A.data, v, i...)
end

Base.convert(::Type{<:Integer}, arr::NArray{1}) = arr.data[1]

function Base.show(io::IO, v::NArray)
	print(v.data)
end

using Random

a = NArray{3}(Random.rand(Int, 3, 3, 3));
a.data[3,3,3]
a[3,3,3]
a.data[3,3,:]
a[3,3,:]
a[3,1:2:3,:] .+ 1
eachindex(a)
map(println, a)

hcat(a[1:2, :, :], Random.rand(Int, 2, 3, 3))