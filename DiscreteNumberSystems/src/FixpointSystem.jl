module FixpointSystem

using Printf

import Base: convert, float, +, -, *, <<, >>

export FixpointScheme, Fixpoint, FixpointArray, CFixpoint, CFixpointArray

export convert, zeros, normalise, cast, sum, *, +, -, power, conj, clamp_wrap, quantise, copy, size, length, show, getindex, setindex!, hcat, lastindex, >>, <<

#########################################################################################
# FixpointScheme Structure
#########################################################################################

"""
Type which holds information about the Fixpoint value and influences its treatment
under arithmetic, logical and conversion operations.

See also: [`Fixpoint`](@ref)
"""
struct FixpointScheme 
    bits :: Integer
    fraction :: Integer
    min :: Integer
    max :: Integer
    unsigned :: Bool
    range :: Integer
    scale :: Integer
    ovflw_behav :: String
    undflw_behav :: RoundingMode
    function FixpointScheme(bits::Integer, fraction::Integer;
         min_int::Union{Integer, Nothing}=nothing,max_int::Union{Integer, Nothing}=nothing,
         unsigned::Union{Bool}=false, ovflw_behav::Union{String}="WRAP",
         undflw_behav::Union{RoundingMode}=RoundNearest)
        scale = 2^fraction;
        range = 2^bits;
        if min_int === nothing
            min = unsigned ? 0 : (-range//2); 
        else
            min = min_int;
        end
        if max_int === nothing
            max = unsigned ? range -1 : (range//2) -1;
        else
            max = max_int; 
        end
        new(bits, fraction, min, max, unsigned, range, scale, ovflw_behav, undflw_behav);
    end
end

#########################################################################################
# Fixpoint Structures
#--------------------
# Each Fixpoint structure has the necessary constructors to convert from Real/Complex
# to Fixpoint. Herein declared are also the Base.convert function overloads
# to convert Fixpoint types back to the Real/Complex domain.
#########################################################################################
"""
Fixpoint type that accepts a single integer accompanied by a FixpointScheme that governs it's handling.

See also: [`FixpointScheme`](@ref)
"""
struct Fixpoint <: Number
    data :: Integer
    scheme :: FixpointScheme
    Fixpoint(fx_data::Integer, scheme::FixpointScheme) = new(fx_data,scheme);
    function Fixpoint(fl_data::Real, scheme::FixpointScheme)
        prod = fl_data * scheme.scale;
        if (scheme.ovflw_behav == "SATURATE")
            data = clamp(round(Integer, prod, scheme.undflw_behav), scheme.min, scheme.max);
        elseif (scheme.ovflw_behav == "WRAP")
            data = clamp_wrap(round(Integer, prod, scheme.undflw_behav), scheme.min, scheme.max);
        else 
            error("No recognisable overflow method specified");
        end
        return new(data,scheme);
    end
end

"""
Overload the Base.convert function to convert from Fixpoint to Real. This process discards
the scheme, only returning the floating point value.
"""
function Base.convert(::Type{<:AbstractFloat},f::Fixpoint) :: AbstractFloat
    return float(f.data)/f.scheme.scale;
end

"""
Convenience function to convert from Fixpoint to a AbstractFloat.
"""
function Base.float(f :: Fixpoint) :: AbstractFloat
    return convert(AbstractFloat, f)
end

"""
FixpointArray type that accepts an integer array accompanied by a FixpointScheme that
governs it's handling.

See also: [`FixpointScheme`](@ref)
"""
struct FixpointArray{N} <: AbstractArray{Fixpoint,N}
    data :: Array{Integer,N}
    scheme :: FixpointScheme
    function FixpointArray{N}(data::Array{<:Integer,N}, scheme::FixpointScheme) where {N}
        new(data, scheme)
    end
    function FixpointArray{N}(fl_data::AbstractArray{<:AbstractFloat,N}, scheme::FixpointScheme) where {N}
        new(map(d -> Fixpoint(d, scheme).data, fl_data), scheme)
    end
end

"""
Overload the Base.convert function to convert from FixpointArray{N} to Array{AbstractFloat,N}. This process discards
the scheme, only returning the floating point array.
"""
function Base.convert(float_arr::Type{<:Array{<:AbstractFloat,N}},f_arr::FixpointArray{N}) where {N}
    return convert.(eltype(float_arr),f_arr.data)./f_arr.scheme.scale
end

"""
Convenience function to convert from FixpointArray to an array of Float64.
"""
function Base.float(f :: FixpointArray{N}) :: Array{<:Float64, N} where {N}
    return convert(Array{<:Float64, N}, f)
end

"""
CFixpoint is the complex extension of Fixpoint that holds two Fixpoint types (real and imag) as its 
Real and Imaginary parts.

See also: [`Fixpoint`](@ref)
"""
struct CFixpoint <: Number
    real :: Fixpoint
    imag :: Fixpoint
    function CFixpoint(real :: Fixpoint, imag :: Fixpoint)
        if real.scheme != imag.scheme
            error("Real and Imag Fixpoint values must have the same scheme.");
        else
            return new(real,imag);
        end
    end
    function CFixpoint(real::Real,imag::Real,scheme::FixpointScheme)
        new(Fixpoint(real,scheme), Fixpoint(imag,scheme))
    end
    CFixpoint(real::Integer,imag::Integer,scheme::FixpointScheme) = new(Fixpoint(real,scheme), Fixpoint(imag,scheme))
end
"""
Overload the Base.convert function to convert from CFixpoint to Complex. This process discards
the scheme, only returning the complex point value.
"""
function Base.convert(::Type{<:Complex{<:AbstractFloat}}, cf::CFixpoint)::Complex{<:AbstractFloat}
    return convert(AbstractFloat, cf.real) +
        convert(AbstractFloat, cf.imag) * 1im;
end

"""
Convenience function to convert from CFixpoint to an array of Complex64.
"""
function Base.float(cf :: CFixpoint) :: ComplexF64
    return convert(ComplexF64, cf)
end

"""
CFixpointArray is the complex extension of Fixpoint that holds two FixpointArray types (real and imag) as its 
Real and Imaginary parts.

See also: [`FixpointArray`](@ref)
"""
struct CFixpointArray{N} <:AbstractArray{CFixpoint,N}
    real :: FixpointArray{N}
    imag :: FixpointArray{N}
    function CFixpointArray{N}(real :: Array{<:Integer,N}, imag :: Array{<:Integer,N}, scheme :: FixpointScheme) where {N}
        return new(FixpointArray{N}(real,scheme),FixpointArray{N}(imag,scheme))
    end
    function CFixpointArray{N}(real :: FixpointArray{N}, imag :: FixpointArray{N}) where {N}
        if real.scheme != imag.scheme
            error("Real and Imag Fixpoint values must have the same scheme.")
        end
       new(real,imag)
    end
end

"""
Overload the Base.convert function to convert from CFixpointArray{N} to Array{Real,N}. This process discards
the scheme, only returning the complex floating point array.
"""
function Base.convert(complex_arr::Type{AbstractArray{<:Complex{<:AbstractFloat},N}},cf_arr::CFixpointArray{N}) ::
                                                         AbstractArray{<:Complex{<:AbstractFloat},N} where {N}                                                    
    return convert(real(eltype(complex_arr)), cf_arr.real) .+ 
    convert(imag(eltype(complex_arr), cf_arr.imag)) .* 1im;
end

"""
Convenience function to convert from CFixpointArray to an array of ComplexF64.
"""
function Base.float(cf :: CFixpointArray{N}) :: Array{ComplexF64,N} where {N}
    return float(cf.real) + 1im * float(cf.imag)
end

#########################################################################################
# Float parsing funtions
#########################################################################################
# """
# Creates a Fixpoint populated with zero values. Creates a CFixpoint if complex is set to true.
# """
# function Base.zeros(fx_scheme :: FixpointScheme, dims :: Tuple; complex :: Bool = false) :: Union{Fixpoint,CFixpoint}
#     if complex
#         return fromComplex(zeros(Float64, dims),zeros(Float64, dims),fx_scheme);
#     else
#         return fromFloat(zeros(Float64, dims), fx_scheme);
#     end
# end

# """
# Fit all data within the min/max values for Fixpoint
# """
# function normalise(f_val :: Fixpoint) :: Fixpoint
#     return Fixpoint(clamp.(f_val.data,f_val.scheme.min,f_val.scheme.max),f_val.scheme);
# end

# """
# Fit all data within the min/max values for CFixpoint
# """
# function normalise(cf_val :: CFixpoint) :: CFixpoint
#     return CFixpoint(normalise(cf_val.real),normalise(cf_val.imag));
# end

# """
# Cast Fixpoint type to Fixpoint type of new scheme
# """
# function cast(f_val :: Fixpoint, f_scheme :: FixpointScheme) :: Fixpoint
#     return Fixpoint(f_val.data,f_scheme);
# end

# """
# Cast CFixpoint type to CFixpoint type of new scheme
# """
# function cast(cf_val :: CFixpoint, f_scheme :: FixpointScheme) :: CFixpoint
#     return CFixpoint(Fixpoint(cf_val.real.data,f_scheme),Fixpoint(cf_val.imag.data,f_scheme));
# end


# #######################################################################################
# # Arithmetic functions
# #######################################################################################

"""
Overload sum function to take a FixpointArray argument.

See also: [`sum`](@ref)
"""
function Base.sum(f :: FixpointArray{N}; dims :: Union{Integer,Colon}=:) where {N}
    sum_val = sum(f.data, dims=dims);
    bits = f.scheme.bits + ceil.(Integer,log2.(length(f.data)/length(sum_val)));
    scheme = FixpointScheme(bits, f.scheme.fraction, min_int=f.scheme.min,
    max_int=f.scheme.max, unsigned=f.scheme.unsigned, ovflw_behav=f.scheme.ovflw_behav,
    undflw_behav=f.scheme.undflw_behav);
    return length(sum_val) == 1 ? Fixpoint(sum_val,scheme) : FixpointArray{N}(sum_val,scheme)
end

"""
Overload sum function to take a CFixpointArray argument.
See also: [`sum`](@ref)
"""
function Base.sum(cf :: CFixpointArray{N}; dims :: Union{Integer,Colon}=:) where {N}
    r_sum_val = sum(cf.real,dims=dims);
    i_sum_val = sum(cf.imag,dims=dims);
    return length(r_sum_val) == 1 ? CFixpoint(r_sum_val,i_sum_val) : CFixpointArray{N}(r_sum_val, i_sum_val)
end

"""
Overload * function to take Fixpoint arguments.

See also: [`*`](@ref)
"""
function *(a :: Fixpoint, b :: Fixpoint) :: Fixpoint
    prod_val = a.data * b.data;
    bits = a.scheme.bits + b.scheme.bits;
    fraction = a.scheme.fraction + b.scheme.fraction;
    unsigned = a.scheme.unsigned & b.scheme.unsigned;
    scheme = FixpointScheme(bits, fraction, unsigned=unsigned, 
    ovflw_behav=a.scheme.ovflw_behav, undflw_behav=a.scheme.undflw_behav);
    return Fixpoint(prod_val,scheme);
end

"""
Overload * function to take FixpointArray arguments.

See also: [`*`](@ref)
"""
function *(a :: FixpointArray{N}, b :: FixpointArray{N}) :: FixpointArray{N} where {N}
    prod_val = a.data .* b.data;
    bits = a.scheme.bits + b.scheme.bits;
    fraction = a.scheme.fraction + b.scheme.fraction;
    unsigned = a.scheme.unsigned & b.scheme.unsigned;
    scheme = FixpointScheme(bits, fraction, unsigned=unsigned, 
    ovflw_behav=a.scheme.ovflw_behav, undflw_behav=a.scheme.undflw_behav);
    return FixpointArray{N}(prod_val,scheme);
end

"""
Overload * function to take CFixpoint arguments.

See also: [`*`](@ref)
"""
function *(a :: CFixpoint, b :: CFixpoint) :: CFixpoint
    function cmult(a, b, c, d)
        # Real part x = a*c - b*d
        x = (a*c)-(b*d);
        # Imaginary part y = a*d + b*c
        y = (a*d)+(b*c);
        return x, y;
    end
    out_real, out_imag = cmult(a.real, a.imag, b.real, b.imag);
    return CFixpoint(out_real, out_imag); 
end

"""
Overload * function to take CFixpointArray arguments.

See also: [`*`](@ref)
"""
function *(a :: CFixpointArray{N}, b :: CFixpointArray{N}) :: CFixpointArray{N} where {N}
    function cmult(a, b, c, d)
        # Real part x = a*c - b*d
        x = (a*c)-(b*d);
        # Imaginary part y = a*d + b*c
        y = (a*d)+(b*c);
        return x, y;
    end
    out_real, out_imag = cmult(a.real, a.imag, b.real, b.imag);
    return CFixpointArray{N}(out_real, out_imag); 
end

"""
Overload + function to take Fixpoint arguments.

See also: [`+`](@ref)
"""
function +(a :: Fixpoint, b :: Fixpoint) :: Fixpoint
    if (a.scheme.scale != b.scheme.scale)
        error("Addition performed between two FixpointArray values of differing scales.");
    end
    add_val = a.data + b.data;
    bits = max(a.scheme.bits,b.scheme.bits) + 1;
    unsigned = a.scheme.unsigned & b.scheme.unsigned;
    scheme = FixpointScheme(bits, a.scheme.fraction, unsigned=unsigned, 
    ovflw_behav=a.scheme.ovflw_behav, undflw_behav=a.scheme.undflw_behav);
    return Fixpoint(add_val,scheme);
end

"""
Overload + function to take FixpointArray arguments.

See also: [`+`](@ref)
"""
function +(a :: FixpointArray{N}, b :: FixpointArray{N}) :: FixpointArray{N} where {N}
    if (a.scheme.scale != b.scheme.scale)
        error("Addition performed between two FixpointArray values of differing scales.");
    end
    add_val = a.data .+ b.data;
    bits = max(a.scheme.bits,b.scheme.bits) + 1;
    unsigned = a.scheme.unsigned & b.scheme.unsigned;
    scheme = FixpointScheme(bits, a.scheme.fraction, unsigned=unsigned, 
    ovflw_behav=a.scheme.ovflw_behav, undflw_behav=a.scheme.undflw_behav);
    return FixpointArray{N}(add_val,scheme);
end

"""
Overload + function to take CFixpoint arguments.

See also: [`+`](@ref)
"""
function +(a :: CFixpoint, b :: CFixpoint) :: CFixpoint
    r_sum = a.real + b.real;
    i_sum = a.imag + b.imag;
    return CFixpoint(r_sum, i_sum);
end

"""
Overload + function to take CFixpointArray arguments.

See also: [`+`](@ref)
"""
function +(a :: CFixpointArray{N}, b :: CFixpointArray{N}) :: CFixpointArray{N} where {N}
    r_sum = a.real + b.real;
    i_sum = a.imag + b.imag;
    return CFixpointArray{N}(r_sum, i_sum);
end

"""
Overload - function to take Fixpoint arguments.
        
See also: [`-`](@ref)
"""
function -(a :: Fixpoint, b :: Fixpoint) :: Fixpoint
    if (a.scheme.scale != b.scheme.scale)
        error("Subtraction performed between two Fixpoint values of differing scales.");
    end 
    sub_val = a.data - b.data;
    bits = max(a.scheme.bits,b.scheme.bits) + 1;
    unsigned = a.scheme.unsigned & b.scheme.unsigned;
    scheme = FixpointScheme(bits, a.scheme.fraction, unsigned=unsigned, 
    ovflw_behav=a.scheme.ovflw_behav, undflw_behav=a.scheme.undflw_behav);
    return Fixpoint(sub_val,scheme);
end

"""
Overload - function to take FixpointArray arguments.
        
See also: [`-`](@ref)
"""
function -(a :: FixpointArray{N}, b :: FixpointArray{N}) :: FixpointArray{N} where {N}
    if (a.scheme.scale != b.scheme.scale)
        error("Subtraction performed between two Fixpoint values of differing scales.");
    end 
    sub_val = a.data .- b.data;
    bits = max(a.scheme.bits,b.scheme.bits) + 1;
    unsigned = a.scheme.unsigned & b.scheme.unsigned;
    scheme = FixpointScheme(bits, a.scheme.fraction, unsigned=unsigned, 
    ovflw_behav=a.scheme.ovflw_behav, undflw_behav=a.scheme.undflw_behav);
    return FixpointArray{N}(sub_val,scheme);
end

"""
Overload - function to take CFixpoint arguments.
        
See also: [`-`](@ref)
"""
function -(a :: CFixpoint, b :: CFixpoint) :: CFixpoint
    r_sub = a.real - b.real;
    i_sub = a.imag - b.imag;
    return CFixpoint(r_sub, i_sub);
end

"""
Overload - function to take CFixpointArray arguments.
        
See also: [`-`](@ref)
"""
function -(a :: CFixpointArray{N}, b :: CFixpointArray{N}) :: CFixpointArray{N} where {N}
    r_sub = a.real - b.real;
    i_sub = a.imag - b.imag;
    return CFixpointArray{N}(r_sub, i_sub);
end

# """
# Returns power of the Fixpoint value given = f.data * f.data.
# """
# function power(f :: Fixpoint) :: Array{Integer}
#     return f.data .* f.data;
# end

# """
# Returns power of the CFixpoint value given = cf * conj(cf).
# See also: [`conj`](@ref)
# """
# function power(cf :: CFixpoint) :: Array{Integer}
#     res = copy(cf) * conj(cf);
#     return res.real;
# end

# """
# Returns conjuage of the CFixpoint value given.
# """
# function Base.conj(cf :: CFixpoint) :: CFixpoint
#     i_res = copy(cf.imag);
#     i_res.data = - copy(cf.imag.data);
#     return CFixpoint(cf.real, i_res);
# end

# """
# Returns conjuage of the CFixpoint value given.
# ! implies inline operation.
# """
# function Base.conj!(cf :: CFixpoint) :: CFixpoint
#     cf.imag.data = -cf.imag.data;
# end
# #######################################################################################
# # Misc Fixpoint type handling functions
# #######################################################################################

"""
Does a clamp operation but wraps the value to min/max rather than saturate 
the value like standard clamp.

See also: [`clamp`](@ref)
"""
function clamp_wrap(i :: Integer, min :: Integer, max :: Integer)
    return ((i - min) % (max - min)) + min;
end

"""
Does a clamp operation but wraps the value to min/max rather than saturate 
the value like standard clamp. Takes a Fixpoint type in this instance.

See also: [`clamp_wrap`](@ref)
"""
function clamp_wrap(f :: Fixpoint, min :: Integer, max :: Integer) ::Fixpoint
    clamp_val = ((f.data .- min) .% (min - max)) .+ min;
    scheme = FixpointScheme(f.scheme.bits,f.scheme.fraction,unsigned=f.scheme.unsigned,
    max_int=max, min_int=min,
    ovflw_behav=f.scheme.ovflw_behav, undflw_behav=f.scheme.undflw_behav);
    return Fixpoint(clamp_val,scheme);        
end

"""
Requantise the data contained in fxpt according to the new scheme provided.
"""
function quantise(fxpt :: FixpointArray{N}, scheme :: FixpointScheme) :: FixpointArray{N} where {N}
    return FixpointArray{N}(float(fxpt), scheme);
end

"""
Requantise the data contained in cfxpt according to the new scheme provided.
"""
function quantise(cfxpt :: CFixpointArray{N}, scheme :: FixpointScheme) :: CFixpointArray{N} where {N}
    return CFixpointArray{N}(float(cfxpt),scheme);
end

# """
# Overload copy() function to copy Fixpoint by value as opposed to reference.

# See also: [`copy`](@ref)
# """
# function Base.copy(f :: Fixpoint) :: Fixpoint
#     tmpscheme = FixpointScheme(f.scheme.bits, f.scheme.fraction, min_int=f.scheme.min,
#     max_int=f.scheme.max, unsigned=f.scheme.unsigned, ovflw_behav=f.scheme.ovflw_behav,
#     undflw_behav=f.scheme.undflw_behav);
#     return Fixpoint(copy(f.data),tmpscheme);
# end

# """
# Overload copy() function to copy CFixpoint by value as opposed to reference.
# """
# function Base.copy(cf :: CFixpoint) :: CFixpoint
#     return CFixpoint(copy(cf.real),copy(cf.imag));
# end

"""
Overload size() function to accept Fixpoint.
"""
function Base.size(f::Fixpoint)
    return size(f.data);
end

"""
Overload size() function to accept FixpointArray.
"""
function Base.size(f::FixpointArray{N}) where {N}
    return size(f.data);
end

"""
Overload size() function to accept CFixpoint.
"""
function Base.size(cf::CFixpoint)
    return size(cf.data);
end

"""
Overload size() function to accept CFixpointArray.
"""
function Base.size(cf::CFixpointArray{N}) where {N}
    return size(cf.real);
end

"""
Overload length() function to accept Fixpoint
"""
function Base.length(f :: Fixpoint) :: Integer 
    return length(f.data)
end

"""
Overload length() function to accept CFixpoint
"""
function Base.length(cf :: CFixpoint) :: Integer 
    return length(cf.data)
end

# """
# Overload show function for printing out Fixpoint summary.
# See also: [`show`](@ref)
# """
# function Base.show(io::IO, f :: Fixpoint)
#     @printf(io,"Fixpoint real %s (%d, %d), shape %s", f.scheme.unsigned ? "unsigned" : "signed",f.scheme.bits, f.scheme.fraction, size(f.data));
# end

# """
# Overload show function for printing out Fixpoint summary.
# See also: [`show`](@ref)
# """
# function Base.show(io::IO, cf :: CFixpoint)
#     @printf(io,"CFixpoint complex %s (%d, %d), shape %s", cf.real.scheme.unsigned ? "unsigned" : "signed",cf.real.scheme.bits, cf.real.scheme.fraction, size(cf.real.data))
# end

# #######################################################################################
# # Dimension mangling Fixpoint type handling functions
# # Note that indexing functions depend on LinearIndexing see:
# # https://docs.julialang.org/en/v1/manual/interfaces/
# #######################################################################################

"""
Overload non-scalar-indexing getindex function for accessing data elements out Fixpoint type.
"""
Base.@inline function Base.getindex(f :: FixpointArray{N}, i...) :: Union{Fixpoint, FixpointArray} where {N}
    @boundscheck checkbounds(f.data, i...);
	@inbounds data = getindex(f.data, i...);
	length(data) == 1 ? Fixpoint(data, f.scheme) : FixpointArray{ndims(data)}(data, f.scheme);
end

"""
Overload non-scalar-indexing getindex function for accessing data elements out CFixpoint type.
"""
Base.@inline function Base.getindex(cf :: CFixpointArray{N}, i...) :: Union{CFixpoint, CFixpointArray} where {N}
    @boundscheck checkbounds(cf.real, i...);
	@inbounds rdata, idata = getindex(cf.real, i...), getindex(cf.imag, i...)
    length(rdata) == 1 ? CFixpoint(rdata, idata) : CFixpointArray{ndims(rdata)}(rdata, idata);
end

"""
Overload scalar setindex function for setting integer data elements of FixpointArray type.
"""
Base.@inline function Base.setindex!(f::FixpointArray{N}, v::Integer, i::Vararg{Int, N}) where {N}
	@boundscheck checkbounds(f.data, i...)
	@inbounds setindex!(f.data, v, i...)
end

"""
Overload scalar setindex function for setting Float data elements of FixpointArray type.
"""
Base.@inline function Base.setindex!(f::FixpointArray{N}, v::AbstractFloat, i::Vararg{Int, N}) where {N}
    @boundscheck checkbounds(f.data, i...)
    @inbounds setindex!(f.data, Fixpoint(v, f.scheme).data, i...)
end

"""
Overload scalar setindex function for setting Fixpoint data elements of FixpointArray type.
"""
Base.@inline function Base.setindex!(f::FixpointArray{N}, v::Fixpoint, i::Vararg{Int, N}) where {N}
    @boundscheck checkbounds(f.data, i...)
    if f.scheme == v.scheme
        @inbounds setindex!(f.data, v.data, i...)
    else
        @inbounds setindex!(f.data, Fixpoint(float(v), f.scheme).data, i...)
    end
end

"""
Overload non-scalar setindex! function for setting integer data elements of FixpointArray type.
"""
Base.@inline function Base.setindex!(f::FixpointArray{N}, v::AbstractArray{<:Integer, M}, i...) where {N, M}
	@boundscheck checkbounds(f.data, i...)
	@inbounds setindex!(f.data, v, i...)
end

"""
Overload non-scalar setindex! function for setting real data elements of FixpointArray type.
"""
Base.@inline function Base.setindex!(f::FixpointArray{N}, v::AbstractArray{<:Real, M}, i...) where {N, M}
	@boundscheck checkbounds(f.data, i...)
	@inbounds setindex!(f.data, FixpointArray{ndims(v)}(v, f.scheme).data, i...)
end

"""
Overload non-scalar setindex! function for setting Fixpoint elements of FixpointArray type.
"""
Base.@inline function Base.setindex!(f::FixpointArray{N}, v::FixpointArray{M}, i...) where {N, M}
	@boundscheck checkbounds(f.data, i...)
	@boundscheck checkbounds(v, i...)
    if f.scheme == v.scheme
	    @inbounds setindex!(f.data, v.data, i...)
    else
	    @inbounds setindex!(f.data, FixpointArray{M}(float(v), f.scheme).data, i...)
    end
end

"""
Overload scalar setindex function for setting integer data elements of CFixpointArray type.
"""
Base.@inline function Base.setindex!(cf::CFixpointArray{N}, v::Complex{<:Integer}, i::Vararg{Int, N}) where {N}
	@boundscheck checkbounds(cf.real.data, i...)
	@inbounds setindex!(cf.real.data, real(v), i...)
	@inbounds setindex!(cf.imag.data, imag(v), i...)
end

"""
Overload scalar setindex function for setting Float data elements of CFixpointArray type.
"""
Base.@inline function Base.setindex!(cf::CFixpointArray{N}, v::Complex{<:AbstractFloat}, i::Vararg{Int, N}) where {N}
    @boundscheck checkbounds(cf.real.data, i...)
    @inbounds setindex!(cf.real.data, Fixpoint(real(v), cf.real.scheme).data, i...)
    @inbounds setindex!(cf.imag.data, Fixpoint(imag(v), cf.imag.scheme).data, i...)
end

"""
Overload scalar setindex function for setting Fixpoint data elements of CFixpointArray type.
"""
Base.@inline function Base.setindex!(cf::CFixpointArray{N}, v::CFixpoint, i::Vararg{Int, N}) where {N}
    @boundscheck checkbounds(cf.real.data, i...)
    if cf.real.scheme == v.real.scheme
        @inbounds setindex!(cf.real.data, v.real.data, i...)
        @inbounds setindex!(cf.imag.data, v.imag.data, i...)
    else
        @inbounds setindex!(cf.real.data, Fixpoint(float(v.real), cf.real.scheme).data, i...)
        @inbounds setindex!(cf.imag.data, Fixpoint(float(v.imag), cf.imag.scheme).data, i...)
    end
end

"""
Overload non-scalar setindex! function for setting complex integer data elements of CCFixpointArray type.
"""
Base.@inline function Base.setindex!(cf::CFixpointArray{N}, v::AbstractArray{<:Complex{<:Integer}, M}, i...) where {N, M}
	@boundscheck checkbounds(cf.real, i...)
	@inbounds setindex!(cf.real, real(v), i...)
	@inbounds setindex!(cf.imag, imag(v), i...)
end

"""
Overload non-scalar setindex! function for setting Complex data elements of CFixpointArray type.
"""
Base.@inline function Base.setindex!(cf::CFixpointArray{N}, v::AbstractArray{<:Complex{<:AbstractFloat}, M}, i...) where {N, M}
	@boundscheck checkbounds(cf.real, i...)
	@inbounds setindex!(cf.real, FixpointArray{ndims(v)}(real(v), cf.real.scheme), i...)
	@inbounds setindex!(cf.imag, FixpointArray{ndims(v)}(imag(v), cf.real.scheme), i...)
end

"""
Overload non-scalar setindex! function for setting CFixpoint elements of CFixpointArray type.
"""
Base.@inline function Base.setindex!(cf::CFixpointArray{N}, v::CFixpointArray{M}, i...) where {N, M}
    if cf.real.scheme == v.real.scheme
        @boundscheck checkbounds(cf.real, i...)
	    @inbounds setindex!(cf.real, v.real, i...)
	    @inbounds setindex!(cf.imag, v.imag, i...)
    else
        # invoke setindex!(c::CFixpointArray{N}, ::AbstractArray{<:Complex{<:AbstractFloat}, M}, i...)
	    @inbounds setindex!(cf, float(v), i...)
    end
end

"""
Overload hcat function to handle horizontal concatenation of FixpointArray.
Requires that schemes match.
"""
function Base.hcat(f_1 :: FixpointArray{N}, f_2 :: FixpointArray{N}) :: FixpointArray where {N}
    #Check schemes match:
    if f_1.scheme == f_2.scheme
        data = hcat(f_1.data,f_2.data)
        return FixpointArray{ndims(data)}(data,f_1.scheme);
    else
        error("FixpointArray args don't share the same scheme.");
    end
end

"""
Overload hcat function to handle horizontal concatenation of CFixpointArray.
Requires that schemes match.
"""
function Base.hcat(cf_1 :: CFixpointArray{N}, cf_2 :: CFixpointArray{N}) :: CFixpointArray{N} where {N}
    #Check real schemes match - imag will match:
    if cf_1.real.scheme == cf_2.real.scheme
        real = hcat(cf_1.real,cf_2.real)
        return CFixpointArray{ndims(real)}(real, hcat(cf_1.imag,cf_2.imag));
    else
        error("CFixpointArray args don't share the same scheme.");
    end
end

# """
# Overload lastindex function to handle slicing of Fixpoint with end
# """
# function Base.lastindex(f :: Fixpoint) :: Int
#     return length(f.data)
# end

# """
# Overload lastindex function to handle slicing of CFixpoint with end
# """
# function Base.lastindex(cf :: CFixpoint) :: Int
#     return length(cf.real.data)
# end

# """
# Overload axes function for Fixpoint
# """
# function Base.axes(f :: Fixpoint) :: AbstractUnitRange{<:Integer}
#     return map(OneTo, size(f))
# end

# """
# Overload axes function for CFixpoint
# """
# function Base.axes(cf :: Fixpoint) :: AbstractUnitRange{<:Integer}
#     return map(OneTo, size(cf))
# end

# # """
# # Overload similar function for Fixpoint
# # """
# # function Base.similar(f :: Fixpoint) 
# #     return 

# #######################################################################################
# # Logical operator functions
# #######################################################################################
"""
Overload >> for Fixpoint.
Apply 'steps' (>=0) right shifts to fxpt. Cannot use >> operator here since we must control rounding.

See also: [`>>`](@ref)
"""
function >>(fxpt :: Fixpoint, steps :: Integer) :: Fixpoint
    if (steps < 0)
        error("Integer value for steps must be greater than or equal to zero.");
    else
        return Fixpoint(round(Integer, fxpt.data/(2^steps), fxpt.scheme.undflw_behav), fxpt.scheme);
    end
end

"""
Overload >> for FixpointArray.
Apply 'steps' (>=0) right shifts to fxpt. Cannot use >> operator here since we must control rounding.

See also: [`>>`](@ref)
"""
function >>(fxpt :: FixpointArray{N}, steps :: Integer) :: FixpointArray{N} where {N}
    if (steps < 0)
        error("Integer value for steps must be greater than or equal to zero.");
    else
        return FixpointArray{N}(round.(Integer, fxpt.data/(2^steps), fxpt.scheme.undflw_behav), fxpt.scheme);
    end
end

"""
Overload >> for CFixpoint.
Apply 'steps' (>=0) right shifts to cfxpt. Cannot use >> operator here since we must control rounding.

See also: [`>>`](@ref)
"""
function >>(cfxpt :: CFixpoint, steps :: Integer) :: CFixpoint
    t_real = cfxpt.real >> steps;
    t_imag = cfxpt.imag >> steps;
    return CFixpoint(t_real, t_imag);
end

"""
Overload >> for CFixpointArray.
Apply 'steps' (>=0) right shifts to cfxpt. Cannot use >> operator here since we must control rounding.

See also: [`>>`](@ref)
"""
function >>(cfxpt :: CFixpointArray{N}, steps :: Integer) :: CFixpointArray{N} where {N}
    t_real = cfxpt.real >> steps;
    t_imag = cfxpt.imag >> steps;
    return CFixpointArray{N}(t_real, t_imag);
end

"""
Overload << for Fixpoint.
Apply 'steps' (>=0) left shifts to fxpt.
See also: [`<<`](@ref)
"""
function <<(fxpt :: Fixpoint, steps :: Integer) :: Fixpoint
    if (steps < 0)
        error("Integer value for steps must be greater than or equal to zero.");
    else    
        return Fixpoint(fxpt.data << steps, fxpt.scheme);
    end
end

"""
Overload << for FixpointArray.
Apply 'steps' (>=0) left shifts to fxpt.
See also: [`<<`](@ref)
"""
function <<(fxpt :: FixpointArray{N}, steps :: Integer) :: FixpointArray{N} where {N}
    if (steps < 0)
        error("Integer value for steps must be greater than or equal to zero.");
    else    
        return FixpointArray{N}(fxpt.data .<< steps, fxpt.scheme);
    end
end

"""
Overload << for CFixpoint.
Apply 'steps' (>=0) left shifts to fxpt.
See also: [`<<`](@ref)
"""
function <<(cfxpt :: CFixpoint, steps :: Integer) :: CFixpoint
    t_real = cfxpt.real << steps;
    t_imag = cfxpt.imag << steps;
    return CFixpoint(t_real,t_imag);
end

"""
Overload << for CFixpointArray.
Apply 'steps' (>=0) left shifts to fxpt.
See also: [`<<`](@ref)
"""
function <<(cfxpt :: CFixpointArray{N}, steps :: Integer) :: CFixpointArray{N} where {N}
    t_real = cfxpt.real << steps;
    t_imag = cfxpt.imag << steps;
    return CFixpointArray{N}(t_real,t_imag);
end

end # Fixpoint module