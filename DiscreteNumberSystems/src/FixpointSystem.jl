module FixpointSystem

using Printf

import Base: +, -, *, <<, >>

export FixpointScheme, Fixpoint, Fixpoint, CFixpoint, CFixpoint

export fromFloat, toFloat, fromComplex, toComplex, zeros, normalise, cast, sum, *, +, -, power, conj, clamp_wrap, quantise, copy, size, length, show, getindex, setindex!, hcat, lastindex, >>, <<

#########################################################################################
# FixpointScheme Structure
#########################################################################################

"""
```
struct FixpointScheme
    bits :: Integer
    fraction :: Integer
    min :: Integer
    max :: Integer
    unsigned :: Bool
    range :: Integer
    scale :: Integer
    ovflw_behav :: String
    undflw_behav :: String
```
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
    undflw_behav :: String
    function FixpointScheme(bits::Integer, fraction::Integer;
         min_int::Union{Integer, Nothing}=nothing,max_int::Union{Integer, Nothing}=nothing,
         unsigned::Union{Bool, Nothing}=false, ovflw_behav::Union{String, Nothing}="WRAP",
         undflw_behav::Union{String, Nothing}="ROUND_EVEN")
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
#########################################################################################
# """
# ```
# struct Fixpoint
#     data :: Integer
#     scheme :: FixpointScheme
# ```
# Fixpoint type that accepts a single integer accompanied by a FixpointScheme that governs it's handling.

# See also: [`FixpointScheme`](@ref)
# """
# struct Fixpoint
#     data :: Integer
#     scheme :: FixpointScheme
# end

"""
Fixpoint type that accepts an integer array or single integer accompanied by a FixpointScheme that
governs it's handling.

See also: [`FixpointScheme`](@ref)
"""
struct Fixpoint
    data :: Array{<:Integer}
    scheme :: FixpointScheme
    Fixpoint(fx_data::Array{<:Integer}, scheme::FixpointScheme) = new(fx_data,scheme);
    Fixpoint(fx_data::Integer,scheme::FixpointScheme) = new([fx_data],scheme);
end

"""
Set to linear indexing to avoid needing to set new getindex/setindex methods for large dim Fixpoint arrays
"""
Base.IndexStyle(::Type{<:Fixpoint}) = IndexLinear()

"""
CFixpoint is the complex extension of Fixpoint that holds two Fixpoint types (real and imag) as its 
Real and Imaginary parts.

See also: [`Fixpoint`](@ref)
"""
# struct CFixpoint
#     real :: Fixpoint
#     imag :: Fixpoint
#     function CFixpoint(real :: Fixpoint, imag :: Fixpoint)
#         if real.scheme != imag.scheme
#             error("Real and Imag Fixpoint values must have the same scheme.");
#         else
#             return new(real,imag);
#         end
#     end
# end

"""
CFixpoint is the complex extension of Fixpoint that holds two Fixpoint types (real and imag) as its 
Real and Imaginary parts.

See also: [`Fixpoint`](@ref)
"""
struct CFixpoint
    real :: Fixpoint
    imag :: Fixpoint
    function CFixpoint(real :: Array{<:Integer}, imag :: Array{<:Integer}, scheme :: FixpointScheme)
        return new(Fixpoint(real,scheme),Fixpoint(imag,scheme));
    end
    function CFixpoint(real :: Fixpoint, imag :: Fixpoint)
        if real.scheme != imag.scheme
            error("Real and Imag Fixpoint values must have the same scheme.");
        else
            return new(real,imag);
        end
    end
end

"""
Set to linear indexing to avoid needing to set new getindex/setindex methods for large dim CFixpoint arrays
"""
Base.IndexStyle(::Type{<:CFixpoint}) = IndexLinear()

#########################################################################################
# Float parsing funtions
#########################################################################################
"""
```
fromFloat(fl_data::Real,scheme::FixpointScheme)
```
Converts a floating point value to Fixpoint according to the FixpointScheme presented.

See also: [`toFloat`](@ref)
"""
function fromFloat(fl_data::Real,scheme::FixpointScheme) :: Fixpoint
    return fromFloat([fl_data],scheme);
end
"""
```
fromFloat(fl_data :: Array{<:Real}, scheme :: FixpointScheme)
```
Converts a floating point array to Fixpoint according to the FixpointScheme presented.

See also: [`toFloat`](@ref)
"""
function fromFloat(fl_data::Array{<:Real}, scheme::FixpointScheme) :: Fixpoint
    prod = fl_data .* scheme.scale;
    rnd_behav = RoundNearest;
    if (scheme.undflw_behav == "ROUND_EVEN")
        rnd_behav = RoundNearest;
    elseif (scheme.undflw_behav =="ROUND_AWAY")
        rnd_behav = RoundNearestTiesAway;
    elseif (scheme.undflw_behav =="TRUNCATE")
        rnd_behav = RoundToZero;
    else
        error("No recognisable rounding method specified");
    end
    if (scheme.ovflw_behav == "SATURATE")
        data = clamp.(round.(Integer,prod, rnd_behav),scheme.min,scheme.max);
    elseif (scheme.ovflw_behav == "WRAP")
        data = clamp_wrap.(round.(Integer,prod, rnd_behav),scheme.min,scheme.max);
    else 
        error("No recognisable overflow method specified");
    end
    return Fixpoint(data,scheme);
end

"""
```
toFloat(f :: Fixpoint) :: Array{Float64}
```
Converts a Fixpoint array to floating point according to the Fixpoint's FixpointScheme.

See also: [`fromFloat`](@ref)
"""
function toFloat(f :: Fixpoint) :: Array{Float64}
    fl_val = Float64.(f.data)./f.scheme.scale;
    return fl_val;
end

"""
Converts a Complex value to CFixpoint according to the FixpointScheme presented.

See also: [`toComplex`](@ref)
"""
function fromComplex(c_data::Union{Array{<:Complex},Complex}, scheme::FixpointScheme) :: CFixpoint
    return CFixpoint(fromFloat(real(c_data),scheme), fromFloat(imag(c_data),scheme));
end

"""
Converts a two floating point values to CFixpoint according to the FixpointScheme presented.

See also: [`toComplex`](@ref)
"""
function fromComplex(r_data::Union{Array{<:Real},Real}, i_data::Union{Array{<:Real},Real}, scheme::FixpointScheme) :: CFixpoint
    return CFixpoint(fromFloat(r_data,scheme),fromFloat(i_data,scheme));
end

"""
Converts a single floating point values to CFixpoint according to the FixpointScheme presented.

See also: [`toComplex`](@ref)
"""
function fromComplex(r_data::Union{Array{<:Real},Real}, scheme::FixpointScheme) :: CFixpoint
    return CFixpoint(fromFloat(r_data,scheme),fromFloat(zeros(size(r_data)),scheme));
end


"""
Converts a CFixpoint array to a complex point array according to the CFixpoint's FixpointScheme.
"""
function toComplex(cfix :: CFixpoint) :: Array{ComplexF64}
    return toFloat(cfix.real) + toFloat(cfix.imag)*im;
end

"""
Creates a Fixpoint populated with zero values. Creates a CFixpoint if complex is set to true.
"""
function Base.zeros(fx_scheme :: FixpointScheme, dims :: Tuple; complex :: Bool = false) :: Union{Fixpoint,CFixpoint}
    if complex
        return fromComplex(zeros(Float64, dims),zeros(Float64, dims),fx_scheme);
    else
        return fromFloat(zeros(Float64, dims), fx_scheme); 
    end
end

"""
Fit all data within the min/max values for Fixpoint
"""
function normalise(f_val :: Fixpoint) :: Fixpoint
    return Fixpoint(clamp.(f_val.data,f_val.scheme.min,f_val.scheme.max),f_val.scheme);
end

"""
Fit all data within the min/max values for CFixpoint
"""
function normalise(cf_val :: CFixpoint) :: CFixpoint
    return CFixpoint(normalise(cf_val.real),normalise(cf_val.imag));
end

"""
Cast Fixpoint type to Fixpoint type of new scheme
"""
function cast(f_val :: Fixpoint, f_scheme :: FixpointScheme) :: Fixpoint
    return Fixpoint(f_val.data,f_scheme);
end

"""
Cast CFixpoint type to CFixpoint type of new scheme
"""
function cast(cf_val :: CFixpoint, f_scheme :: FixpointScheme) :: CFixpoint
    return CFixpoint(Fixpoint(cf_val.real.data,f_scheme),Fixpoint(cf_val.imag.data,f_scheme));
end


#######################################################################################
# Arithmetic functions
#######################################################################################

"""
Overload sum function to take Fixpoint type array as argument.

See also: [`sum`](@ref)
"""
function Base.sum(f :: Fixpoint; dims :: Union{Integer,Colon}=:) :: Fixpoint
    sum_val = sum(f.data, dims=dims);
    bits = f.scheme.bits + ceil.(Integer,log2.(length(f.data)/length(sum_val)));
    scheme = FixpointScheme(bits, f.scheme.fraction, min_int=f.scheme.min,
    max_int=f.scheme.max, unsigned=f.scheme.unsigned, ovflw_behav=f.scheme.ovflw_behav,
    undflw_behav=f.scheme.undflw_behav);
    return Fixpoint(sum_val,scheme);
end

"""
Overload sum function to take CFixpoint type array as argument.
See also: [`sum`](@ref)
"""
function Base.sum(cf :: CFixpoint; dims :: Union{Integer,Colon}=:) :: CFixpoint
    r_sum_val = sum(cf.real,dims=dims);
    i_sum_val = sum(cf.imag,dims=dims);
    return CFixpoint(r_sum_val,i_sum_val);
end

"""
Overload * function to take Fixpoint type arrays as arguments.
        
See also: [`*`](@ref)
"""
function *(a :: Fixpoint, b :: Fixpoint) :: Fixpoint
    prod_val = a.data .* b.data;
    bits = a.scheme.bits + b.scheme.bits;
    fraction = a.scheme.fraction + b.scheme.fraction;
    unsigned = a.scheme.unsigned & b.scheme.unsigned;
    scheme = FixpointScheme(bits, fraction, unsigned=unsigned, 
    ovflw_behav=a.scheme.ovflw_behav, undflw_behav=a.scheme.undflw_behav);
    return Fixpoint(prod_val,scheme);
end

"""
Overload * function to take CFixpoint type arrays as arguments.

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
Overload + function to take Fixpoint type arrays as arguments.
        
See also: [`+`](@ref)
"""
function +(a :: Fixpoint, b :: Fixpoint) :: Fixpoint
    if (a.scheme.scale != b.scheme.scale)
        error("Addition performed between two Fixpoint values of differing scales.");
    end
    add_val = a.data .+ b.data;
    bits = max(a.scheme.bits,b.scheme.bits) + 1;
    unsigned = a.scheme.unsigned & b.scheme.unsigned;
    scheme = FixpointScheme(bits, a.scheme.fraction, unsigned=unsigned, 
    ovflw_behav=a.scheme.ovflw_behav, undflw_behav=a.scheme.undflw_behav);
    return Fixpoint(add_val,scheme);
end

"""
Overload + function to take CFixpoint type arrays as arguments.
        
See also: [`+`](@ref)
"""
function +(a :: CFixpoint, b :: CFixpoint) :: CFixpoint
    r_sum = a.real + b.real;
    i_sum = a.imag + b.imag;
    return CFixpoint(r_sum, i_sum);
end

"""
Overload - function to take Fixpoint type arrays as arguments.
        
See also: [`-`](@ref)
"""
function -(a :: Fixpoint, b :: Fixpoint) :: Fixpoint
    if (a.scheme.scale != b.scheme.scale)
        error("Subtraction performed between two Fixpoint values of differing scales.");
    end 
    sub_val = a.data .- b.data;
    bits = max(a.scheme.bits,b.scheme.bits) + 1;
    unsigned = a.scheme.unsigned & b.scheme.unsigned;
    scheme = FixpointScheme(bits, a.scheme.fraction, unsigned=unsigned, 
    ovflw_behav=a.scheme.ovflw_behav, undflw_behav=a.scheme.undflw_behav);
    return Fixpoint(sub_val,scheme);
end

"""
Overload - function to take CFixpoint type arrays as arguments.
        
See also: [`-`](@ref)
"""
function -(a :: CFixpoint, b :: CFixpoint) :: CFixpoint
    r_sub = a.real - b.real;
    i_sub = a.imag - b.imag;
    return CFixpoint(r_sub, i_sub);
end

"""
Returns power of the Fixpoint value given = f.data * f.data.
"""
function power(f :: Fixpoint) :: Array{Integer}
    return f.data .* f.data;
end

"""
Returns power of the CFixpoint value given = cf * conj(cf).
See also: [`conj`](@ref)
"""
function power(cf :: CFixpoint) :: Array{Integer}
    res = copy(cf) * conj(cf);
    return res.real;
end

"""
Returns conjuage of the CFixpoint value given.
"""
function Base.conj(cf :: CFixpoint) :: CFixpoint
    i_res = copy(cf.imag);
    i_res.data = - copy(cf.imag.data);
    return CFixpoint(cf.real, i_res);
end

"""
Returns conjuage of the CFixpoint value given.
! implies inline operation.
"""
function Base.conj!(cf :: CFixpoint) :: CFixpoint
    cf.imag.data = -cf.imag.data;
end
#######################################################################################
# Misc Fixpoint type handling functions
#######################################################################################

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
function quantise(fxpt :: Fixpoint, scheme :: FixpointScheme) :: Fixpoint
    return fromFloat(toFloat(fxpt), scheme);
end

"""
Requantise the data contained in cfxpt according to the new scheme provided.
"""
function quantise(cfxpt :: CFixpoint, scheme :: FixpointScheme) :: CFixpoint
    return fromComplex(toComplex(cfxpt),scheme);
end

"""
Overload copy() function to copy Fixpoint by value as opposed to reference.

See also: [`copy`](@ref)
"""
function Base.copy(f :: Fixpoint) :: Fixpoint
    tmpscheme = FixpointScheme(f.scheme.bits, f.scheme.fraction, min_int=f.scheme.min,
    max_int=f.scheme.max, unsigned=f.scheme.unsigned, ovflw_behav=f.scheme.ovflw_behav,
    undflw_behav=f.scheme.undflw_behav);
    return Fixpoint(copy(f.data),tmpscheme);
end

"""
Overload copy() function to copy CFixpoint by value as opposed to reference.
"""
function Base.copy(cf :: CFixpoint) :: CFixpoint
    return CFixpoint(copy(cf.real),copy(cf.imag));
end

"""
Overload size() function to accept Fixpoint.
"""
function Base.size(f::Fixpoint)
    return size(f.data);
end

"""
Overload size() function to accept CFixpoint.
"""
function Base.size(cf::CFixpoint)
    return size(cf.real);
end

"""
Overload length() function to accept Fixpoint
"""
function Base.length(f :: Fixpoint)
    return prod(size(f))
end

"""
Overload length() function to accept CFixpoint
"""
function Base.length(cf :: CFixpoint)
    return prod(size(cf))
end

"""
Overload show function for printing out Fixpoint summary.
See also: [`show`](@ref)
"""
function Base.show(io::IO, f :: Fixpoint)
    @printf(io,"Fixpoint real %s (%d, %d), shape %s", f.scheme.unsigned ? "unsigned" : "signed",f.scheme.bits, f.scheme.fraction, size(f.data));
end

"""
Overload show function for printing out Fixpoint summary.
See also: [`show`](@ref)
"""
function Base.show(io::IO, cf :: CFixpoint)
    @printf(io,"CFixpoint complex %s (%d, %d), shape %s", cf.real.scheme.unsigned ? "unsigned" : "signed",cf.real.scheme.bits, cf.real.scheme.fraction, size(cf.real.data))
end

#######################################################################################
# Dimension mangling Fixpoint type handling functions
# Note that indexing functions depend on LinearIndexing see:
# https://docs.julialang.org/en/v1/manual/interfaces/
#######################################################################################

"""
Overload getindex function for accessing data elements out Fixpoint type.
"""
function Base.getindex(f :: Fixpoint, i :: Int)
    return Fixpoint(f.data[i],f.scheme);
end

"""
Overload getindex function for accessing data elements out CFixpoint type.
"""
function Base.getindex(cf :: CFixpoint, i :: Int)
    return CFixpoint(cf.real[i],cf.imag[i]);
end

"""
Overload getindex function for accessing data elements out Fixpoint type. This 
overload provides a fall back to the one above in the instance that the Fixpoint is
indexed with a multidimensional set of indices
"""
function Base.getindex(f :: Fixpoint, I...)
    return f[I]
end

"""
Overload getindex function for accessing data elements out CFixpoint type. This 
overload provides a fall back to the one above in the instance that the Fixpoint is
indexed with a multidimensional set of indices
"""
function Base.getindex(cf :: CFixpoint, I...)
    return cf[I]
end

"""
Overload setindex function for setting data elements out Fixpoint type.
"""
function Base.setindex!(f :: Fixpoint, v :: Fixpoint , i :: Int) :: Nothing
    f.data[i] = v.data;
end

"""
Overload setindex function for setting data elements out CFixpoint type.
"""
function Base.setindex!(cf :: CFixpoint, v :: CFixpoint, i :: Int) :: Nothing
    cf.real[i] = v.real;
    cf.imag[i] = v.imag;
end

"""
Overload setindex function for setting data elements out Fixpoint type.
Falls back to earlier setindex! function in the event of multidimensional
indexing.
"""
function Base.setindex!(f :: Fixpoint, v :: Fixpoint , I...) :: Nothing
    f.data[I] = v.data;
end

"""
Overload setindex function for setting data elements out CFixpoint type.
Falls back to earlier setindex! function in the event of multidimensional
indexing.
"""
function Base.setindex!(cf :: CFixpoint, v :: CFixpoint, I...) :: Nothing
    cf.real[I] = v.real;
    cf.imag[I] = v.imag;
end

"""
Overload hcat function to handle horizontal concatenation of Fixpoint types.
Requires that schemes match.
"""
function Base.hcat(f_1 :: Fixpoint, f_2 :: Fixpoint) :: Fixpoint 
    #Check schemes match:
    if f_1.scheme == f_2.scheme
        return Fixpoint(hcat(f_1.data,f_2.data),f_1.scheme);
    else
        error("Fixpoint args don't share the same scheme.");
    end
end

"""
Overload hcat function to handle horizontal concatenation of CFixpoint types.
Requires that schemes match.
"""
function Base.hcat(cf_1 :: CFixpoint, cf_2 :: CFixpoint) :: CFixpoint 
    #Check real schemes match - imag will match:
    if cf_1.real.scheme == cf_2.real.scheme
        return CFixpoint(hcat(cf_1.real,cf_2.real), hcat(cf_1.imag,cf_2.imag));
    else
        error("CFixpoint args don't share the same scheme.");
    end
end

"""
Overload lastindex function to handle slicing of Fixpoint with end
"""
function Base.lastindex(f :: Fixpoint) :: Int
    return length(f.data)
end

"""
Overload lastindex function to handle slicing of CFixpoint with end
"""
function Base.lastindex(cf :: CFixpoint) :: Int
    return length(cf.real.data)
end

"""
Overload axes function for Fixpoint
"""
function Base.axes(f :: Fixpoint) :: AbstractUnitRange{<:Integer}
    return map(OneTo, size(f))
end

"""
Overload axes function for CFixpoint
"""
function Base.axes(cf :: Fixpoint) :: AbstractUnitRange{<:Integer}
    return map(OneTo, size(cf))
end

# """
# Overload similar function for Fixpoint
# """
# function Base.similar(f :: Fixpoint) 
#     return 

#######################################################################################
# Logical operator functions
#######################################################################################
"""
Overload >> function for Fixpoint args.
Apply 'steps' (>=0) right shifts to fxpt. Cannot use >> operator here since we must control rounding.

See also: [`>>`](@ref)
"""
function >>(fxpt :: Fixpoint, steps :: Integer) :: Fixpoint
    if (steps < 0)
        error("Integer value for steps must be greater than or equal to zero.");
    else
        if (fxpt.scheme.undflw_behav == "ROUND_EVEN")
            rnd_behav = RoundNearest;
        elseif (fxpt.scheme.undflw_behav =="ROUND_AWAY")
            rnd_behav = RoundNearestTiesAway;
        elseif (fxpt.scheme.undflw_behav =="TRUNCATE")
            rnd_behav = RoundToZero;
        else
            error("No recognisable rounding method specified");
        end
        return Fixpoint(round.(Integer, fxpt.data/(2^steps),rnd_behav),fxpt.scheme);
    end
end

"""
Overload >> function for CFixpoint args.
Apply 'steps' (>=0) right shifts to cfxpt. Cannot use >> operator here since we must control rounding.

See also: [`>>`](@ref)
"""
function >>(cfxpt :: CFixpoint, steps :: Integer) :: CFixpoint
    t_real = cfxpt.real >> steps;
    t_imag = cfxpt.imag >> steps;
    return CFixpoint(t_real, t_imag);
end

"""
Overload << function for Fixpoint args.
Apply 'steps' (>=0) left shifts to fxpt.
See also: [`<<`](@ref)
"""
function <<(fxpt :: Fixpoint, steps :: Integer) :: Fixpoint
    t_fxpt = copy(fxpt);
    if (steps < 0)
        error("Integer value for steps must be greater than or equal to zero.");
    else    
        t_fxpt.data .<<= steps;
    end
    return t_fxpt;
end

"""
Overload << function for Fixpoint args.
Apply 'steps' (>=0) left shifts to fxpt.
See also: [`<<`](@ref)
"""
function <<(cfxpt :: CFixpoint, steps :: Integer) :: CFixpoint
    t_real = cfxpt.real << steps;
    t_imag = cfxpt.imag << steps;
    return CFixpoint(t_real,t_imag);
end

end # Fixpoint module