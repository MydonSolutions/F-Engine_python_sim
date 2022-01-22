import Base: sum, +, -, *, <<, >>, typemin, typemax, show, copy, getindex, setindex!,
        size, zeros, hcat, axes
using Printf

#########################################################################################
# FixpointArray Structures
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
Type which holds information about the FixpointArray value and influences its treatment
under arithmetic, logical and conversion operations.

See also: [`FixpointArray`](@ref)
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

"""
FixpointArray type that accepts an integer array or single integer accompanied by a FixpointScheme that
governs it's handling.

See also: [`FixpointScheme`](@ref)
"""
struct FixpointArray
    data :: Array{<:Integer}
    scheme :: FixpointScheme
    FixpointArray(fx_data::Array{<:Integer}, scheme::FixpointScheme) = new(fx_data,scheme);
    FixpointArray(fx_data::Integer,scheme::FixpointScheme) = new([fx_data],scheme);
end

Base.IndexStyle(::Type{<:FixpointArray}) = IndexLinear()

"""
CFixpointArray is the complex extension of FixpointArray that holds two FixpointArray types (real and imag) as its 
Real and Imaginary parts.

See also: [`FixpointArray`](@ref)
"""
struct CFixpointArray
    real :: FixpointArray
    imag :: FixpointArray
    function CFixpointArray(real :: Array{<:Integer}, imag :: Array{<:Integer}, scheme :: FixpointScheme)
        return new(FixpointArray(real,scheme),FixpointArray(imag,scheme));
    end
    function CFixpointArray(real :: FixpointArray, imag :: FixpointArray)
        if real.scheme != imag.scheme
            error("Real and Imag FixpointArray values must have the same scheme.");
        else
            return new(real,imag);
        end
    end
end

Base.IndexStyle(::Type{<:CFixpointArray}) = IndexLinear()
#########################################################################################
# Float parsing funtions
#########################################################################################
"""
```
fromFloat(fl_data::Real,scheme::FixpointScheme)
```
Converts a floating point value to FixpointArray according to the FixpointScheme presented.

See also: [`toFloat`](@ref)
"""
function fromFloat(fl_data::Real,scheme::FixpointScheme) :: FixpointArray
    return fromFloat([fl_data],scheme);
end
"""
```
fromFloat(fl_data :: Array{<:Real}, scheme :: FixpointScheme)
```
Converts a floating point array to FixpointArray according to the FixpointScheme presented.

See also: [`toFloat`](@ref)
"""
function fromFloat(fl_data::Array{<:Real}, scheme::FixpointScheme) :: FixpointArray
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
    return FixpointArray(data,scheme);
end

"""
```
toFloat(f :: FixpointArray) :: Array{Float64}
```
Converts a FixpointArray array to floating point according to the FixpointArray's FixpointScheme.

See also: [`fromFloat`](@ref)
"""
function toFloat(f :: FixpointArray) :: Array{Float64}
    fl_val = Float64.(f.data)./f.scheme.scale;
    return fl_val;
end

"""
Converts a Complex value to CFixpointArray according to the FixpointScheme presented.

See also: [`toComplex`](@ref)
"""
function fromComplex(c_data::Union{Array{<:Complex},Complex}, scheme::FixpointScheme) :: CFixpointArray
    return CFixpointArray(fromFloat(real(c_data),scheme), fromFloat(imag(c_data),scheme));
end

"""
Converts a two floating point values to CFixpointArray according to the FixpointScheme presented.

See also: [`toComplex`](@ref)
"""
function fromComplex(r_data::Union{Array{<:Real},Real}, i_data::Union{Array{<:Real},Real}, scheme::FixpointScheme) :: CFixpointArray
    return CFixpointArray(fromFloat(r_data,scheme),fromFloat(i_data,scheme));
end

"""
Converts a single floating point values to CFixpointArray according to the FixpointScheme presented.

See also: [`toComplex`](@ref)
"""
function fromComplex(r_data::Union{Array{<:Real},Real}, scheme::FixpointScheme) :: CFixpointArray
    return CFixpointArray(fromFloat(r_data,scheme),fromFloat(zeros(size(r_data)),scheme));
end


"""
Converts a CFixpointArray array to a complex point array according to the CFixpointArray's FixpointScheme.
"""
function toComplex(cfix :: CFixpointArray) :: Array{ComplexF64}
    return toFloat(cfix.real) + toFloat(cfix.imag)*im;
end

"""
Creates a FixpointArray populated with zero values. Creates a CFixpointArray if complex is set to true.
"""
function zeros(fx_scheme :: FixpointScheme, dims :: Tuple; complex :: Bool = false) :: Union{FixpointArray,CFixpointArray}
    if complex
        return fromComplex(zeros(Float64, dims),zeros(Float64, dims),fx_scheme);
    else
        return fromFloat(zeros(Float64, dims), fx_scheme); 
    end
end

"""
Fit all data within the min/max values for FixpointArray
"""
function normalise(f_val :: FixpointArray) :: FixpointArray
    return FixpointArray(clamp.(f_val.data,f_val.scheme.min,f_val.scheme.max),f_val.scheme);
end

"""
Fit all data within the min/max values for CFixpointArray
"""
function normalise(cf_val :: CFixpointArray) :: CFixpointArray
    return CFixpointArray(normalise(cf_val.real),normalise(cf_val.imag));
end

"""
Cast FixpointArray type to FixpointArray type of new scheme
"""
function cast(f_val :: FixpointArray, f_scheme :: FixpointScheme) :: FixpointArray
    return FixpointArray(f_val.data,f_scheme);
end

"""
Cast CFixpointArray type to CFixpointArray type of new scheme
"""
function cast(cf_val :: CFixpointArray, f_scheme :: FixpointScheme) :: CFixpointArray
    return CFixpointArray(FixpointArray(cf_val.real.data,f_scheme),FixpointArray(cf_val.imag.data,f_scheme));
end


#######################################################################################
# Arithmetic functions
#######################################################################################

"""
```
sum(f :: FixpointArray, dims :: Union{Integer,Colon})
```
Overload sum function to take FixpointArray type array as argument.

See also: [`sum`](@ref)
"""
function sum(f :: FixpointArray; dims :: Union{Integer,Colon}=:) :: FixpointArray
    sum_val = sum(f.data, dims=dims);
    bits = f.scheme.bits + ceil.(Integer,log2.(length(f.data)/length(sum_val)));
    scheme = FixpointScheme(bits, f.scheme.fraction, min_int=f.scheme.min,
    max_int=f.scheme.max, unsigned=f.scheme.unsigned, ovflw_behav=f.scheme.ovflw_behav,
    undflw_behav=f.scheme.undflw_behav);
    return FixpointArray(sum_val,scheme);
end

"""
```
sum(cf :: CFixpointArray, dims :: Union{Integer,Colon}=:)
```
Overload sum function to take CFixpointArray type array as argument.
See also: [`sum`](@ref)
"""
function sum(cf :: CFixpointArray; dims :: Union{Integer,Colon}=:) :: CFixpointArray
    r_sum_val = sum(cf.real,dims=dims);
    i_sum_val = sum(cf.imag,dims=dims);
    return CFixpointArray(r_sum_val,i_sum_val);
end

"""
```
*(a :: FixpointArray, b :: FixpointArray)
```
Overload * function to take FixpointArray type arrays as arguments.
        
See also: [`*`](@ref)
"""
function *(a :: FixpointArray, b :: FixpointArray) :: FixpointArray
    prod_val = a.data .* b.data;
    bits = a.scheme.bits + b.scheme.bits;
    fraction = a.scheme.fraction + b.scheme.fraction;
    unsigned = a.scheme.unsigned & b.scheme.unsigned;
    scheme = FixpointScheme(bits, fraction, unsigned=unsigned, 
    ovflw_behav=a.scheme.ovflw_behav, undflw_behav=a.scheme.undflw_behav);
    return FixpointArray(prod_val,scheme);
end

"""
```
*(a :: CFixpointArray, b :: CFixpointArray) :: CFixpointArray
```
Overload * function to take CFixpointArray type arrays as arguments.

See also: [`*`](@ref)
"""
function *(a :: CFixpointArray, b :: CFixpointArray) :: CFixpointArray
    function cmult(a, b, c, d)
        # Real part x = a*c - b*d
        x = (a*c)-(b*d);
        # Imaginary part y = a*d + b*c
        y = (a*d)+(b*c);
        return x, y;
    end
    out_real, out_imag = cmult(a.real, a.imag, b.real, b.imag);
    return CFixpointArray(out_real, out_imag); 
end

"""
```
+(a :: FixpointArray, b :: FixpointArray)
```
Overload + function to take FixpointArray type arrays as arguments.
        
See also: [`+`](@ref)
"""
function +(a :: FixpointArray, b :: FixpointArray) :: FixpointArray
    if (a.scheme.scale != b.scheme.scale)
        error("Addition performed between two FixpointArray values of differing scales.");
    end
    add_val = a.data .+ b.data;
    bits = max(a.scheme.bits,b.scheme.bits) + 1;
    unsigned = a.scheme.unsigned & b.scheme.unsigned;
    scheme = FixpointScheme(bits, a.scheme.fraction, unsigned=unsigned, 
    ovflw_behav=a.scheme.ovflw_behav, undflw_behav=a.scheme.undflw_behav);
    return FixpointArray(add_val,scheme);
end
"""
```
+(a :: CFixpointArray, b :: CFixpointArray)
```
Overload + function to take CFixpointArray type arrays as arguments.
        
See also: [`+`](@ref)
"""
function +(a :: CFixpointArray, b :: CFixpointArray) :: CFixpointArray
    r_sum = a.real + b.real;
    i_sum = a.imag + b.imag;
    return CFixpointArray(r_sum, i_sum);
end

"""
```
-(a :: FixpointArray, b :: FixpointArray)
```
Overload - function to take FixpointArray type arrays as arguments.
        
See also: [`-`](@ref)
"""
function -(a :: FixpointArray, b :: FixpointArray) :: FixpointArray
    if (a.scheme.scale != b.scheme.scale)
        error("Subtraction performed between two FixpointArray values of differing scales.");
    end 
    sub_val = a.data .- b.data;
    bits = max(a.scheme.bits,b.scheme.bits) + 1;
    unsigned = a.scheme.unsigned & b.scheme.unsigned;
    scheme = FixpointScheme(bits, a.scheme.fraction, unsigned=unsigned, 
    ovflw_behav=a.scheme.ovflw_behav, undflw_behav=a.scheme.undflw_behav);
    return FixpointArray(sub_val,scheme);
end

"""
```
-(a :: CFixpointArray, b :: CFixpointArray)
```
Overload - function to take CFixpointArray type arrays as arguments.
        
See also: [`-`](@ref)
"""
function -(a :: CFixpointArray, b :: CFixpointArray) :: CFixpointArray
    r_sub = a.real - b.real;
    i_sub = a.imag - b.imag;
    return CFixpointArray(r_sub, i_sub);
end

"""
```
power(f :: FixpointArray)
```
Returns power of the FixpointArray value given = f.data * f.data.
"""
function power(f :: FixpointArray) :: Array{Integer}
    return f.data .* f.data;
end

"""
```
power(cf :: CFixpointArray)
```
Returns power of the CFixpointArray value given = cf * conj(cf).
See also: [`conj`](@ref)
"""
function power(cf :: CFixpointArray) :: Array{Integer}
    res = copy(cf) * conj(cf);
    return res.real;
end

"""
```
conj(cf :: CFixpointArray)
```
Returns conjuage of the CFixpointArray value given.
"""
function conj(cf :: CFixpointArray) :: CFixpointArray
    i_res = copy(cf.imag);
    i_res.data = - copy(cf.imag.data);
    return CFixpointArray(cf.real, i_res);
end

"""
```
conj!(cf :: CFixpointArray)
```
Returns conjuage of the CFixpointArray value given.
! implies inline operation.
"""
function conj(cf :: CFixpointArray) :: CFixpointArray
    cf.imag.data = -cf.imag.data;
end
#######################################################################################
# Misc FixpointArray type handling functions
#######################################################################################

"""
```
clamp_wrap(f :: FixpointArray, min :: Integer, max :: Integer)
```
An overload of clamp_wrap to take a FixpointArray array argument instead of an Integer.
        
See also: [`clamp_wrap`](@ref)
"""
function clamp_wrap(f :: FixpointArray, min :: Integer, max :: Integer)
    clamp_val = ((f.data .- min) .% (min - max)) .+ min;
    scheme = FixpointScheme(f.scheme.bits,f.scheme.fraction,unsigned=f.scheme.unsigned,
    max_int=max, min_int=min,
    ovflw_behav=f.scheme.ovflw_behav, undflw_behav=f.scheme.undflw_behav);
    return FixpointArray(clamp_val,scheme);        
end

"""
```
clamp_wrap(f :: Integer, min :: Integer, max :: Integer)
```
Does a clamp operation but wraps the value to min/max rather than saturate 
the value like standard clamp.
        
See also: [`clamp`](@ref)
"""
function clamp_wrap(i :: Integer, min :: Integer, max :: Integer)
    return ((i - min) % (max - min)) + min;
end

"""
Requantise the data contained in fxpt according to the new scheme provided.
"""
function quantise(fxpt :: FixpointArray, scheme :: FixpointScheme) :: FixpointArray
    return fromFloat(toFloat(fxpt), scheme);
end

"""
Requantise the data contained in cfxpt according to the new scheme provided.
"""
function quantise(cfxpt :: CFixpointArray, scheme :: FixpointScheme) :: CFixpointArray
    return fromComplex(toComplex(cfxpt),scheme);
end

"""
Overload copy() function to copy FixpointArray by value as opposed to reference.

See also: [`copy`](@ref)
"""
function copy(f :: FixpointArray)
    tmpscheme = FixpointScheme(f.scheme.bits, f.scheme.fraction, min_int=f.scheme.min,
    max_int=f.scheme.max, unsigned=f.scheme.unsigned, ovflw_behav=f.scheme.ovflw_behav,
    undflw_behav=f.scheme.undflw_behav);
    return FixpointArray(copy(f.data),tmpscheme);
end

"""
```
function copy(cf :: CFixpointArray)
```
Overload copy() function to copy CFixpointArray by value as opposed to reference.
"""
function copy(cf :: CFixpointArray)
    return CFixpointArray(copy(cf.real),copy(cf.imag));
end

"""
```
function size(f :: FixpointArray)
```
Overload size() function to accept FixpointArray.
"""
function size(f::FixpointArray)::Integer
    return size(f.data);
end

"""
```
function size(cf :: CFixpointArray)
```
Overload size() function to accept CFixpointArray.
"""
function size(cf::CFixpointArray)::Integer
    return size(cf.real);
end


"""
```
show(io :: IO, f :: FixpointArray)
```
Overload show function for printing out FixpointArray summary.
See also: [`show`](@ref)
"""
function show(io::IO, f :: FixpointArray)
    @printf(io,"FixpointArray real %s (%d, %d), shape %s", f.scheme.unsigned ? "unsigned" : "signed",f.scheme.bits, f.scheme.fraction, size(f.data));
end

"""
```
show(io :: IO, f :: FixpointArray)
```
Overload show function for printing out FixpointArray summary.
See also: [`show`](@ref)
"""
function show(io::IO, cf :: CFixpointArray)
    @printf(io,"CFixpointArray complex %s (%d, %d), shape %s", cf.real.scheme.unsigned ? "unsigned" : "signed",cf.real.scheme.bits, cf.real.scheme.fraction, size(cf.real.data))
end

"""
```
getindex(f :: FixpointArray, i :: Int)
```
Overload getindex function for accessing data elements out FixpointArray type.
"""
function getindex(f :: FixpointArray, i :: Int) :: FixpointArray
    return FixpointArray(f.data[i],f.scheme);
end

"""
```
getindex(cf :: CFixpointArray, i :: Int)
```
Overload getindex function for accessing data elements out CFixpointArray type.
"""
function getindex(cf :: CFixpointArray, i :: Int) :: CFixpointArray
    return CFixpointArray(cf.real[i],cf.imag[i]);
end

"""
```
getindex(f :: FixpointArray, i :: UnitRange{Int64})
```
Overload getindex function for accessing data elements out FixpointArray type.
"""
function getindex(f :: FixpointArray, i :: Vararg{UnitRange{Int64},N}) :: FixpointArray where {N}
    return FixpointArray(f.data[i...],f.scheme);
end

"""
```
getindex(cf :: CFixpointArray, i :: UnitRange{Int64})
```
Overload getindex function for accessing data elements out CFixpointArray type.
"""
function getindex(cf :: CFixpointArray, i :: Vararg{UnitRange{Int64},N}) :: CFixpointArray where {N}
    return CFixpointArray(cf.real[i...],cf.imag[i...]);
end

"""
```
getindex(f :: FixpointArray, i :: Vector{Int})
```
Overload getindex function for accessing data elements out FixpointArray type.
"""
function getindex(f :: FixpointArray, i :: Vector{Int}) :: FixpointArray
    return FixpointArray(f.data[i],f.scheme);
end

"""
```
getindex(cf :: CFixpointArray, i ::Vector{Int})
```
Overload getindex function for accessing data elements out CFixpointArray type.
"""
function getindex(cf :: CFixpointArray, i :: Vector{Int}) :: CFixpointArray
    return CFixpointArray(cf.real[i],cf.imag[i]);
end

"""
```
setindex!(f :: FixpointArray, i ::Vector{Int})
```
Overload setindex function for accessing data elements out FixpointArray type.
"""
function setindex!(f :: FixpointArray, val :: FixpointArray, i :: Vector{Int}) :: Nothing
    f.data[i] = val.data;
end

"""
```
setindex!(cf :: CFixpointArray, i ::Vector{Int})
```
Overload setindex function for accessing data elements out CFixpointArray type.
"""
function setindex!(cf :: CFixpointArray, val :: CFixpointArray, i :: Vector{Int}) :: Nothing
    cf.real[i] = val.real;
    cf.imag[i] = val.imag;
end

"""
```
setindex!(f :: FixpointArray, i ::UnitRange{Int})
```
Overload setindex function for accessing data elements out FixpointArray type.
"""
function setindex!(f :: FixpointArray, val :: FixpointArray, i :: UnitRange{Int}) :: Nothing
    f.data[i] = val.data;
    return;
end

"""
```
setindex!(cf :: CFixpointArray, i ::UnitRange{Int})
```
Overload setindex function for accessing data elements out CFixpointArray type.
"""
function setindex!(cf :: CFixpointArray, val :: CFixpointArray, i :: UnitRange{Int}) :: Nothing
    cf.real[i] = val.real;
    cf.imag[i] = val.imag;
    return;
end

"""
````
setindex!(f :: FixpointArray, i :: Int)
```
Overload setindex! function for accessing data elements out FixpointArray type.
"""
function setindex!(f :: FixpointArray, val :: FixpointArray, i :: Int) :: Nothing
    f.data[i] = val.data[i];
end

"""
```
setindex!(cf :: CFixpointArray, i :: Int)
```
Overload setindex! function for accessing data elements out CFixpointArray type.
"""
function setindex!(cf :: CFixpointArray, val :: CFixpointArray, i :: Int) :: CFixpointArray
    cf.real[i] = val.real;
    cf.imag[i] = val.imag;
end

"""
Overload axes function to handle FixpointArray types.
"""
function axes(f :: FixpointArray, i :: Int64) :: AbstractUnitRange
    return axes(f.data,i);
end

"""
Overload axes function to handle CFixpointArray types.
"""
function axes(cf :: CFixpointArray, i :: Int64) :: AbstractUnitRange
    return axes(cf.real.data,i);
end

"""
Overload hcat function to handle horizontal concatenation of FixpointArray types.
Requires that schemes match.
"""
function hcat(f_1 :: FixpointArray, f_2 :: FixpointArray) :: FixpointArray 
    #Check schemes match:
    if f_1.scheme == f_2.scheme
        return FixpointArray(hcat(f_1.data,f_2.data),f_1.scheme);
    else
        error("FixpointArray args don't share the same scheme.");
    end
end

"""
Overload hcat function to handle horizontal concatenation of CFixpointArray types.
Requires that schemes match.
"""
function hcat(cf_1 :: CFixpointArray, cf_2 :: CFixpointArray) :: CFixpointArray 
    #Check real schemes match - imag will match:
    if cf_1.real.scheme == cf_2.real.scheme
        return CFixpointArray(hcat(cf_1.real,cf_2.real), hcat(cf_1.imag,cf_2.imag));
    else
        error("CFixpointArray args don't share the same scheme.");
    end
end

#######################################################################################
# Logical operator functions
#######################################################################################
"""
```
>>(fxpt :: FixpointArray, steps :: Integer)
```
Overload >> function for FixpointArray args.
Apply 'steps' (>=0) right shifts to fxpt. Cannot use >> operator here since we must control rounding.

See also: [`>>`](@ref)
"""
function >>(fxpt :: FixpointArray, steps :: Integer) :: FixpointArray
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
        return FixpointArray(round.(Integer, fxpt.data/(2^steps),rnd_behav),fxpt.scheme);
    end
end

"""
```
>>(cfxpt :: CFixpointArray, steps :: Integer)
```
Overload >> function for CFixpointArray args.
Apply 'steps' (>=0) right shifts to cfxpt. Cannot use >> operator here since we must control rounding.

See also: [`>>`](@ref)
"""
function >>(cfxpt :: CFixpointArray, steps :: Integer) :: CFixpointArray
    t_real = cfxpt.real >> steps;
    t_imag = cfxpt.imag >> steps;
    return CFixpointArray(t_real, t_imag);
end

"""
```
<<(fxpt :: FixpointArray, steps :: Integer)
```
Overload << function for FixpointArray args.
Apply 'steps' (>=0) left shifts to fxpt.
See also: [`<<`](@ref)
"""
function <<(fxpt :: FixpointArray, steps :: Integer) :: FixpointArray
    t_fxpt = copy(fxpt);
    if (steps < 0)
        error("Integer value for steps must be greater than or equal to zero.");
    else    
        t_fxpt.data .<<= steps;
    end
    return t_fxpt;
end

"""
```
<<(fxpt :: FixpointArray, steps :: Integer)
```
Overload << function for FixpointArray args.
Apply 'steps' (>=0) left shifts to fxpt.
See also: [`<<`](@ref)
"""
function <<(cfxpt :: CFixpointArray, steps :: Integer) :: CFixpointArray
    t_real = cfxpt.real << steps;
    t_imag = cfxpt.imag << steps;
    return CFixpointArray(t_real,t_imag);
end