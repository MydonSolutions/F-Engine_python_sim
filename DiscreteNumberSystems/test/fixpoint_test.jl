using Test
using Base: abs
using Pkg
Pkg.activate("./DiscreteNumberSystems")

using DiscreteNumberSystems.FixpointSystem
"""
FixpointArray testing. 
"""
#Create our scheme:
f_scheme_1 = FixpointScheme(18,17);
f_scheme_2 = FixpointScheme(13,12);

#Create our floating point values:
scalar = 0.14159265
v1 = [0.15566998822, 0.0001, 0.45];
v2 = [0.333222111991, 0.3344888, 0.779]
val1 = [0.15566998822 0.0001 0.45; 0.119 0.55 0.37711];
val2 = [0.333222111991 0.3344888 0.779; 0.123 0.6622 0.4621];

#Cast them to fixed point using the above scheme
f_scalar = Fixpoint(scalar, f_scheme_1)
f_val1 = FixpointArray{ndims(val1)}(val1,f_scheme_1);
f_val2 = FixpointArray{ndims(val2)}(val2,f_scheme_1);

#Test Fixpoint conversion
@test isapprox(float(f_scalar), scalar, atol=0.0001)
@test isapprox(float(f_val1), val1, atol=0.0001)
@test isapprox(float(f_val2), val2, atol=0.0001)

#Test addition
@test isapprox(float(f_scalar + f_scalar), scalar + scalar, atol=0.0001)
@test isapprox(float(f_val1 + f_val2), val1 .+ val2, atol=0.0001)

#Test multiplication
@test isapprox(float(f_scalar * f_scalar), scalar * scalar, atol=0.0001)
@test isapprox(float(f_val1 * f_val2), val1 .* val2, atol=0.0001)

#Test summing the vectors
@test isapprox(float(sum(f_val1,dims=1)), sum(val1,dims=1), atol=0.0001)

#Test subtracting the vectors
@test isapprox(float(f_scalar - f_scalar), scalar - scalar, atol=0.0001)
@test isapprox(float(f_val1 - f_val2), val1 .- val2, atol=0.0001)

#Test right shifting the vector values
ideal_rshift_scalar = f_scalar.data .>> 1;
f_rshift_scalar = f_scalar >> 1;
@test isapprox(f_rshift_scalar.data, ideal_rshift_scalar, atol=1.0)
ideal_rshift = f_val1.data .>> 1;
f_rshift = f_val1 >> 1;
@test isapprox(f_rshift.data, ideal_rshift, atol=0.0001)

#Test left shifting the vector values
ideal_lshift_scalar = f_scalar.data .>> 1;
f_lshift_scalar = f_scalar >> 1;
@test isapprox(f_lshift_scalar.data, ideal_lshift_scalar, atol=1.0)
ideal_lshift = f_val1.data .<< 1;
f_lshift = f_val1 << 1;
@test isapprox(f_lshift.data, ideal_lshift, atol=0.0001)

#Test hcat of Fixpoint values
ideal_hcat = hcat(val1,val2)
f_hcat = hcat(f_val1,f_val2)
@test isapprox(float(f_hcat), ideal_hcat, atol=0.0001)

#Test quantise of Fixpoint values
@test isapprox(float(quantise(f_val1,f_scheme_2)), val1, atol=0.001)

#Test slicing
sliced_f_val1 = f_val1[2:end-1]
@test isa(sliced_f_val1, FixpointArray)
@test isapprox(float(sliced_f_val1), val1[2:end-1], atol=0.0001)
sliced_f_val1[:] = f_val1[1:end-2].data
@test isapprox(float(sliced_f_val1), val1[1:end-2], atol=0.0001)
sliced_f_val1[:] = f_val1[2:end-1]
@test isapprox(float(sliced_f_val1), val1[2:end-1], atol=0.0001)
sliced_f_val1[:] = val2[1:end-2]
@test isapprox(float(sliced_f_val1), val2[1:end-2], atol=0.0001)

@test isa(sliced_f_val1[1], Fixpoint)
sliced_f_val1[1] = f_scalar.data
@test isapprox(float(sliced_f_val1[1]), scalar, atol=0.0001)
sliced_f_val1[1] = f_scalar
@test isapprox(float(sliced_f_val1[1]), scalar, atol=0.0001)
sliced_f_val1[1] = scalar
@test isapprox(float(sliced_f_val1[1]), scalar, atol=0.0001)
 
"""
CFixpointArray testing
"""
scalar_complex = scalar + 2im*scalar
cf_scalar = CFixpoint(real(scalar_complex), imag(scalar_complex), f_scheme_1)
cf_val1 = CFixpointArray{ndims(f_val1)}(f_val1.data, f_val2.data, f_scheme_1);
cf_val2 = CFixpointArray{ndims(f_val1)}(f_val2, f_val1);
c_val1 = val1 + val2*im;
c_val2 = val2 + val1*im;

#Test conversion complex -> CFixpoint works
@test isapprox(float(cf_val2), c_val2, atol=0.0001);
@test isapprox(float(cf_val1), c_val1, atol=0.0001);

#Test addition of CFixpoint types
c_add = cf_val1 + cf_val2;
ideal_cadd = c_val1 .+ c_val2;
@test isapprox(float(c_add), ideal_cadd, atol=0.0001);

# #Test multiplication of CFixpoint types
c_mul = cf_val1 * cf_val2;
ideal_cmul = c_val1 .* c_val2;
@test isapprox(float(c_mul), ideal_cmul, atol= 0.0001)

#Test subtraction of CFixpoint types
c_min = cf_val1 - cf_val2;
ideal_cmin = c_val1 .- c_val2;
@test isapprox(float(c_min), ideal_cmin, atol=0.0001)

#Test summation of CFixpoint types
c_sum = sum(cf_val1,dims=1);
ideal_csum = sum(c_val1,dims=1);
@test isapprox(float(c_sum), ideal_csum, atol = 0.0001)

#Test rightshift of CFixpoint types
c_rshift = cf_val1 >> 1;
ideal_rshift = (val1 + 1im .* val2)/2;
@test isapprox(float(c_rshift.real), real(ideal_rshift), atol = 0.0001)
@test isapprox(float(c_rshift.imag), imag(ideal_rshift), atol = 0.0001)

#Test leftshift of CFixpointArray
c_lshift = cf_val1 << 2;
ideal_lshift = (val1 + 1im .* val2)*4;
@test isapprox(float(c_lshift.real), real(ideal_lshift), atol = 0.0001)
@test isapprox(float(c_lshift.imag), imag(ideal_lshift), atol = 0.0001)

#Test slicing
sliced_cf_val1 = cf_val1[2:end-1]
@test isa(sliced_cf_val1, CFixpointArray)
@test isapprox(float(sliced_cf_val1), c_val1[2:end-1], atol=0.0001)
sliced_cf_val1[:] = cf_val1[1:end-2].real.data + 1im * cf_val1[1:end-2].imag.data
@test isapprox(float(sliced_cf_val1), c_val1[1:end-2], atol=0.0001)
sliced_cf_val1[:] = cf_val1[2:end-1]
@test isapprox(float(sliced_cf_val1), c_val1[2:end-1], atol=0.0001)
sliced_cf_val1[:] = c_val2[1:end-2]
@test isapprox(float(sliced_cf_val1), c_val2[1:end-2], atol=0.0001)

@test isa(sliced_cf_val1[1], CFixpoint)
sliced_cf_val1[1] = cf_scalar.real.data + 1im*cf_scalar.imag.data
@test isapprox(float(sliced_cf_val1[1]), scalar_complex, atol=0.0001)
sliced_cf_val1[1] = cf_scalar
@test isapprox(float(sliced_cf_val1[1]), scalar_complex, atol=0.0001)
sliced_cf_val1[1] = scalar_complex
@test isapprox(float(sliced_cf_val1[1]), scalar_complex, atol=0.0001)