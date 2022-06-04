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
sliced_f_val1[:] = f_val1[1:end-2]
 
"""
CFixpointArray testing
"""
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

# #Test subtraction of CFixpoint types
# c_min = cf_val1 - cf_val2;
# ideal_cmin = c_val1 .- c_val2;
# @test isapprox(float(c_min), ideal_cmin, atol=0.0001)

# #Test summation of CFixpoint types
# c_sum = sum(cf_val1,dims=1);
# ideal_csum = sum(c_val1,dims=1);
# @test any(abs.(toComplex(c_sum) .- ideal_csum) .<0.0001)

# #Test rightshift of CFixpoint types
# c_rshift = cf_val1 >> 1;
# ideal_rshift_re = cf_val1.real.data .>> 1
# ideal_rshift_im = cf_val1.imag.data .>> 1
# @test any((abs.(c_rshift.real.data .- ideal_rshift_re).-abs.(c_rshift.imag.data .- ideal_rshift_im)) .<0.0001)

