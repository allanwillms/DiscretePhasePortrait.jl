using DiscretePhasePortrait
using Test

@testset "DiscretePhasePortrait.jl" begin
    # Write your tests here.
	function ricker(parms)
		K,a,L,b = parms
		rickerbase(u,v,KK,aa) = exp(KK - u - aa*v)
		F(x,y) = x * rickerbase(x,y,K,a)
		G(x,y) = y * rickerbase(y,x,L,b)
		function Jac(x,y)
			J = Array{Float64}(undef,2,2)
			J[1,1] = rickerbase(x,y,K,a) - F(x,y)
			J[1,2] = -a * F(x,y)
			J[2,1] = -b * G(x,y)
			J[2,2] = rickerbase(y,x,L,b) - G(x,y)
			return J
		end
		return (F,G,Jac)
	end
	@testset "output tests" begin
		lims = [1.5 2.5; 0.0 0.1]
		K = 1.5
		a = 0.4
		L = 1.9
		b = 0.3  
		F,G,Jac = ricker([K,a,L,b])
		x0,y0 = 1.3,1.3
		p = phaseportrait(F,G,lims)
		@test p == nothing
		p = phaseportrait(F,G,lims,fixpt=[x0,y0],Jac=Jac)
		@test p[1][1] ≈ 0.8409090 atol = 1.0e-6
		@test p[1][2] ≈ 1.6477272 atol = 1.0e-6
		@test p[2][1] ≈-0.8179117 atol = 1.0e-6
		@test p[2][2] ≈ 0.3292753 atol = 1.0e-6
		@test p[4][1] ≈ 1.2392252 atol = 1.0e-6
		@test p[4][2] ≈ 2.6731933 atol = 1.0e-6
	end


end
