#################################
#########PROGRAM PROPPER#########
#################################

push!(LOAD_PATH,"./")
using LinearAlgebra
import DustyShearingWave 
const DS=DustyShearingWave


global a=1.5
Σ₀(x)=10.0.*DS.smoothG(a.*(x.-0.35).+0.5).*DS.smoothG(0.5.-a.*(x.-0.65)).+0.46
global P=0.45
global κ₂=1.0
global cₛ=1.0
global tₛ=0.01

W(x,ω)=Σ₀(x)./(P./(ω^2 -κ₂) .+im.*(cₛ^2*Σ₀(x) .- P)./ω)

global ω_guess=2.101-0.01im

#guess=DS.ω_find(8,ω_guess,W)
#Calculate a guess for ω
function getΩ!(guess)
	for i∈[8 10 12 14 15 16]
		guess=DS.ω_find(i,guess,W)
		println("N=",i)
		println(guess)
		guess=Complex(vec(guess.zero)[1],vec(guess.zero)[2])
	end
	return guess
end

global ω=getΩ!(ω_guess)
println((ω))
#global λ,v̄=DS.Eig(64,guess,W)

#println(λ)
