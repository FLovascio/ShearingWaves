module DustyShearingWave

using AbstractFFTs
using FFTW
using LinearAlgebra
using NLsolve
using ForwardDiff
using Roots
using IntervalArithmetic, IntervalRootFinding

export smoothG, ω_find, B̂, Eig, Â

function smoothF(t)
	retArray=0.0.*t
	retArray[t.>0.0]=ℯ.^(-1.0./t[t.>0.0])
	return retArray
end

smoothG(t)=smoothF(t)./(smoothF(t).+smoothF(1.0.-t))

function B̂(w)
	N=length(w)
	halfN=Int(floor(N/2))
	q=fft(w)
	B=zeros(Complex{Float64},N,N)
	for i=(-halfN+1):(halfN)
		for j=(-halfN+1):(halfN)
			B[halfN+i,halfN+j]=-2.0*π*i*q[abs(j-i)]/float(N)
		end
	end
	return B
end

function Â(w)
	N=length(w)
	halfN=Int(floor(N/2))
	q=fft(1.0 ./ w)
	#println(N)
	A=zeros(Complex{Float64}, N, N)
	#println(A)
	for i=(1-halfN):(halfN)
		for j=(1-halfN):(halfN)
			#println(i,j)
			A[halfN+i,halfN+j]=4.0*(π^2)*(j^2)*q[abs(j-i)+1]./float(N)
		end
	end
	return A
end

function Eig(N,ω,W::Function)
	x=collect(1.0:1.0:N)
	w=W(x,ω)
	return eigen(Â(w))
end

Determinant(w)=det(Â(w)-I)

function ω_find(N,ω₀,W)
	x=collect(1.0:1.0:N)
  w(ω₁)=W(x,ω₁)
	function f!(F,z)
		Z=Complex(z[1],z[2])
		F[1] = real(Determinant(w(Z)))
		F[2] = imag(Determinant(w(Z)))
	end
	ωguess=[real(ω₀) imag(ω₀)]
	return nlsolve(f!,ωguess)
end 

end









