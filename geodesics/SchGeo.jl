
using Plots

function velocityVerlet(ODE, lmbda0, q0, lmbdaf, dlmbda)
	N = floor(Int,(lmbdaf - lmbda0)/dlmbda )#Time steps
	lmbda = zeros(N)
	q = zeros(N, 6)
	lmbda[1] = lmbda0
	q[1,:] = q0
	for i in 1:N-1 
		lmbda[i+1] = lmbda[i] + dlmbda
		v_half = q[i,4:6] + 0.5*ODE(q[i,:])[4:6]*dlmbda
		q[i+1,1:3] = q[i,1:3] + v_half*dlmbda
		q[i+1,4:6] = v_half + 0.5*ODE(q[i,:])[4:6]*dlmbda
	end
	return q
end


function geo(q0)
	G = 1.5607939e-13 # [lyr^3 sunMass^-1 yr^-1]
	M = G*1. # BH Mass in geometric units
	sfun = 2*M - q0[2]
	f = zeros(6)
	f[1:3] = q0[4:6]
	f[4] = 2*M*q0[1]*q0[2]/(sfun*q0[2])
	f[5] = M*sfun*q0[4]^2/(q0[2]^3) - M*q0[5]^2/(sfun*q0[2]) - sfun*q0[6]^2
	f[6] = -2*q0[5]*q0[6]/q0[2]
	return f
end

function flat(q0)
	f = zeros(6)
	f[1:3] = q0[4:6]
	f[4] = 0
	f[5] = 0
	f[6] = 0
	return f
end


t0 = 0.
r0 = 500.
ph0 = 0.

vr0 = -sqrt(2)/2
vph0 = sqrt(2)/2
vt0 = sqrt(vr0^2 + vph0^2)
q0 = [t0, r0, ph0, vt0, vr0, vph0]

Q = velocityVerlet(flat,0,q0,10,1)

x = Q[:,2].*cos.(Q[:,3])
y = Q[:,2].*sin.(Q[:,3])

println(Q[:,3])

display(plot(x,y))
readline()



