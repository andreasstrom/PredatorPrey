import matplotlib
matplotlib.use('Agg')
from dolfin import *
import numpy as np
from matplotlib import pyplot as plt
import random

# Create mesh and define function space
mesh = Mesh("circle_fine.xml.gz")

# Construct the finite element space
W = VectorFunctionSpace(mesh, 'P', 1)

# Define parameters :
T = 1000
dt = 0.5
alpha = 0.4
beta = 2.0
gamma = 0.8
task = 2

# Class representing the intial conditions
if task == 1:
	class InitialConditions(UserExpression):
		def eval(self, values, x):
			values[0] = 4.0/15.0-(2.0e-7)*(x[0]-0.1*x[1]-225.0)*(x[0]-0.1*x[1]-675.0)
			values[1] = 22.0/45.0-(3.0e-5)*(x[0]-450.0)-(1.2e-4)*(x[1]-150.0)

		def value_shape(self):
			return(2,)
elif task == 2:
	class InitialConditions(UserExpression):
		def eval(self, values, x):
			values[0] = 0.5*(1.0-random.uniform(0,1))
			values[1] = 0.25+0.5*random.uniform(0,1)

		def value_shape(self):
			return(2,)

# Define initial condition
indata = InitialConditions(degree=2)
u0 = Function(W)
u0 = interpolate(indata, W)

u,w = TrialFunction(W), TestFunction(W)

a0 = 2*u[0]*w[0]*dx - dt*u[0]*w[0]*dx \
	+ dt*inner(grad(u[0]),grad(w[0]))*dx

a1 = 2*u[1]*w[1]*dx + gamma*dt*u[1]*w[1]*dx \
	+dt*inner(grad(u[1]),grad(w[1]))*dx

L0 = 2*u0[0]*w[0]*dx + dt*u0[0]*w[0]*dx - dt*inner(grad(u0[0]),grad(w[0]))*dx \
	-2*dt*((u0[0]**2)+(u0[0]*u0[1])/(u0[0]+alpha))*w[0]*dx

L1 = 2*u0[1]*w[1]*dx - gamma*dt*u0[1]*w[1]*dx - dt*inner(grad(u0[1]),grad(w[1]))*dx \
	+2*dt*beta*(u0[0]*u0[1])/(u0[0]+alpha)*w[1]*dx

a = a0+a1
L = L0+L1
out_file = File("results/"+str(task)+"sln.pvd","compressed")

t = 0.0
u = Function(W)
u.assign(u0)
time = np.arange(0,T+dt,dt).tolist()
pop_u = list()
pop_v = list()
out_file << (u,t)

while t<=T:
	u0.assign(u)
	A = assemble(a)
	b = assemble(L)
	solve(A,u.vector(),b,"bicgstab","default")
	pop_u.append(assemble(u[0]*dx))
	pop_v.append(assemble(u[1]*dx))
	if t in [50,100,150,1000]:
		print("Saving sln")
		out_file << (u,t)
	t += dt

plt.plot(time,pop_u,label="Prey")
plt.plot(time,pop_v,label="Predators")
plt.xlabel("Time (s)")
plt.ylabel("Population Rate")
plt.legend()
plt.title("Population Rates")
plt.savefig(str(task)+"pop_rates.eps")
