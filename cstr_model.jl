using JuMP
using Ipopt

ntfe = 10
time = zeros(ntfe)

m = Model(Ipopt.Optimizer)

print(m)
