
using Drifters
import Drifters.CyclicArrays as CyclicArrays
import Drifters.OrdinaryDiffEq as OrdinaryDiffEq 

function cyclicarray_example()
    𝑃=cyclicarray_setup()

    x0=collect(0:20:490).%100 .+ 170
    y0=reshape(repeat(collect(0:20:80),1,5)',25) .+ 160
    uInit=Float64.([x0';y0'])

    N=1500
    tspan = (0.0,N*𝑃.dt)

    prob = OrdinaryDiffEq.ODEProblem(dxy_dt_CyclicArray,uInit,tspan,𝑃)
    #sol = OrdinaryDiffEq.solve(prob,Tsit5(),reltol=1e-5,abstol=1e-5)
    sol = OrdinaryDiffEq.solve(prob,Drifters.Euler(),dt=𝑃.dt)

    x=[sol[1,i,end] for i in 1:length(x0)]
    y=[sol[2,i,end] for i in 1:length(x0)]

    return x,y
end

function cyclicarray_setup()
    nx=ny=100
    nfaces=1;

    faces=zeros(nfaces,2,2,4);
    faces[1,1,1,:]=[1,1,2,0];
    faces[1,1,2,:]=[1,1,1,0];
    faces[1,2,1,:]=[1,2,2,0];
    faces[1,2,2,:]=[1,2,1,0];

    g=CyclicArrays.CyclicArray(faces);

    ind=range(1, 360, length = 360)
    x = y = ind/180*pi
    f(y,x) = sin(x) + sin(y)

    X = repeat(reshape(x, 1, :), length(y), 1)
    Y = repeat(y, 1, length(x))

    F=f.(Y,X)

    phi=CyclicArrays.CyclicArray(F,g)
    xg=CyclicArrays.CyclicArray(repeat(ind, 1, length(x))',g)
    yg=CyclicArrays.CyclicArray(repeat(ind, 1, length(x)),g)
    u=diff(phi,dims=2)
    v=-diff(phi,dims=1)

    (u = u, v = v, xg = xg, yg = yg, dt = 100)
end
