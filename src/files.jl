import CSV

import Base: write

"""
    write(I::Individuals; file=tempname()*".csv")

Write record to a csv file, `fil`, which can later be read using `Drifters.csv_read(fil)`.

"""
write(I::Individuals; file=tempname()*".csv") = CSV.write(file, I.🔴)

"""
    csv_read(fil::String; D=missing, time=missing)

```
Drifters.csv_read(fil,D=D)
```
"""
function csv_read(fil::String; D=missing, time=missing)
	df0=CSV.read(fil,DataFrame)

	df=ismissing(time) ? df0 : DataFrame(groupby(df0,:t)[time])
	in("fid",names(df)) ? nothing : (df.fid=Int.(df.f))

	for col in names(df)
	    df[!, col] = Vector(df[!, col])
	end
	
	ismissing(D) ? df : add_lonlat!(df,D.XC,D.YC)
end

## special case 

"""
    csv_read(ii=14,jj=1; D=missing, output=false, path=pwd(), time=missing)

(special case)

```
Drifters.csv_read(14,1,D=D,path=p0)
Drifters.csv_read_all(D=D,path=p0);
```
"""
function csv_read(ii=14,jj=1; D=missing, output=false, path=pwd(), time=missing)
	if output
		fil=joinpath(path,"global_ocean_circulation_outputs","initial_$(ii)_$(jj)_▶▶.csv")
	else
		fil=joinpath(path,"global_ocean_circulation_inputs","initial_$(ii)_$(jj).csv")
	end
	df=csv_read(fil,D=D,time=time)
	df.ii.=ii; df.jj.=jj;
	df
end

"""
    csv_read_all(; D=missing, output=false, path=pwd(), time=missing)

(special case)

```
Drifters.csv_read_all(D=D,path=p0);
```
"""
function csv_read_all(; D=missing, output=false, path=pwd(), time=missing)
    all_positions=DataFrame()
    for ii in 1:14
        for jj in 1:6
        tmp1=csv_read(ii,jj,D=D,output=output,time=time,path=path)
        append!(all_positions,tmp1)
        end
    end
    all_positions
end

## earlier code 

"""
    read_initial_positions(np ::Int; filename="global_ocean_circulation.csv")

"""
function read_initial_positions(np ::Int; filename="global_ocean_circulation.csv")
    if filename=="global_ocean_circulation.csv"
        p=dirname(pathof(Drifters))
        fil=joinpath(p,"../examples/worldwide/global_ocean_circulation.csv")
    else
        fil=filename
    end
    return DataFrame(CSV.File(fil))[1:np,:]
end
