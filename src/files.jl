import CSV

function csv_read(fil::String; D=missing, time=missing)
	df0=CSV.read(fil,DataFrame)

	df=ismissing(time) ? df0 : DataFrame(groupby(df0,:t)[time])
	in("fid",names(df)) ? (df.fid=Int.(df.f)) : nothing

	for col in names(df)
	    df[!, col] = Vector(df[!, col])
	end
	
	ismissing(D) ? df : add_lonlat!(df,D.XC,D.YC)
end

## special case 

function csv_read(ii=14,jj=1; D=missing, input=false, path=pwd(), time=missing)
	if input
		fil=joinpath(path,"global_ocean_circulation_inputs","initial_$(ii)_$(jj).csv")
	else
		fil=joinpath(path,"global_ocean_circulation_outputs","initial_$(ii)_$(jj)_▶▶.csv")
	end
	df=csv_read(fil,D=D)
	df.ii.=ii; df.jj.=jj;
	df
end

function csv_read_all(; D=missing, input=false, path=pwd(), time=missing)
    all_positions=DataFrame()
    for ii in 1:14
        for jj in 1:6
        tmp1=csv_read(ii,jj,D=D,input=input,time=time,path=path)
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
