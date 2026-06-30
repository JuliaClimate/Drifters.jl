
function Drifters_read_csv(ii=14,jj=1,D=missing; input=false, path=pwd(), time=1)
    	if input
		fil=joinpath(path,"global_ocean_circulation_inputs","initial_$(ii)_$(jj).csv")
	else
		fil=joinpath(path,"global_ocean_circulation_outputs","initial_$(ii)_$(jj)_▶▶.csv")
	end
	df0=CSV.read(fil,DataFrame)

	df=input ? df0 : DataFrame(groupby(df0,:t)[time])
	input ? (df.fid=Int.(df.f)) : nothing

	for col in names(df)
	    df[!, col] = Vector(df[!, col])
	end
	
    df.ii.=ii; df.jj.=jj;
	add_lonlat!(df,D.XC,D.YC)
end

function Drifters_read_csv(D=missing; input=false, path=pwd(), time=1)
    all_positions=DataFrame()
    for ii in 1:14
        for jj in 1:6
        tmp1=Drifters_read_csv(ii,jj,D,input=input, time=time, path=path)
        append!(all_positions,tmp1)
        end
    end
    all_positions
end

