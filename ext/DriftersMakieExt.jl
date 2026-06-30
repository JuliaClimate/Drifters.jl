module DriftersMakieExt

	using Makie, Drifters
	import Drifters: Individuals, DriftersDataset, DataFrame, demo, gcdist, MeshArrays
	import Makie: plot

	function plot(I::Individuals,plot_type=:plot_start_end)
		plot(DriftersDataset(data=(I=I,), options=(plot_type=plot_type,)))
	end

	function plot(x::DriftersDataset)
		if !isempty(x.options)
			o=x.options
			if string(o.plot_type)=="simple_plot1"
				simple_plot1(x.data[:I],x.data[:ϕ])
			elseif string(o.plot_type)=="simple_plot2"
				simple_plot2(x.data[:I])
			elseif string(o.plot_type)=="global_plot1"
				global_plot1(x.data[:I],x.data[:df])
			elseif string(o.plot_type)=="plot_start_end"
				plot_start_end(x)
			elseif string(o.plot_type)=="jcon_drifters"
				plot_drifters_jcon(x.data.gdf;x.options...)
			elseif string(o.plot_type)=="Oscar_plot"
				Oscar_plot(x.data.df;x.options...)
			elseif string(o.plot_type)=="density_z"
				Drifters_density1d_plot(x.data.df)
			elseif string(o.plot_type)=="density_xy"
				Drifters_density2d_plot(x.data.df)
			elseif string(o.plot_type)=="initial_positions"
				plot_initial_positions(x.data.D,x.data.df00,x.data.df0)
			else
				println("unknown option (b)")	
			end
		else
			println("unknown option (a)")
		end
	end

##

include("Makie/StochasticDiffEq.jl")
include("Makie/simple.jl")
include("Makie/global.jl")
include("Makie/jcon_paper.jl")
include("Makie/density_plots.jl")

end
