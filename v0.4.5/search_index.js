var documenterSearchIndex = {"docs":
[{"location":"particle_cloud/#Particle-Cloud","page":"Particle Cloud","title":"Particle Cloud","text":"","category":"section"},{"location":"particle_cloud/","page":"Particle Cloud","title":"Particle Cloud","text":"(Image: ) (Image: )","category":"page"},{"location":"particle_cloud/","page":"Particle Cloud","title":"Particle Cloud","text":"Using the same setup as detailed_look.jl or example2(), here we simulate a point cloud getting advected by the flow field. For additional documentation e.g. see : 1, 2, 3, 4","category":"page"},{"location":"particle_cloud/#1.-Import-Software","page":"Particle Cloud","title":"1. Import Software","text":"","category":"section"},{"location":"particle_cloud/","page":"Particle Cloud","title":"Particle Cloud","text":"using IndividualDisplacements, Statistics\n\nimport IndividualDisplacements.OrdinaryDiffEq: solve, Tsit5\nimport IndividualDisplacements.DataFrames: DataFrame\n\np=dirname(pathof(IndividualDisplacements))\ninclude(joinpath(p,\"../examples/more/example123.jl\"));\ninclude(joinpath(p,\"../examples/more/recipes_plots.jl\"))","category":"page"},{"location":"particle_cloud/#2.-Setup-Problem","page":"Particle Cloud","title":"2. Setup Problem","text":"","category":"section"},{"location":"particle_cloud/","page":"Particle Cloud","title":"Particle Cloud","text":"𝑃,Γ=example2_setup()\n\nii1=5:5:40; ii2=5:5:25\nx=vec([x-0.5 for x in ii1, y in ii2])\ny=vec([y-0.5 for x in ii1, y in ii2])\nxy = permutedims([[x[i];y[i];1.0] for i in eachindex(x)])\n\nsolv(prob) = IndividualDisplacements.ensemble_solver(prob,solver=Tsit5(),reltol=1e-6,abstol=1e-6)\ntr = DataFrame(ID=Int[], x=Float64[], y=Float64[], t=Float64[])\n\n#𝐼 = Individuals{Float64,2}(📌=xy[:,:], 🔴=tr, 🆔=collect(1:size(xy,2)),","category":"page"},{"location":"particle_cloud/","page":"Particle Cloud","title":"Particle Cloud","text":"                    🚄 = dxdt!, ∫ = solv, 🔧 = postprocess_xy, 𝑃=𝑃);","category":"page"},{"location":"particle_cloud/","page":"Particle Cloud","title":"Particle Cloud","text":"I=(position=xy,record=deepcopy(tr),velocity=dxdt!,\n   integration=solv,postprocessing=postprocess_xy,parameters=𝑃)\n𝐼=Individuals(I)","category":"page"},{"location":"particle_cloud/#3.-Compute-Trajectories","page":"Particle Cloud","title":"3. Compute Trajectories","text":"","category":"section"},{"location":"particle_cloud/","page":"Particle Cloud","title":"Particle Cloud","text":"𝑇 = (0.0,2998*3600.0)\n∫!(𝐼,𝑇)","category":"page"},{"location":"particle_cloud/#4.-Display-results","page":"Particle Cloud","title":"4. Display results","text":"","category":"section"},{"location":"particle_cloud/","page":"Particle Cloud","title":"Particle Cloud","text":"𝐼.🔴.lon=5000*𝐼.🔴.x\n𝐼.🔴.lat=5000*𝐼.🔴.y\nplt=plot_paths(𝐼.🔴,size(xy,2),100000.0)","category":"page"},{"location":"particle_cloud/","page":"Particle Cloud","title":"Particle Cloud","text":"Compare with trajectory output from MITgcm","category":"page"},{"location":"particle_cloud/","page":"Particle Cloud","title":"Particle Cloud","text":"#IndividualDisplacements.flt_example_download()\n#df=read_flt(IndividualDisplacements.flt_example_path,Float32)\n#ref=plot_paths(df,size(xy,2),100000.0)","category":"page"},{"location":"particle_cloud/","page":"Particle Cloud","title":"Particle Cloud","text":"","category":"page"},{"location":"particle_cloud/","page":"Particle Cloud","title":"Particle Cloud","text":"This page was generated using Literate.jl.","category":"page"},{"location":"examples/#Examples","page":"Examples","title":"Examples","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"The four examples outlined below form a tutorial of sorts that complements the User Guide. They hopefully provide a useful jumping off point in order to configure IndividualDisplacements.jl for new problems.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"Output generated by IndividualDisplacements.jl is in DataFrames.jl tabular format which comes with powerful and convenient analysis methods. They can  readily be plotted in space and time using, e.g. , Plots.jl or Makie.jl.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"To rerun an example yourself, the recommended method is to copy the corresponding notebook (code) link, paste it into the Pluto.jl prompt, and click open.","category":"page"},{"location":"examples/#Simple-Two-Dimensional-Flow","page":"Examples","title":"Simple Two-Dimensional Flow","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"notebook (html) ➭ notebook (code)","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"Simulate an ensemble of displacements (and trajectories) in a simple 2D configuration. ","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"The FlowFields constructor readily wraps a flow field provided in the standard Array format, adds a time range, and returns a FlowFields data structure 𝐹. ","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"All that is left to do at this stage is to define initial positions for the individuals, put them together with 𝐹 into the Individuals data structure 𝐼, and call ∫!(𝐼).","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"Exercises include the non-periodic domain case, statistics made easy via DataFrames.jl, and replacing the flow field with your own.","category":"page"},{"location":"examples/#Simple-Three-Dimensional-Flow","page":"Examples","title":"Simple Three-Dimensional Flow","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"notebook (html) ➭ notebook (code)","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"Set up a three-dimensional flow field u,v,w, initialize a single particle at position 📌, and wrap everything up within an Individuals data structure 𝐼.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"𝐼 is displaced by integrating the individual velocity, moving along through space, over time 𝑇.  This is the main computation done in this package – interpolating u,v,w to individual positions 𝐼.📌 on the fly, using 𝐼.🚄, and integrating through time, using 𝐼.∫.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"The flow field consists of rigid body rotation, plus a convergent term, plus a sinking term in the vertical direction. This flow field generates a downward, converging spiral – a idealized version of a relevant case in the Ocean.","category":"page"},{"location":"examples/#Global-Ocean-Circulation","page":"Examples","title":"Global Ocean Circulation","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"notebook (html) ➭ notebook (code)","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"A simulation of floating particles over the Global Ocean which illustrates (1) using time variable velocity fields, (2) global connections, (3) particle re-seeding, and (4) output statistics. ","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"The flow field is based on a data-constrained ocean model solution. The problem is configured in a way to mimic, albeit very crudely, the near-surface tranport of plastics or planktons.","category":"page"},{"location":"examples/#Three-Dimensional-Pathways","page":"Examples","title":"Three Dimensional Pathways","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"notebook (html) ➭ notebook (code)","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"A simulation of particles that follow the three-dimensional ocean circulation. This example illustrates (1) the 3D case in a relatistic configuration, (2) tracking the advent or origin of a water patch, and (3) multifacted visualizations in 3D.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"The flow field is based on a data-constrained, realistic, ocean model. The problem configuration mimics, albeit very approximately, ocean tracers / coumpounds transported by water masses.","category":"page"},{"location":"examples/#Additional-Examples","page":"Examples","title":"Additional Examples","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"Interactivity : notebook (html) ➭ notebook (code)\nAtmosphere : notebook (html) ➭ notebook (code)\nMITgcm : Particle cloud, Detailed look \nPlotting : Plots.jl, Makie.jl, PyPlot.jl","category":"page"},{"location":"examples/#Animations","page":"Examples","title":"Animations","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"Below are animations of results generated in Global Ocean Circulation and Three Dimensional Pathways.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"(Image: simulated particle movie (5m))","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"(Image: simulated particle movie (3D))","category":"page"},{"location":"workflow/#User-Guide","page":"User Guide","title":"User Guide","text":"","category":"section"},{"location":"workflow/","page":"User Guide","title":"User Guide","text":"As shown in the Examples section, the typical worflow is:","category":"page"},{"location":"workflow/","page":"User Guide","title":"User Guide","text":"set up a FlowFields data structure (𝐹)\nset up Individuals (𝐼) with initial position 📌 and 𝐹\ndisplace 𝐼.📌 by\t∫!(𝐼,𝑇) following 𝐼.𝐹 over 𝑇 \npost-process by 𝐼.🔧 and record information in 𝐼.🔴\ngo back to step 2 and continue if needed","category":"page"},{"location":"workflow/","page":"User Guide","title":"User Guide","text":"The data structures for steps 1 and 2 are documented below. Steps 3 and 4 normally take place as part of ∫! which post-processes results, using 🔧, records them in 🔴, and finally updates the positions of individuals in 📌. Since 🔴 is a DataFrame, it is easily manipulated, plotted, or saved in step 4 or after the fact.","category":"page"},{"location":"workflow/#Overview","page":"User Guide","title":"Overview","text":"","category":"section"},{"location":"workflow/","page":"User Guide","title":"User Guide","text":"A central goal of this package is to support scientific analysis of model simulations and observations of materials and particles within the Climate System. The package scope thus includes, for example, drifting plastics in the Ocean or chemical compounds in the Atmosphere. ","category":"page"},{"location":"workflow/","page":"User Guide","title":"User Guide","text":"As a starting point, the package supports all types of gridded model output (on Arakawa C-grids) from the MIT General Circulation Model via the MeshArrays.jl package (docs found here). ","category":"page"},{"location":"workflow/","page":"User Guide","title":"User Guide","text":"By convention, IndividualDisplacements.jl expects input flow fields to be provided in a uniform fashion (see FlowFields) summarized below: ","category":"page"},{"location":"workflow/","page":"User Guide","title":"User Guide","text":"normalized to grid index units (i.e. in 1/s rather than m/s units)\npositive towards increasing indices\nusing the Arakawa C-grid, with u (resp v) staggered by -0.5 point in direction 1 (resp 2) from grid cell centers. ","category":"page"},{"location":"workflow/","page":"User Guide","title":"User Guide","text":"The Examples section documents various, simple methods to prepare and ingest such flow fields (time varying or not; in 2D or 3D) and derive individual displacements / trajectories from them. They cover simple grids often used in process studies, Global Ocean simulations normally done on more complex grids, plotting tools, and data formats. ","category":"page"},{"location":"workflow/","page":"User Guide","title":"User Guide","text":"For an overview of the examples, please refer to the example guide. The rest of this section is focused on the package's data structures and core functions.","category":"page"},{"location":"workflow/#Data-Structures","page":"User Guide","title":"Data Structures","text":"","category":"section"},{"location":"workflow/","page":"User Guide","title":"User Guide","text":"The Individuals struct contains a FlowFields struct (incl. e.g. arrays), initial positions for the individuals, and the other elements (e.g. functions) involved in ∫!(𝐼,𝑇) as documented hereafter.","category":"page"},{"location":"workflow/","page":"User Guide","title":"User Guide","text":"Modules = [IndividualDisplacements]\nOrder   = [:type]","category":"page"},{"location":"workflow/#IndividualDisplacements.FlowFields","page":"User Guide","title":"IndividualDisplacements.FlowFields","text":"abstract type FlowFields\n\nData structure that provide access to flow fields (gridded as arrays) which will be  used to interpolate velocities to individual locations later on (once embedded in an Individuals struct). \n\nFollowing the C-grid convention also used in MITgcm (https://mitgcm.readthedocs.io)  flow fields are expected to be staggered as follows: grid cell i,j has its center located at i-1/2,j-1/2 while the corresponding u[i,j] (resp. `v[i,j]) is located at i-1,j-1/2 (resp. i-1/2,j-1). \n\nAlso by convention, velocity fields are expected to have been normalized to grid units (e.g. 1/s rather than m/s) before sending them to one of the supported FlowFields constructors (using either Array or MeshArray):\n\n𝐹_Array2D(u0,u1,v0,v1,𝑇)\n𝐹_Array3D(u0,u1,v0,v1,w0,w1,𝑇)\n𝐹_MeshArray2D(u0,u1,v0,v1,𝑇,update_location!)\n𝐹_MeshArray3D(u0,u1,v0,v1,w0,w1,𝑇,update_location!)\n\nUsing the FlowFields constructor which gets selected by the type of u0 etc. For example :\n\n𝐹=FlowFields(u,u,v,v,0*w,1*w,[0.0,10.0])\n𝐹=FlowFields(u,u,v,v,[0.0,10.0],func)\n\nas shown in the online documentation examples.\n\n\n\n\n\n","category":"type"},{"location":"workflow/#IndividualDisplacements.Individuals","page":"User Guide","title":"IndividualDisplacements.Individuals","text":"struct Individuals{T,N}\n\nData:           📌 (position),   🔴(record), 🆔 (ID), 𝑃 (FlowFields)\nFunctions:      🚄 (velocity),   ∫ (integration), 🔧(postprocessing)\nNamedTuples:    𝐷 (diagnostics),      𝑀 (metadata)\n\nThe velocity function 🚄 typically computes velocity at individual positions (📌 to start) within the  specified space-time domain by interpolating gridded variables (provided via 𝑃). Individual trajectories  are computed by integrating (∫) interpolated velocities through time. Normally, integration is done by  calling ∫! which updates 📌 at the end and records results in 🔴 via 🔧. Ancillary data, for use in  🔧 for example, can be provided in 𝐷 and metadata stored in 𝑀.\n\nUnicode cheatsheet:\n\n📌=\\:pushpin:<tab>,          🔴=\\:red_circle:<tab>, 🆔=\\:id:<tab>\n🚄=\\:bullettrain_side:<tab>, ∫=\\int<tab>,          🔧=\\:wrench:<tab>\n𝑃=\\itP<tab>,                 𝐷=\\itD<tab>,           𝑀=\\itM<tab>\n\nSimple constructors that use FlowFields to choose adequate defaults:\n\nIndividuals(𝐹::𝐹_Array2D,x,y)\nIndividuals(𝐹::𝐹_Array3D,x,y,z)\nIndividuals(𝐹::𝐹_MeshArray2D,x,y,fid)\nIndividuals(𝐹::𝐹_MeshArray3D,x,y,z,fid)\n\nFurther customization is achievable via keyword constructors:\n\ndf=DataFrame( ID=[], x=[], y=[], z=[], t = [])\n𝐼=Individuals{Float64,2}(📌=zeros(3,10),🆔=1:10,🔴=deepcopy(df))\n𝐼=Individuals(📌=zeros(3,2),🆔=collect(1:2),🔴=deepcopy(df))\n\nOr via the plain text (or no-unicode) constructors:\n\ndf=DataFrame( ID=[], x=[], y=[], z=[], t = [])\nI=(position=zeros(3,2),ID=1:2,record=deepcopy(df))\nI=Individuals(I)\n\n\n\n\n\n","category":"type"},{"location":"workflow/#Main-Functions","page":"User Guide","title":"Main Functions","text":"","category":"section"},{"location":"workflow/","page":"User Guide","title":"User Guide","text":"∫!(𝐼,𝑇) displaces individuals 𝐼 continuously over time period 𝑇 according to velocity function 🚄, temporal integration method ∫, and post-processor 🔧 (all embedded within 𝐼).","category":"page"},{"location":"workflow/","page":"User Guide","title":"User Guide","text":"∫!","category":"page"},{"location":"workflow/#IndividualDisplacements.∫!","page":"User Guide","title":"IndividualDisplacements.∫!","text":"∫!(𝐼::Individuals,𝑇::Tuple)\n\nDisplace simulated individuals continuously through space over time period 𝑇 starting from position 📌. \n\nThis is typically achieved by computing the cumulative integral of velocity experienced by each individual along its trajectory (∫ 🚄 dt).\nThe current default is solve(prob,Tsit5(),reltol=1e-8,abstol=1e-8) but all solver options from the OrdinaryDiffEq.jl package are available.\nAfter this, ∫! is also equipped to postprocess results recorded into 🔴 via the 🔧 workflow, and the last step in ∫! consists in updating 📌 to be ready for continuing in a subsequent call to ∫!.\n\n\n\n\n\n∫!(𝐼::Individuals)\n\nCall ∫!(𝐼::Individuals,𝐼.𝑃.𝑇)\n\n\n\n\n\n","category":"function"},{"location":"workflow/","page":"User Guide","title":"User Guide","text":"The velocity interpolation functions (🚄) carry out the central computation of this package – interpolating gridded flow fields to individual positions 📌. It is normally called via ∫! to integrate velocity 🚄 over a chosen time period. ","category":"page"},{"location":"workflow/","page":"User Guide","title":"User Guide","text":"Velocity interpolation for several array and grid types.\nPreprocessing and postprocessing methods.\nI/O routines to read (write) results from (to) file.","category":"page"},{"location":"workflow/","page":"User Guide","title":"User Guide","text":"and other functionalities provided in src/compute.jl and src/data_wrangling.jl are further documented in the Tool Box section. Ingestion of trajectory data which have been collected by the Ocean Drifting Buoy Program (movie) is also supported.","category":"page"},{"location":"#IndividualDisplacements.jl","page":"Introduction","title":"IndividualDisplacements.jl","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"IndividualDisplacements.jl computes elementary point displacements over gridded domains. Inter-operability with global climate model output (on Arakawa C-grids in general) that can be represented via MeshArrays.jl (see these docs) is a central goal of IndividualDisplacements.jl. ","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"The initial example suite relies on gridded ocean transports estimated in OCCA (Forget 2010), ECCOv4 (Forget et al. 2015), and CBIOMES (Forget, 2018). For the Atmosphere, an additional example based on an idealized model simulation is provided in MITgcmTools.jl. ","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"Typical applications include the simulation and analysis of materials moving through atmospheric flows (e.g. dust or chemicals) or oceanic flows (e.g. plastics or planktons). An illustration of the type of observations that IndividualDisplacements.jl can simulate is provided below. These data collected by Ocean Drifting Buoy can be accessed via OceanRobots.jl.","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"(Image: Global Drifter Program data)","category":"page"},{"location":"detailed_look/#Detailed-Look","page":"Detailed Look","title":"Detailed Look","text":"","category":"section"},{"location":"detailed_look/","page":"Detailed Look","title":"Detailed Look","text":"(Image: ) (Image: )","category":"page"},{"location":"detailed_look/","page":"Detailed Look","title":"Detailed Look","text":"A more detailed look at spatial interpolation, integration through time, and I/O. For additional documentation e.g. see 1, 2, 3, 4. Here we illustrate a few things in more detail:","category":"page"},{"location":"detailed_look/","page":"Detailed Look","title":"Detailed Look","text":"reading velocities from file.\ngridded velocity output (Udata, Vdata)\npre-computed trajectory output (float_traj*data)\ninterpolating U,V from gridded output to individual locations\ncompared with u,v from float_traj*data\ncomputing trajectories (location v time) using OrdinaryDiffEq.jl\ncompared with x(t),y(t) from float_traj*data","category":"page"},{"location":"detailed_look/#1.-Import-Software","page":"Detailed Look","title":"1. Import Software","text":"","category":"section"},{"location":"detailed_look/","page":"Detailed Look","title":"Detailed Look","text":"using IndividualDisplacements, MITgcmTools\nimport IndividualDisplacements.OrdinaryDiffEq as OrdinaryDiffEq\nimport IndividualDisplacements.DataFrames as DataFrames\np=dirname(pathof(IndividualDisplacements))\ninclude(joinpath(p,\"../examples/more/example123.jl\"))\ninclude(joinpath(p,\"../examples/more/recipes_plots.jl\"))","category":"page"},{"location":"detailed_look/#2.-Read-Trajectory-Output","page":"Detailed Look","title":"2. Read Trajectory Output","text":"","category":"section"},{"location":"detailed_look/","page":"Detailed Look","title":"Detailed Look","text":"from MITgcm/pkg/flt","category":"page"},{"location":"detailed_look/","page":"Detailed Look","title":"Detailed Look","text":"IndividualDisplacements.flt_example_download()\ndirIn=IndividualDisplacements.flt_example_path\nprec=Float32\ndf=read_flt(dirIn,prec);\n\nplt=plot_paths(df,300,100000.0)","category":"page"},{"location":"detailed_look/#3.-Read-Gridded-Variables","page":"Detailed Look","title":"3. Read Gridded Variables","text":"","category":"section"},{"location":"detailed_look/","page":"Detailed Look","title":"Detailed Look","text":"using MeshArrays.jl and e.g. a NamedTyple","category":"page"},{"location":"detailed_look/","page":"Detailed Look","title":"Detailed Look","text":"𝑃,Γ=example2_setup();","category":"page"},{"location":"detailed_look/#4.-Visualize-Velocity-Fields","page":"Detailed Look","title":"4. Visualize Velocity Fields","text":"","category":"section"},{"location":"detailed_look/","page":"Detailed Look","title":"Detailed Look","text":"plt=heatmap(Γ.mskW[1,1].*𝑃.u0,title=\"U at the start\")\nplt=heatmap(Γ.mskW[1,1].*𝑃.u1-𝑃.u0,title=\"U end - U start\")","category":"page"},{"location":"detailed_look/#5.-Visualize-Trajectories","page":"Detailed Look","title":"5. Visualize Trajectories","text":"","category":"section"},{"location":"detailed_look/","page":"Detailed Look","title":"Detailed Look","text":"(select one trajectory)","category":"page"},{"location":"detailed_look/","page":"Detailed Look","title":"Detailed Look","text":"tmp=df[df.ID .== 200, :]\ntmp[1:4,:]","category":"page"},{"location":"detailed_look/","page":"Detailed Look","title":"Detailed Look","text":"Super-impose trajectory over velocity field (first for u ...)","category":"page"},{"location":"detailed_look/","page":"Detailed Look","title":"Detailed Look","text":"x=Γ.XG.f[1][:,1]\ny=Γ.YC.f[1][1,:]\nz=transpose(Γ.mskW[1].*𝑃.u0);\n\nplt=contourf(x,y,z,c=:delta)\nplot!(tmp[:,:lon],tmp[:,:lat],c=:red,w=4,leg=false)","category":"page"},{"location":"detailed_look/","page":"Detailed Look","title":"Detailed Look","text":"Super-impose trajectory over velocity field (... then for v)","category":"page"},{"location":"detailed_look/","page":"Detailed Look","title":"Detailed Look","text":"x=Γ.XC.f[1][:,1]\ny=Γ.YG.f[1][1,:]\nz=transpose(Γ.mskW[1].*𝑃.v0);\n\nplt=contourf(x,y,z,c=:delta)\nplot!(tmp[:,:lon],tmp[:,:lat],c=:red,w=4,leg=false)","category":"page"},{"location":"detailed_look/#6.-Interpolate-Velocities","page":"Detailed Look","title":"6. Interpolate Velocities","text":"","category":"section"},{"location":"detailed_look/","page":"Detailed Look","title":"Detailed Look","text":"dx=Γ.dx\nuInit=[tmp[1,:lon];tmp[1,:lat]]./dx\nnSteps=Int32(tmp[end,:time]/3600)-2\ndu=fill(0.0,2);","category":"page"},{"location":"detailed_look/","page":"Detailed Look","title":"Detailed Look","text":"Visualize and compare with actual grid point values – jumps on the tangential component are expected with linear scheme:","category":"page"},{"location":"detailed_look/","page":"Detailed Look","title":"Detailed Look","text":"tmpu=fill(0.0,100)\ntmpv=fill(0.0,100)\ntmpx=fill(0.0,100)\nfor i=1:100\n    tmpx[i]=500.0 *i./dx\n    dxdt!(du,[tmpx[i];0.499./dx],𝑃,0.0)\n    tmpu[i]=du[1]\n    tmpv[i]=du[2]\nend\n\nplt=plot(tmpx,tmpu,label=\"u (interp)\")\nplot!(Γ.XG.f[1][1:10,1]./dx,𝑃.u0[1:10,1],marker=:o,label=\"u (C-grid)\")\nplot!(tmpx,tmpv,label=\"v (interp)\")\nplot!(Γ.XG.f[1][1:10,1]./dx,𝑃.v0[1:10,1],marker=:o,label=\"v (C-grid)\")","category":"page"},{"location":"detailed_look/","page":"Detailed Look","title":"Detailed Look","text":"And similarly in the other direction","category":"page"},{"location":"detailed_look/","page":"Detailed Look","title":"Detailed Look","text":"tmpu=fill(0.0,100)\ntmpv=fill(0.0,100)\ntmpy=fill(0.0,100)\nfor i=1:100\n    tmpy[i]=500.0 *i./dx\n    dxdt!(du,[0.499./dx;tmpy[i]],𝑃,0.0)\n    tmpu[i]=du[1]\n    tmpv[i]=du[2]\nend\n\nplt=plot(tmpx,tmpu,label=\"u (interp)\")\nplot!(Γ.YG.f[1][1,1:10]./dx,𝑃.u0[1,1:10],marker=:o,label=\"u (C-grid)\")\nplot!(tmpx,tmpv,label=\"v (interp)\")\nplot!(Γ.YG.f[1][1,1:10]./dx,𝑃.v0[1,1:10],marker=:o,label=\"v (C-grid)\")","category":"page"},{"location":"detailed_look/","page":"Detailed Look","title":"Detailed Look","text":"Compare recomputed velocities with those from pkg/flt","category":"page"},{"location":"detailed_look/","page":"Detailed Look","title":"Detailed Look","text":"nSteps=2998\ntmpu=fill(0.0,nSteps); tmpv=fill(0.0,nSteps);\ntmpx=fill(0.0,nSteps); tmpy=fill(0.0,nSteps);\nrefu=fill(0.0,nSteps); refv=fill(0.0,nSteps);\nfor i=1:nSteps\n    dxy_dt_replay(du,[tmp[i,:lon],tmp[i,:lat]],tmp,tmp[i,:time])\n    refu[i]=du[1]./dx\n    refv[i]=du[2]./dx\n    dxdt!(du,[tmp[i,:lon],tmp[i,:lat]]./dx,𝑃,Float64(tmp[i,:time]))\n    tmpu[i]=du[1]\n    tmpv[i]=du[2]\nend\n\nplt=plot(tmpu,label=\"u\")\nplot!(tmpv,label=\"v\")\nplot!(refu,label=\"u (ref)\")\nplot!(refv,label=\"v (ref)\")","category":"page"},{"location":"detailed_look/#6.-Compute-Trajectories","page":"Detailed Look","title":"6. Compute Trajectories","text":"","category":"section"},{"location":"detailed_look/","page":"Detailed Look","title":"Detailed Look","text":"Solve through time using OrdinaryDiffEq.jl with","category":"page"},{"location":"detailed_look/","page":"Detailed Look","title":"Detailed Look","text":"dxdt! is the function computing d(position)/dt\nuInit is the initial condition u @ tspan[1]\ntspan is the time interval\n𝑃 are parameters for dxdt!\nTsit5 is the time-stepping scheme\nreltol and abstol are tolerance parameters","category":"page"},{"location":"detailed_look/","page":"Detailed Look","title":"Detailed Look","text":"tspan = (0.0,nSteps*3600.0)\n#prob = OrdinaryDiffEq.ODEProblem(dxy_dt_replay,uInit,tspan,tmp)\nprob = OrdinaryDiffEq.ODEProblem(dxdt!,uInit,tspan,𝑃)\nsol = OrdinaryDiffEq.solve(prob,OrdinaryDiffEq.Tsit5(),reltol=1e-8,abstol=1e-8)\nsol[1:4]","category":"page"},{"location":"detailed_look/","page":"Detailed Look","title":"Detailed Look","text":"Compare recomputed trajectories with originals from MITgcm/pkg/flt","category":"page"},{"location":"detailed_look/","page":"Detailed Look","title":"Detailed Look","text":"ref=transpose([tmp[1:nSteps,:lon] tmp[1:nSteps,:lat]])\nmaxLon=80*5.e3\nmaxLat=42*5.e3\n#show(size(ref))\nfor i=1:nSteps-1\n    ref[1,i+1]-ref[1,i]>maxLon/2 ? ref[1,i+1:end]-=fill(maxLon,(nSteps-i)) : nothing\n    ref[1,i+1]-ref[1,i]<-maxLon/2 ? ref[1,i+1:end]+=fill(maxLon,(nSteps-i)) : nothing\n    ref[2,i+1]-ref[2,i]>maxLat/2 ? ref[2,i+1:end]-=fill(maxLat,(nSteps-i)) : nothing\n    ref[2,i+1]-ref[2,i]<-maxLat/2 ? ref[2,i+1:end]+=fill(maxLat,(nSteps-i)) : nothing\nend\nref=ref./dx;\n\nplt=plot(sol[1,:],sol[2,:],linewidth=5,title=\"Using Recomputed Velocities\",\n     xaxis=\"lon\",yaxis=\"lat\",label=\"Julia Solution\") # legend=false\nplot!(ref[1,:],ref[2,:],lw=3,ls=:dash,label=\"MITgcm Solution\")","category":"page"},{"location":"detailed_look/","page":"Detailed Look","title":"Detailed Look","text":"","category":"page"},{"location":"detailed_look/","page":"Detailed Look","title":"Detailed Look","text":"This page was generated using Literate.jl.","category":"page"},{"location":"API/#Velocity-Interpolation","page":"Tool Box","title":"Velocity Interpolation","text":"","category":"section"},{"location":"API/","page":"Tool Box","title":"Tool Box","text":"The dxdt! etc functions compute the tracked individual velocity. ","category":"page"},{"location":"API/","page":"Tool Box","title":"Tool Box","text":"dxdt!\ndxy_dt_replay\ndxy_dt_CyclicArray","category":"page"},{"location":"API/#IndividualDisplacements.dxdt!","page":"Tool Box","title":"IndividualDisplacements.dxdt!","text":"dxdt!(du,u,p::𝐹_MeshArray3D,tim)\n\nInterpolate velocity from gridded fields (3D; with halos) to position u (x,y,z,fIndex) to compute the derivative of position v time  du_dt.\n\nExtended help\n\nusing IndividualDisplacements\nu,v,w,pos,func=vortex_flow_field(format=:MeshArray)\n𝐹=FlowFields(u,u,v,v,0*w,1*w,[0,3*pi],func)\n𝐼=Individuals(𝐹,pos...)\n∫!(𝐼)\n\nref=[6.16, 7.23, 1.29, 1.0] # hide\nprod(isapprox.(𝐼.📌,ref,atol=1.0)) # hide\n\n\n\n\n\ndxdt!(du,u,p::𝐹_MeshArray2D,tim)\n\nInterpolate velocity from gridded fields (2D; with halos) to position u (x,y,fIndex) to compute the derivative of position v time  du_dt.\n\nExtended help\n\nusing IndividualDisplacements\nu,v,w,pos,func=random_flow_field(format=:MeshArray);\n𝐹=FlowFields(u,u,v,v,[0,1.0],func);\n𝐼=Individuals(𝐹,pos...);\n∫!(𝐼);\n\nisa(𝐼.📌,Vector)\n\n\n\n\n\ndxdt!(du,u,𝑃::𝐹_Array3D,tim)\n\nInterpolate velocity from gridded fields (3D; NO halos) to position u (x,y,z) to compute the derivative of position v time  du_dt.\n\nExtended help\n\nusing IndividualDisplacements\nu,v,w,pos=vortex_flow_field(format=:Array)\n𝐹=FlowFields(u,u,v,v,0*w,1*w,[0,3*pi])\n𝐼=Individuals(𝐹,pos...)\n∫!(𝐼)\n\nref=[9.35, 7.93, 1.28] # hide\nprod(isapprox.(𝐼.📌,ref,atol=1.0)) # hide\n\n\n\n\n\ndxdt!(du,u,𝑃::𝐹_Array2D,tim)\n\nInterpolate velocity from gridded fields (2D; NO halos) to position u (x,y) to compute the derivative of position v time  du_dt.\n\nExtended help\n\nusing IndividualDisplacements\nu,v,w,pos=random_flow_field(format=:Array);\n𝐹=FlowFields(u,u,v,v,[0,1.0]);\n𝐼=Individuals(𝐹,pos...);\n∫!(𝐼);\n\nisa(𝐼.📌,Vector)\n\n\n\n\n\n","category":"function"},{"location":"API/#IndividualDisplacements.dxy_dt_replay","page":"Tool Box","title":"IndividualDisplacements.dxy_dt_replay","text":"dxy_dt_replay(du,u,p::DataFrame,t)\n\nInterpolate velocity from MITgcm float_trajectories output and return position increment du.\n\nExtended help\n\nusing IndividualDisplacements\np=dirname(pathof(IndividualDisplacements))\ninclude(joinpath(p,\"../examples/more/detailed_look.jl\"))\nprod(isapprox.(sol[:,end],ref[:,end],atol=1.0))\n\n\n\n\n\n","category":"function"},{"location":"API/#IndividualDisplacements.dxy_dt_CyclicArray","page":"Tool Box","title":"IndividualDisplacements.dxy_dt_CyclicArray","text":"dxy_dt_CyclicArray(du,u,𝑃::NamedTuple,tim)\n\nNearest neighbor (?) velocity from gridded fields (2D; NO halos but not needed when CyclicArrays is used to extend valid indice ranges).\n\nExtended help\n\nnotes: spatial interpolation & temporal interpolation are lacking\n\nusing IndividualDisplacements\np=dirname(pathof(IndividualDisplacements))\ninclude(joinpath(p,\"../examples/more/example_CyclicArray.jl\"))\n(x,y)=cyclicarray_example()\n\nref=[330.5 290.5]\nprod(isapprox.([x[end] y[end]],ref,atol=1.0))\n\n\n\n\n\n","category":"function"},{"location":"API/#Setup-And-Postprocessing","page":"Tool Box","title":"Setup And Postprocessing","text":"","category":"section"},{"location":"API/","page":"Tool Box","title":"Tool Box","text":"Convenience functions to initialize a simulation and post-process the output are provided. ","category":"page"},{"location":"API/","page":"Tool Box","title":"Tool Box","text":"postprocess_xy\npostprocess_MeshArray\nadd_lonlat!","category":"page"},{"location":"API/#IndividualDisplacements.postprocess_xy","page":"Tool Box","title":"IndividualDisplacements.postprocess_xy","text":"postprocess_xy(sol,𝑃::FlowFields,𝐷::NamedTuple; id=missing, 𝑇=missing)\n\nCopy sol to a DataFrame & map position to x,y coordinates, and define time axis for a simple doubly periodic domain\n\n\n\n\n\n","category":"function"},{"location":"API/#IndividualDisplacements.postprocess_MeshArray","page":"Tool Box","title":"IndividualDisplacements.postprocess_MeshArray","text":"postprocess_MeshArray(sol,𝑃::FlowFields,𝐷::NamedTuple; id=missing, 𝑇=missing)\n\nCopy sol to a DataFrame & map position to lon,lat coordinates using \"exchanged\" 𝐷.XC, 𝐷.YC via add_lonlat!\n\n\n\n\n\n","category":"function"},{"location":"API/#IndividualDisplacements.add_lonlat!","page":"Tool Box","title":"IndividualDisplacements.add_lonlat!","text":"add_lonlat!(df::DataFrame,XC,YC)\n\nAdd lon & lat to dataframe using \"exchanged\" XC, YC\n\n\n\n\n\nadd_lonlat!(df::DataFrame,XC,YC,func::Function)\n\nAdd lon & lat to dataframe using \"exchanged\" XC, YC after updating  subdomain indices (via func) if needed (MeshArrays.locationisout)\n\n\n\n\n\n","category":"function"},{"location":"API/","page":"Tool Box","title":"Tool Box","text":"Basic geography support:","category":"page"},{"location":"API/","page":"Tool Box","title":"Tool Box","text":"gcdist\ndiff\nstproj\nstproj_inv\nrandn_lonlat\ninterp_to_lonlat\ninterp_to_xy\nnearest_to_xy","category":"page"},{"location":"API/#IndividualDisplacements.gcdist","page":"Tool Box","title":"IndividualDisplacements.gcdist","text":"gcdist(𝐼::Individuals)\n\nGreat circle distance (gcd in radians) between final and initial positions.\n\n\n\n\n\n","category":"function"},{"location":"API/#Base.diff","page":"Tool Box","title":"Base.diff","text":"Base.diff(𝐼::Individuals)\n\nDifference in grid unit coordinates (dx,dy) between final and initial positions.\n\n\n\n\n\n","category":"function"},{"location":"API/#IndividualDisplacements.stproj","page":"Tool Box","title":"IndividualDisplacements.stproj","text":"stproj(XC,YC;XC0=0.0,YC0=90.0)\n\nApply to XC,YC (longitude, latitude) the stereographic projection which would put XC0,YC0 (longitude, latitude) at x,y=0,0\n\n\n\n\n\n","category":"function"},{"location":"API/#IndividualDisplacements.stproj_inv","page":"Tool Box","title":"IndividualDisplacements.stproj_inv","text":"stproj_inv(xx,yy;XC0=0.0,YC0=90.0)\n\nApply to xx,yy (stereographic projection coordinates) the reverse  of the stereographic projection which puts XC0,YC0 (longitude,  latitude) at x,y=0,0\n\n\n\n\n\n","category":"function"},{"location":"API/#IndividualDisplacements.randn_lonlat","page":"Tool Box","title":"IndividualDisplacements.randn_lonlat","text":"randn_lonlat(nn=1,seed=missing)\n\nRandomly distributed longitude, latitude positions on the sphere.\n\n\n\n\n\n","category":"function"},{"location":"API/#IndividualDisplacements.interp_to_lonlat","page":"Tool Box","title":"IndividualDisplacements.interp_to_lonlat","text":"interp_to_lonlat(X::MeshArray,Γ::NamedTuple,lon,lat)\n\nUse MeshArrays.Interpolate() to interpolate to arbitrary positions (e.g., a regular grid for plotting).\n\nExtended help\n\nusing IndividualDisplacements\nimport IndividualDisplacements: MeshArrays\nγ=MeshArrays.GridSpec(\"LatLonCap\",MeshArrays.GRID_LLC90)\nΓ=MeshArrays.GridLoad(γ,option=\"full\")\n\nlon=[i for i=20.:20.0:380., j=-70.:10.0:70.]\nlat=[j for i=20.:20.0:380., j=-70.:10.0:70.]\ntmp1=interp_to_lonlat(Γ.Depth,Γ,lon,lat)\n\n(f,i,j,w,_,_,_)=MeshArrays.InterpolationFactors(Γ,vec(lon),vec(lat))\nIntFac=(lon=lon,lat=lat,f=f,i=i,j=j,w=w)\ntmp1=interp_to_lonlat(Γ.Depth,IntFac)\n    \nprod(isapprox(maximum(tmp1),5896.,atol=1.0))\n\n# output\n\ntrue\n\n\n\n\n\ninterp_to_lonlat(X::MeshArray,IntFac::NamedTuple)\n\nUse MeshArrays.Interpolate() to interpolate to arbitrary positions (e.g., a regular grid for plotting).\n\n\n\n\n\n","category":"function"},{"location":"API/#IndividualDisplacements.interp_to_xy","page":"Tool Box","title":"IndividualDisplacements.interp_to_xy","text":"interp_to_xy(df::DataFrame,Zin::MeshArray)\n\nInterpolate \"exchanged\" / \"hallo-included\" Zin to df[!,:x], df[!,:y] on df[!,:fid]\n\n\n\n\n\n","category":"function"},{"location":"API/#IndividualDisplacements.nearest_to_xy","page":"Tool Box","title":"IndividualDisplacements.nearest_to_xy","text":"nearest_to_xy(α::MeshArray,x,y,f)\n\nValue of α at eachindex of the grid cell center nearest to x,y on subdomain array / facet f\n\n\n\n\n\nnearest_to_xy(α::Array,x,y)\n\nValue of α at eachindex of the grid cell center nearest to x,y\n\n\n\n\n\n","category":"function"},{"location":"API/#Toy-Problems","page":"Tool Box","title":"Toy Problems","text":"","category":"section"},{"location":"API/","page":"Tool Box","title":"Tool Box","text":"These are used to demonstrate and test the package functionalities:","category":"page"},{"location":"API/","page":"Tool Box","title":"Tool Box","text":"random_flow_field\nvortex_flow_field","category":"page"},{"location":"API/#IndividualDisplacements.random_flow_field","page":"Tool Box","title":"IndividualDisplacements.random_flow_field","text":"random_flow_field(;component=:Rotational,np=12,nq=18,format=:Array)\n\nGenerate random flow fields on a grid of np x nq points for use in simple examples.\n\nThe :Rotational component option is most similar to what is done  in the standard example.\n\nThe :Divergent component option generates a purely divergent flow field instead.\n\n(U,V,Φ)=random_flow_field(component=:Rotational)\nF=convert_to_FlowFields(U,V,10.0)\nI=Individuals(F,x,y,fill(1,length(x)))\n\n\n\n\n\n","category":"function"},{"location":"API/#IndividualDisplacements.vortex_flow_field","page":"Tool Box","title":"IndividualDisplacements.vortex_flow_field","text":"vortex_flow_field(; np=12,nz=4,format=:Array)\n\nSet up an idealized flow field which consists of  rigid body rotation,  plus a convergent term, plus a sinking term.\n\nu,v,w,func=vortex_flow_field(format=:MeshArray)\n\n\n\n\n\n","category":"function"},{"location":"API/#Read-External-Files","page":"Tool Box","title":"Read External Files","text":"","category":"section"},{"location":"API/","page":"Tool Box","title":"Tool Box","text":"Trajectory simulated by the MITgcm or observed by the global drifter program can be read from file using, respectively MITgcmTools.read_flt (from MITgcmTools.jl) or  OceanRobots.drifters_hourly_read (from OceanRobots.jl).","category":"page"}]
}
