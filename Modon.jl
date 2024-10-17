# This script combines Modon_1_Layer.jl and Modon_2_Layer.jl into one file and allows for looping through
# several cases for different values of β and different numbers of layers.

using GeophysicalFlows, NetCDF, Random, QGDipoles
using LinearAlgebra: mul!, ldiv!
using Random: seed!
Random.seed!(1)

# Define Case parameters:

β₁, β₂ = [0.5, 0], [0, 0]			# β in upper (1) and lower (2) layers in each case
Layers = [1, 1]					# number of layers for each case
case_name = ["1L0.5", "1L0.5_0"]
IC_type = [0, 1001]				# 0 - use modon IC, N > 0 - use Nth time from IC_file
IC_file = ["0", "1L0.5_all"] 			# file containing IC field(s)

#β₁, β₂ = [0, 1, 0.5, 0, 1, 0, 1], [0, 0, 0, 0, 0, 1, 1]	# β in upper (1) and lower (2) layers in each case
#Layers = [1, 1, 1, 2, 2, 2, 2]				# number of layers for each case
#case_name = ["1L0", "1L1", "1L0.5", "2L0:0", "2L1:0", "2L0:1", "2L1:1"]
#IC_type = [0, 0, 0, 0, 0, 0, 0]				# 0 - use modon IC, N > 0 - use Nth time from IC_file
#IC_file = ["_1", "_2", "_3", "_4", "_5", "_6", "_7"] 	# file containing IC field(s)

# Define other physical parameters:

U, a = 1, 1					# vortex speed and radius for initial condition
R = [1, 1]					# Rossby radius in each layer, R = sqrt(H)
ActiveLayers = [1, 0]				# active/passive layers for initial vortex (2-layer)
init_β = 1					# initiate solution as 1: β-plane modon, 0: f-plane modon

# Define grid parameters:

Nx, Ny = 2048, 2048				# number of gridpoints in x and y directions
Lx, Ly = 20.48, 20.48				# domain size in x and y directions

# Define numerical scheme parameters:

T = 1000					# simulation stop time
C, M₀ = 1, 0					# frame speed and sponge layer damping rate
r_c, σ = Lx/2-2, 0.5				# cutting parameters; radius, width of smoothing region
dev = GPU()					# device, CPU() or GPU() (GPU is much faster)
stepper = "FilteredRK4"				# timestepping method, e.g. "RK4", "LSRK54" or "FilteredRK4"
aliased_fraction = 0				# fraction of wavenumbers zeroed out in dealiasing
κ₁, κ₂ = 2.5π, 5π				# range of wavenumbers where random noise is added to IC

# Define save parameters:

savename = "Test_data/Case"			# Save filename prefix for NetCDF data file
IC_pref	= "Test_data/Case_"			# IC filename prefix, common to all IC files
save_all, save_window = 1, 1			# flags for saving all fields and windowed fields
Ns, Nw = 1000, 5000				# number of full field and windowed saves
Lx₁, Lx₂, Ly₁, Ly₂ = -6, 2, -2, 2		# x and y ranges of save window

# Calculate derived parameters:

Δt = 0.625*((Lx/Nx)+(Ly/Ny))/(5*U)		# timestep, set using approx CFL
Nt = Int(ceil(T/Δt))				# total number of iterations
H = R.^2					# layer depth (required for 2-layer only)
Nc = length(β₁)					# number of cases to run
t_all = LinRange(0,T,Ns+1)			# time vector for all field saves
t_window = LinRange(0,T,Nw+1)			# time vector for windowed saves
ix = Int.(Nx/2+1 .+ (Lx₁*Nx/Lx:Lx₂*Nx/Lx))	# x index for windowed save
iy = Int.(Ny/2+1 .+ (Ly₁*Ny/Ly:Ly₂*Ny/Ly))	# y index for windowed save
Nx_w, Ny_w = length(ix), length(iy)		# size of windowed arrays

# Create sponge layer mask:

grid = TwoDGrid(dev; nx=Nx, ny=Ny, Lx, Ly)	   # create dummy grid to define sponge layer function M
x, y = gridpoints(grid)
M = @. -M₀ * (exp(-(x+Lx/2)^2/(0.01*Lx)^2) + exp(-(x-Lx/2)^2/(0.01*Lx)^2)) # predefine M to reduce function evaluations

# Create cutting function:

function f_cut(x, y, x_c, y_c)
	r = @. sqrt((x - x_c)^2 + (y - y_c)^2)
	mask = r .< r_c
	return @. mask + exp(-(r - r_c)^2 / σ^2) * !mask
end

# Define functions:

fstring(num) = string(round(num,sigdigits=8))
istring(num) = string(Int(num))
to_CPU(f) = device_array(CPU())(f)

function calcF1!(Fh, sol, t, clock, vars, params, grid)	# 1-layer sponge layer forcing
	@. vars.qh = sol
	ldiv!(vars.q, grid.rfftplan, vars.qh)
	mul!(Fh, grid.rfftplan, vars.q .* M )
	return nothing
end

function calcF2!(Fh, sol, t, clock, vars, params, grid) # 2-layer sponge layer forcing
	@. vars.qh = sol
	invtransform!(vars.q, vars.qh, params)
	fwdtransform!(Fh, vars.q .* M, params)
	return nothing
end

function create_file_all(grid, Nl, filename, β, K)
	if Nl == 1
		nccreate(filename, "psi", "x", grid.x, "y", grid.y, "t", t_all)
  		nccreate(filename, "Q", "x", grid.x, "y", grid.y, "t", t_all)
  		nccreate(filename, "E", "t")
  		nccreate(filename, "Z", "t")
	else
		nccreate(filename, "psi_1", "x", grid.x, "y", grid.y, "t", t_all)
  		nccreate(filename, "psi_2", "x", grid.x, "y", grid.y, "t", t_all)
  		nccreate(filename, "Q_1", "x", grid.x, "y", grid.y, "t", t_all)
  		nccreate(filename, "Q_2", "x", grid.x, "y", grid.y, "t", t_all)
  		nccreate(filename, "KE_1", "t")
  		nccreate(filename, "KE_2", "t")
  		nccreate(filename, "PE", "t")
	end
	ncputatt(filename," ", Dict("R" => R[1:Nl], "U" => U, "a" => a, "b" => β, "K" => K, "C" => C))
	return nothing
end

function create_file_window(grid, Nl, filename, β, K)
	if Nl == 1
		nccreate(filename, "Q", "x", grid.x[ix], "y", grid.y[iy], "t", t_window)
	else
		nccreate(filename, "Q_1", "x", grid.x[ix], "y", grid.y[iy], "t", t_window)
  		nccreate(filename, "Q_2", "x", grid.x[ix], "y", grid.y[iy], "t", t_window)
	end
	ncputatt(filename," ", Dict("R" => R[1:Nl], "U" => U, "a" => a, "b" => β, "K" => K, "C" => C))
	return nothing
end

function get_centre(problem, grid, Nl)
	if Nl == 1
		_, imax = findmax(problem.vars.q); _, imin = findmin(problem.vars.q);
	else
		_, imax = findmax(problem.vars.q[:,:,1]); _, imin = findmin(problem.vars.q[:,:,1]);
	end
	return (grid.x[imax[1]]+grid.x[imin[1]])/2
end

function save_field_data(problem, grid, diagnostics, Nl, filename, i, iter)
	if Nl == 1
		ψ, q = reshape(to_CPU(problem.vars.ψ),(Nx, Ny, 1)), reshape(to_CPU(problem.vars.q),(Nx, Ny, 1))
		ncwrite(ψ, filename, "psi", start = [1, 1, i+1], count = [Nx, Ny, 1])
		ncwrite(q, filename, "Q", start = [1, 1, i+1], count = [Nx, Ny, 1])

		E, Z = diagnostics
		ncwrite([Lx*Ly*E.data[i+1]],filename, "E", start = [i+1], count = [1])
		ncwrite([Lx*Ly*Z.data[i+1]],filename, "Z", start = [i+1], count = [1])

		qmax, imax = findmax(problem.vars.q); qmin, imin = findmin(problem.vars.q); Δq = qmax + qmin
		print("Iteration: " * istring(iter) * ", t = " * fstring(problem.clock.t) * ", (E, Z) = (" 
    			* fstring(Lx*Ly*E.data[i+1]) * ", " * fstring(Lx*Ly*Z.data[i+1]) * "), (x₀, Δq) = (" *
			fstring((grid.x[imax[1]]+grid.x[imin[1]])/2) * ", " * fstring(Δq) * ") \n")
	else
		ψ₁, q₁ = reshape(to_CPU(problem.vars.ψ[:,:,1]),(Nx, Ny, 1)), reshape(to_CPU(problem.vars.q[:,:,1]),(Nx, Ny, 1))
		ψ₂, q₂ = reshape(to_CPU(problem.vars.ψ[:,:,2]),(Nx, Ny, 1)), reshape(to_CPU(problem.vars.q[:,:,2]),(Nx, Ny, 1))
		ncwrite(ψ₁, filename, "psi_1", start = [1, 1, i+1], count = [Nx, Ny, 1])
		ncwrite(q₁, filename, "Q_1", start = [1, 1, i+1], count = [Nx, Ny, 1])
		ncwrite(ψ₂, filename, "psi_2", start = [1, 1, i+1], count = [Nx, Ny, 1])
		ncwrite(q₂, filename, "Q_2", start = [1, 1, i+1], count = [Nx, Ny, 1])

		E = diagnostics[1]
		ncwrite([Lx*Ly*E.data[i+1][1][1]],filename, "KE_1", start = [i+1], count = [1])
		ncwrite([Lx*Ly*E.data[i+1][1][2]],filename, "KE_2", start = [i+1], count = [1])
		ncwrite([Lx*Ly*E.data[i+1][2]],filename, "PE", start = [i+1], count = [1])
		
		qmax, imin = findmax(problem.vars.q[:,:,1]); qmin, imax = findmin(problem.vars.q[:,:,1]); Δq = qmax + qmin
		print("Iteration: " * istring(iter) * ", t = " * fstring(problem.clock.t) * ", E = " 
			* fstring(Lx*Ly*(E.data[i+1][1][1]+E.data[i+1][1][2]+E.data[i+1][2])) * ", (x₀, Δq) = (" *
			fstring((grid.x[imax[1]]+grid.x[imin[1]])/2) * ", " * fstring(Δq) * ") \n")
	end
	return nothing
end

function save_window_data(problem, Nl, filename, i, iter)
	Q = to_CPU(problem.vars.q)
	S_w = (Nx_w, Ny_w, 1)

	if Nl == 1
		_, i_m = findmax(Q)
  		ix_c = Int.(mod.((-(Nx/2-1):(Nx/2)) .+ i_m[1] .+ 1, Nx))
  		ix_c[ix_c .== 0] .= Nx

  		Q_s = reshape(Q[ix_c, :][ix, iy], S_w)
		ncwrite(Q_s, filename, "Q", start = [1, 1, i+1], count = [Nx_w, Ny_w, 1])
	else
		Q_1, Q_2 = Q[:, :, 1], Q[:, :, 2]

		_, i_m = findmax(Q_1)
		ix_c = Int.(mod.((-(Nx/2-1):(Nx/2)) .+ i_m[1] .+ 1, Nx))
		ix_c[ix_c .== 0] .= Nx

  		Q1_s, Q2_s = reshape(Q_1[ix_c, :][ix, iy], S_w), reshape(Q_2[ix_c, :][ix, iy], S_w)
		ncwrite(Q1_s, filename, "Q_1", start = [1, 1, i+1], count = [Nx_w, Ny_w, 1])
		ncwrite(Q2_s, filename, "Q_2", start = [1, 1, i+1], count = [Nx_w, Ny_w, 1])
	end
	print("Iteration: " * istring(iter) * ", t = " * fstring(problem.clock.t) * ", "
		* fstring(problem.clock.t/T * 100) * "% complete \n")
	return nothing
end
	

# Loop through all cases:

for ic in 1:Nc

	# Create problem:

	if Layers[ic] == 1
		prob = SingleLayerQG.Problem(dev; nx=Nx, ny=Ny, Lx, Ly, β=β₁[ic], U = -C, deformation_radius=R[1],
			dt=Δt, stepper, calcF=calcF1!, aliased_fraction)
	else
		prob = MultiLayerQG.Problem(2, dev; nx=Nx, ny=Ny, Lx, Ly, β=β₁[ic], U = -[C, C], H, b = [1, 0], dt=Δt,
			topographic_pv_gradient = (0, β₂[ic]-β₁[ic]), stepper, calcFq=calcF2!, aliased_fraction)
	end

	# Set initial condition (modon plus random noise):

	κ = @.sqrt(prob.grid.Krsq)
	
	if IC_type[ic] == 0

		if Layers[ic] == 1
			β, AL, R₁ = β₁[ic], ActiveLayers[1], R[1]
		else
			β, AL, R₁ = [β₁[ic], β₂[ic]], ActiveLayers, R
		end
		_, q₀, K = CreateModonLQG(prob.grid, 10, U, a, R₁, init_β*β, AL; K₀=4*ones(sum(ActiveLayers)), tol=1e-10)
		q₀ .+= 1e-6*Nx*irfft(device_array(dev)(exp.(im*2π*randn(Int(Nx/2+1), Int(Nx)))).*(κ.>κ₁).*(κ.<κ₂), Nx)

	else
		
		IC_dat, IC_i = IC_pref * IC_file[ic] * ".nc", IC_type[ic]

		if Layers[ic] == 1
			Q = ncread(IC_dat, "Q")[:, :, IC_i]
			β, q₀, K = β₁[ic], device_array(dev)(Q), 0
		else
			Q₁, Q₂ = ncread(IC_dat, "Q_1")[:, :, IC_i], ncread(IC_dat, "Q_2")[:, :, IC_i]
			β, q₀, K = [β₁[ic], β₂[ic]], device_array(dev)(cat(Q₁, Q₂; dims = 3)), [0, 0]
		end

	end

	if Layers[ic] == 1
		q₀ = reshape(q₀, Nx, Ny)
		SingleLayerQG.set_q!(prob, q₀)
	else
		MultiLayerQG.set_q!(prob, q₀)
	end

	# Define diagnostics:

	if Layers[ic] == 1
		E = Diagnostic(SingleLayerQG.energy, prob; nsteps=Nt, freq=Int(Nt/Ns))
		Z = Diagnostic(SingleLayerQG.enstrophy, prob; nsteps=Nt, freq=Int(Nt/Ns))
		diags = [E, Z]
	else
		E = Diagnostic(MultiLayerQG.energies, prob; nsteps=Nt, freq=Int(Nt/Ns))
		diags = [E]
	end

	# Create save files and perform initial saves:

	if save_all == 1

		filename_all = savename * "_" * case_name[ic] * "_all.nc"
		if isfile(filename_all); rm(filename_all); end
		create_file_all(prob.grid, Layers[ic], filename_all, β, K)

		nccreate(filename_all, "M", "x", "y"); ncwrite(to_CPU(M),filename_all, "M")
		
		if Layers[ic] == 1
			nccreate(filename_all, "q0", "x", "y"); ncwrite(to_CPU(q₀),filename_all, "q0")
		else
			nccreate(filename_all, "q0_1", "x", "y"); ncwrite(to_CPU(q₀[:,:,1]),filename_all, "q0_1")
			nccreate(filename_all, "q0_2", "x", "y"); ncwrite(to_CPU(q₀[:,:,2]),filename_all, "q0_2")
		end

		save_field_data(prob, prob.grid, diags, Layers[ic], filename_all, 0, 0)

	end

	if save_window == 1

		filename_window = savename * "_" * case_name[ic] * "_window.nc"
		if isfile(filename_window); rm(filename_window); end
		create_file_window(prob.grid, Layers[ic], filename_window, β, K)
	
		save_window_data(prob, Layers[ic], filename_window, 0, 0)

	end

	# Timestep:

	I = gcdx(Int(Nt/Ns), Int(Nt/Nw))[1]	# Number of iterations per loop

	for i1 in 1:ceil(Nt/I)

		stepforward!(prob, diags, I)		# evolve problem in time
		
		if Layers[ic] == 1			# update solution
			SingleLayerQG.updatevars!(prob)
		else
			MultiLayerQG.updatevars!(prob)
		end

		if i1*I % Int(Nt/Ns) == 0			# save all at specified frequency
			
			x_c = get_centre(prob, prob.grid, Layers[ic])
			q₁ = prob.vars.q .* (f_cut(x .- Lx, y, x_c, 0) .+ f_cut(x, y, x_c, 0) .+ f_cut(x .+ Lx, y, x_c, 0))

			if Layers[ic] == 1
				SingleLayerQG.set_q!(prob, q₁)
			else
				MultiLayerQG.set_q!(prob, q₁)
			end
			
			if save_all == 1
				save_field_data(prob, prob.grid, diags, Layers[ic], filename_all, Int(i1*I*Ns/Nt), i1*I)
			end
		end

		if i1*I % Int(Nt/Nw) == 0			# save window at specified frequency
			if save_window == 1
				save_window_data(prob, Layers[ic], filename_window, Int(i1*I*Nw/Nt), i1*I)
			end
		end
		
	end

end