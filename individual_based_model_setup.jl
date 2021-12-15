using Pkg
Pkg.activate(".")

Pkg.add(name = "PlanktonIndividuals", version = "0.1.9")
Pkg.add(name = "Oceananigans", version = "0.55.0")
pkg"add Random, Statistics, Printf, YAML, Serialization, SeawaterPolynomials"


using Random, Statistics, Printf, YAML, Serialization
using Oceananigans
using Oceananigans.Units: day
using SeawaterPolynomials
using PlanktonIndividuals

#params_path = "params_autotroph.yml" # parameters for autotrophs
params_path = "params_mixotroph.yml" # parameters for mixotrophs
paths_path = "paths_nut_init.yml"
par_path = "Input/par_doc.bin"
nday = 30*9
println(params_path)

Ogrid = RegularRectilinearGrid(size = (25, 1, 50), extent = (100, 4, 200),topology = (Periodic, Periodic, Bounded),)

buoyancy = SeawaterBuoyancy(
    equation_of_state = SeawaterPolynomials.TEOS10EquationOfState(),
    constant_salinity = 35.0 # psu
)

N² = 1e-5 # s⁻²
α = SeawaterPolynomials.thermal_expansion(20, 35, 0, buoyancy.equation_of_state)
g = buoyancy.gravitational_acceleration
ρᵣ = buoyancy.equation_of_state.reference_density
∂T∂z = ρᵣ * N² / (α * g)
bottom_temperature_boundary_condition = GradientBoundaryCondition(∂T∂z)

peak_outgoing_radiation = 200 # Watts / m²
heat_capacity = 3991 # J / kg / ᵒC
reference_density = buoyancy.equation_of_state.reference_density # kg m⁻³
peak_outgoing_flux = peak_outgoing_radiation / (reference_density * heat_capacity)
Qᵇ = α * g * peak_outgoing_flux

@inline diurnal_cycle(t, day) = max(0, - cos(2π * t / day))
@inline nocturnal_cycle(t, day) = max(0, cos(2π * t / day))
@inline outgoing_flux(x, y, t, p) = p.peak * nocturnal_cycle(t, p.day)

surface_temperature_boundary_condition = FluxBoundaryCondition(outgoing_flux, parameters=(day=day, peak=peak_outgoing_flux))

T_bcs = TracerBoundaryConditions(Ogrid, bottom = bottom_temperature_boundary_condition,
                                 top = surface_temperature_boundary_condition)

light_attenuation_scale = 20 # m
surface_solar_insolation = 400 # Watts / m²
surface_solar_temperature_flux = surface_solar_insolation / (reference_density * heat_capacity)

@inline daylight(z, t, λ, day) = exp(z / λ) * diurnal_cycle(t, day)
@inline solar_flux_divergence(z, t, Qᴵ, λ, day) = Qᴵ / λ * daylight(z, t, λ, day)
@inline diurnal_solar_flux_divergence(x, y, z, t, p) =
    max(0, solar_flux_divergence(z, t, p.surface_flux, p.attenuation, p.day))

interior_heating = Forcing(diurnal_solar_flux_divergence,
                                 parameters = (surface_flux = surface_solar_temperature_flux,
                                               attenuation = light_attenuation_scale,
                                               day = day))

Omodel = IncompressibleModel(architecture = CPU(),
                                    grid = Ogrid,
                                 closure = AnisotropicMinimumDissipation(),
                                coriolis = FPlane(f=-1e-4),
                                 tracers = :T, #tracer_names,
                                buoyancy = buoyancy,
                     boundary_conditions = (T=T_bcs,),
                                 forcing = (T=interior_heating,),
)

initial_temperature(x, y, z) = (20 + ∂T∂z * z + ∂T∂z * Ogrid.Lz * 1e-4 * randn() * exp(z / (8 * Ogrid.Δz)))

set!(Omodel, T = initial_temperature,)

wizard = TimeStepWizard(cfl=0.2, Δt=1.0, max_change=1.1, max_Δt=20.0)

simulation = Simulation(Omodel, iteration_interval = 1, stop_time=86400, Δt = wizard)

run!(simulation)

resultpath = PrepRunDir()
phy_grid = read_Ogrids(Omodel.grid, resultpath);
RunParam.nTime = 1440*nday
RunParam.ΔT = 60

# update parameters
phy_params = YAML.load(open(params_path))
update_params!(RunParam.params,phy_params)

paths = YAML.load(open(paths_path))
phy_model = PI_Model(phy_grid, RunParam,
                     PAR = read_IR_input(RunParam.ΔT, phy_grid, par_path),
                     nutrients = load_nut_initials(paths, phy_grid));

Nsimulation = Simulation(Omodel, Δt=20.0, stop_time=86400, iteration_interval=1)

for i in 1:RunParam.nTime
    vel_field =[]
    for j in 1:3
        Nsimulation.stop_time += 20
        run!(Nsimulation)
        u = Omodel.velocities.u.data.parent
        v = Omodel.velocities.v.data.parent
        w = Omodel.velocities.w.data.parent
        vel = PlanktonIndividuals.velocity(u, v, w)
        push!(vel_field,vel)
    end
    vel_itps = (generate_vel_itp(phy_model.grid, vel_field[1]),
                generate_vel_itp(phy_model.grid, vel_field[2]),
                generate_vel_itp(phy_model.grid, vel_field[3]))
    PI_advectRK4!(phy_model, RunParam.ΔT, vel_itps)
    PI_TimeStep!(phy_model, RunParam.ΔT, vel_field[end], resultpath)

    phyts_sp = sort_species(phy_model.individuals.phytos, phy_model.params["P_Nsp"])
    write_species_dynamics(phy_model.t, phyts_sp, resultpath)

end

open(resultpath*"diags_spcs.bin", "w") do io
    serialize(io, phy_model.diags.spcs)
end
open(resultpath*"diags_pop.bin", "w") do io
    serialize(io, phy_model.diags.pop)
end
open(resultpath*"diags_tr.bin", "w") do io
    serialize(io, phy_model.diags.tr)
end
