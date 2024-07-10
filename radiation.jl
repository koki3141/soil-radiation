using Plots
using Unitful
using Unitful.DefaultSymbols
using Dates
using Statistics
using DataFrames


function solar_radiation(time::Unitful.Time)
    hours = uconvert(u"hr", time) % 24u"hr"
    if 6u"hr" <= hours <= 18u"hr"
        θ = 2π * (hours + 12u"hr") / 24u"hr"
        R_s = 800u"W/m^2" * cos(θ)
    else
        R_s = 0u"W/m^2"
    end
    return R_s
end

function antoine_equation(
    soil_temperature::Unitful.Temperature,
    A = 8.07131 + log10(101325/760),
    B = 1730.63,
    C = 233.425)
    T_s = ustrip(uconvert(u"°C", soil_temperature))
    p = (10^(A - B / (C + T_s)))u"Pa"
    return p
end


function atmospheric_emissivity(
    atmospheric_temperature::Unitful.Temperature,
    soil_temperature::Unitful.Temperature)
    p_w = antoine_equation(soil_temperature)
    ϵ_a = 1.24 * (ustrip(uconvert(u"hPa",p_w)) /  ustrip(atmospheric_temperature))^(1/7)
    return ϵ_a
end

# 0.97は水の割合が0.3以上の時の値
function quantity_of_heat_by_radiation(
    solar_radiation, 
    soil_temperature::Unitful.Temperature,
    atmospheric_temperature::Unitful.Temperature,
    side,
    area = side^2,
    α_s = 0.075,
    ϵ_s = 0.97,
    σ = 5.67e-8 * u"W/(m^2 * K^4)")
    ϵ_a = atmospheric_emissivity(atmospheric_temperature, soil_temperature)
    atmospheric_temperature_kelvin = uconvert(u"K", atmospheric_temperature)
    soil_temperature_kelvin = uconvert(u"K", soil_temperature)
    quantity_of_heat = solar_radiation * (1 - α_s) * area + σ * ϵ_s * (ϵ_a * (atmospheric_temperature_kelvin^4) - soil_temperature_kelvin^4) * area * 5
    return quantity_of_heat
end



function soil_heat_capacity(
    θ_w, # soil volumetric water content
    ρ_b::Unitful.Quantity = 950u"kg/m^3",  # soil bulk density
    c_o::Unitful.Quantity = 1920u"J/kg/K",  # specific heat of soil organic matter
    c_m::Unitful.Quantity = 900u"J/kg/K",  # specific heat of soil mineral
    C_w::Unitful.Quantity = 4.186e6u"J/m^3/K",  # volumetric heat capacity of water
    f_o_m::Float64 = 0.0745  # soil organic matter fraction by mass
    )::Unitful.Quantity

    soil_volumetric_heat_capacity = ρ_b * (c_o * f_o_m+ c_m * (1 - f_o_m))

    C_s = soil_volumetric_heat_capacity + C_w * θ_w

    return C_s
end

function delta_temperature(
    Q,
    θ_w,
    Δt,
    side,
    volume::Unitful.Quantity = side^3 # soil volume
    )::Unitful.Temperature
    C_s = soil_heat_capacity(θ_w)
    ΔT = Q / (C_s * volume )  * Δt
    return ΔT
end

function calculate_soil_temperature(
    time::Unitful.Time, 
    T_s::Unitful.Temperature,
    θ_w, 
    Δt,
    side = 0.01u"m")
    R_s = solar_radiation(time)
    atmospheric_temperature = atmoshperic_temperature(time)
    quantity_of_heat = quantity_of_heat_by_radiation(R_s, T_s, atmospheric_temperature, side)
    ΔT = delta_temperature(quantity_of_heat, θ_w, Δt, side)
    T_s += ΔT
    return T_s
end

function calculate_radiation_budget(times, θ_w, Δt, first_soil_temperature = 303.15u"K")
    soil_temperatures = []
    for time in times
        first_soil_temperature = calculate_soil_temperature(time, first_soil_temperature, θ_w, Δt)
        push!(soil_temperatures, first_soil_temperature)
    end
    return soil_temperatures
end


function atmoshperic_temperature(time::Unitful.Time)
    hours = uconvert(u"hr", time)
    θ = 2π * (hours + 10u"hr") / 24u"hr"
    atmoshperic_temperature = uconvert(u"K", 30u"°C") + 10u"K"*cos(θ)
    return atmoshperic_temperature
end


function calculate_daily_average_temperature(T_s_array, Δt)
    daily_seconds = Int(ustrip(uconvert(u"s", Day(1)) / Δt))  # 1日の秒数
    T_s_celsius = [ustrip(uconvert(u"°C", T_s)) for T_s in T_s_array]
    daily_stats = DataFrame(Day = Int[], T_Average_°C = Float64[], T_Min_Hour = Float64[], T_min_°C = Float64[], T_Max_Hour = Float64[], T_Max_°C = Float64[])
    
    day_counter = 1
    for i in 1:daily_seconds:(length(T_s_array)-1)
        maxid  = argmax(T_s_celsius[i:i+daily_seconds])
        max = maximum(T_s_celsius[i:i+daily_seconds])
        minid = argmin(T_s_celsius[i:i+daily_seconds])
        min = minimum(T_s_celsius[i:i+daily_seconds])
        average = mean(T_s_celsius[i:i+daily_seconds])
        push!(daily_stats, (day_counter, average, (minid%daily_seconds)/3600, min, (maxid%daily_seconds)/3600, max))
        day_counter += 1
    end
    show(stdout, "text/plain", daily_stats)
end

function plot_radiation_budget(T_s_array, times)
    T_s_celsius = [ustrip(uconvert(u"°C", T_s)) for T_s in T_s_array]
    times_day = [ustrip(uconvert(u"hr", time)) for time in times] / 24
    plot(times_day,
        T_s_celsius,
        label = "Soil Temperature",
        xlabel = "Day",
        ylabel = "Temperature (°C)",
        title = "Soil Temperature over $(maximum(times_day)) day",
        xticks = 0:1:maximum(times_day),
        legend = true,
        grid = true)
end

function plot_atmospheric_temperature(times)
    atmospheric_temperatures = [atmoshperic_temperature(time) for time in times]
    atmospheric_temperatures_celsius = [ustrip(uconvert(u"°C", temp)) for temp in atmospheric_temperatures]
    times_day = [ustrip(uconvert(u"hr", time)) for time in times] / 24

    plot(times_day,
        atmospheric_temperatures_celsius,
        label = "Atmospheric Temperature",
        xlabel = "day",
        ylabel = "Temperature (°C)",
        title = "Atmospheric Temperature over $(maximum(times_day)) day",
        xticks = 0:1:maximum(times_day),
        legend = true,
        grid = true)
end

function plot_radiation(times)
    radiation_values = [ustrip(solar_radiation(time)) for time in times]
    times_day = [ustrip(uconvert(u"hr", time)) for time in times] / 24

    plot(times_day,
        radiation_values,
        label = "Radiation",
        xlabel = "day",
        ylabel = "Radiation (W/m^2)",
        title = "Radiation over $(maximum(times_day)) day",
        xticks = 0:1:maximum(times_day),
        legend = true,
        grid = true,
        left_margin = 10Plots.mm)
end


function main()
    final_time = Day(20)
    Δt = 1u"s"
    times =  0u"s":Δt:uconvert(u"s", final_time)
    soil_volumetric_water_content = 0.3
    is_plot = false
    soil_temperature = calculate_radiation_budget(times, soil_volumetric_water_content, Δt)
    calculate_daily_average_temperature(soil_temperature, Δt)

    if (is_plot)
        # 各プロットを生成
        p1 = plot_radiation_budget(soil_temperature, times)
        p2 = plot_atmospheric_temperature(times)
        p3 = plot_radiation(times)

        # レイアウトにまとめる
        plot(p1, p2, p3, layout = (3,1), size = (800, 1200))
    end
end

main()
