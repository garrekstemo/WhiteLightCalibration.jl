using DelimitedFiles
using CSV
using DataFrames
using DataInterpolations
using CairoMakie
using GLMakie

poly(λ, p) = @. -p[1] + p[2] * λ - p[3] * λ^2 + p[4] * λ^3

"""
    calibrate(λs, time, raw, shift_func)

Calibrate the time axis of the raw data by shifting each row in the 2D matrix `raw` by the
corresponding value in `shift_func`.
The function `shift_func` is a function of the wavelength `λ` and returns the shift in fs.
"""
function calibrate(λs, time, raw, interp_func)
    data = copy(raw)

    for (i, λ) in enumerate(λs)
        shift_by = interp_func(λ)  # Get the time shift for this wavelength

        # Iterate over each time step
        for (j, t) in enumerate(time)
            new_time = t + shift_by

            # If the new time is outside the range of the raw data, skip it
            # set the value to 0 for clarity in the plot
            if new_time < first(time) || new_time > last(time)
                data[j, i] = 0.0
                continue
            end

            # Interpolate the raw data to find the value at the new time point
            itp = LinearInterpolation(raw[:, i], time)
            data[j, i] = itp(new_time)
        end
    end

    return data
end


# Load raw data
# truncate raw data to account for shorter wavelength range in both λs_raw and calibration data
raw = readdlm("raw/CCDABS_230609_190026.lvm", skipstart=2)[:, 1:end-2]
λs_raw = readdlm("raw/CCDspec_230510_173207.txt", skipstart=1)[1:end-1, 1]
λs600 = DataFrame(CSV.File("calibration/600_300.txt", header=true))[!, 1]
λs650 = DataFrame(CSV.File("calibration/650_300.txt", header=true))[!, 1]
fs = range(-4000, stop=40000, length = size(raw, 1))

# Create interpolation functions for the calibration data
poly600 = poly(λs600, [30780, 112.12, 0.1383, 5.9476e-5])
poly650 = poly(λs650, [13720, 52.688, 0.067453, 2.9517e-5])

interp600 = LinearInterpolation(poly600, λs600)
interp650 = LinearInterpolation(poly650, λs650)

# Check that the interpolation functions are correct
λs_raw[1]
interp600(λs_raw[1])


# Make test data
test = fill(0.8, 329, 2046)
test[50:51, :] .= -0.5
test[end-10:end-9, :] .= -0.5
test[1:2, :] .= -0.5

# Calibrate the raw data
cal = calibrate(λs_raw, collect(fs), raw, interp600);


# Plots
fig = Figure()
DataInspector()
ax1, hm1 = heatmap(fig[1, 1][1, 1],
        fs ./ 1000, λs_raw, raw,
        colormap = :RdBu, colorrange = (-1, 1),
        )

Colorbar(fig[1, 1][1, 2], hm1, label = "Intensity (a.u.)")
ax1.ylabel = "Wavelength (nm)"

ax2, hm2 = heatmap(fig[2, 1][1, 1],
        fs ./ 1000, λs_raw, cal,
        colormap = :RdBu, colorrange = (-1, 1),
        )

Colorbar(fig[2, 1][1, 2], hm2, label = "Intensity (a.u.)")
ax2.xlabel = "Pump delay (ps)"
ax2.ylabel = "Wavelength (nm)"
# lines!(poly600, λs600)

fig
# save("figures/test.pdf", fig, backend=CairoMakie)

##

# Test the `circshift` built-in Julia function

function shift_rows(data)
    newdata = copy(data)
    shift = 1
    for row in eachrow(newdata)
        idx = round(Int, shift)
        row[:] = circshift(row, idx)
        shift = shift * 1.2
    end
    newdata
end


A = rand(1:10, 8, 10)
A[:, 1] = [999 for i in 1:size(A)[1]]
shift_rows(A)