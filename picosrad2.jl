using CSV
using DataFrames
using Plots
using LsqFit

# Modelo lineal para la regresión
model(x, p) = p[1] * x 

global m
global fit


# Función para analizar todos los picos de un espectro
function analyze_all_peaks(filepath::String; channel_min::Int = 0, channel_max::Int = 1024, sensitivity::Float64 = 0.1)
    df = CSV.read(filepath, DataFrame)
    channels = df.Channel
    counts = df.Counts
    
    # Filtrar el rango de canales
    range_indices = findall(channels .>= channel_min .&& channels .<= channel_max)
    if isempty(range_indices)
        println("No se encontraron datos en el rango especificado para $filepath")
        return DataFrame(PeakChannel=Int[], PeakCounts=Int[], FWHM=Float64[], Area=Float64[])
    end
    
    channels_range = channels[range_indices]
    counts_range = counts[range_indices]
    
    # Sensibilidad: umbral mínimo para considerar un pico
    threshold = maximum(counts_range) * sensitivity

    # Encontrar máximos locales por encima del umbral
    peaks = findall(i -> (i > 1 && i < length(counts_range) && 
                          counts_range[i] > counts_range[i-1] && 
                          counts_range[i] > counts_range[i+1] && 
                          counts_range[i] >= threshold), 1:length(counts_range))
    
    # Inicializar resultados
    results = DataFrame(PeakChannel=Int[], PeakCounts=Int[], FWHM=Float64[], Area=Float64[])
    
    # Analizar cada pico
    for peak_index in peaks
        peak_channel = channels_range[peak_index]
        peak_counts = counts_range[peak_index]
        
        # Calcular el half_max
        half_max = peak_counts / 2

        # Encontrar el punto más cercano a la izquierda del pico
        left_idx = findlast(x -> x <= half_max, counts_range[1:peak_index])
        
        # Encontrar el punto más cercano a la derecha del pico
        right_idx = findfirst(x -> x <= half_max, counts_range[peak_index:end])
        
        # Ajustar el índice derecho para el rango completo
        if !isnothing(left_idx) && !isnothing(right_idx)
            right_idx += peak_index - 1
            # Calcular el FWHM
            fwhm = abs(channels_range[right_idx] - channels_range[left_idx])
            
            # Calcular el área bajo la curva (cuentas totales en el rango FWHM)
            area_under_curve = sum(counts_range[left_idx:right_idx])
            
            # Agregar resultados
            push!(results, (PeakChannel=peak_channel, PeakCounts=peak_counts, FWHM=fwhm, Area=area_under_curve))
        end
    end
    
    return results, channels, counts
end

# Función para realizar la regresión lineal
function perform_linear_regression(peaks::DataFrame, energy_values::Vector{Float64})
    x = peaks.PeakChannel
    y = energy_values

    if length(x) != length(y)
        throw(DimensionMismatch("El número de elementos en x y y debe coincidir"))
    end

    # Ajustar un modelo lineal: y = m * x 
    initial_params = [1.0]
    fit = curve_fit(model, x, y, initial_params)

    m = fit.param
    println("Pendiente (m): $m")

    return m, fit
end

# Analizar archivos CSV
files = ["Det BGO cal14_11_fuentes2_Ba133.csv"]
all_results = Dict()

for file in files
    results, channels, counts = analyze_all_peaks(file; channel_min=0, channel_max=800, sensitivity=0.2)
    all_results[file] = results

    println("Picos detectados en $file: ")
    println(results)
    
    println("Selecciona los canales de interés para los picos (por ejemplo: [524, 526, 531]):")
    selected_channels = eval(Meta.parse(readline()))

    if all(channel in results.PeakChannel for channel in selected_channels)
        println("Canales seleccionados correctamente.")
    else
        println("Algunos de los canales seleccionados no están en los picos detectados.")
        continue
    end

    println("Introduce los valores de energía correspondientes a los canales seleccionados:")
    energy_values = eval(Meta.parse(readline()))

    if length(selected_channels) != length(energy_values)
        println("El número de canales seleccionados no coincide con el número de valores de energía proporcionados.")
        continue
    end

    selected_peaks = filter(row -> row.PeakChannel in selected_channels, eachrow(results))
    global m, fit = perform_linear_regression(DataFrame(selected_peaks), energy_values)

    ees = model(channels, fit.param)
    p = plot(ees, counts, label="Calibración", xlabel="Nivel Energético [keV]", ylabel="Cuentas", title="Espectro Ba-133", linewidth=2, color=:pink,legend=:outerbottom, dpi = 320)
    
    ens = model(selected_peaks.PeakChannel, fit.param)
    for (idx, value) in enumerate(ens)
        scatter!(p, [value], [selected_peaks.PeakCounts[idx]], label="$(round(value, digits = 3)) keV", marker=:triangle, color=:red)
    end

    display(p)
    savefig(p, "$(basename(file))_selected_peaks_and_fit.png")
    println("Gráfica guardada: $(basename(file))_selected_peaks_and_fit.png")
    
    # Crear DataFrame con los datos seleccionados, incluyendo PeakCounts
    output_df = DataFrame(
        PeakChannel=selected_peaks.PeakChannel,
        EnergyLevel=ens,
        PeakCounts=selected_peaks.PeakCounts,
        FWHM=selected_peaks.FWHM,
        Area=selected_peaks.Area
    )
    
    # Guardar en CSV
    output_csv = "$(basename(file))_selected_peaks.csv"
    CSV.write(output_csv, output_df)
    println("Datos guardados en: $output_csv")
end

#=
# Mostrar resultados en consola
for (file, results) in all_results
    println("Archivo: $file")
    println(results)
end
=#

# Archivos CSV a procesar
files2 = [ "Det BGO Am241 300s 0cm.csv" , "Det BGO Am241 300s 5cm.csv"]
names2 = [0.0, 0.05]
# Procesar todos los archivos
for (idx, file) in enumerate(files2)

    results, channels, counts = analyze_all_peaks(file; channel_min=0, channel_max=200, sensitivity=0.1)
    all_results[file] = results

    println("Picos detectados en $file: ")
    println(results)
    
    println("Selecciona los canales de interés para los picos (por ejemplo: [524, 526, 531]):")
    selected_channels = eval(Meta.parse(readline()))

    if all(channel in results.PeakChannel for channel in selected_channels)
        println("Canales seleccionados correctamente.")
    else
        println("Algunos de los canales seleccionados no están en los picos detectados.")
        continue
    end

    selected_peaks = filter(row -> row.PeakChannel in selected_channels, eachrow(results))

    ees = model(channels, fit.param)
    p = plot(ees, counts, label="d = $(names2[idx]) m", xlabel="Nivel Energético [MeV]", ylabel="Cuentas", title="Espectro Am-241", linewidth=2, color=:pink, legend=:outerbottom, dpi = 320)
    
    ens = model(selected_peaks.PeakChannel, fit.param)
    for (idx, value) in enumerate(ens)
        scatter!(p, [value], [selected_peaks.PeakCounts[idx]], label="$(round(value, digits = 3)) keV", marker=:triangle, color=:red)
    end

    display(p)
    savefig(p, "$(basename(file))_selected_peaks_and_fit.png")
    println("Gráfica guardada: $(basename(file))_selected_peaks_and_fit.png")
    
    # Crear DataFrame con los datos seleccionados, incluyendo PeakCounts
    output_df = DataFrame(
        PeakChannel=selected_peaks.PeakChannel,
        EnergyLevel=ens,
        PeakCounts=selected_peaks.PeakCounts,
        FWHM=selected_peaks.FWHM,
        Area=selected_peaks.Area
    )
    
    # Guardar en CSV
    output_csv = "$(basename(file))_selected_peaks.csv"
    CSV.write(output_csv, output_df)
    println("Datos guardados en: $output_csv")

end
