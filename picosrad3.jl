using CSV
using DataFrames
using Plots
using GLM  # Para realizar regresión lineal

# Función para analizar todos los picos de un espectro
function analyze_all_peaks(filepath::String; channel_min::Int = 0, channel_max::Int = 1024, sensitivity::Float64 = 0.1)
    # Cargar archivo CSV
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
        right_idx = findfirst(x -> x <= half_max, counts_range[peak_index:end]) + peak_index - 1
        
        # Validar índices y calcular el FWHM
        if isnothing(left_idx) || isnothing(right_idx)
            continue  # Si no se puede calcular el FWHM, pasar al siguiente pico
        end
        
        fwhm = abs(channels_range[right_idx] - channels_range[left_idx])
        
        # Calcular el área bajo la curva (cuentas totales en el rango FWHM)
        area_under_curve = sum(counts_range[left_idx:right_idx])
        
        # Agregar resultados
        push!(results, (PeakChannel=peak_channel, PeakCounts=peak_counts, FWHM=fwhm, Area=area_under_curve))
    end
    
    return results
end

# Regresión lineal entre canales y energías conocidas
function calibrate_energy(channels::Vector{Int}, energies::Vector{Float64})
    df = DataFrame(Channel=channels, Energy=energies)
    model = lm(@formula(Energy ~ Channel), df)
    println("Modelo de calibración:\n", model)
    return model
end

# Aplicar el modelo para predecir energías
function predict_energies(results::DataFrame, model)
    results.Energy = predict(model, DataFrame(Channel=results.PeakChannel))
    return results
end

# Energías conocidas para calibración
known_channels = [200, 250, 300]  # Ejemplo de canales conocidos
known_energies = [661.7, 1173.2, 1332.5]  # Ejemplo de energías conocidas (en keV)

# Calibrar modelo de energía
energy_model = calibrate_energy(known_channels, known_energies)

# Analizar archivos CSV y aplicar calibración
files = ["Det_BGO_Co60_300s_5cm.csv", "Det_CsI_Co60_300s_5cm.csv", "Det_NaI_Co60_300s_5cm.csv"]
all_results = Dict()

for file in files
    # Analizar espectro y obtener picos
    results = analyze_all_peaks(file; channel_min=150, channel_max=350, sensitivity=0.5)
    
    # Predecir energías usando el modelo calibrado
    results = predict_energies(results, energy_model)
    all_results[file] = results

    # Graficar espectro con picos identificados y energías
    df = CSV.read(file, DataFrame)
    p = plot(df.Channel, df.Counts, label="Espectro $(basename(file))", xlabel="Canal", ylabel="Cuentas", title="Espectro con Energías Identificadas", linewidth=2)
    
    for row in eachrow(results)
        scatter!(p, [row.PeakChannel], [row.PeakCounts], color=:red, marker=:circle, label="E = $(round(row.Energy, digits=1)) keV")
    end
    
    display(p)
    savefig(p, "$(basename(file))_energies.png")
    println("Gráfica guardada: $(basename(file))_energies.png")
end

# Mostrar resultados en consola
for (file, results) in all_results
    println("Archivo: $file")
    println(results)
end
