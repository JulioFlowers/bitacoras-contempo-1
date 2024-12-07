
using CSV
using DataFrames
using StatsBase
using Plots

pk_min = 180
pk_max = 300

# Función para analizar un espectro y calcular el FWHM y el área bajo la curva en un rango específico
function analyze_spectrum(filepath::String, peak_min::Int, peak_max::Int, channel_min::Int = 0, channel_max::Int = 1024)
    # Cargar archivo CSV
    df = CSV.read(filepath, DataFrame)
    channels = df.Channel
    counts = df.Counts
    
    # Filtrar el rango de canales
    range_indices = findall(channels .>= channel_min .&& channels .<= channel_max)
    if isempty(range_indices)
        println("No se encontraron datos en el rango especificado para $filepath")
        return (filepath, NaN, NaN, NaN, NaN, NaN)
    end
    
    channels_range = channels[range_indices]
    counts_range = counts[range_indices]
    
    # Filtrar el rango donde se buscará el pico
    peak_indices = findall(channels_range .>= peak_min .&& channels_range .<= peak_max)
    if isempty(peak_indices)
        println("No se encontró un pico en el rango especificado para $filepath")
        return (filepath, NaN, NaN, NaN, NaN, NaN)
    end
    
    # Encontrar el pico en el rango especificado
    peak_index_local = argmax(counts_range[peak_indices])
    peak_index = peak_indices[peak_index_local]
    peak_channel = channels_range[peak_index]
    peak_counts = counts_range[peak_index]
    
    # Calcular el half_max
    # Calcular el half_max y encontrar los puntos más cercanos
    half_max = peak_counts / 2

    # Encontrar el punto más cercano a la izquierda del pico en el rango
    left_idx = findmax([(half_max - counts_range[i]) <= 0 ? abs(half_max - counts_range[i]) : Inf for i in 1:peak_index])[2]

    # Encontrar el punto más cercano a la derecha del pico en el rango
    right_idx = findmin([(half_max - counts_range[i]) <= 0 ? abs(half_max - counts_range[i]) : Inf for i in peak_index:length(counts_range)])[2] + peak_index - 1
    
    # Validar índices y calcular el FWHM en valor absoluto
    if isnothing(left_idx) || isnothing(right_idx)
        println("No se encontró el FWHM correctamente para $filepath")
        return (filepath, peak_channel, peak_counts, NaN, NaN, peak_index)
    end
    fwhm = abs(channels_range[right_idx] - channels_range[left_idx])
    
    # Calcular el área bajo la curva (cuentas totales en el rango FWHM)
    area_under_curve = sum(counts_range[left_idx:right_idx])
    
    # Retornar los resultados
    return (filepath, peak_channel, peak_counts, fwhm, area_under_curve, peak_index)
end

# Analizar los tres archivos CSV en el rango de canales 150 a 350, buscando el pico en un rango específico
#files = ["Det BGO Co60 300s 0cm.csv", "Det CsI Co60 300s 0cm.csv", "Det NaI Co60 300s 0cm.csv"] #distancia de 0c m
files = ["Det BGO Co60 300s 5cm.csv", "Det CsI Co60 300s 5cm.csv", "Det NaI Co60 300s 5cm.csv"] #distancia de 5 cm
results = [analyze_spectrum(file, 0, 1024, pk_min, pk_max) for file in files]

# Mostrar resultados para cada archivo
println("Resultados por archivo:")
for (filepath, peak_channel, peak_counts, fwhm, area_under_curve, _) in results
    println("Archivo: $filepath")
    println("  Canal del pico: $peak_channel")
    println("  Cuentas en el pico: $peak_counts")
    println("  FWHM: $fwhm")
    println("  Área bajo la curva: $area_under_curve\n")
end

# Comparar eficiencia y resolución
best_resolution = argmin(map(x -> x[4], results))  # índice del menor FWHM
best_efficiency = argmax(map(x -> x[5], results))  # índice del mayor área

println("Comparación final:")
println("Archivo con mejor resolución (menor FWHM): $(results[best_resolution][1])")
println("Archivo con mayor eficiencia (mayor área bajo la curva): $(results[best_efficiency][1])")

img = plot()
# Visualizar todos los espectros con indicadores de pico en el rango
for (filepath, peak_channel, peak_counts, fwhm, area, peak_index) in results
    
    

    df = CSV.read(filepath, DataFrame)
    plot!(img, df.Channel, df.Counts, label=basename(filepath))
    
    # Añadir un marcador en el pico dentro del rango
    if !isnan(peak_channel) && !isnan(peak_counts)
        scatter!(img, [peak_channel], [peak_counts], color=:red, marker=:circle, label="Pico - $(basename(filepath))")
    end

    
end

savefig(img, "identific2.png")
display(img)

xlabel!("Canal")
ylabel!("Cuentas")
title!("Espectros Comparativos (Canales 150 a 350)")
