using CSV
using DataFrames
using Plots
using LsqFit
using Statistics

# Cargar el archivo CSV
file_path = "caesium/inversocuadradocesio.csv"  # Cambia esto al nombre de tu archivo
data = CSV.read(file_path, DataFrame)

# Extraer las columnas
x = data.distancia[2:end]
y = data.area[2:end]

# Modelo para el ajuste
model(x, p) = p[1] ./ (x .^ 2)  # y = 1/(x^2)

# Parámetros iniciales para el ajuste
p0 = [0.00000000000000001]

# Realizar el ajuste
fit = curve_fit(model, x, y, p0)

# Coeficiente de determinación R^2
y_fit = model(x, fit.param)
residuals = y - y_fit
ss_res = sum(residuals .^ 2)
ss_tot = sum((y .- mean(y)) .^ 2)
r2 = 1 - ss_res / ss_tot

error_x = fill(0.0005, length(x))
# Gráfica
p = plot(
    x,
    y,
    xerror = error_x,
    color = :pink,
    label = "Datos experimentales",
    xlabel = "Distancia [m]",
    ylabel = "Área",
    title = "Total de cuentas del fotopico del Cs-137 en función de la distancia.",
    dpi=320
)
plot!(p,
    x,
    y_fit,
    label = "Ajuste: R² = $(round(r2, digits=3))",
    lw = 2,
    size = (800, 600),
	margin = 5Plots.mm,
)

display(p)

savefig(p, "inversocesio.png")
# Agregar incertidumbre en x como barras de error

