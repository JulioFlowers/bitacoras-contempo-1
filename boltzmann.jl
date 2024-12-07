using CSV # para cargar archivos csv
using Plots # para gráficar
using DataFrames # para convertir un CSV en DataFrame
using LsqFit # Ajustes no lineales
using Statistics #usar promedio
using LaTeXStrings # Para escribir mates
using LinearAlgebra # manipular matrices

# Tipografía para coincidir las gráficas con el documento
Plots.default(fontfamily=("Computer Modern"))

# Se lee el archivo CSV y se almacena en un DataFrame
data = CSV.read("boltzamannesmalt.csv", DataFrame)

at4 = data.AT4
errat4 = data.errat4
pot = data.Potencia
errpot = data.errpot

# Se define la función de modelo (lineal)
function model(x, a)
    return a[1] .* x
end

    a0 = [0.0]
    # Se realiza el ajuste utilizando el algoritmo de Levenberg-Marquardt
    fit_result = curve_fit(model, at4, pot, a0)

    # Se calcula la función de ajuste para los datos
    y_fit = model(at4, fit_result.param)

    # Se extrae la matriz de covarianza
    cov_matrix = estimate_covar(fit_result)

    # Se extraen los elementos diagonales (varianzas) de la matriz de covarianza
    variances = diag(cov_matrix)

    # Se calculan los errores estándar como la raíz cuadrada de las varianzas
    std_errors = sqrt.(variances)

    # Se imprimen los errores estándar de los coeficientes
    for (j, std_err) in enumerate(std_errors)
        println("Error estándar del coeficiente $(j): ", std_err)
    end

    # Se imprimen los valores de los coeficientes del ajuste
    a_fit = fit_result.param[1]
    println("Coeficiente ajustado a: ", a_fit)

    # Calcular el coeficiente de determinación R^2
    ss_res = sum((pot.- y_fit) .^ 2)
    ss_tot = sum((pot.- mean(pot)) .^ 2)
    r_squared = round(1 - (ss_res / ss_tot), digits=3)
    println("Coeficiente de determinación R^2 : ", r_squared)

    # Crear la gráfica
    img = plot(
        at4, pot,
        yerror=hcat(errpot...),
        xerror=hcat(errat4...),
        label="medición",
        seriestype=:scatter,
        xlabel=L"$AT^{4}\ [m^{2}K^{4}]$",
        ylabel=L"$ P\ [W] $",
        title=L"Determinación $\epsilon$, esmalte acrílico negro.",
        linecolor=:blue,
        dpi=620,
        legend=:outerbottom
    )

    plot!(img, at4, model(at4, fit_result.param),
        linecolor=:purple,
        label=L"ajuste lineal $R^2\ =\ %$(r_squared)$ 
                $\epsilon \sigma^{*}  = %$(round(a_fit, digits =10)),\ \delta\ \epsilon \sigma^{*}  = %$(round(std_errors[1], digits =10)) $"
    )

    # Mostrar la gráfica
    display(img)

    # Exportar la gráfica
    savefig(img, "epsilon.png")
