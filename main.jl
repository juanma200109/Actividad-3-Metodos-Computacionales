# Importar las bibliotecas necesarias
using LinearAlgebra
using DataFrames
using CSV
using Plots
using Dates
using Clustering
using Statistics

# Cargando los dataframe

lines = DataFrame(CSV.File("Clase_3/lines.csv"))
nodes = DataFrame(CSV.File("Clase_3/nodes.csv"))
demandas = DataFrame(CSV.File("Clase_3/Demandas.csv"))
generacion = DataFrame(CSV.File("Clase_3/solardata.csv"))

function  calcular_ynn_yns(lines,nodes)
    """
    Entradas: 
        - lineas: DataFrame
        - nodes: DataFrame
    Salida: 
        - Ynn: Matriz de admitacias diferentes al nodo slack
        - yns: Vector fila de admitacias respecto al nodo slack
    """

    num_nodes = nrow(nodes)
    num_lines = nrow(lines)

    # Creando una matriz NXN de ceros, que sea de complejos
    Ybus = zeros(num_nodes,num_nodes)*1im

    s = string(nodes[nodes.TYPE .== 3, "NAME"][1][2])
    s = parse(Int64, s)

    for k = 1:num_lines
        # nodo de envio
        n1 = string(lines.FROM[k])
        n1= parse(Int64,n1[2:end])
        # nodo de recibo
        n2 = string(lines.TO[k])
        n2= parse(Int64,n2[2:end])
        # Admitancia
        yl = 1/(lines.R[k]+lines.X[k]*1im)
        # Susceptancia
        Bs = lines.B[k]*1im/2
        # Extrayendo el valor del tap
        t = lines.TAP[k]
        if lines.TAP[k] == 0
            Ybus[n1,n1] += yl + Bs # Dentro de la diagonal suma
            Ybus[n1,n2] -= yl # Fuera de la diagonal resta
            Ybus[n2,n1] -= yl # Fuera de la diagonal resta
            Ybus[n2,n2] += yl +  Bs # Dentro de la diagonal suma; el efecto del elemento shunt se tiene en cuenta para añadir el nodo a tierra, sino la matriz posiblemente sera singular
        else
            Ybus[n1,n1] += ((t)^2-t)*yl # Dentro de la diagonal suma
            Ybus[n1,n2] -= (t)*yl # Fuera de la diagonal resta
            Ybus[n2,n1] -= (t)*yl 
            Ybus[n2,n2] += (1-t)*yl
        end
    end
    #  Extrayendo la matriz de admitancias Y_nn (sin el nodo slack)
    Y_nn = Ybus[setdiff(1:num_nodes, s), setdiff(1:num_nodes, s)]
    # Extrayendo la fila de admitancias Y_ns (con respecto al nodo slack)
    Y_ns = Ybus[s+1:end, s]
    return Y_nn, Y_ns, Ybus
end

Y_nn, Y_ns, Ybus = calcular_ynn_yns(lines,nodes)
# display(Y_nn)
# display(Y_ns)
# display(Ybus)

function flujo_punto_fijo(lines, nodes)
    """
    Entradas: 
        - lineas: DataFrame
        - nodes: DataFrame
    Salida: 
        - 
    """

    # Número de nodos
    num_nodes = nrow(nodes)
    # Número de líneas
    num_lines = nrow(lines)

    # Calculando la Ybus
    Y_nn, Y_ns, Ybus = calcular_ynn_yns(lines, nodes)

    # Inicializando los valores de flujo
    Vn = ones(num_nodes-1) + zeros(num_nodes-1)*1im
    # Definiendo el nodo slack
    Vs = 1 + 0*1im

    # Calculando la inversa de Y_nn
    Y_nn_inv = inv(Y_nn)

    # Definiendo las potencias nodales (sin el nodo slack)
    Sn = (nodes.PGEN[2:end] - nodes.PLOAD[2:end]) + (nodes.QGEN[2:end] - nodes.QLOAD[2:end])*1im

    # Tamaño del array tanto para errores como para iteracciones
    max_iter = 4

    # Empezando el proceso iterativo
    # Vector de errores
    errores = zeros(max_iter)
    # Vector de iteraciones
    iteraciones = zeros(max_iter)
    for k in 1:max_iter
        # Almacenando las tensiones de la iteración anterior
        Vn_ant = Vn
        # Expresión de vn para punto fijo
        Vn = Y_nn_inv*(conj.(Sn./Vn) .- Y_ns*Vs)
        # Calculando el error
        error = norm(Vn - Vn_ant)
        errores[k] = error
        # Calculando el número de iteraciones
        iteraciones[k] = k
        # Si el error es menor a 1e-6, se detiene la iteración
        if error < 1e-6
            println("Convergencia alcanzada en la iteración $k")
            break
        end
    end

    display(errores)

    # Creando un grafico de convergencia logaritmico
    P = plot(iteraciones, errores, yscale=:log10, xlabel = "iteraciones", ylabel = "error")
    display(P)
    return Vn
end

################################################################################################################################################

# --- Procesamiento de datos de generación ---
# Suponiendo que 'generacion' es un DataFrame con columnas 'Fecha' y 'Potencia'
# Agrupar los datos por fecha
SolarData_df2 = DataFrame()
for i = 1:288:nrow(generacion)
    SolarData_df2[!, generacion.Fecha[i]]= generacion.Potencia[i:i+287]
end
# --- Exportar datos a CSV ---
# Se guarda el dataframe en un CSV
CSV.write("Clase_3/solardata_matriz.csv", SolarData_df2)

# Convertir SolarData_df2 a una matriz numérica (288 filas × 365 columnas)
matriz_potencia_solar = Matrix(SolarData_df2)

# Se gráfica las curvas de generación por día para los 365 días.
theme(:dark)
p = plot(size=(800, 600), legend=false)
for i in 1:365
    plot!(p, matriz_potencia_solar[:, i], label = false, alpha=0.5)
end
title!("Potencia Generada por día para los 365 días")
xlabel!("Intervalos de 5 min")
ylabel!("Potencia Generada")

# --- Clustering con K-means ---
# Elegir el número de clusters (3: alta, baja e intermedia generación)
k = 3
resultado = kmeans(matriz_potencia_solar, k)  # Aplicar K-means
asignaciones = assignments(resultado)   # Asignación de cada día a un cluster
centroides = resultado.centers          # Centroides de los clusters

# Creando un DataFrame para los centroides
dias_tipicos_df = DataFrame(centroides, Symbol.("DT " .* string.(1:k)))

# Se guarda el dataframe de dias_tipicos en un CSV
CSV.write("Clase_3/dias_tipicos.csv", dias_tipicos_df)

# Crear un gráfico de los centroides
p1 = plot(size=(800, 600), legend=false)
theme(:dark)
for i in 1:k
    plot!(p1, centroides[:, i], label=false, alpha=0.5)  # Graficar cada centroide
end
title!("Centroides de los clusters")
xlabel!("Intervalo (5 min)")
ylabel!("Potencia")

# Combinar ambos gráficos en un layout vertical
p3 = plot(p, p1, layout=(2, 1), size=(800, 600))
display(p3)  # Mostrar el gráfico combinado


# --- Procesamiento de datos de demanda ---
# Suponiendo que 'demandas' es un DataFrame con perfiles de demanda por minuto
# Crear una columna 'grupo' para agrupar cada 5 minutos
demandas.grupo = div.(0:(nrow(demandas)-1), 5) .+ 1
# Calcular el promedio de las columnas (excepto 'grupo') por cada grupo de 5 minutos
demandas_prom = combine(groupby(demandas, :grupo), names(demandas, Not(:grupo)) .=> mean .=> names(demandas, Not(:grupo)))
select!(demandas_prom, Not(:grupo))  # Eliminar la columna auxiliar 'grupo'

# --- Renombrar columnas de demandas_prom según nodos PQ ---
# Suponiendo que 'nodes' es un DataFrame con información de los nodos, incluyendo 'PLOAD' (carga activa)
# Crear una lista de encabezados solo para los nodos con carga (PLOAD != 0)
encabezados = String[]
for i in 1:nrow(nodes)
    if nodes.PLOAD[i] != 0
        push!(encabezados, "N$i")  # Convertir el índice del nodo a string
    end
end

# Renombrar las columnas de 'demandas_prom' usando los encabezados de los nodos PQ
rename!(demandas_prom, names(demandas_prom) .=> encabezados)

# --- Exportar el DataFrame renombrado a CSV ---
# Guardar el DataFrame en un archivo CSV en la carpeta "Clase_3"
CSV.write("Clase_3/perfiles_de_demanda_promediados.csv", demandas_prom)

# --- Graficar las potencias demandadas ---
# Definir la potencia base (100 MVA) para escalar los valores
p_base = 100  # Potencia base de 100 MVA

# Crear un gráfico base con tamaño 800x600, sin leyenda y tema oscuro
p4 = plot(size=(800, 600), legend=false)
theme(:dark)

# Graficar cada columna de 'demandas_prom' (excepto la última en el código original)
for i in 1:ncol(demandas_prom)-1
    plot!(p4, demandas_prom[:, i] .* p_base, label=false, alpha=0.5)  # Escalar por p_base
end

# Añadir título y etiquetas a los ejes
title!("Potencia demandada")
xlabel!("Intervalo (5 min)")
ylabel!("Potencia")

# Mostrar el gráfico
display(p4)

########################### FLUJO CUASI-DINAMICO ###########################

function flujo_cuasi_dinamico(lines, nodes, Demandas_prom, dias_tipicos)

    """
    Entradas: 
        - lineas: DataFrame
        - nodes: DataFrame
        - demandas_prom: DataFrame
        - dias_tipicos: DataFrame
    Salida: 
        - Archivos CSV con tensiones Vn para cada caso e intervalo de tiempo.
    """
    # Número de nodos
    num_nodes = nrow(nodes)
    # Número de líneas
    num_lines = nrow(lines)

    # Calculando la Ybus
    Y_nn, Y_ns, Ybus = calcular_ynn_yns(lines, nodes)

    # Inicializando los valores de flujo
    Vn = ones(num_nodes-1) + zeros(num_nodes-1)*1im
    # Definiendo el nodo slack
    Vs = 1 + 0*1im

    # Calculando la inversa de Y_nn
    Y_nn_inv = inv(Y_nn)

    # Calculando el flujo cuasi DINAMICO
    # Recorriendo los casos de estudio (dias_tipicos)
    for i in 1:ncol(dias_tipicos)
        # Inicializando el DataFrame para las tensiones (Vn)
        Vns = DataFrame()
        # Iterando sobre los intervalos de tiempos de 5 minutos
        for j in 1:nrow(dias_tipicos)
            # Inicializando las tensiones y potencias
            Vn = ones(num_nodes-1) + zeros(num_nodes-1)*1im
            Vs = 1 + 0*1im
            # Potencias netas inyectadas en los nodos no-Slack
            Sn = zeros(num_nodes-1)
            # Calculando potencias netas inyectadas en los nodos no-Slack
            for k in 1:ncol(Demandas_prom)
                nodo = parse(Int64, names(Demandas_prom)[k][2:end]) - 1
                if (nodo + 1) == 34
                    Sn[nodo] = dias_tipicos[j, i] - Demandas_prom[j, k]
                else
                    Sn[nodo] = - Demandas_prom[j, k]
                end
            end
            # Método de punto fijo para calcular tensiones
            for r = 1:4
                Vn = Y_nn_inv * (conj.(Sn ./ Vn) .- Y_ns * Vs)
            end
            # Almacenamiento de resultados en el DataFrame
            nombre_col = "int $(generacion.Hora[j])"
            Vns[!, nombre_col[1:9]] = Vn
        end
        nombre_archivo = "Clase_3/Vn_caso_$i".*".csv"
        CSV.write(nombre_archivo, Vns)
    end
end

flujo_cuasi_dinamico(lines, nodes, demandas_prom, dias_tipicos_df)