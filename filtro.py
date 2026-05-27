# filtrar_simple.py - Sin pandas, solo Python estándar

VOLUMEN_MIN = 102.0
VOLUMEN_MAX = 113.0
AREA_MIN = 132.0
AREA_MAX = 142.0

contador_total = 0
contador_filtrado = 0

with open('resultados.txt', 'r') as f_in:
    with open('filtro.txt', 'w') as f_out:
        # Escribir encabezado
        f_out.write('A\tB\tC\tb_limite\tVolumen\tArea\n')
        
        for linea in f_in:
            # Saltar líneas de comentario
            if linea.startswith('#') or linea.startswith('=') or linea.strip() == '':
                continue
            
            # Leer valores (formato: columnas separadas por espacios)
            datos = linea.split()
            if len(datos) >= 6:
                try:
                    A = float(datos[0])
                    B = float(datos[1])
                    C = float(datos[2])
                    b_limite = float(datos[3])
                    Volumen = float(datos[4])
                    Area = float(datos[5])
                    
                    contador_total += 1
                    DIFERENCIA_MINIMA = 30.0
                    DIFERENCIA_MAXIMA = 40.0
                    diferencia_area_volumen = abs(Area - Volumen)
                    # Aplicar filtro
                    if VOLUMEN_MIN <= Volumen <= VOLUMEN_MAX and \
                       AREA_MIN <= Area <= AREA_MAX and DIFERENCIA_MINIMA <= diferencia_area_volumen <= DIFERENCIA_MAXIMA:   #se le agrega un filtro para que la diferencia entre el area superficiall y el volumen sea acorde a lo que se visualiza experimentalmente, esto no s ayuda a filtrar otras figuras que no estaría acordes con el modelo
                        f_out.write(f"{A:.6f}\t{B:.6f}\t{C:.6f}\t{b_limite:.6f}\t{Volumen:.6f}\t{Area:.6f}\n")
                        contador_filtrado += 1
                        
                except ValueError:
                    pass  # Ignorar líneas con formato incorrecto

print(f"Total de combinaciones procesadas: {contador_total}")
print(f"Combinaciones que cumplen el filtro: {contador_filtrado}")
print("Resultados guardados en 'filtro_simple.txt'")