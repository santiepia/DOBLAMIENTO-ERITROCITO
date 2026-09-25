Recordar que se le bajó la precision a las integrales, lo mejor es tenerlo en     h = 0.01 

Orden con el que se deben correr los archivos:

1. main.f90
2. filtro.py
3. generar_energias.py
4. grafica.py

Tenemos experiemntalmente los valores de eritrocitos de area y volumen:

thin cell area: 168, volumen:128 ; 

average cell area: 137,volumen: 107 ;

thick cell area: 110, volumen: 80


Se detecta que los valores para obtener los valores reales de un eritrocito( o sea, 110 de area y 90 de volumen) es hacer intervalos de aproxiamdamente 2.0 y 3.0 para los valores de A, B y C.

PROCEDIMIENTO CORRER EL CODIGO: 

1. main.f90
2. filtro.py
3. grafica.py


Tambien se hizo un filtro de las figuras para que la diferencia entre el area superficial y el volumen esté entre 30 y 40 que es lo que se tiene registrado experimentalmente.


Para pdoer obtener la grafica en Figura 6 del articulo en donde muestra la figura que adopta el eritrocito teniendo los valores dr area superficial y voluemn fijo y permitir variar al parametro B, podemos organizar los valores obtenidos ne filtro en donde se tome la B de menor a mayor y de esa fora podemos calcular la energia de dobamiento de cada valor y hacer la grafica de E vs B que se tiene en el artiulo.

La energia de doblamiento es: E = (D / 2) · ∫ (1/R₁ + 1/R₂)² dA

D es la rigidez de doblamiento, en este caso, es la misma para todas las figuras, por esta razón podemos descartar esta constante para todos debido a que tendran la misma, entonces relativamente tenemos: E = ∫ (1/R₁ + 1/R₂)² dA
 Ahora R1 y R2 las hallamos de la forma:


 R₁ = [1 + (dY/dx)²]^(3/2) / |d²Y/dx²|


 R₂ = Y(x) · √[1 + (dY/dx)²] / |dY/dx|


 Ahora la energia de doblamiento en coordenadas polares es de la forma:

 E = 4π ∫₀^{x_max} [ (1/R₁ + 1/R₂)² · Y(x) · √(1 + (dY/dx)²) ] dx

 -2π viene de integrar en el ángulo azimutal (rotación)

 -Otro 2 viene de integrar de -x_max a x_max (simetría par)


 En el archivo fortran, no es necesario iterar para valores de A y C negativos, debido a que en la ecuacion tomada todos los valores A y B tienen potencia par, así que si se toma negativo pues dará lo mismo, por otro lado, para el parametro B tampoco tomamoss negativos debido a que B solo es la amplitud de la ecuacion 





 