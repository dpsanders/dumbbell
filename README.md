<script type="text/javascript" src="http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=default"></script>

dumbbell
========

En el *repo* están los archivos para obtener el espectro de Lyapunov de la mancuerna, tanto de Python como de Julia. 
La verdad es que el código con Julia quedó mucho mejor, pero aún tengo que agregar algunos comentarios para hacero más legible.

Sólo uasando Julia, hay tres *scripts* necesarios para todos los cálculos:
 - `dumbbell_type.jl` que contiene la definición de un `type` llamado `dumbbell`, así como las funciones necesarias para calcular los tiempos de colisión y las transformaciones en $\vec{v}$ y  $\omega$ en cada choque.
 - `Gram_Schimidt.jl` que es un script muy corto que implementa el proceso de ortogonalización de Gram--Schimdt para un número arbitrario de vectores, de cualquier dimensión. En ese mismo archivo se define la función `Gram_Schimidt` que devuelve una tupla cuya primera entrada son los vectores ortogonales y la segunda los vectores orto*normales*.
 - `dumbbell_lyapunov_functions.jl` contiene las funciones para manipular los vectores de desplazamiento (los $\delta \Gamma$ del espacio tangente) y así poder calcular el espectro de Lyapunov. Probablemente este es el archivo más críptico, porque no he añadido comentarios sobre qué hace cada función, pero básicamente o es el mapeo de la colisión, o su derivada aplicada a los vectores rellevantes.

Finalmente, el archivo `dumbbell_lyapunov_spectrum.jl` es el cálculo del espetro de Lyapunov, para 10,000 condiciones. Al terminar de ejecutarse el programa grafica el espectro como función del tiempo. Un ejemplo de lo que se obtiene es 

![espectro-lyapunov](https://copy.com/x8Vjqhb9xc0faPWW)

Las líneas continuas son los 6 exponentes calculados usando el número de colisiones y las líneas continuas son los exponentes usando el tiempo transucurrido. Para verificar, la suma de exponentes da algo del orden de $10^{-13}$, como debería ocurrir para el sistema conservativo. Como puede verse, hay tres exponentes positivos y tres negativos, esencialmente simétricos.
