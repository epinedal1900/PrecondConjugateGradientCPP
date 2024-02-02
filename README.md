# Método gradiente conjugado precondicionado en c++

Este documento proporciona instrucciones para configurar tu entorno de desarrollo y ejecutar un programa en C que implementa el método de Gradiente Conjugado Precondicionado para resolver sistemas de ecuaciones lineales.

## Instalación de VSCode y MinGW

- **VSCode**: Sigue este [video tutorial](https://www.youtube.com/watch?v=yhTw84ivFUk) para instalar Visual Studio Code y configurarlo para el desarrollo en C/C++.
- **MinGW**: Descarga e instala MinGW desde [SourceForge](https://sourceforge.net/projects/mingw/) para compilar y ejecutar tu código en C.

### Notas:

- **Verifica la Ruta de MinGW**: si instalaste MinGW en una ruta diferente (mi ruta fue `C:\MinGW`), necesitas ajustar los archivos de configuración preexistentes en el repositorio (`launch.json`, `tasks.json`, `c_cpp_properties.json`) para reflejar la ruta correcta.
- **Reinicio del Sistema**: reinicia tu PC después de la instalación y la actualización del PATH para que se apliquen los cambios.

## Configuración y ejecución del programa

Los archivos de configuración necesarios ya están en este repositorio. Asegúrate de que la ruta especificada para MinGW en `c_cpp_properties.json` coincide con tu configuración local.

Para **ejecutar el programa**, abre el archivo `main.cpp` en VSCode, selecciona `Ejecutar > Ejecutar sin depuración` o presiona `Ctrl+F5`. Esto compilará y ejecutará el programa.

### Exploración del código

El código de este repositorio es un intento de aislar la función `ConjugateGradientIncompleteCholesky` que se encuentra en `source/Tools/MatSolver/main.cpp`,del repositorio original de FEMT. Puedes encontrar el [repositorio original aquí](http://personal.cimat.mx:8181/~miguelvargas/FEMT/FEMT-beta36.tar.xz), y revisar la [documentación aquí](http://personal.cimat.mx:8181/~miguelvargas/FEMT/#Documentation). Para examinar las funciones auxiliares, coloca el cursor sobre una función, presiona `Ctrl` y haz clic para acceder a su definición.

Buena suerte!
