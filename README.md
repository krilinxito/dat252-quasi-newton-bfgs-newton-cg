# Métodos Cuasi-Newtonianos, Algoritmo BFGS y Método de Newton-CG

**Universidad Mayor de San Andrés (UMSA)** — Facultad de Ciencias Puras y Naturales — Carrera de Informática  
**Materia:** Métodos Numéricos II (`DAT-252`)  
**Docente:** M.Sc. Carlos Mullisaca Choque  
**Estudiantes:** Gómez Mallo, Maximiliano · Medina Balboa, Iver Samuel  
**Fecha de entrega:** La Paz, Bolivia — 25 de Noviembre de 2026

---

## Contenido del Repositorio

| Archivo | Descripción |
| :--- | :--- |
| `Informe_QuasiNewton_BFGS_NewtonCG_DAT252.docx` | Informe académico completo en formato oficial `DAT-252` (ecuaciones OMML nativas, tablas experimentales, pseudocódigos y figuras). |
| `Presentacion_QuasiNewton_BFGS_NewtonCG.html` | Presentación interactiva HTML5 (14 diapositivas con KaTeX, notas del expositor y 3 simuladores interactivos en vivo). |
| `Guia_Para_Dummies_BFGS_NewtonCG.md` | Guía intuitiva de estudio con analogías, diccionario de símbolos, ejemplo de BFGS a mano y preguntas de defensa. |
| `Exposicion_QuasiNewton_BFGS_NewtonCG.ipynb` | Notebook de Jupyter reproducible con implementaciones desde cero (GD, DFP, BFGS, L-BFGS, Newton-CG exacto y *Hessian-Free*) y validación con `scipy.optimize.minimize`. |
| `requirements.txt` | Dependencias Python fijadas para reproducir todos los experimentos y figuras. |
| `figuras/` | Figuras en alta resolución (`fig1_rosenbrock.png`, `fig2_curvatura_negativa.png`, `fig3_cuadratica_escalado.png`) y resumen JSON generados por el notebook. |

## Ejecución del Notebook (Entorno Virtual)

```powershell
python -m venv venv
.\venv\Scripts\Activate.ps1
pip install -r requirements.txt
jupyter lab Exposicion_QuasiNewton_BFGS_NewtonCG.ipynb
```

## Presentación HTML

Abre `Presentacion_QuasiNewton_BFGS_NewtonCG.html` directamente en cualquier navegador moderno:
- `←` / `→` o `Espacio`: cambiar de diapositiva
- `N`: mostrar/ocultar notas del expositor
- `O`: vista general de diapositivas
- `D`: alternar tema claro/oscuro
- `F`: pantalla completa
