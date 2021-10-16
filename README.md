# nongrav_yarkovsky

Fortran file including nongravitational Yarkovsky effect for SWIFT N-body simulator. See [Hal Levison's](https://www.boulder.swri.edu/~hal/swift.html) 
website for download and documentation on SWIFT. Yarkovsky calculation follows [Vokrouhlicky 2000](https://doi.org/10.1006/icar.2000.6469) for diurnal
acceleration vector.

\begin{equation}
    \mathbf{a}_d = \frac{4\alpha}{9}\frac{\Phi(r)}{1+\lambda}G[\sin\delta + \cos\delta\mathbf{s}\times]\frac{\mathbf{r}\times\mathbf{s}}{r}
\end{equation}
