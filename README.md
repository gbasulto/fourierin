# fourierin
Guillermo Basulto-Elias  

This is a package in `R` to numerically calculate Fourier-type integrals of multivariate functions with compact support. Specifically, integrals of the type

<!--
$$
I \left[f(t), \boldsymbol{a}, \boldsymbol{b};r, s \right]
  = \left[ \frac{|s|}{(2\pi)^{1 - r}}\right]^{n/2}
  \int_{a_1}^{b_1}\int_{a_2}^{b_2}\cdots\int_{a_n}^{b_n}
  f(\boldsymbol{t})e^{\imath s 
  \langle \boldsymbol{w}, \boldsymbol{t}\rangle} \text{d}\boldsymbol{t},
$$
-->
<img src="eq01_main_eq.png" width="600px" />

where,

<!--
$\boldsymbol{a} = (a_1, \ldots, a_n)$, 
$\boldsymbol{b} = (b_1, \ldots, b_n)$, 
$\boldsymbol{t} = (t^{(1)}, \ldots, t^{(n)})$, 
$\boldsymbol{w} = (w^{(1)}, \ldots, w^{(n)})$, 
$a_l \leq t^{(l)} \leq b_l$,
$c_l \leq w^{(l)} \leq c_l$.
-->
<img src="eq02_details.png" width="200px" />

Common values for _r_ are -1, 0 and -1, while common values for _s_ are _-2\pi_, -1, 1 and $2\pi$. For example, if _f_ is a density function, _s = 1_ and _r = 1_ could be used to obtain the characteristic function of _f_. Conversely, if _f_ is the characteristic function of a probability density function, then _r = -1_ and _s = -1_ could be used to recover the density.

