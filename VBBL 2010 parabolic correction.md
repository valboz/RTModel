## Orignial integral
$$
\begin{aligned}

dA_I &= \frac{1}{2}\int_{\theta_i}^{\theta_i+\Delta\theta} \vec{x} \wedge \vec{x'} d\theta \\
\vec{x} &= \vec{x}_{\theta_0} + \vec{x'}_{\theta_0}d\theta+\frac{1}{2}\vec{x''_{\theta_0}}d\theta^2+\frac{1}{6}\vec{x'''_{\theta_0}}d\theta^3\\
\vec{x'} &= \vec{x'}_{\theta_0} + \vec{x''}_{\theta_0}d\theta+\frac{1}{2}\vec{x'''_{\theta_0}}d\theta^2+\frac{1}{6}\vec{x''''_{\theta_0}}d\theta^3\\
\end{aligned}
$$
so we have the Taylor expansion of this integration: we only consider the first four order (fifth order after integration)
$$
\begin{aligned}
x\wedge x' &= x_{\theta_0} \wedge x'_{\theta_0} +x_{\theta_0} \wedge x''_{\theta_0} d \theta + \frac{1}{2}x_{\theta_0} \wedge x'''_{\theta_0} d\theta^2 + \frac{1}{6} x_{\theta_0} \wedge x''''_{\theta_0}d\theta^3 \\ 

&+ x'_{\theta_0} \wedge x''_{\theta_0} d \theta^2 + \frac{1}{2}x'_{\theta_0} \wedge x'''_{\theta_0} d\theta^3 + \frac{1}{6} x'_{\theta_0} \wedge x''''_{\theta_0}d\theta^4 \\

&+ \frac{1}{2}x''_{\theta_0} \wedge x'_{\theta_0} d\theta^2 + \frac{1}{4}x''_{\theta_0} \wedge x'''_{\theta_0} d\theta^4 \\ 

&+ \frac{1}{6}x'''_{\theta_0} \wedge x'_{\theta_0} d\theta^3 + \frac{1}{6}x'''_{\theta_0} \wedge x''_{\theta_0} d\theta^4 
\end{aligned}
$$
simplify this equation, we get:
$$
\begin{aligned}
x\wedge x' &= x_{\theta_0} \wedge x'_{\theta_0} +x_{\theta_0} \wedge x''_{\theta_0} d \theta + \frac{1}{2}x_{\theta_0} \wedge x'''_{\theta_0} d\theta^2 + \frac{1}{6} x_{\theta_0} \wedge x''''_{\theta_0}d\theta^3 \\ 

&+ \frac{1}{2}x'_{\theta_0} \wedge x''_{\theta_0} d \theta^2 + \frac{1}{3}x'_{\theta_0} \wedge x'''_{\theta_0} d\theta^3 + \frac{1}{6} x'_{\theta_0} \wedge x''''_{\theta_0}d\theta^4 \\

& + \frac{1}{12}x''_{\theta_0} \wedge x'''_{\theta_0} d\theta^4 \\ 
\end{aligned}
$$
so $\int_{\theta_i}^{\theta_i+\delta \theta} x\wedge x' d\theta$ equals:
$$
\begin{aligned}
\int x\wedge x' d\theta &= x_{\theta_0} \wedge x'_{\theta_0} d\theta +\frac{1}{2}x_{\theta_0} \wedge x''_{\theta_0} d \theta^2 + \frac{1}{6}x_{\theta_0} \wedge x'''_{\theta_0} d\theta^3 + \frac{1}{24} x_{\theta_0} \wedge x''''_{\theta_0}d\theta^4 \\ 

&+ \frac{1}{6}x'_{\theta_0} \wedge x''_{\theta_0} d \theta^3 + \frac{1}{12}x'_{\theta_0} \wedge x'''_{\theta_0} d\theta^4 + \frac{1}{30} x'_{\theta_0} \wedge x''''_{\theta_0}d\theta^5 \\

& + \frac{1}{60}x''_{\theta_0} \wedge x'''_{\theta_0} d\theta^5 \\ 
\end{aligned}
$$

## Trapezoidal formula 
trapezoidal formula is given by:
$$
\begin{aligned}
dA_t &= \frac{1}{2}x_\theta \wedge x_{\theta + d\theta} \\
x_\theta \wedge x_{\theta + d\theta} &= x_\theta \wedge x'_{\theta_0}d\theta+ \frac{1}{2} x_\theta \wedge x''_{\theta_0}d\theta^2+ \frac{1}{6} x_\theta \wedge x'''_{\theta_0}d\theta^3 + \frac{1}{24} x_\theta \wedge x'''' d\theta^4 + \frac{1}{120} x_\theta \wedge x''''' d\theta^5\\
\end{aligned}
$$
## Error of trapezoidal formula
$$
\begin{aligned}
T = dA_I - dA_t &= \frac{1}{2}(\frac{1}{6}x'_{\theta_0} \wedge x''_{\theta_0} d \theta^3 + \frac{1}{12}x'_{\theta_0} \wedge x'''_{\theta_0} d\theta^4 + \frac{1}{30} x'_{\theta_0} \wedge x''''_{\theta_0}d\theta^5 + \frac{1}{60}x''_{\theta_0} \wedge x'''_{\theta_0} d\theta^5 \\
&-\frac{1}{120} x_\theta \wedge x''''' d\theta^5) \\
& =  \frac{1}{12}x'_{\theta_0} \wedge x''_{\theta_0} d \theta^3 + \frac{1}{24}x'_{\theta_0} \wedge x'''_{\theta_0} d\theta^4 + \frac{1}{60} x'_{\theta_0} \wedge x''''_{\theta_0}d\theta^5 + \frac{1}{120}x''_{\theta_0} \wedge x'''_{\theta_0} d\theta^5 \\
&-\frac{1}{240} x_\theta \wedge x''''' d\theta^5
\end{aligned}
$$
## Parabolic correction
Parabolic correction is shown as:
$$
\begin{aligned}
dA_p &= \frac{1}{24}[(x'\wedge x'')|_\theta + (x'\wedge x'')|_{\theta+d\theta}] d\theta^3 \\
&= \frac{1}{24}[2*(x'\wedge x'')|_\theta + (x'\wedge x''') |_{\theta} d\theta+ \frac{1}{2} (x''\wedge x'''+x'\wedge x'''')|_\theta d\theta^2 ] d\theta^3 \\
&= \frac{1}{12}(x'\wedge x'')|_\theta d\theta^3+ \frac{1}{24}(x'\wedge x''') |_{\theta} d\theta^4+ \frac{1}{48} (x''\wedge x'''+x'\wedge x'''')|_\theta d\theta^5
\end{aligned}
$$
so the error after the parabolic correction is
$$
\begin{aligned}
T_p = dA_I - dA_t-dA_p &= \frac{1}{60} x'_{\theta_0} \wedge x''''_{\theta_0}d\theta^5 + \frac{1}{120}x''_{\theta_0} \wedge x'''_{\theta_0} d\theta^5 -\frac{1}{240} x_\theta \wedge x''''' d\theta^5 \\
& -\frac{1}{48} (x''\wedge x'''+x'\wedge x'''')|_\theta d\theta^5
\end{aligned}
$$