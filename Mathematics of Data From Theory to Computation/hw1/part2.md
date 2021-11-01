The solution before Part 2.4 is hand-written in anther PDF file.
### 2.1
I implemet all the credited method in algorithm.py
### 2.2/2.3
It is contained in the written part.
### 2.4
Go to folder ./question2. 
#### a
The observed convergence rates is not entirely consistent with the theoretical ones. The theoretical ones provide an upper bound. I adopt a tighter bound by noticing $\sigma(t)(1-\sigma(t))\leq \frac{1}{4}$, since we could have
$$
\|\nabla f_\mu(x)\|_2 \leq \frac{1}{4}\|A^\top A\|_2+\mu
$$
The results are shown in Line 33 in question_2_4.py and *_tight.pdf files. We observe that it is a tighter estimate than before.

### b
The observed convergence rates of GD and GDstr is linear. Since our objective is strongly convex, setting $\alpha=\frac{1}{L}$ and $\alpha=\frac{2}{L+\mu}$ will result in linear convergence. The regression result confirms such idea.

### d
The empirical speed of convergence is faster than the theoretical one. Since the theoretical result is an upper bound. The bash output for ./question_2_4.py is

GD results: a_GD=-0.000180, b_GD=1.196370
GD theoretical rate: a_GD=-0.000053, b_GD=1.282497
GD rate diff =-0.00012623199407315507
GD_str results: a_GD=-0.000289, b_GD=1.119041
GD_str theoretical rate: a_GD=-0.000107, b_GD=1.282084436227831
GD_str rate diff =-0.000182

and the resulting figures are stored in ./figs/fig_ex2_4_convergence_rate,pdf and ./figs/fig_ex2_4_convergence_rate_full.pdf. If we use a tighter upper bound as discussed in 2.4.a, the resulting figures are stored in ./figs/fig_ex2_4_convergence_rate_tight,pdf and ./figs/fig_ex2_4_convergence_rate_full_tighter.pdf

All the method is tested in ./log_reg.py, the classification error for ./log_reg.py is reported below:

GD : 0.1386861313868613

GDstr : 0.0948905109489051

AGD : 0.058394160583941604

AGDstr : 0.058394160583941604

AGDR : 0.072992700729927

AdaGrad : 0.058394160583941604

SGD : 0.072992700729927

SAG : 0.051094890510948905

SVR : 0.058394160583941604

SubG : 0.12408759124087591

L1-prox

ISTA : 0.1386861313868613

FISTA : 0.145985401459854

FISTAR : 0.145985401459854

PROXSG : 0.145985401459854

L2-prox

ISTA : 0.11678832116788321

FISTA : 0.06569343065693431

FISTAR : 0.058394160583941604

PROXSG : 0.12408759124087591

we can see that FISTA and FISTAR, AdaGrad, Stochastic method (SGD, SAG, SVR) would be better choice.


### 3.3.1
#### (a) 
$$
\nabla (f_{\ell_1} + g_{\ell_1})(\alpha) = \mathbf{W}\mathbf{P}_\Omega^\top(\mathbf{P}_\Omega \mathbf{W}^\top \alpha-\mathbf{b})+\lambda_{\ell_1} \mathrm{sign}(\alpha)
$$
where 
$$
\mathrm{sign}(x) = \begin{cases}
    1 & \text{if } x >0 \\
    -1 & \text{if } x<0 \\
    [-1,1]&\text{if }x=0
\end{cases}
$$
and applied element-wisely. Then
$$
\nabla(f_{\bf TV}+g_{\bf TV})(x)=\mathbf{P}_\Omega^\top(\mathbf{P}_\Omega x - b) + \lambda_{\bf TV}\nabla g_{\bf TV} (x)
$$

- For the isotropic case.
Define
$$
D_1 = \left(\begin{matrix}
    -1 & 1\\
     & -1 & 1\\
     &&\ddots &\ddots \\
     &&&-1 & 1\\
     0&&\cdots&&0
    \end{matrix}\right) \in \mathbb{R}^{m\times m}
$$
(1) We start with the differencing matrix by column. It is just the differencing matrix of a TV term of dimension $m$, namely, $D_1$, i.e., $\nabla x^1_{:,j} =  D_1 \cdot x_{:,j}$.

(2) We can then proceed by column, i.e., $\nabla x_{i,:} =  x_{i,:}\cdot D_1^\top $.

(3) If we stacking the matrix into a vector, define $\text{vec}(X)_{i+(j-1)m}=X_{i,j}$, we could have $\text{vec}(\nabla x^1) = (I\otimes D_1) \cdot\text{vec}(x)$, where $I$ is a $m\times m$ identity matrix.

(4) Since two consecutive entries of the same row are separated
by $m$ positions after vectorization, we have $\text{vec}(\nabla x^2) = (D_1\otimes I) \cdot\text{vec}(x)$.

(5) In total, we can write
 $$\|x\|_{\bf TV,\ell_1}=\|(I\otimes D_1) \cdot\text{vec}(x)\|_1+\| (D_1\otimes I) \cdot\text{vec}(x)\|_1$$
hence the gradient is
$$\nabla\|x\|_{\bf TV,\ell_1}=(I\otimes D_1^\top){\rm sign}((I\otimes D_1) \cdot\text{vec}(x))+(D_1^\top\otimes I){\rm sign}( (D_1\otimes I) \cdot\text{vec}(x))$$
by the chain rule.

- For the anisotropic case, the matrix form does not simplify the calculation. Denote $x_g = \nabla \|x\|_{\bf TV,\ell_2}$. Notice $x_{i,j}$ only appears in the term $\nabla x^1_{i,j}$, $\nabla x^2_{i,j}$, $\nabla x^1_{i-1,j}$, $\nabla x^2_{i,j-1}$. We have
$$
(\nabla \|x\|_{\bf TV,\ell_2})_{i,j}=\frac{\partial (\|\nabla x^1_{i,j}\|_2+\|\nabla x^2_{i,j}\|_2+\|\nabla x^1_{i-1,j}\|_2+\|\nabla x^2_{i,j-1}\|_2)}{\partial x_{i,j}}= -\frac{\nabla x^1_{i,j}+\nabla x^2_{i,j}}{\|\nabla x_{i,j}\|_2} + \frac{\nabla x^1_{i-1,j}}{\|\nabla x_{i-1,j}\|_2}+\frac{\nabla x^2_{i,j-1}}{\|\nabla x_{i,j-1}\|_2}
$$
which is the gradient.
#### (b)
By the definition of $\|A\|_2=\sup_{x\neq 0}\frac{\|Ax\|_2}{\|x\|_2}$. And by SVD decomposition, $A=U\Sigma V^\top$ we have $A^\top A = V \Sigma^2 V^\top $. Hence, the operator norm of $A^\top A$ (equals the largest singular value) is the square of the operator norm of $A$.



Since 
$$
\nabla f_{\ell_1}(\alpha)=\mathbf{W}\mathbf{P}_\Omega^\top(\mathbf{P}_\Omega \mathbf{W}^\top \alpha-\mathbf{b})
$$
Hence, the Lipschitz constant would be
$$
\|\mathbf{W}\mathbf{P}_\Omega^\top\mathbf{P}_\Omega\mathbf{W}^\top\|_2=\|\mathbf{P}_\Omega\mathbf{W}^\top\|_2^2 = \|\mathbf{P}_\Omega\|_2^4=1
$$
since $\mathbf{W}$ is unitary matrix, and $\mathbf{P}_\Omega \mathbf{P}_\Omega^\top$ is a diagnoal matrix which only has 1 as its nonzero element on its diagonal. Similarly, the Lipschitz constant for $\nabla f_{\bf TV}$ is
$$
\|\mathbf{P}_\Omega^\top\mathbf{P}_\Omega\|_2=\|\mathbf{P}_\Omega\|_2^2=1
$$



3.3.2 
Go to folder ./question3. 
The  reconstructed images are plotted in ./results/log_lambda_{lambda}.png, when setting lambda = $10^{\{-3,-2.5,-2,\cdots,0.5,1\}}$, and the PSNR-lambda relation is plotted in ./results/PSNR-lambda.png. We can see that setting lambda around 0.01-0.1 would be the best choice. Setting $\lambda$ too large will make $\ell_1$ penalized method output zero.

In general, from the recovered image, we find that $\ell_1$ norm tends to make the reconstructed image dimmer than the original ones, and it cannot make clear distinctions between different patch of color, but it can recover some details of the original image. While the reconstruction image by TV norm will lose detailed information, but has almost the same brightness, and can make clear distinction between dark and light colors.
