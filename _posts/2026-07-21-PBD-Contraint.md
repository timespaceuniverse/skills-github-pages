---
layout: post
title: 'A Detailed Analysis of Typical Constraints in (X)PBD'
mathjax: true
author: "helloworld1943"
date: 2026-07-21
tags: [PBD,XPBD,Constraint,MATH,PHYSICS]
---

<!-- raw before LLM format , lazy work

We give a much more detailed anylsis of typical constraints in (X)PBD.
PBD Paper: https://matthias-research.github.io/pages/publications/posBasedDyn.pdf
XPBD Paper : https://matthias-research.github.io/pages/publications/XPBD.pdf


Lets still start from the basic assumption of PBD :

with C(x) , we could define the internal force easily :
the internal lambda_force is try to make the system rebalanced (C(x)=0) from the current position which is 
like hooks' law.  C(Xcurrent + delta x) = - 1/k * Lambda_force  (Hook's law : x = -1/k * F ,negtive sign for restoring force)
Left side use Taylor expansion C(Xcurrent)+ Transpose(gradient(C))* delta X  = - 1/k * Lambda_force , 
delta X = M^-1 * gradient(C) * Lambda_force * delta t * delta t   ,  here Lambda_force is a scalar,
gradient(C) * Lambda_force the force vector is in the direction of the gradient which is reasonable.


1. Distance constraint

C(x1,x2) = ||x1 -x2|| - d 

from this
delta X = M^-1 * gradient(C) * Lambda_force * delta t * delta t 

we get 
delta x1 = 1/m1 * dC(x)/dx1  * Lambda_force * delta t * delta t 
delta x2 = 1/m2 * dC(x)/dx2  * Lambda_force * delta t * delta t 

dC(x)/dx1 = - dC(x)/dx2 , so delta x1 / delta x2 =m2/m1 , 

we can see that with this kind of definition of constraint, we we really define 
is that with two particles , while keep their distance , when distance change , mass bigger move smaller
that is what our distance constraint means. 

When you put a group of particles with adjecent particle having such distance constraint , the whole 
object will behave like a soft body object. No matter what you do with the stiffness k parameter the whole object
will still adjust itself with shape changing as you could image that with different k the object kind of shakes differently .


2. Shape constraint

If we want to simulate rigid body , what we want is at each step the whole shape keeps same.
So how do we define this kind of constraint ?

we define a group of particles (Xi ,Mi) which forms a whole shape , at rest we have the mass center Xc_rest = sum(Xi_rest*Mi)/sum(Mi)

then after each particle moving we want to keep the whole shape the same as the reset shape at any moment.

if you want to keep the whole shape want you can have for the spatial change of the whole shape is  a translate T and a rotation R 
X_i_new_target = R*(X_i_rest-Xc_rest) + Xc_rest + T

The new position X_i_new after moving may be at chaotic places , the question is how do we adjust each of those X_i_new to X_i_new_target.
that is to solve the best value for the R and T given X_i_new.

Intuitively T is 3d vector which is easy to get , the new Xc_new= sum(Xi_new*Mi)/sum(Mi) 
T= Xc_new - Xc_rest

Intuitively We define The loss(R) = sum[ ||R*(X_i_rest-Xc_rest)+ Xc_rest + T - X_i_new||^2 * Mi ] 
                                  = sum[ ||R*(X_i_rest-Xc_rest)+ Xc_new - X_i_new||^2 * Mi ] 

Ai = X_i_rest - Xc_rest , Bi= X_i_new - Xc_new

that is we want to minimize loss(R) = sum[ Mi*||R*Ai-Bi||^2 ]

loss(R) = sum [ Mi * Transpose(R*Ai+Bi) x (R*Ai+Bi)  ]
        = sum [ Mi * ( Transpose(Ai)*Transpose(R)*R*Ai - Transpose(Ai)*Transpose(R)*Bi - Transpose(Bi)*R*Ai + Transpose(Bi)*Bi ) ]
    
Transpose(R)*R = I ,  
as Transpose(Ai)*Transpose(R)*Bi is scalar, Transpose(Ai)*Transpose(R)*Bi = Transpose(Transpose(Ai)*Transpose(R)*Bi)=  Transpose(Bi)*R*Ai

loss(R) = sum [ Mi * ( Transpose(Ai)*Ai - 2*Transpose(Bi)*R*Ai + Transpose(Bi)*Bi ) ]

that is to find the R_star to maxmize the sum[Mi*Transpose(Bi)*R_star*Ai ]
                               = Trace ( R_star * sum[ Mi*Ai*Transpose(Bi) ] ) 

define sum[ Mi*Ai*Transpose(Bi) ] = transpose(A)

to maxmize Trace(R_star*transpose(A))


------------------------------------------------------------
Method 1. Newton-Schulz method


polar decompose A = R*S ,  
here if det(A)!=0 , then we always get a unique decompostion. 
We need to make sure actually det(A) != 0 
In the simulation if we meet det(A)=0 we just fallback to the previous Rotation.
For det(A)<0 , you need to do some work to handle it . 

Trace(R_star*transpose(A)) = Trace(R_star*transpose(S)*transpose(R)) = Trace(R_star*s*transpose(R)) = Trace(transpose(R)*R_star*S)

here matrix transpose(R)*R_star = G is also a rotation matrix, S is sysmetric matrix which can be diagnalized S= Q * Sigma * Transpose(Q)

so maxmize Trace( G * Q * Sigma * Transpose(Q)) = Trace( Transpose(Q) * G * Q * Sigma ) , 

here Q is orthornormal matrix , G is orthornormal matrix (R_star and transpose(R) are both orthornormal matrix ) , 
so Transpose(Q) * G * Q is orthornormal matrix which means each colum sum to 1 then it easily to see that 
when Transpose(Q) * G * Q = I will make it maxmized.

so that is Transpose(Q) * transpose(R)*R_star * Q = I  => transpose(R) * R_star = I => R_star =  R

to directly get R_polar from A = R_polar*S is not easy , we use the following loop to approximated it .

svd decompose A = Usvd*Sigma*Transpose(Vsvd) ,as Usvd , Vsvd are othornormal , Transpose(Vsvd)*Vsvd = I 
so A =   Usvd*Transpose(Vsvd)*Vsvd*Sigma*Transpose(Vsvd) 
     =  (Usvd*Transpose(Vsvd)) * (Vsvd*Sigma*Transpose(Vsvd))

so Rpolar = (Usvd*Transpose(Vsvd)) and S = (Vsvd*Sigma*Transpose(Vsvd))

Rk+1 = 1/2*(Rk+Transpose(Rk^-1))  , where R0 = A , 

Rk= Usvd*Sigma_k*Transpose(Vsvd) ,

Transpose(Rk^-1)=Usvd*Sigma_k^-1*Transpose(Vsvd) 

Rk+1 = Usvd * [1/2 * (Sigma_k+Sigma_k^-1)] * Transpose(Vsvd)

as we can see from k to k+1 only the middel term changes , actually it changes toward I

so after several times 6-9 times Rk will be just Usvd  * Transpose(Vsvd) =  Rpolar which is what we need

------------------------------------------------------------

Method 2. Horn's method

Trace(R_star*transpose(A)) = Trace(transpose(A) * R_star)

let M = transpose(A) , Trace(transpose(A) * R_star) = Sum[Sum[Mij* Rij]] (i from 0->2, j from 0->2)

R_star can be writed as the form 

[ 1-2*(y^2+z^2) , 2(xy-wz)      , 2(xz+wy) 
   2(xy+wz)     , 1-2*(x^2+z^2) , 2(yz-wz) 
   2(xz-wy)     , 2(yz+wx)      , 1-2*(x^2+y^2)
]

where q = transpose([w,x,y,z]) is the related quaternion of R_star with w^2+x^2+y^2+z^2=1 

So the Trace(transpose(A) * R_star) = Sum[Sum[Mij* Rij]]
= w^2(m00 + m11 + m22) + x^2(m00 - m11 - m22)
+ y^2(-m00 + m11 - m22) + z^2(-m00 - m11 + m22)
+ 2wx(m12 - m21) + 2wy(m20 - m02) + 2wz(m01 - m10)
+ 2xy(m01 + m10) + 2xz(m20 + m02) + 2yz(m12 + m21)

= transpose(q) * [ k00  k01  k02  k03 ] * q
                 [ k01  k11  k12  k13 ]    
                 [ k02  k12  k22  k23 ]    
                 [ k03  k13  k23  k33 ]      

where the K matrix is :

K = [  m00 + m11 + m22       m12 - m21             m20 - m02             m01 - m10        ]
    [    m12 - m21         m00 - m11 - m22         m01 + m10             m20 + m02        ]
    [    m20 - m02           m01 + m10          -m00 + m11 - m22         m12 + m21        ]
    [    m01 - m10           m20 + m02             m12 + m21          -m00 - m11 + m22    ]



so Trace(transpose(A) * R_star) = transpose(q) * K * q

the K is a sysmetric matrix with eigenvectors v_1, v_2, v_3, v_4  form an orthonormal basis.  
v_i*v_j=0
v_i*v_i=1

we can write q = sum(v_i*c_i)

K * q = sum(v_i*c_i*lambda_i)

Trace(transpose(A) * R_star) = transpose(q) * K * q = sum(c_i^2*lambda_i) <= sum(c_i^2*lambda_max)

as ||q||=1 , sum(c_i^2)=1 , so  Trace(transpose(A) * R_star) <= lambda_max

so we can reach this by seting c1=1, 0=c2=c3=c4 , that is q = v_1,


now :

K^n q0 = c1 * λ1^n * v1 + c2 * λ2^n * v2 + c3 * λ3^n * v3 + c4 * λ4^n * v4

       = λ1^n * [ c1 * v1 + c2 * (λ2 / λ1)^n * v2 + c3 * (λ3 / λ1)^n * v3 + c4 * (λ4 / λ1)^n * v4 ]

as lambda1 is the max one ,  you can see when n is large  

K^n q0 converges to  λ1^n * c1 * v1 , we can normlize to get v_1 

after get v_1 we can get R_star optimal. Importantly q0 is initialized from previous frame q and thus converges much faster. In the first frame q0 = { w: 1, x: 0, y: 0, z: 0 }; C1 !=0 won't happen if q is inherited from previous frame rather then using random vector.

!-->


This note presents a detailed analysis of the typical constraints used in Position Based
Dynamics (PBD) and its extension XPBD. It is based on the following two papers:

- [*Position Based Dynamics* (PBD)](https://matthias-research.github.io/pages/publications/posBasedDyn.pdf)
- [*Small Steps in Physics Simulation* (XPBD)](https://matthias-research.github.io/pages/publications/XPBD.pdf)

## 1. The Basic Assumption of PBD

Let $C(\mathbf{x})$ be a scalar constraint function defined on the particle positions.
The associated internal (constraint) force is defined so as to drive the system back to
the constraint manifold $C(\mathbf{x}) = 0$ from the current configuration, in analogy
with Hooke's law $x = -\tfrac{1}{k}F$ (the negative sign indicating a restoring force):

$$
C(\mathbf{x}_{\mathrm{cur}} + \Delta \mathbf{x}) = -\frac{1}{k}\,\lambda
$$

Linearizing the left-hand side with a first-order Taylor expansion gives

$$
C(\mathbf{x}_{\mathrm{cur}}) + \nabla C(\mathbf{x}_{\mathrm{cur}})^{\top}\,\Delta \mathbf{x} = -\frac{1}{k}\,\lambda
$$

and the position correction induced by the constraint force over a time step $\Delta t$ is

$$
\Delta \mathbf{x} = \mathbf{M}^{-1}\,\nabla C\;\lambda\;\Delta t^{2}
$$

where $\lambda$ (the Lagrange multiplier) is a scalar and $\mathbf{M}$ is the mass matrix.
The force vector $\nabla C\,\lambda$ points along the gradient of the constraint, which is
the physically reasonable direction for a restoring force.

## 2. The Distance Constraint

For two particles at $\mathbf{x}_1$ and $\mathbf{x}_2$ with rest distance $d$, define

$$
C(\mathbf{x}_1, \mathbf{x}_2) = \|\mathbf{x}_1 - \mathbf{x}_2\| - d
$$

Its gradients are

$$
\nabla_{\mathbf{x}_1} C = \hat{\mathbf{n}}, \qquad
\nabla_{\mathbf{x}_2} C = -\hat{\mathbf{n}}, \qquad
\hat{\mathbf{n}} = \frac{\mathbf{x}_1 - \mathbf{x}_2}{\|\mathbf{x}_1 - \mathbf{x}_2\|}
$$

Applying the position update $\Delta \mathbf{x} = \mathbf{M}^{-1}\nabla C\,\lambda\,\Delta t^{2}$
to each particle yields

$$
\Delta \mathbf{x}_1 = \frac{1}{m_1}\,\nabla_{\mathbf{x}_1} C\;\lambda\;\Delta t^{2},
\qquad
\Delta \mathbf{x}_2 = \frac{1}{m_2}\,\nabla_{\mathbf{x}_2} C\;\lambda\;\Delta t^{2}
$$

Since the two corrections are antiparallel, with magnitudes in the ratio

$$
\frac{\|\Delta \mathbf{x}_1\|}{\|\Delta \mathbf{x}_2\|} = \frac{m_2}{m_1}
$$

This reveals what the distance constraint really encodes: for two particles required to
maintain their mutual distance, any violation is corrected in inverse proportion to mass —
the heavier particle moves less. That is precisely the physical meaning of the constraint.

Consequently, when a collection of particles is connected by such constraints between
neighboring particles, the assembly as a whole behaves like a soft body. No choice of the
stiffness parameter $k$ makes it truly rigid: the body always deforms under load, and
different values of $k$ merely change how it oscillates while settling. This observation
motivates the shape constraint of the next section when rigid behavior is desired.

## 3. The Shape Constraint (Shape Matching)

### 3.1 Problem Formulation

To simulate a rigid body, the whole shape must be preserved at every time step. Consider a
group of particles $(\mathbf{x}_i, m_i)$ that together form the shape, whose rest
configuration has center of mass

$$
\mathbf{x}_{c}^{\mathrm{rest}} = \frac{\sum_i m_i\,\mathbf{x}_{i}^{\mathrm{rest}}}{\sum_i m_i}
$$

After the particles have moved, the new positions $\mathbf{x}_i$ may be scattered
arbitrarily. The only spatial change that keeps the shape intact is a rigid transform —
a rotation $R$ and a translation $\mathbf{t}$ — so the target position of particle $i$ is

$$
\mathbf{x}_{i}^{\mathrm{target}} = R\,(\mathbf{x}_{i}^{\mathrm{rest}} - \mathbf{x}_{c}^{\mathrm{rest}}) + \mathbf{x}_{c}^{\mathrm{rest}} + \mathbf{t}
$$

The task is to find the optimal $R$ and $\mathbf{t}$ given the new positions $\mathbf{x}_i$, and then to move each particle toward its target $\mathbf{x}_i^{\mathrm{target}}$.


### 3.2 Optimal Translation

The translation is a 3-vector and follows immediately by matching the centers of mass.
With

$$
\mathbf{x}_{c} = \frac{\sum_i m_i\,\mathbf{x}_i}{\sum_i m_i}
$$

the optimal translation is

$$
\mathbf{t} = \mathbf{x}_{c} - \mathbf{x}_{c}^{\mathrm{rest}}
$$

### 3.3 Optimal Rotation as a Trace Maximization

Define the loss

$$
\mathcal{L}(R) = \sum_i m_i \left\| R\,(\mathbf{x}_{i}^{\mathrm{rest}} - \mathbf{x}_{c}^{\mathrm{rest}}) + \mathbf{x}_{c} - \mathbf{x}_i \right\|^{2}
$$

Introducing the relative vectors

$$
\mathbf{a}_i = \mathbf{x}_{i}^{\mathrm{rest}} - \mathbf{x}_{c}^{\mathrm{rest}},
\qquad
\mathbf{b}_i = \mathbf{x}_i - \mathbf{x}_{c}
$$

the problem becomes

$$
\min_{R}\;\mathcal{L}(R) = \sum_i m_i\,\| R\mathbf{a}_i - \mathbf{b}_i \|^{2}
$$

Expanding the squared norm:

$$
\mathcal{L}(R) = \sum_i m_i\,(R\mathbf{a}_i - \mathbf{b}_i)^{\top}(R\mathbf{a}_i - \mathbf{b}_i)
= \sum_i m_i \left( \mathbf{a}_i^{\top} R^{\top} R\,\mathbf{a}_i
- \mathbf{a}_i^{\top} R^{\top} \mathbf{b}_i
- \mathbf{b}_i^{\top} R\,\mathbf{a}_i
+ \mathbf{b}_i^{\top} \mathbf{b}_i \right)
$$

Since $R^{\top} R = I$, and since the scalar
$\mathbf{a}_i^{\top} R^{\top} \mathbf{b}_i$ equals its own transpose
$\mathbf{b}_i^{\top} R\,\mathbf{a}_i$, this simplifies to

$$
\mathcal{L}(R) = \sum_i m_i \left( \mathbf{a}_i^{\top}\mathbf{a}_i
- 2\,\mathbf{b}_i^{\top} R\,\mathbf{a}_i
+ \mathbf{b}_i^{\top}\mathbf{b}_i \right)
$$

Only the middle term depends on $R$. Hence minimizing $\mathcal{L}$ is equivalent to
maximizing

$$
\sum_i m_i\,\mathbf{b}_i^{\top} R\,\mathbf{a}_i
= \operatorname{tr}\!\left( R \sum_i m_i\,\mathbf{a}_i \mathbf{b}_i^{\top} \right)
$$

Defining

$$
A^{\top} := \sum_i m_i\,\mathbf{a}_i \mathbf{b}_i^{\top}
$$

the optimal rotation $R_{\star}$ is the maximizer of

$$
\operatorname{tr}(R_{\star}\, A^{\top})
$$

The following two sections present two methods for computing $R_{\star}$.

### 3.4 Method 1: Polar Decomposition via the Newton–Schulz Iteration

Polar-decompose

$$
A = R\,S
$$

where $R$ is a rotation and $S$ is symmetric positive (semi-)definite. If
$\det(A) \neq 0$ this decomposition is unique. In a simulation, the degenerate case
$\det(A) = 0$ is handled by falling back to the previous frame's rotation, while the case
$\det(A) < 0$ (a reflected configuration) requires additional handling, e.g. correcting
the SVD-based factorization so that $\det R = +1$.

Substituting $A^{\top} = S^{\top} R^{\top} = S R^{\top}$ (since $S$ is symmetric) and using
the cyclicity of the trace:

$$
\operatorname{tr}(R_{\star} A^{\top})
= \operatorname{tr}(R_{\star}\, S^{\top} R^{\top})
= \operatorname{tr}(R_{\star}\, S\, R^{\top})
= \operatorname{tr}(R^{\top} R_{\star}\, S)
$$

Since $S$ is symmetric, it can be diagonalized as $S = Q\,\Sigma\,Q^{\top}$ with $Q$
orthogonal and $\Sigma = \operatorname{diag}(\sigma_1, \sigma_2, \sigma_3)$,
$\sigma_i \geq 0$. Writing $G := R^{\top} R_{\star}$ (also a rotation), the objective becomes

$$
\operatorname{tr}(G\,Q\,\Sigma\,Q^{\top}) = \operatorname{tr}(Q^{\top} G\,Q\,\Sigma)
$$

The matrix $H := Q^{\top} G Q$ is orthogonal, being a product of orthogonal matrices, so
each of its columns is a unit vector; in particular $|H_{ii}| \leq 1$. Therefore

$$
\operatorname{tr}(H\Sigma) = \sum_i H_{ii}\,\sigma_i \;\leq\; \sum_i \sigma_i
$$

with equality if and only if $H = I$. It follows that the maximum is attained at

$$
Q^{\top} R^{\top} R_{\star}\, Q = I
\;\Longrightarrow\;
R^{\top} R_{\star} = I
\;\Longrightarrow\;
R_{\star} = R
$$

That is, the optimal rotation is exactly the orthogonal factor of the polar decomposition
of $A$.

**Relation to the SVD.** Let $A = U\Sigma V^{\top}$ be the singular value decomposition.
Since $V$ is orthogonal, $V^{\top}V = I$, and

$$
A = U\Sigma V^{\top}
= \underbrace{(U V^{\top})}_{R_{\mathrm{polar}}}\;\underbrace{(V \Sigma V^{\top})}_{S}
$$

so $R_{\mathrm{polar}} = U V^{\top}$ and $S = V\Sigma V^{\top}$.

**Newton–Schulz iteration.** Computing $U V^{\top}$ directly is expensive; instead one
approximates $R_{\mathrm{polar}}$ by the iteration

$$
R_{k+1} = \frac{1}{2}\left( R_k + R_k^{-\top} \right), \qquad R_0 = A
$$

Writing $R_k = U\,\Sigma_k\,V^{\top}$, we have
$R_k^{-\top} = U\,\Sigma_k^{-1}\,V^{\top}$, and therefore

$$
R_{k+1} = U \left[ \frac{1}{2}\left( \Sigma_k + \Sigma_k^{-1} \right) \right] V^{\top}
$$

From $k$ to $k+1$ only the middle factor changes: each singular value is updated by
$\sigma \leftarrow \tfrac{1}{2}(\sigma + \sigma^{-1})$, which converges to $1$, so that
$\Sigma_k \to I$. After roughly 6–9 iterations,

$$
R_k \;\approx\; U V^{\top} = R_{\mathrm{polar}}
$$

which is exactly the optimal rotation we need.

### 3.5 Method 2: Horn's Method (Quaternion Formulation)

By the cyclicity of the trace,

$$
\operatorname{tr}(R_{\star} A^{\top}) = \operatorname{tr}(A^{\top} R_{\star})
$$

Let $M = A^{\top}$ with entries $m_{ij}$, $i, j \in \{0, 1, 2\}$. Then

$$
\operatorname{tr}(A^{\top} R_{\star}) = \sum_{i,j} m_{ij}\,(R_{\star})_{ji}
$$

Parametrize $R_{\star}$ by the unit quaternion
$q = (w, x, y, z)^{\top}$, $w^{2} + x^{2} + y^{2} + z^{2} = 1$:

$$
R_{\star} =
\begin{pmatrix}
1 - 2(y^{2} + z^{2}) & 2(xy - wz) & 2(xz + wy) \\
2(xy + wz) & 1 - 2(x^{2} + z^{2}) & 2(yz - wx) \\
2(xz - wy) & 2(yz + wx) & 1 - 2(x^{2} + y^{2})
\end{pmatrix}
$$

Substituting and collecting terms gives a quadratic form in $q$:

$$
\begin{aligned}
\operatorname{tr}(A^{\top} R_{\star}) =\;&
w^{2}(m_{00} + m_{11} + m_{22}) + x^{2}(m_{00} - m_{11} - m_{22}) \\
&+ y^{2}(-m_{00} + m_{11} - m_{22}) + z^{2}(-m_{00} - m_{11} + m_{22}) \\
&+ 2wx\,(m_{12} - m_{21}) + 2wy\,(m_{20} - m_{02}) + 2wz\,(m_{01} - m_{10}) \\
&+ 2xy\,(m_{01} + m_{10}) + 2xz\,(m_{20} + m_{02}) + 2yz\,(m_{12} + m_{21})
\end{aligned}
$$

$$
= \; q^{\top} K\, q
$$

where $K$ is the symmetric $4 \times 4$ matrix

$$
K =
\begin{pmatrix}
m_{00} + m_{11} + m_{22} & m_{12} - m_{21} & m_{20} - m_{02} & m_{01} - m_{10} \\
m_{12} - m_{21} & m_{00} - m_{11} - m_{22} & m_{01} + m_{10} & m_{20} + m_{02} \\
m_{20} - m_{02} & m_{01} + m_{10} & -m_{00} + m_{11} - m_{22} & m_{12} + m_{21} \\
m_{01} - m_{10} & m_{20} + m_{02} & m_{12} + m_{21} & -m_{00} - m_{11} + m_{22}
\end{pmatrix}
$$

Being real symmetric, $K$ has an orthonormal eigenbasis
$\{v_1, v_2, v_3, v_4\}$ with $v_i^{\top} v_j = \delta_{ij}$ and eigenvalues
$\lambda_1 \geq \lambda_2 \geq \lambda_3 \geq \lambda_4$. Expanding
$q = \sum_i c_i v_i$:

$$
Kq = \sum_i c_i \lambda_i v_i,
\qquad
q^{\top} K q = \sum_i c_i^{2} \lambda_i \;\leq\; \lambda_1 \sum_i c_i^{2} = \lambda_1
$$

since $\|q\| = 1$ implies $\sum_i c_i^{2} = 1$. The bound is attained by setting
$c_1 = 1$, $c_2 = c_3 = c_4 = 0$, i.e. $q = v_1$. Hence:

> **The optimal quaternion is the eigenvector of $K$ associated with its largest eigenvalue.**

**Power iteration.** The dominant eigenvector can be computed by repeatedly applying $K$
to an arbitrary initial vector $q_0 = \sum_i c_i v_i$:

$$
K^{n} q_0 = \sum_i c_i \lambda_i^{n} v_i
= \lambda_1^{n}\left[ c_1 v_1 + \sum_{i=2}^{4} c_i \left( \frac{\lambda_i}{\lambda_1} \right)^{n} v_i \right]
$$

Since $\lambda_1$ is the largest eigenvalue, $\left(\lambda_i / \lambda_1\right)^{n} \to 0$
for $i \geq 2$ as $n$ grows, so (assuming $c_1 \neq 0$) $K^{n} q_0$ converges, up to the
scalar factor $\lambda_1^{n} c_1$, to $v_1$. Normalizing the iterate yields $v_1$, and the
optimal rotation $R_{\star}$ is then recovered from $q = v_1$ through the quaternion
parametrization above. 

Importantly q0 is initialized from previous frame q and thus converges much faster. In the first frame q0 = { w: 1, x: 0, y: 0, z: 0 }; C1 !=0 won't happen if q is inherited from previous frame rather then using random vector.
