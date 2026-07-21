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


svd decompose A = Usvd*Sigma*Transpose(Vsvd) ,as Usvd , Vsvd are othornormal , Transpose(Vsvd)*Vsvd = I 
so A =   Usvd*Transpose(Vsvd)*Vsvd*Sigma*Transpose(Vsvd) 
     =  (Usvd*Transpose(Vsvd)) * (Vsvd*Sigma*Transpose(Vsvd))

so Rpolar = (Usvd*Transpose(Vsvd)) and S = (Vsvd*Sigma*Transpose(Vsvd))

Rather than using SVD decomposition which is hard we use Newton-Schulz alogrithm to calculate the Rpolar which is fast and easy to implement.

Rk+1 = 1/2*(Rk+Transpose(Rk^-1))  , where R0 = A , 

Rk= Usvd*Sigma_k*Transpose(Vsvd) ,

Transpose(Rk^-1)=Usvd*Sigma_k^-1*Transpose(Vsvd) 

Rk+1 = Usvd * [1/2 * (Sigma_k+Sigma_k^-1)] * Transpose(Vsvd)

as we can see from k to k+1 only the middel term changes , actually it changes toward I

so after several times 6-9 times Rk will be just Usvd  * Transpose(Vsvd) which is what we need



!-->


This note presents a detailed analysis of the typical constraints used in (X)PBD.

**References**

- PBD paper: [Position Based Dynamics](https://matthias-research.github.io/pages/publications/posBasedDyn.pdf)
- XPBD paper: [XPBD: Position-Based Simulation of Compliant Constrained Dynamics](https://matthias-research.github.io/pages/publications/XPBD.pdf)

## 1. Basic Assumption of PBD

Given a constraint function $C(\mathbf{x})$, the internal force is defined such that it drives the system back to equilibrium, i.e. $C(\mathbf{x}) = 0$, from the current configuration. This is analogous to Hooke's law, $x = -\frac{1}{k}F$, where the negative sign indicates a restoring force:

$$
C(\mathbf{x}_{\text{current}} + \Delta\mathbf{x}) = -\frac{1}{k}\,\lambda_{\text{force}}
$$

Applying a first-order Taylor expansion to the left-hand side yields

$$
C(\mathbf{x}_{\text{current}}) + \nabla C(\mathbf{x})^{\top}\,\Delta\mathbf{x} = -\frac{1}{k}\,\lambda_{\text{force}}
$$

which leads to the position update

$$
\Delta\mathbf{x} = M^{-1}\,\nabla C\;\lambda_{\text{force}}\;\Delta t^{2}
$$

Here $\lambda_{\text{force}}$ is a scalar, and the force vector $\nabla C\,\lambda_{\text{force}}$ acts along the gradient of the constraint, which is physically reasonable.

## 2. Distance Constraint

The distance constraint between two particles is defined as

$$
C(\mathbf{x}_1, \mathbf{x}_2) = \|\mathbf{x}_1 - \mathbf{x}_2\| - d
$$

Substituting this into the position update rule gives

$$
\Delta\mathbf{x}_1 = \frac{1}{m_1}\frac{\partial C(\mathbf{x})}{\partial \mathbf{x}_1}\,\lambda_{\text{force}}\,\Delta t^{2},
\qquad
\Delta\mathbf{x}_2 = \frac{1}{m_2}\frac{\partial C(\mathbf{x})}{\partial \mathbf{x}_2}\,\lambda_{\text{force}}\,\Delta t^{2}
$$

Since $\frac{\partial C(\mathbf{x})}{\partial \mathbf{x}_1} = -\frac{\partial C(\mathbf{x})}{\partial \mathbf{x}_2}$, the two displacements point in opposite directions and their magnitudes satisfy

$$
\frac{\|\Delta\mathbf{x}_1\|}{\|\Delta\mathbf{x}_2\|} = \frac{m_2}{m_1}
$$

**Interpretation.** This constraint enforces the distance between two particles while distributing the positional correction in inverse proportion to their masses: the heavier particle moves less. This is precisely what the distance constraint means.

When a group of particles is connected by such distance constraints between adjacent particles, the whole object behaves like a soft body. No matter how the stiffness parameter $k$ is chosen, the object will still deform under external influence; the value of $k$ only changes the way the object oscillates while it adjusts itself.

## 3. Shape Constraint

To simulate a rigid body, the overall shape must remain unchanged at every time step. The question is how to define such a constraint.

Consider a group of particles $(\mathbf{x}_i, m_i)$ that together form a shape. At rest, the center of mass is

$$
\mathbf{x}_c^{\text{rest}} = \frac{\sum_i m_i\,\mathbf{x}_i^{\text{rest}}}{\sum_i m_i}
$$

After the particles have moved, we want the shape to coincide with the rest shape at all times. The most general spatial change that preserves the shape is a translation $\mathbf{t}$ together with a rotation $R$.

The translation is determined by the center of mass of the new configuration,

$$
\mathbf{x}_c^{\text{new}} = \frac{\sum_i m_i\,\mathbf{x}_i^{\text{new}}}{\sum_i m_i},
\qquad
\mathbf{t} = \mathbf{x}_c^{\text{new}} - \mathbf{x}_c^{\text{rest}}
$$

so that the target position of particle $i$ is

$$
\mathbf{x}_i^{\text{target}} = R\,(\mathbf{x}_i^{\text{rest}} - \mathbf{x}_c^{\text{rest}}) + \mathbf{x}_c^{\text{rest}} + \mathbf{t}
= R\,(\mathbf{x}_i^{\text{rest}} - \mathbf{x}_c^{\text{rest}}) + \mathbf{x}_c^{\text{new}}
$$

The new positions $\mathbf{x}_i^{\text{new}}$ may be scattered arbitrarily; the problem is to find the optimal $R$ and $\mathbf{t}$ that map the rest shape onto them.

### 3.1 The Least-Squares Problem

We define the loss

$$
\mathcal{L}(R) = \sum_i m_i \,\big\| R\,(\mathbf{x}_i^{\text{rest}} - \mathbf{x}_c^{\text{rest}}) + \mathbf{x}_c^{\text{new}} - \mathbf{x}_i^{\text{new}} \big\|^{2}
$$

Introducing the relative vectors

$$
A_i = \mathbf{x}_i^{\text{rest}} - \mathbf{x}_c^{\text{rest}},
\qquad
B_i = \mathbf{x}_i^{\text{new}} - \mathbf{x}_c^{\text{new}}
$$

the problem becomes

$$
\min_{R}\ \mathcal{L}(R) = \sum_i m_i\,\|R A_i - B_i\|^{2}
$$

*Remark:* adding any constant vector to every $B_i$ (for instance, defining $B_i = \mathbf{x}_i^{\text{new}} - \mathbf{t}$ instead) does not affect the optimal rotation, because $\sum_i m_i A_i = \mathbf{0}$.

### 3.2 Reduction to a Trace Maximization

Expanding the squared norm,

$$
\mathcal{L}(R) = \sum_i m_i\,(R A_i - B_i)^{\top}(R A_i - B_i)
= \sum_i m_i\left( A_i^{\top} R^{\top} R A_i - A_i^{\top} R^{\top} B_i - B_i^{\top} R A_i + B_i^{\top} B_i \right)
$$

Since $R^{\top}R = I$, and since $A_i^{\top}R^{\top}B_i$ is a scalar and therefore satisfies
$A_i^{\top}R^{\top}B_i = \left(A_i^{\top}R^{\top}B_i\right)^{\top} = B_i^{\top}R A_i$, the loss simplifies to

$$
\mathcal{L}(R) = \sum_i m_i\left( A_i^{\top}A_i - 2\,B_i^{\top}R A_i + B_i^{\top}B_i \right)
$$

Minimizing $\mathcal{L}$ over $R$ is thus equivalent to finding the rotation $R^{\star}$ that maximizes

$$
\sum_i m_i\,B_i^{\top} R^{\star} A_i
= \operatorname{tr}\!\left( R^{\star} \sum_i m_i\,A_i B_i^{\top} \right)
$$

Defining

$$
A := \left(\sum_i m_i\,A_i B_i^{\top}\right)^{\top}
$$

the objective is to maximize $\operatorname{tr}(R^{\star} A^{\top})$.

### 3.3 Optimal Rotation via Polar Decomposition

Polar decomposition factorizes $A$ as

$$
A = R\,S
$$

with $R$ a rotation matrix and $S$ a symmetric matrix. If $\det(A) \neq 0$, the decomposition is unique.

**Practical considerations.** One must ensure $\det(A) \neq 0$. If $\det(A) = 0$ occurs during simulation, a common fallback is to reuse the rotation from the previous step. The case $\det(A) < 0$ requires special handling (a standard remedy is to negate the column of $U$ associated with the smallest singular value in the SVD of $A$).

Using $A = RS$ and the cyclicity of the trace,

$$
\operatorname{tr}(R^{\star}A^{\top})
= \operatorname{tr}(R^{\star}S^{\top}R^{\top})
= \operatorname{tr}(R^{\star}S R^{\top})
= \operatorname{tr}(R^{\top}R^{\star}S)
$$

The product $G := R^{\top}R^{\star}$ is itself a rotation matrix. Since $S$ is symmetric, it can be diagonalized as $S = Q\Sigma Q^{\top}$ with $Q$ orthonormal, so we must maximize

$$
\operatorname{tr}(G\,Q\Sigma Q^{\top}) = \operatorname{tr}(Q^{\top} G\,Q\,\Sigma)
$$

Both $Q$ and $G$ are orthonormal (as $R^{\star}$ and $R^{\top}$ are orthonormal), hence $Q^{\top}GQ$ is orthonormal as well. Its columns therefore have unit Euclidean norm, so each of its diagonal entries is at most $1$. Since the diagonal entries of $\Sigma$ are non-negative, the trace is maximized exactly when

$$
Q^{\top}GQ = I
\;\Longleftrightarrow\;
Q^{\top}R^{\top}R^{\star}Q = I
\;\Longleftrightarrow\;
R^{\top}R^{\star} = I
\;\Longleftrightarrow\;
R^{\star} = R
$$

That is, the optimal rotation is precisely the rotational factor of the polar decomposition of $A$.

### 3.4 Constructing the Polar Factor via SVD

Let the singular value decomposition of $A$ be $A = U\Sigma V^{\top}$, where $U$ and $V$ are orthonormal, so that $V^{\top}V = I$. Then

$$
A = U\,\underbrace{V^{\top}V}_{I}\,\Sigma V^{\top}
= \underbrace{\left(U V^{\top}\right)}_{R_{\text{polar}}}\;\underbrace{\left(V\Sigma V^{\top}\right)}_{S}
$$

Hence

$$
R_{\text{polar}} = U V^{\top},
\qquad
S = V\Sigma V^{\top}
$$

### 3.5 Newton–Schulz Iteration

Rather than computing an SVD, which is expensive, the polar factor can be obtained with the Newton–Schulz iteration, which is fast and easy to implement:

$$
R_{k+1} = \frac{1}{2}\left(R_k + (R_k^{-1})^{\top}\right),
\qquad R_0 = A
$$

Writing $R_k = U\Sigma_k V^{\top}$, we have $(R_k^{-1})^{\top} = U\Sigma_k^{-1}V^{\top}$, and therefore

$$
R_{k+1} = U\left[\frac{1}{2}\left(\Sigma_k + \Sigma_k^{-1}\right)\right]V^{\top}
$$

From iteration $k$ to $k+1$, only the middle factor changes, and it converges toward the identity. After roughly 6–9 iterations, $R_k \approx UV^{\top}$, which is exactly the polar rotation we need.

