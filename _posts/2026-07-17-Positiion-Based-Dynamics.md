<!-- raw before LLM format , lazy work

We give a much more logic clear way to represent the math of PBD and XPBD in this blog.
PBD Paper: https://matthias-research.github.io/pages/publications/posBasedDyn.pdf
XPBD Paper : https://matthias-research.github.io/pages/publications/XPBD.pdf

ABBREviate: PBD (Position based dynamics) and XPBD (Extended PBD) can simulate world physics(rigid body ,soft body, fluid)
using a unified particles physics (Newton's Law no quantum physics involved). 
XPBD framework can be in the future integrated into AI as a basic building block to construct a world model.
The math in the XPBD paper is fussy, we here give a much clear step by step deduction.


the math of (X)PBD is easy :
Given a n of partiles with initial position pi(x) ,massi , 
the newton's moment law gives: Delta Position = Delta T * V = Delta T * F / M  
Simulation : Delta T , Mass are knowns in computer simulation environment , Force is not known.
Force can be split to two categroies : 1. active force (external) 2. inactive force .(internal)
Examples: Active force like people push someobject which can be defined as sytem external force .
Inactive force like friction force , like spring internal force.  
In Reality it is the way how we could define it easily . For External force which can be directly 
defined and for internal force we don't define it directly but deduced from Constraint(x). 
There are serval types of constraint which used to describe different type of object (fluid,rigid,soft)
So for the given many paritiles we can define many constraints with different constraints types.

Example :  
          given two particles in one dimension space . (x1,m1) , (x2,m2) , a distance constraint 
          can be defined as C(x1,x2)= ||x1-x2|| - L , you can think of it as we always want the two 
          particles keep a distance L , if they are too close C < 0 , too far C> 0 , distant L when C=0.
          So we can always change x1 , x2 to let it get balanced C=0 .  Imagine a group of particles with 
          each has a distance constraint to close particles . This kind of structure could represent 
          a softbody. For other type of constraints(rigid body ,fluid ) i will not explain here , 
          you can refer to the paper for details. 

with C(x) , we could define the internal force easily :
the internal lambda_force is try to make the system rebalanced (C(x)=0) from the current position which is 
like hooks' law.  C(Xcurrent + delta x) = - 1/k * Lambda_force  (Hook's law : x = -1/k * F ,negtive sign for restoring force)
Left side use Taylor expansion C(Xcurrent)+ Transpose(gradient(C))* delta X  = 1/k * Lambda_force , 
delta X = M^-1 * gradient(C) * Lambda_force * delta t * delta t   ,  here Lambda_force is a scalar,
gradient(C) * Lambda_force the force vector is in the direction of the gradient which is reasonable.
So put together replace delta x gives :
C(x) + Transpose(gradient(C)) * (M^-1 * gradient(C) * Lambda_force * delta t * delta t ) =  - 1/k * Lambda_force ,
for the paper in original XPBD they define lambda = lambda_force * delta t * delta t , and 1/k = alpha  , 1/(k * delta t * delta t) = alpha~
so we can write :
C(x) + Transpose(gradient(C)) * (M^-1 * gradient(C) * Lambda_force * delta t * delta t ) =  - 1/(k * delta t * delta t) * Lambda_force *(delta t * delta t)
which is :
C(x) + Transpose(gradient(C)) * (M^-1 * gradient(C) * Lambda ) =  - alpha~ * Lambda , 

In the Simulation the Lambda is calculating in K steps loop , for each loop the Position x get updated and 
the Lambda is updated accordingly. The actually logic is you get deltaX as you solve constraint Cj(X) one by one ,
After you have the delta x you are going to solve a delta lambda 

for the Step K you have :
C(Xk+deltax) = C(Xk) + Transpose ( gradient C(x)) *(M^-1* gradient C(Xk)*lambda k *delta t ^ 2) 

solve the lambda k , solve the delta x , get X k+1 = X k 

for the Step K+1 you have : 

C(Xk+deltax + small delta x ) = C(Xk) + Transpose ( Gradient C(Xk))*(M^-1 *Gradient C(Xk)(Lamdba k + small delta Lambda) * delta t ^ 2)

= C(Xk+1) +  Transpose ( Gradient C(Xk))*(M^-1 *Gradient C(Xk)(small delta Lambda) * delta t ^ 2)  =  - alpha * ( Lambda k + small delta Lambda ) * delta t ^ 2

Intersting here :

You can assume : Gradient C(Xk) = Gradient C(Xk+1) or you stick to C(X0) for better performance.

from the last equation :
M^-1 is just a diagnal matrix , redefine the variables with delta t ^2 abbsorbed we finally have (same with the paper) : 

small delta Lambda  =  (-C(Xk+1)-Lambda k * ~alpha)/( Transpose ( Gradient C(Xk)) * M^-1 *Gradient C(Xk) + ~alpha)

--> 
---
layout: post
title: "A Rigorous Mathematical Formulation of PBD and XPBD"
mathjax: true
---


## 1. Fundamentals of Particle Dynamics and Constraints

Consider a system composed of $n$ particles. Each particle $i$ is characterized by its position $\mathbf{x}_i \in \mathbb{R}^3$ and its mass $m_i \in \mathbb{R}$. Let $\mathbf{x} = [\mathbf{x}_1^T, \mathbf{x}_2^T, \dots, \mathbf{x}_n^T]^T \in \mathbb{R}^{3n}$ denote the generalized system coordinate vector. 

From Newton's second law of motion, the relationship governing the incremental displacement $\Delta \mathbf{x}$ over a discrete time step $\Delta t$ is expressed as:
$$\Delta \mathbf{x} = \Delta t \mathbf{v} = \Delta t^2 \mathbf{M}^{-1} \mathbf{F}$$
