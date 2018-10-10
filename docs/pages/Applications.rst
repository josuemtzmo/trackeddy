============
Applications
============

This code was developed with the goal of analyse and extract the energy
contained by Eddies in the ocean. However it can be used for a variety
of processes. To know more about these processes, please se below:

It may be implemented to track eddies and some other applications, and it could
be capable of track eddies in pressure levels, isopycnals and vorticity.

Kinetic Energy decomposition
----------------------------

In many parts of the ocean, transient processes contain more kinetic energy (commonly known as eddy kinetic energy, or :math:`EKE`) than the men kinetic energy (:math:`MKE`). Variations in the oceanic energy balance can be diagnosed through :math:`EKE`, which allows the analysis  of positive or negative feedbacks on climate change. However, one difficulty in determining the role of eddies in the ocean transient adjustment to climate change is the lack of a process-based definition of :math:`EKE`.

The aim of this algorithm is to decompose and analyse :math:`EKE` according to different ocean processes. Specifically, the separation of kinetic energy will be recustructed using a 2D gaussian fitting for each closedd eddy detected (:math:`EKE_{eddy}`) from the eddy kinetic energy due to meandering jets (:math:`EKE_{jets}`) and the background eddy kinetic energy (:math:`EKE_{background}`):

.. math::
   EKE = EKE_{eddy} + \underbrace{EKE_{jets} + EKE_{background}}_{EKE_{residual}}
..

However, this decomposition represents several challenges like:

- Second order terms which maybe important in the energy balance.

.. math::
   KE = MKE + EKE
..

Expanding this equation we obtain:

.. math::
   KE = MKE + EKE_{eddy} + \underbrace{EKE_{jets} + EKE_{background}}_{EKE_{residual}} + EKE'_{O^1}
..

Replacing the kinetic energy definition:

.. math::
   \hspace{-3cm}\frac{1}{2}\rho_0 (u^2+v^2) = \frac{1}{2}\rho_0 (\bar{u}^2 + \bar{v}^2) + \frac{1}{2}\rho_0 (u_{eddy}^2 + v_{eddy}^2) +
..

.. math::
   \hspace{2.7cm}\frac{1}{2}\rho_0 (u_{jets}^2 + v_{jets}^2) + \frac{1}{2}\rho_0 (u_{background}^2 + v_{background}^2) +
..

.. math::
   \hspace{0.6cm}\rho_0 (\bar{u}u_{eddy} + \bar{v}v_{eddy}) + \rho_0 (\bar{u}u_{jets} + \bar{v}v_{jets}) +
..

.. math::
   \hspace{4cm}\rho_0 (\bar{u}u_{background} + \bar{v}v_{background}) + \rho_0 (u_{eddy}u_{jets} + v_{eddy}v_{jets}) +
..

.. math::
   \hspace{0cm} \rho_0 (u_{eddy}u_{background} + v_{eddy}v_{background}) +
..

.. math::
   \hspace{-0.6cm} \rho_0 (u_{jets}u_{background} + v_{jets}v_{background})
..

where :math:`u = \bar{u} + u_{eddy} + u_{jet} + u_{background}`. Assuming :math:`\iff` :math:`\bar{EKE'_{O^1}} \rightarrow 0` :math:`\implies` we can ingore those terms. However, this implications is really hard to prove unless we define an exact way to extract the velocity field for each eddy.

Because of this, the first approach to this problem will be the decomposition of the components.


Heatwaves
---------



Oceanic Tracers
---------------


Regional Studies
----------------
