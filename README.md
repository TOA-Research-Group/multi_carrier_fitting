## Usage
It accepts PPMS output files as input. Please specify the channels (x) you used for measuring the longitudinal resistance and the transverse resistance. The program read the columns of the magnetic filed and Bridge x Resistances.

```python
popt, pcov = HRMR_PPMS.multi_carrier_fit(h_, sigma_xx, sigma_yx , p0=[1e26, 1e-2, -1e26, -1e-2])
```
h_: np.array of magnetic fields , p0: Initial values ​​of carrier parameters. Here, carriers are defined in the order of electrons and holes, with a carrier concentration of 1e26\[/m^3\] and a mobility of 1e-2\[m^2/Vs\].

See example.py for further usage.

## Physics
When multiple bands contribute to conduction, it is difficult to estimate the mobility and carrier concentration solely from the Hall voltage. Based on the standard multi-carrier model, the longitudinal conductivity $\sigma_{xx}$ and transverse conductivity $\sigma_{xy}$ are expressed by the following equations []:

$$\sigma_{xx}(\mathbf{B}) = \sum_i \frac{n_i e \mu_i}{1 + (\mu_i B)^2}$$
$$\sigma_{yx}(\mathbf{B}) = \sum_i \frac{s_i n_i e \mu_i^2 B}{1 + (\mu_i B)^2}$$

Here, $i$ is the carrier index, $e$ is the elementary charge, $\mu_i$ is the mobility, and $s_i$ takes the value $-1$ for holes and $+1$ for electrons. By using conductivity, the total conductivity can be expressed as the sum of the conductivities for each carrier type. Since the values obtained from standard measurements are resistances, such as longitudinal resistance and Hall resistance, they must be converted into conductivities. @schematicHallBar shows a typical setup for Hall measurements.

![](./HallDiagram.png)

Here, the $x, y, z$ axes are defined as shown in the figure. The longitudinal voltage is $V_x$ and the Hall voltage is $V_y$. Let $I$ be the current, $W$ the width, $L$ the length, and $t$ the thickness of the Hall bar. The resistivity is expressed as $E_i = \rho_{ij}(\mathbf{B}) J_j$, with $\mathbf{E}$ as the electric field and $\mathbf{J}$ as the current density. In the following equations, the $(\mathbf{B})$ notation is omitted. Using these variables, the longitudinal resistivity and Hall resistivity can be calculated from the following equations. A geometry with $L/W = 4$ is commonly used.
$$E_x = V_x / L = \rho_{xx} J_x = \rho_{xx} \frac{I}{W t}$$
$$E_y = V_y / W = \rho_{yx} J_x = \rho_{yx} \frac{I}{W t}$$
Conductivity is expressed as $J_i = \sigma_{ij} E_j$, and the relation $\sigma = \rho^{-1}$ holds. We expand $\rho$ as:
$$\rho = \begin{pmatrix} \rho_{xx} & \rho_{xy} & \rho_{xz} \\ \rho_{yx} & \rho_{yy} & \rho_{yz} \\ \rho_{zx} & \rho_{zy} & \rho_{zz} \end{pmatrix}$$
Here, we consider the case where the magnetic field is applied in the $z$-direction, as in the setup of @schematicHallBar. In this situation, $\rho$ can be expressed as:

$$\rho = \begin{pmatrix} \rho_{xx} & \rho_{xy} & 0 \\ \rho_{yx} & \rho_{yy} & 0 \\ 0 & 0 & \rho_{zz} \end{pmatrix}$$

Applying the formula for the inverse of a 3x3 matrix:

$$\rho^{-1} = \frac{1}{(\rho_{xx}\rho_{yy}-\rho_{xy}\rho_{yx})\rho_{zz}} \begin{pmatrix} \rho_{yy}\rho_{zz} & -\rho_{xy}\rho_{zz} & 0 \\ -\rho_{yx}\rho_{zz} & \rho_{xx}\rho_{zz} & 0 \\ 0 & 0 & \rho_{xx}\rho_{yy}-\rho_{xy}\rho_{yx} \end{pmatrix}$$

Therefore, the longitudinal and transverse components of the conductivity $\sigma = \rho^{-1}$ are expressed in terms of resistivity as:
$$\sigma_{xx} = \frac{\rho_{yy}}{\rho_{xx}\rho_{yy}-\rho_{xy}\rho_{yx}}$$
$$\sigma_{xy} = \frac{-\rho_{xy}}{\rho_{xx}\rho_{yy}-\rho_{xy}\rho_{yx}}$$
$$\sigma_{yx} = \frac{-\rho_{yx}}{\rho_{xx}\rho_{yy}-\rho_{xy}\rho_{yx}}$$
$$\sigma_{yy} = \frac{\rho_{xx}}{\rho_{xx}\rho_{yy}-\rho_{xy}\rho_{yx}}$$

These equations require four resistivity components ($\rho_{xx}, \rho_{yy}, \rho_{xy}, \rho_{yx}$). However, in the Hall bar configuration shown in @schematicHallBar, only $\rho_{xx}$ and $\rho_{yx}$ can be measured. While one method is to perform measurements on a Hall bar rotated by 90 degrees, this is often unnecessary due to the following relations.
### Onsagar's reciprocal relations
The Onsager theorem states that the off-diagonal transport coefficients $L_{ij}$ and $L_{ji}$ are equal ($L_{ij} = L_{ji}$) when time-reversal symmetry is present. The Onsager relations must be modified in systems where time-reversal symmetry is broken, such as when a magnetic field is applied. In the case with $\mathbf{B} = (0, 0, B_z)^T$, the Onsager reciprocal relation is given by []:

$$\sigma_{yx}(\mathbf{B}) = \sigma_{xy}(-\mathbf{B}) = -\sigma_{xy}(\mathbf{B})$$

### Isotropic resistivity
This is not strictly true in some cases depending on the symmetry.
$$\rho_{xx} = \rho_{yy}$$

By applying these relations, the conductivity can be expressed solely in terms of $\rho_{xx}$ and $\rho_{yx}$ as:

$$\sigma_{xx} = \frac{\rho_{xx}}{\rho_{xx}^2 + \rho_{yx}^2}$$
$$\sigma_{yx} = \frac{-\rho_{yx}}{\rho_{xx}^2 + \rho_{yx}^2}$$