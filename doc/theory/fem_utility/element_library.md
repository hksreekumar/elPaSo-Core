(Element_library)=
# Element library
```{important}
Refer [Feature Overview](../feature_overview.md) for availability.
```

## Beam elements
| **Beam**                            | **Nodes** | **Dofs**               | **Material**         |  |
|--------------------------------------|:------:|:--------------------|:--------------|:------:|
| BeamBernoulli | 2 | $x_{1,2}; w_3$      | $E, A, I_y, t$        | <img src = "../../images/elements/line2d.png" width="250px"> |
| BeamTimoshenko           | 2 | $x_{1,2}; w_3$  | $E$, $A$, $I_y$, $\rho$, $\nu$ | <img src = "../../images/elements/line2d.png" width="250px"> |
| BeamBernoulli10           | 2 | $ x_{1,2,3}$; $ w_{2,3}$  | $E$, $A$, $Iy$, $Iz$, $\rho$ | <img src = "../../images/elements/line3d.png" width="250px"> |
| BeamTimoshenko10           | 2 | $ x_{1,2,3}$; $ w_{2,3}$  | $E$, $A$, $Iy$, $Iz$, $\rho$, $\nu$ | <img src = "../../images/elements/line3d.png" width="250px"> |
| BeamBernoulli12           | 2 | $ x_{1,2,3}$; $ w_{1,2,3}$  | $E$, $A$, $Ix$, $Iy$, $Iz$, $\rho$ | <img src = "../../images/elements/line3d.png" width="250px"> |
| BeamTimoshenko12           | 2 | $ x_{1,2,3}$; $ w_{1,2,3}$  | $E$, $A$, $Ix$, $Iy$, $Iz$, $\rho$, $\nu$  | <img src = "../../images/elements/line3d.png" width="250px"> |

## Brick elements
| **Brick**                            | **Nodes** | **Dofs**               | **Material**         |  |
|--------------------------------------|:------:|:--------------------|:--------------|:------:|
| Brick8 | 8 | $ x_{1,2,3}$      | $E$, $\nu$, $\rho$        | <img src = "../../images/elements/hex8.png" width="150px"> |
| Brick20 | 20 | $ x_{1,2,3}$      | $E$, $\nu$, $\rho$        | <img src = "../../images/elements/hex20.png" width="150px"> |
| Brick27 | 27 | $ x_{1,2,3}$      | $E$, $\nu$, $\rho$        | <img src = "../../images/elements/hex27.png" width="150px"> |

## Cable or truss elements

| **Cable**                            | **Nodes** | **Dofs**               | **Material**         |  |
|--------------------------------------|:------:|:--------------------|:--------------|:------:|
| Cable2d | 2 |$ x_{1,2}$      | $E$, $A$, $\rho$, $Fi$        | <img src = "../../images/elements/line2d.png" width="250px"> |
| Cable3d | 2 |$ x_{1,2,3}$      | $E$, $A$, $\rho$, $Fi$        | <img src = "../../images/elements/line3d.png" width="250px"> |

## Plate elements
| **Plate**                            | **Nodes** | **Dofs**               | **Material**         |  |
|--------------------------------------|:------:|:--------------------|:--------------|:------:|
| Kirch4 | 4 | $x_{3}; \frac{dw}{dx}; \frac{dw}{dy}; \frac{dw}{dxy}$ | $E, ν, ρ,t$ | <img src = "../../images/elements/quad4.png" width="150px">  |
| DSG3   | 4 | $x_{3}; w_{1,2}$                                      | $E, ν, ρ,t$ | <img src = "../../images/elements/dummy.png" width="150px">  |
| DSG4   | 4 | $x_{3}; w_{1,2}$                                      | $E, ν, ρ,t$ | <img src = "../../images/elements/quad4.png" width="150px">  |
| DSG9   | 9 | $x_{3}; w_{1,2}$                                      | $E, ν, ρ,t$ | <img src = "../../images/elements/quad9.png" width="150px">  |

## Disc elements
| **Disc** | **Nodes** | **Dofs**               | **Material**         |  |
|--------------------------------------|:------:|:--------------------|:--------------|:------:|
| Disc3                                | 3    | $x_{1,2}$          | $E, ν, ρ,t$  |   <img src = "../../images/elements/dummy.png" width="150px">   |
| Disc4                                | 4    | $x_{1,2}$          | $E, ν, ρ,t$  |   <img src = "../../images/elements/quad4.png" width="150px">   |
| Disc4Dr                              | 4    | $x_{1,2}$ ;$w_{3}$ | $ E, ν, ρ,t$ |   <img src = "../../images/elements/quad4.png" width="150px">   |
| Disc9                                | 9    | $x_{1,2}$          | $E, ν, ρ,t$  |   <img src = "../../images/elements/quad9.png" width="150px">   |
| Disc9s                               | 9    | $x_{1,2}$          | $E, ν, ρ,t$  |   <img src = "../../images/elements/quad9.png" width="150px">   |

## Shell elements

| **Shell** | **Nodes** | **Dofs**               | **Material**         |  |
|--------------------------------------|:------:|:--------------------|:--------------|:------:|
| PlShell3                                | 3    | $ x_{1,2,3}$; $ w_{1,2,3}$          |    |   <img src = "../../images/elements/dummy.png" width="150px">   |
| PlShell4                                | 4    | $ x_{1,2,3}$; $ w_{1,2,3}$          |    |   <img src = "../../images/elements/quad4.png" width="150px">   |
| PlShell4Dr                                | 4    | $ x_{1,2,3}$; $ w_{1,2,3}$          |    |   <img src = "../../images/elements/quad4.png" width="150px">   |
| PlShell9                                | 9    | $ x_{1,2,3}$; $ w_{1,2,3}$          |    |   <img src = "../../images/elements/quad9.png" width="150px">   |

## Fluid elements
| **Fluid2d** | **Nodes** | **Dofs**               | **Material**         |  |
|--------------------------------------|:------:|:--------------------|:--------------|:------:|
| Fluid2d4                              | 4  | fluid |  |  <img src = "../../images/elements/quad4.png" width="150px"> |
| Fluid2d9                              | 9  | fluid |  |  <img src = "../../images/elements/quad9.png" width="150px"> |
| **Fluid3d** | **Nodes** | **Dofs**               | **Material**         |  |
| Fluid4                                | 4  | fluid |  |  <img src = "../../images/elements/quad4.png" width="150px"> |
| Fluid8                                | 8  | fluid |  |  <img src = "../../images/elements/hex8.png" width="150px"> |
| Fluid27                               | 27 | fluid |  |  <img src = "../../images/elements/hex27.png" width="150px"> |

## Fluid flow elements
| **Fluidflow** | **Nodes** | **Dofs**               | **Material**         |  |
|--------------------------------------|:------:|:--------------------|:--------------|:------:|
| FF4                                          | 4 | $p, v_{1}, v_{2}, v_{3}$ | not necessary | <img src = "../../images/elements/quad4.png" width="150px">  |
| FF9                                          | 9 | $p, v_{1}, v_{2}, v_{3}$ | not necessary | <img src = "../../images/elements/quad9.png" width="150px">  |

## Porous elements
| **Poro3dUP** | **Nodes** | **Dofs**               | **Material**         |  |
|--------------------------------------|:------:|:--------------------|:--------------|:------:|
| Poro3dUP8                                   | 8  | $ x_{1,2,3}; fluid$                         |  |  <img src = "../../images/elements/dummy.png" width="50px"> |
| Poro3dUP27                                  | 27 | $x_{1,2,3}; fluid$                          |  |  <img src = "../../images/elements/dummy.png" width="50px"> |
| **Poro3dUU** | **Nodes** | **Dofs**               | **Material**         |  |
| Poro3dUU8                                   | 8  | $ x_{1,2,3}; w_{1,2,3}$                     |  |  <img src = "../../images/elements/dummy.png" width="50px"> |
| **PoroDisc** | **Nodes** | **Dofs**               | **Material**         |  |
| DiscQ2P1                                    |   | $ x_{1,2}; w_{3}; pore_{0}$                 |  |  <img src = "../../images/elements/dummy.png" width="50px"> |
| **PoroPlate** | **Nodes** | **Dofs**               | **Material**         |  |
| PlateQ2P1                                   |   | $ x_{3}; w_{1,2}; pore_{1}$                 |  |  <img src = "../../images/elements/dummy.png" width="50px"> |
| DiscQ2P1                                    |   | $ x_{1,2}; w_{3}; pore_{0}$                 |  | <img src = "../../images/elements/dummy.png" width="50px">  |
| **PoroPlate** | **Nodes** | **Dofs**               | **Material**         |  |
| PlateQ2P1                                   |   | $ x_{3}; w_{1,2}; pore_{1}$                 |  |  <img src = "../../images/elements/dummy.png" width="50px"> |
| PlateQ2P2                                   |   | $ x_{3}; w_{1,2}; pore_{1}; xd3$            |  | <img src = "../../images/elements/dummy.png" width="50px">  |
| PlateQ2P3                                   |   | $ x_{3}; w_{1,2}; pore_{1,3};xd3; wd_{1,2}$ |  | <img src = "../../images/elements/dummy.png" width="50px">  |
| PlateQ3P1                                   |   | $ x_{3}; w_{1,2}; pore_{1}$                 |  | <img src = "../../images/elements/dummy.png" width="50px">  |
| PlateQ3P3                                   |   | $ x_{3}; w_{1,2}; pore_{1,3};xd3; wd_{1,2}$ |  |  <img src = "../../images/elements/dummy.png" width="50px"> |
| **PoroShell** | **Nodes** | **Dofs**               | **Material**         |  |
| ShellQ2P1                                   |   | $ x_{3}; w_{1,2}; pore_{0,1}$               |  | <img src = "../../images/elements/dummy.png" width="50px">  |
