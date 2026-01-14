# SPG (Spin Point Groups) 库

> [English Version](README.md)

这是一个基于 Lean 4 的物理计算库。它的核心目标非常明确：**用绝对精确的代数方法，解决凝聚态物理中的磁群对称性分析问题。**

不同于常见的数值计算代码（Python/MATLAB），这里没有浮点数，没有 `1.0000001` 这种误差。所有的矩阵、向量、系数都建立在有理数域（`ℚ`）上。这意味着当你问程序“这个项是不是严格禁止的？”时，它给出的 `true` 是数学上可证明的真。

## 这个库解决什么问题？

在研究磁性材料（如反铁磁、Altermagnet）时，我们经常需要回答以下问题：

1.  **对称性枚举**：给定了晶格结构和磁序（比如磁矩指向），这个体系到底还剩什么对称性？（求磁点群 MPG）
2.  **不变量分析**：在这个对称性下，哈密顿量里允许出现哪些项？（例如：是否允许 $k_x \sigma_y$ 这种自旋分裂项？是否允许 $d$-wave 形式的能带分裂？）
3.  **宏观物性**：这个磁结构是否允许产生净磁矩（铁磁）？是否允许产生电极化（铁电）？

SPG 库提供了一套自动化的工作流来回答这些问题。

---

## 核心设计：它是如何工作的？

要使用这个库，你只需要理解它是如何对物理世界建模的。

### 1. 什么是“磁群操作”？

在代码中，一个群元素 (`SPGElement`) 被定义为一对矩阵：

$$g = (R, S)$$

*   **`spatial` (R)**: 描述原子位置如何移动的 $3 \times 3$ 矩阵（旋转、镜面、反演）。
*   **`spin` (S)**: 描述时间反演操作。
    *   在本库目前的设计中，我们只关心“是否包含时间反演”。
    *   因此 $S$ 只有两种取值：`I` (单位阵，无时间反演) 或 `-I` (含时间反演)。

这种设计对应了磁群理论中的 Corepresentation 思想。群乘法就是简单的矩阵分量相乘。

### 2. 物理量怎么变换？（严查事实）

当你对一个物理系统施加对称操作 $g$ 时，不同的物理量有不同的变法。这个库严格区分了以下几类对象，并在源码中硬编码了它们的变换规则：

*   **极矢量 (Polar Vector)**，如位置 $r$、电场 $E$、动量 $k$：
    *   变换规则：$v' = R v$。
    *   *注意*：对于动量 $k$，如果操作包含时间反演（$T$），动量会反向。所以在处理 $k \cdot p$ 哈密顿量时，库会自动处理 $k \to -k$ 的变号。
*   **轴矢量 (Axial Vector)**，如磁矩 $M$、自旋 $\sigma$：
    *   变换规则：$v' = (\det R) \cdot (R v)$。
    *   因子 $\det R$ 保证了它们在空间反演（$R=-I$）下是不变的。
    *   *时间反演*：磁矩和自旋在时间反演下都会变号。代码中通过判断 $S$ 是否为 `-I` 来自动施加这个负号。

### 3. 如何寻找不变量？

这是本库最强大的功能。寻找“允许的哈密顿量项”本质上是在解一个巨大的线性方程组。

比如，你想知道是否允许 $c \cdot k_x \sigma_y$ 这样的项。库会这样做：
1.  构建一个包含所有可能的 $k$ 二次项和 $\sigma$ 分量的基底空间。
2.  对群里的每一个元素 $g$，计算它把基底变成了什么（比如 $k_x \to -k_y, \sigma_y \to \sigma_x$）。
3.  建立方程：$g \cdot H \cdot g^{-1} = H$。
4.  在有理数域上解这个线性方程组（求 Nullspace），得到所有线性独立的解。

---

## 上手指南：从零分析一个磁性体系

假设你正在研究一种四方晶系的反铁磁体，你想知道它是不是最近很火的 **Altermagnet**（交替磁体）。

### 第一步：引入库

```lean
import SPG

-- 引入所有必要的命名空间
open SPG
open SPG.Interface         -- 也就是 Op[...] 这种记号
open SPG.Algebra           -- 群运算
open SPG.Physics.Hamiltonian -- k·p 理论相关
open SPG.Geometry.SpatialOps -- 预定义的矩阵（如 mat_4_z）
```

### 第二步：定义你的对称性

你需要告诉程序这个晶体有哪些对称操作。对于 D4h 这种群，通常只需要写出生成元。

```lean
-- 定义一个包含时间反演的操作：C4z * T
-- Op[ 空间矩阵, 自旋部分 ]
-- ^-1 代表包含时间反演，^1 代表不包含
def g_C4z_T : SPGElement := Op[mat_4_z, ^-1]

-- 定义普通的 C2x 旋转（不含时间反演）
def g_C2x   : SPGElement := Op[mat_2_x, ^1]

-- 定义空间反演（不含时间反演）
def g_Inv   : SPGElement := Op[mat_inv, ^1]
```

### 第三步：生成群

程序会自动算出这几个生成元能生出的所有群元素（群闭包）。

```lean
def my_group : List SPGElement := generate_group [g_C4z_T, g_C2x, g_Inv]
#eval my_group.length  -- 看看群里有多少个元素
```

### 第四步：问物理问题

现在有了群，我们可以问它任何问题。

**问题 1：这个对称性允许 $d$-wave 的自旋分裂吗？**

Altermagnet 的标志性特征是允许像 $(k_x^2 - k_y^2)\sigma_z$ 这样的项。我们可以直接构建这个项并测试它是不是不变量。

```lean
-- 定义多项式 (kx^2 - ky^2)
def p_dwave : Poly := (kx * kx) - (ky * ky)

-- 构造哈密顿量项：多项式 * σz
-- singleTerm 接受一个多项式和一个自旋分量 (.I, .x, .y, .z)
def term_dwave : KPHam := singleTerm p_dwave .z

-- 判定！
#eval isInvariantHam my_group term_dwave
-- 如果输出 true，恭喜你，这是一个 Altermagnet。
```

**问题 2：把所有允许的项都列给我看看？**

你可以让程序自动求解所有 2 阶以下的允许项：

```lean
-- 求解 2 阶以内的所有矢量不变量 (Vector Invariants)
-- 返回结果是按阶数分类的列表
def allowed_terms := invariants_vector_by_degree_solve my_group 2

-- 打印出来看
#eval allowed_terms
```

### 第五步：验证宏观物性（对称性破缺）

假设你给这个系统加了一个沿着 $(1,1,0)$ 方向的磁场（或者磁矩指向这里），对称性会降低成什么样？

```lean
def neel_vector : Vec3 := ![1, 1, 0]

-- 计算磁点群 (MPG)：即原群中能保持这个磁矩不变的子群
def mpg_110 : List SPGElement := get_mpg my_group neel_vector

-- 在这个新群下，允许产生 z 方向的铁电极化吗？
-- 判据：所有元素必须保持 z 方向的极矢量不变
#eval allows_z_polarization mpg_110
```

---

## 源码导读

如果你想深入了解实现细节，可以按以下路径阅读源码（括号内是关键概念）：

1.  **[SPG/Algebra/Basic.lean](file:///home/yizhou/SPG/SPG/Algebra/Basic.lean)**
    *   定义了 `SPGElement` 结构体。
    *   你会看到矩阵是用 `Mathlib` 的 `Matrix` 定义的，元素类型锁死为 `ℚ`。

2.  **[SPG/Algebra/Actions.lean](file:///home/yizhou/SPG/SPG/Algebra/Actions.lean)**
    *   **核心物理逻辑都在这**。
    *   `magnetic_action`: 仔细看那个 `detR`，那就是轴矢量的灵魂。
    *   `electric_action`: 极矢量的标准变换。

3.  **[SPG/Physics/Hamiltonian/Ham.lean](file:///home/yizhou/SPG/SPG/Physics/Hamiltonian/Ham.lean)**
    *   定义了 `KPHam`（k·p 哈密顿量）。
    *   它是如何存储的？其实就是 4 个多项式：一个标量部分（乘单位阵 $I$），三个向量部分（乘 $\sigma_x, \sigma_y, \sigma_z$）。

4.  **[SPG/Physics/Hamiltonian/Invariants.lean](file:///home/yizhou/SPG/SPG/Physics/Hamiltonian/Invariants.lean)**
    *   这是最硬核的部分。
    *   `invariants_vector_by_degree_solve`: 这里实现了把群论问题转化为线性代数问题（求 Nullspace）的算法。

## 局限性

为了保持精确性和可判定性，本库目前做了一些取舍：

*   **只支持点群**：目前不处理平移操作（Space Groups）。
*   **不处理一般自旋旋转**：目前的自旋操作被锁定为“共线磁性”模型（即自旋只有翻转/不翻转，或者简单的轴向变换）。对于需要处理 Spin-Orbit Coupling 极强且涉及复杂自旋旋转的体系，可能需要扩展 `spin` 矩阵的定义。
