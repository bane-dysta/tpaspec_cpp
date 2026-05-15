# tpaspec

一个 C++17 单文件工具，用于从 Multiwfn 风格的 `SOS.txt` 文件计算有限态求和（SOS）近似下的简并双光子吸收（degenerate TPA）谱。

## 功能概览

`tpaspec` 可以完成两种工作流：

1. 直接读取已有的 Multiwfn 风格 `SOS.txt`：

```bash
tpaspec --sos SOS.txt --out mol_tpa
```

2. 从 Gaussian `.fch/.fchk` 和 TDDFT `.log/.out` 自动调用 Multiwfn 生成 `SOS.txt`，然后计算 TPA：

```bash
tpaspec mol.fch mol.log --multiwfn /path/to/Multiwfn --out mol_tpa
```

也支持自动寻找同名文件：

```bash
tpaspec mol.fch
tpaspec mol.log
```

> 注意：笔者对TPA是初学者，不能保证本程序实现的正确性。

## 编译

需要支持 C++17 的编译器。

```bash
g++ -std=c++17 -O2 -Wall -Wextra -pedantic tpaspec.cpp -o tpaspec
```

如果只用 `--sos SOS.txt` 模式，不需要 Multiwfn。但你自己生成SOS.txt还是要Multiwfn，所以需要安装 Multiwfn。
如果希望自动生成 `SOS.txt`，需要安装 Multiwfn，并且 Gaussian 的 `.fch/.fchk` 和 TDDFT 输出文件可用。

## 输入文件格式

程序期望的 `SOS.txt` 为 Multiwfn 进行SOS计算后输出的文本文件：

```text
N
state_index energy_eV
state_index energy_eV
...
i j mux muy muz
i j mux muy muz
...
```

格式：

- 第一行 `N`：激发态数量。
- 随后 `N` 行：每个激发态的编号和激发能，单位为 eV。
- 后续行：偶极矩矩阵元，单位为 atomic unit。
- `0` 表示基态，`1...N` 表示激发态。
- `0 j` 是基态到激发态的跃迁偶极矩。
- `i j` 且 `i,j > 0` 是激发态之间的偶极矩。
- `i i` 是第 `i` 个态的永久偶极矩。
- `0 0` 是基态永久偶极矩。

示例：

```text
40
1    4.998400
2    5.220500
...
0    0    5.266078E-01  -3.939612E+00  -1.607734E+00
0    1    2.517155E-01   4.283881E-01  -3.553494E-01
1    1    6.419194E+00  -2.714914E+00  -3.413476E+00
1    2    3.830755E-02   1.357284E-01  -1.991823E-01
...
```

> 程序假设偶极矩矩阵为实数且对称，读入 `i j` 后会自动令 `mu(i,j) = mu(j,i)`。这对普通非相对论实波函数计算通常成立；如果涉及 SOC、复波函数、磁场或非厄米响应等，这个假设未必适用。

## 原理

### 简并双光子吸收

简并 TPA 指两个相同频率的光子同时被吸收：

$$
E_f - E_0 = 2\omega
$$

程序中默认基态能量为零，因此目标态 `f` 的单光子能量为：

$$
\omega_f = \frac{E_f}{2}
$$

这里 $E_f$ 是总激发能，单位在输入中为 eV，程序内部会转换成 Hartree。

### 2. SOS TPA transition tensor

对每个目标激发态 `f`，程序计算一个二阶张量：

$$
S_{ab}^{0\to f}(\omega_f)
=
\sum_{n=1}^{N}
\frac{
\mu_{0n}^{a}\,\bar{\mu}_{nf}^{b}
+
\mu_{0n}^{b}\,\bar{\mu}_{nf}^{a}
}
{E_n-\omega_f-i\eta}
$$

其中：

- $a,b$ 是笛卡尔方向 $x,y,z$。
- $n$ 是 SOS 中间态，遍历输入文件中的所有激发态。
- $N$ 是输入 `SOS.txt` 中包含的激发态数。
- $E_n$ 是中间态能量。
- $\eta$ 是可选的 SOS 分母阻尼，对应命令行参数 `--eta-sos-eV`。

代码中没有显式遍历基态 $n=0$。这是因为基态项已经和 $n=f$ 项合并成了 difference-dipole 修正。

具体地说，当 $n != f$ 时：

$$
\bar{\mu}_{nf} = \mu_{nf}
$$

当 $n = f$ 时：

$$
\bar{\mu}_{ff} = \mu_{ff} - \mu_{00}
$$

因此 $n=f$ 项为：

$$
\frac{
\mu_{0f}^{a}(\mu_{ff}^{b}-\mu_{00}^{b})
+
\mu_{0f}^{b}(\mu_{ff}^{a}-\mu_{00}^{a})
}
{E_f-\omega_f-i\eta}
$$

在简并 TPA 中，$E_f - omega_f = omega_f$。这个写法等价于把普通 SOS 公式里的 $n=0$ 基态永久偶极项和 $n=f$ 终态永久偶极项合并，能够显式体现 difference dipole。

> 注意：不同文献和程序可能在 $S$ 的定义里放入额外的 $1/2$ 对称化因子。`tpaspec` 使用的是代码中这种 “full SOS tensor” 约定，并通过 `--n-factor` 控制后续 GM 转换的整体因子。和其他程序比较绝对值时，一定要检查这个约定，否则可能差一个 2 或 4 的倍数。

### 4. 旋转平均 TPA strength

得到张量 $S_ab$ 后，程序计算三个旋转不变量：

$$
A = \left|S_{xx}+S_{yy}+S_{zz}\right|^2
$$

$$
B = \sum_{a,b} S_{ab} S_{ab}^{*}
$$

$$
C = \sum_{a,b} S_{ab} S_{ba}^{*}
$$

然后根据光偏振类型计算 TPA strength：

$$
\delta = \frac{F A + G B + H C}{30}
$$

程序支持三种偏振平均：

| `--pol` | 含义 | $(F, G, H)$ |
|---|---|---:|
| `linear` | 平行线偏振，默认 | `(2, 2, 2)` |
| `perpendicular` | 垂直线偏振 | `(-1, 4, -1)` |
| `circular` | 圆偏振 | `(-2, 3, 3)` |

输出文件中的 `delta_au` 就是这个 $\delta$，单位为 atomic unit。

### 5. 从 TPA strength 到 GM

程序使用下面的转换形式：

$$
\sigma_{\mathrm{GM}}(E)
=
N_{\mathrm{factor}}
\frac{\pi^3\alpha a_0^5}{c}
\omega_f^2
\,g(E-E_f)\,
\delta_f
\times 10^{50}
$$

其中：

- $\alpha$ 是精细结构常数。
- $a_0$ 是 Bohr 半径，程序中以 cm 表示。
- $c$ 是光速，程序中以 cm/s 表示。
- $\omega_f = E_f/2$，必须是单个光子的能量，不是总激发能。
- $g(E-E_f)$ 是归一化线形函数，单位为 $Hartree^{-1}$。
- $\delta_f$ 是目标态 `f` 的 TPA strength。
- $N_{\mathrm{factor}}$ 由 `--n-factor` 指定，默认值为 `4.0`。
- $1 GM = 10^{-50} cm^4 s photon^{-1} molecule^{-1}$。

默认 $N_{\mathrm{factor}}=4$ 对应代码采用的 full SOS tensor convention 和常见 single-beam GM 转换写法。如果你要和其他程序或文献严格比较，必须确认它们使用的 `S` 定义、`delta` 定义、光场约定和 GM 转换因子是否相同。

### 6. 谱线展宽

程序先对每个 final state 在共振位置 $\omega = E_f/2$ 计算一次 $S_{ab}$ 和 $\delta_f$，得到 stick spectrum；然后用归一化线形函数把 sticks 展宽成连续谱。

支持两种线形：

#### Lorentzian

$$
g(x) = \frac{1}{\pi}\frac{\gamma}{x^2+\gamma^2}
$$

#### Gaussian

$$
g(x) = \frac{\sqrt{\ln 2}}{\sqrt{\pi}\gamma}
\exp\left[-\ln 2\frac{x^2}{\gamma^2}\right]
$$

这里 `gamma` 是 HWHM。如果命令行给的是 FWHM，程序会自动除以 2：

```bash
--width-eV 0.20 --width-kind fwhm
```

等价于 HWHM = 0.10 eV。

> 注意：`peak_sigma_GM` 是展宽后的峰高，强烈依赖线形和展宽宽度。报告 GM 峰值时必须同时报告 `shape`、`width-eV`、`width-kind` 和 `n-factor`。

## 命令行参数

| 参数 | 默认值 | 说明 |
|---|---:|---|
| `--sos FILE` | 无 | 读取已有 `SOS.txt` |
| `--fch FILE` | 无 | Gaussian formatted checkpoint |
| `--log FILE` | 无 | Gaussian TDDFT log/out 文件 |
| `--multiwfn EXE` | `Multiwfn` | Multiwfn 可执行文件路径 |
| `--out PREFIX` | 输入文件名 + `_tpa` | 输出文件前缀 |
| `--pol linear\|perpendicular\|circular` | `linear` | 偏振平均类型 |
| `--shape lorentzian\|gaussian` | `lorentzian` | 谱线函数 |
| `--width-eV VALUE` | `0.10` | 展宽宽度 |
| `--width-kind hwhm\|fwhm` | `hwhm` | `--width-eV` 的含义 |
| `--eta-sos-eV VALUE` | `0.0` | SOS 分母阻尼 |
| `--n-factor VALUE` | `4.0` | GM 转换整体因子 |
| `--x nm\|exc_ev\|photon_ev` | `nm` | gnuplot 横轴 |
| `--points N` | `2000` | 连续谱网格点数 |
| `--states N` | 全部 | 只输出前 N 个 final states |
| `--grid-min-eV VALUE` | 自动 | 总激发能网格下限 |
| `--grid-max-eV VALUE` | 自动 | 总激发能网格上限 |
| `--dry-run` | 关闭 | 只打印 Multiwfn 脚本和命令，不实际运行 |
| `--keep-work` | 关闭 | 保留 Multiwfn 临时目录 |
| `--no-gnuplot` | 关闭 | 不输出 gnuplot 脚本 |

## 输出文件

假设输出前缀为 `mol_tpa`，程序会生成：

### `mol_tpa.sticks.csv`

每个 final state 一行，包含 stick 信息。

主要列：

| 列名 | 含义 |
|---|---|
| `state` | final state 编号 |
| `E_exc_eV` | 总激发能 `E_f` |
| `photon_eV` | 单光子能量 `E_f/2` |
| `tpa_nm` | 单光子波长 |
| `delta_au` | TPA strength，a.u. |
| `A_au`, `B_au`, `C_au` | 旋转平均不变量 |
| `peak_sigma_GM` | 该态展宽后线中心峰高，GM |
| `nearest_intermediate_state` | 最接近半频共振的中间态 |
| `min_abs_En_minus_omega_eV` | `min |E_n - omega_f|` |
| `polarization` | 偏振类型 |
| `shape` | 线形 |
| `width_eV` | 展宽参数 |
| `width_kind` | HWHM/FWHM |
| `eta_sos_eV` | SOS 分母阻尼 |
| `n_factor` | GM 转换因子 |

### `mol_tpa.curve.csv`

连续谱数据。

主要列：

| 列名 | 含义 |
|---|---|
| `x` | 当前绘图横轴数值 |
| `nm` / `exc_ev` / `photon_ev` | 横轴类型标签 |
| `E_exc_eV` | 总激发能 |
| `photon_eV` | 单光子能量 |
| `tpa_nm` | 单光子波长 |
| `sigma_GM` | 展宽后的 TPA cross section |

### `mol_tpa.gnuplot`

简单的 gnuplot 脚本，可用于画图：

```bash
gnuplot mol_tpa.gnuplot
```

会生成 `mol_tpa.png`。

### `mol_tpa.SOS.txt`

实际使用的 `SOS.txt` 副本。

### `mol_tpa.report.json`

记录本次计算的参数、态数和常数，方便追踪。

### `mol_tpa.multiwfn.out`

只有自动调用 Multiwfn 时生成，保存 Multiwfn 标准输出。

## 技巧
由于SOS会随着计算的态数增加不断收敛，不像响应理论一样想算哪个态算哪个就行，推荐多算一些激发态来确认感兴趣的态是否收敛。

### 1. 先用少量 final states 检查流程

```bash
./tpaspec --sos SOS.txt --states 5 --out test_s1_s5
```

确认输出、单位、波长和峰位置合理。

### 2. 做 SOS 态数收敛测试

例如分别生成包含不同数量激发态的 `SOS.txt`：

```text
20 states
40 states
60 states
80 states
100 states
```

然后比较同一个 final state 的 `delta_au`，例如 S1、S2、S3：

```bash
./tpaspec --sos SOS_20.txt  --states 3 --out tpa_20
./tpaspec --sos SOS_40.txt  --states 3 --out tpa_40
./tpaspec --sos SOS_60.txt  --states 3 --out tpa_60
./tpaspec --sos SOS_80.txt  --states 3 --out tpa_80
./tpaspec --sos SOS_100.txt --states 3 --out tpa_100
```

如果 `delta_au` 随态数变化很大，说明有限 SOS 还没收敛，不能把绝对 GM 当作可靠结果。

### 3. 检查半频共振

如果输出中出现 warning：

```text
WARNING: near-resonant SOS denominators detected; consider --eta-sos-eV.
```

说明某个中间态 `E_n` 非常接近 `E_f/2`。这种情况下 SOS 分母很小，结果会对阻尼和态能非常敏感。

可以尝试：

```bash
./tpaspec --sos SOS.txt --eta-sos-eV 0.02 --out eta002
./tpaspec --sos SOS.txt --eta-sos-eV 0.05 --out eta005
./tpaspec --sos SOS.txt --eta-sos-eV 0.10 --out eta010
```

并报告 `eta-sos-eV`。

### 4. 同时报告 stick strength 和展宽谱

建议论文或报告中至少写清楚：

```text
SOS states included: 例如 100 excited states
TPA tensor convention: full SOS tensor, no extra 1/2 in S_ab
Polarization average: linear / circular / perpendicular
N factor: 4.0
Lineshape: Lorentzian or Gaussian
Width: 例如 HWHM = 0.10 eV
eta_sos: 例如 0.00 eV 或 0.05 eV
```

其中 `delta_au` 比 `peak_sigma_GM` 更适合做方法比较，因为 `peak_sigma_GM` 会随展宽宽度改变。
