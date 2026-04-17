# 手算带权 POD

> 参考资料：
> 
> [康斯坦茨大学数学与统计系 POD 讲义](https://www.math.uni-konstanz.de/numerik/personen/volkwein/teaching/POD-Book.pdf)
> 
> [Stanford University CME345 课程讲义](https://web.stanford.edu/group/frg/course_work/CME345/CA-AA216-CME345-Ch4.pdf)

在带权 POD 里，内积写成

$$
\langle u,v\rangle_W=u^\top Wv,
$$

其中 $W$ 是对称正定矩阵。累计能量通常写成

$$
E(\ell)=\frac{\sum_{i=1}^{\ell}\lambda_i}{\sum_i \lambda_i}.
$$

离散 POD 既可以写成特征值问题，也可以和 SVD 对应起来。

---

## 1. 原始快照

我们取 4 个二维快照：

$$
x_1=
\begin{bmatrix}
2 \\
1
\end{bmatrix},
\quad
x_2=
\begin{bmatrix}
0 \\
1
\end{bmatrix},
\quad
x_3=
\begin{bmatrix}
1 \\
2
\end{bmatrix},
\quad
x_4=
\begin{bmatrix}
1 \\
0
\end{bmatrix}.
$$

把它们排成快照矩阵：

$$
X=
\begin{bmatrix}
2 & 0 & 1 & 1 \\
1 & 1 & 2 & 0
\end{bmatrix}.
$$

这里每一列是一个快照，这正是 POD 里最常见的离散组织方式。

---

## 2. 先减均值

先求均值快照：

$$
\bar x=\frac{1}{4}(x_1+x_2+x_3+x_4)=
\frac{1}{4}
\begin{bmatrix}
2+0+1+1 \\
1+1+2+0
\end{bmatrix}=
\begin{bmatrix}
1 \\
1
\end{bmatrix}.
$$

于是中心化后的快照为

$$
y_j=x_j-\bar x.
$$

逐个算得：

$$
y_1=
\begin{bmatrix}
1 \\
0
\end{bmatrix},
\quad
y_2=
\begin{bmatrix}
-1 \\
0
\end{bmatrix},
\quad
y_3=
\begin{bmatrix}
0 \\
1
\end{bmatrix},
\quad
y_4=
\begin{bmatrix}
0 \\
-1
\end{bmatrix}.
$$

所以中心化快照矩阵是

$$
Y=
\begin{bmatrix}
1 & -1 & 0 & 0 \\
0 & 0 & 1 & -1
\end{bmatrix}.
$$

流体问题里常常先减均值，再对波动部分做 POD。数学上就是把快照矩阵从 $X$ 换成中心化后的 $Y$。

---

## 3. 加权内积

现在我们不使用普通 Euclidean 内积，而是选

$$
W=
\begin{bmatrix}
2 & 0 \\
0 & 1
\end{bmatrix}.
$$

于是带权内积为

$$
\langle u,v\rangle_W=u^\top Wv,
$$

带权范数为

$$
\|u\|_W^2=u^\top Wu.
$$

这表示：**第一个分量的权重是第二个分量的 2 倍**。你可以把它理解为第一个自由度更“重要”，或者对应更大的单元体积、积分权重或物理能量权重。

例如对任意

$$
u=
\begin{bmatrix}
a \\
b
\end{bmatrix},
$$

有

$$
\|u\|_W^2=
\begin{bmatrix}
a & b
\end{bmatrix}
\begin{bmatrix}
2 & 0 \\
0 & 1
\end{bmatrix}
\begin{bmatrix}
a \\
b
\end{bmatrix}
=2a^2+b^2.
$$

---

## 4. 第一 POD 模态：直接从优化问题手算

带权 POD 的第一模态，是解下面这个问题：

$$
\max_{\phi^\top W\phi=1}
\sum_{j=1}^4 \langle y_j,\phi\rangle_W^2.
$$

这就是“在带权意义下，让总投影能量最大”的方向。它和普通 POD 完全同构，只是把内积换成了 $W$-内积。

设

$$
\phi=
\begin{bmatrix}
a \\
b
\end{bmatrix}.
$$

### 4.1 约束条件

因为要求 $\phi$ 是 $W$-单位向量，所以

$$
\phi^\top W\phi=1
\quad\Longrightarrow\quad
2a^2+b^2=1.
$$

### 4.2 目标函数

依次计算 4 个快照在 $\phi$ 上的带权投影：

$$
\langle y_1,\phi\rangle_W=
\phi^\top Wy_1=
\begin{bmatrix}
a & b
\end{bmatrix}
\begin{bmatrix}
2 \\
0
\end{bmatrix}
=2a,
$$

$$
\langle y_2,\phi\rangle_W=-2a,
\qquad
\langle y_3,\phi\rangle_W=b,
\qquad
\langle y_4,\phi\rangle_W=-b.
$$

所以目标函数为

$$
J(a,b)
=(2a)^2+(-2a)^2+b^2+(-b)^2
=8a^2+2b^2.
$$

利用约束 $2a^2+b^2=1$，可写成

$$
b^2=1-2a^2,
$$

因此

$$
J(a,b)=8a^2+2(1-2a^2)=4a^2+2.
$$

而由 $2a^2+b^2=1$ 可知

$$
0\le a^2\le \frac{1}{2}.
$$

所以 $J$ 在 $a^2=\frac{1}{2}$ 时最大。

于是

$$
a=\pm\frac{1}{\sqrt{2}},
\qquad
b=0.
$$

模态正负号等价，我们取正号：

$$
\boxed{
\phi_1=
\begin{bmatrix}
\frac{1}{\sqrt{2}} \\
0
\end{bmatrix}
}
$$

对应最大值

$$
\lambda_1=J_{\max}=4.
$$

---

## 5. 第二 POD 模态

第二模态要满足两个条件：

1. $W$-单位长度：

$$
\phi_2^\top W\phi_2=1.
$$

2. 与第一模态 $W$-正交：

$$
\phi_2^\top W\phi_1=0.
$$

设

$$
\phi_2=
\begin{bmatrix}
a \\
b
\end{bmatrix}.
$$

先看正交条件：

$$
\phi_2^\top W\phi_1=
\begin{bmatrix}
a & b
\end{bmatrix}
\begin{bmatrix}
2 & 0 \\
0 & 1
\end{bmatrix}
\begin{bmatrix}
\frac{1}{\sqrt{2}} \\
0
\end{bmatrix}=
\sqrt{2}\,a
=0.
$$

所以

$$
a=0.
$$

再用单位长度条件：

$$
2a^2+b^2=1
\quad\Longrightarrow\quad
b^2=1.
$$

取正号：

$$
\boxed{
\phi_2=
\begin{bmatrix}
0 \\
1
\end{bmatrix}
}
$$

对应能量为

$$
\lambda_2=
\sum_{j=1}^4 \langle y_j,\phi_2\rangle_W^2
=0^2+0^2+1^2+(-1)^2
=2.
$$

---

## 6. 检查 $W$-正交归一

现在验证一下：

$$
\phi_1^\top W\phi_1=
\begin{bmatrix}
\frac{1}{\sqrt{2}} & 0
\end{bmatrix}
\begin{bmatrix}
2 & 0 \\
0 & 1
\end{bmatrix}
\begin{bmatrix}
\frac{1}{\sqrt{2}} \\
0
\end{bmatrix}
=1,
$$

$$
\phi_2^\top W\phi_2=
\begin{bmatrix}
0 & 1
\end{bmatrix}
\begin{bmatrix}
2 & 0 \\
0 & 1
\end{bmatrix}
\begin{bmatrix}
0 \\
1
\end{bmatrix}
=1,
$$

$$
\phi_1^\top W\phi_2=0.
$$

所以这两个模态确实构成一组 **$W$-正交归一基**。

---

## 7. 用矩阵方法复核一遍

把刚才的结果用矩阵写法再验一次。

先算

$$
YY^\top=
\begin{bmatrix}
1 & -1 & 0 & 0 \\
0 & 0 & 1 & -1
\end{bmatrix}
\begin{bmatrix}
1 & 0 \\
-1 & 0 \\
0 & 1 \\
0 & -1
\end{bmatrix}=
\begin{bmatrix}
2 & 0 \\
0 & 2
\end{bmatrix}.
$$

对于带权 POD，把普通推导中的内积换成 $W$-内积后，可得到广义特征值问题

$$
WYY^\top W\phi=\lambda W\phi,
$$

若 $W$ 可逆，则等价于

$$
YY^\top W\phi=\lambda \phi.
$$

代入本例：

$$
YY^\top W=
\begin{bmatrix}
2 & 0 \\
0 & 2
\end{bmatrix}
\begin{bmatrix}
2 & 0 \\
0 & 1
\end{bmatrix}=
\begin{bmatrix}
4 & 0 \\
0 & 2
\end{bmatrix}.
$$

所以特征值问题就是

$$
\begin{bmatrix}
4 & 0 \\
0 & 2
\end{bmatrix}
\phi=
\lambda \phi.
$$

显然有：

- $\lambda_1=4$，对应方向 $[1,0]^\top$，再做 $W$-归一化后得到

$$
\phi_1=
\begin{bmatrix}
\frac{1}{\sqrt{2}} \\
0
\end{bmatrix}.
$$

- $\lambda_2=2$，对应方向 $[0,1]^\top$，它本身已经是 $W$-归一化的。

和前面从优化问题直接算出的结果完全一致。

---

## 8. 模态系数

带权 POD 里，第 $j$ 个快照在第 $i$ 个模态上的系数是

$$
a_{ij}=\langle y_j,\phi_i\rangle_W=\phi_i^\top Wy_j.
$$

矩阵写法是

$$
A=\Phi^\top W Y.
$$

### 第一模态的系数

$$
a_{11}=\langle y_1,\phi_1\rangle_W=\sqrt{2},
\quad
a_{12}=\langle y_2,\phi_1\rangle_W=-\sqrt{2},
\quad
a_{13}=0,
\quad
a_{14}=0.
$$

### 第二模态的系数

$$
a_{21}=0,
\quad
a_{22}=0,
\quad
a_{23}=1,
\quad
a_{24}=-1.
$$

所以系数矩阵是

$$
A=
\begin{bmatrix}
\sqrt{2} & -\sqrt{2} & 0 & 0 \\
0 & 0 & 1 & -1
\end{bmatrix}.
$$

---

## 9. 用模态把快照重构出来

POD 展开是

$$
y_j=\sum_{i=1}^{2} a_{ij}\phi_i.
$$

逐个验证。

### 对 $y_1$

$$
y_1=\sqrt{2}\,\phi_1+0\cdot \phi_2=
\sqrt{2}
\begin{bmatrix}
\frac{1}{\sqrt{2}} \\
0
\end{bmatrix}=
\begin{bmatrix}
1 \\
0
\end{bmatrix}.
$$

### 对 $y_2$

$$
y_2=-\sqrt{2}\,\phi_1=
\begin{bmatrix}
-1 \\
0
\end{bmatrix}.
$$

### 对 $y_3$

$$
y_3=\phi_2=
\begin{bmatrix}
0 \\
1
\end{bmatrix}.
$$

### 对 $y_4$

$$
y_4=-\phi_2=
\begin{bmatrix}
0 \\
-1
\end{bmatrix}.
$$

所以用两个模态可以 **精确重构** 全部中心化快照。

再加回均值，就得到原始快照：

$$
x_j=\bar x+y_j.
$$

例如

$$
x_1=
\begin{bmatrix}
1 \\
1
\end{bmatrix}
+
\begin{bmatrix}
1 \\
0
\end{bmatrix}=
\begin{bmatrix}
2 \\
1
\end{bmatrix},
\qquad
x_3=
\begin{bmatrix}
1 \\
1
\end{bmatrix}
+
\begin{bmatrix}
0 \\
1
\end{bmatrix}=
\begin{bmatrix}
1 \\
2
\end{bmatrix}.
$$

---

## 10. 只保留第一模态时会怎样

如果只取第一模态，近似变成

$$
\hat y_j^{(1)}=a_{1j}\phi_1.
$$

于是

$$
\hat y_1^{(1)}=
\sqrt{2}\,\phi_1=
\begin{bmatrix}
1 \\
0
\end{bmatrix},
\qquad
\hat y_2^{(1)}=
-\sqrt{2}\,\phi_1=
\begin{bmatrix}
-1 \\
0
\end{bmatrix},
$$

$$
\hat y_3^{(1)}=0,
\qquad
\hat y_4^{(1)}=0.
$$

加回均值后，对原始快照的一模态重构为

$$
\hat x_1^{(1)}=
\begin{bmatrix}
2 \\
1
\end{bmatrix},
\quad
\hat x_2^{(1)}=
\begin{bmatrix}
0 \\
1
\end{bmatrix},
\quad
\hat x_3^{(1)}=
\begin{bmatrix}
1 \\
1
\end{bmatrix},
\quad
\hat x_4^{(1)}=
\begin{bmatrix}
1 \\
1
\end{bmatrix}.
$$

也就是说，只保留第一模态时，模型完全保留了“第一个方向上的波动”，而把“第二个方向上的波动”全部抹掉了。

---

## 11. 能量占比和截断误差

本例中两个特征值是

$$
\lambda_1=4,
\qquad
\lambda_2=2.
$$

总能量为

$$
\lambda_1+\lambda_2=6.
$$

所以累计能量占比为

$$
E(1)=\frac{4}{6}=\frac{2}{3},
\qquad
E(2)=1.
$$

如果只保留第一模态，那么丢掉的能量就是

$$
\lambda_2=2.
$$

你也可以直接用误差验证这一点。对

$$
y_3=
\begin{bmatrix}
0 \\
1
\end{bmatrix},
$$

一模态近似是 $0$，所以误差为

$$
\|y_3\|_W^2=
\begin{bmatrix}
0 & 1
\end{bmatrix}
\begin{bmatrix}
2 & 0 \\
0 & 1
\end{bmatrix}
\begin{bmatrix}
0 \\
1
\end{bmatrix}
=1.
$$

同理，

$$
\|y_4\|_W^2=1.
$$

而 $y_1,y_2$ 被第一模态精确表示，误差为 0。

所以总截断误差是

$$
0+0+1+1=2=\lambda_2.
$$

---

## 12. 和 SVD 对起来

定义加权后的快照矩阵

$$
\widetilde Y=W^{1/2}Y.
$$

因为

$$
W^{1/2}=
\begin{bmatrix}
\sqrt{2} & 0 \\
0 & 1
\end{bmatrix},
$$

所以

$$
\widetilde Y=
\begin{bmatrix}
\sqrt{2} & -\sqrt{2} & 0 & 0 \\
0 & 0 & 1 & -1
\end{bmatrix}.
$$

于是

$$
\widetilde Y\widetilde Y^\top=
\begin{bmatrix}
4 & 0 \\
0 & 2
\end{bmatrix}.
$$

因此 $\widetilde Y$ 的奇异值是

$$
\sigma_1=2,
\qquad
\sigma_2=\sqrt{2},
$$

因为奇异值的平方等于特征值。

$\widetilde Y$ 的左奇异向量显然是

$$
u_1=
\begin{bmatrix}
1 \\
0
\end{bmatrix},
\qquad
u_2=
\begin{bmatrix}
0 \\
1
\end{bmatrix}.
$$

再变回原变量：

$$
\phi_i=W^{-1/2}u_i.
$$

因为

$$
W^{-1/2}=
\begin{bmatrix}
\frac{1}{\sqrt{2}} & 0 \\
0 & 1
\end{bmatrix},
$$

所以

$$
\phi_1=W^{-1/2}u_1=
\begin{bmatrix}
\frac{1}{\sqrt{2}} \\
0
\end{bmatrix},
\qquad
\phi_2=W^{-1/2}u_2=
\begin{bmatrix}
0 \\
1
\end{bmatrix}.
$$

这和前面直接手算出来的带权 POD 模态完全一致。

---

## 13. 这个例子里，你应该真正看懂什么

你可以把完整流程记成这 6 步：

1. **原始快照**：$x_j$
2. **减均值**：$y_j=x_j-\bar x$
3. **选内积**：$\langle u,v\rangle_W=u^\top Wv$
4. **求模态**：最大化带权投影能量
5. **算系数**：$a_{ij}=\langle y_j,\phi_i\rangle_W$
6. **重构原场**：

$$
x_j\approx \bar x+\sum_{i=1}^r a_{ij}\phi_i
$$

在这个具体例子里：

- 均值是

$$
\bar x=
\begin{bmatrix}
1 \\
1
\end{bmatrix}
$$

- 第一模态是

$$
\phi_1=
\begin{bmatrix}
\frac{1}{\sqrt{2}} \\
0
\end{bmatrix}
$$

- 第二模态是

$$
\phi_2=
\begin{bmatrix}
0 \\
1
\end{bmatrix}
$$

- 能量分别是

$$
\lambda_1=4,
\qquad
\lambda_2=2
$$

- 第一模态保留了

$$
66.7\%
$$

的加权能量。

---

## 14. 一句话解释这个例子的物理意义

因为我们给了第一个自由度更大的权重 $W_{11}=2$，所以 POD 会优先保留“第一个方向上的波动”。这就是为什么第一模态沿着第一个坐标方向，而不是像无权情形那样对两个方向一视同仁。这里的“最优”永远是相对于你选定的内积而言的。

---

## 15. 代码验证

```python
import numpy as np

np.set_printoptions(precision=6, suppress=True)

# ============================================================
# 1) 原始快照：每一列是一个快照
#    x1 = [2,1]^T, x2 = [0,1]^T, x3 = [1,2]^T, x4 = [1,0]^T
# ============================================================
X = np.array([
    [2, 0, 1, 1],
    [1, 1, 2, 0]
], dtype=float)

# ============================================================
# 2) 带权内积矩阵 W
#    <u,v>_W = u^T W v
# ============================================================
W = np.diag([2.0, 1.0])

# ============================================================
# 3) 先减均值
# ============================================================
x_bar = X.mean(axis=1, keepdims=True)   # 形状 (2,1)
Y = X - x_bar

print("=== 原始快照矩阵 X ===")
print(X)

print("\n=== 均值快照 x_bar ===")
print(x_bar)

print("\n=== 中心化后的快照矩阵 Y = X - x_bar ===")
print(Y)

# ============================================================
# 4) 方法一：直接解带权 POD 的特征值问题
#
#    对应我们手算里的：
#    (Y Y^T W) phi = lambda phi
#
#    注意：
#    - 这里不是普通正交，而是 W-正交
#    - 求出来的特征向量需要做 W-归一化
# ============================================================
M = Y @ Y.T @ W

eigvals, eigvecs = np.linalg.eig(M)

# 由于数值计算可能顺序乱掉，这里按特征值从大到小排序
idx = np.argsort(eigvals)[::-1]
eigvals = np.real(eigvals[idx])
eigvecs = np.real(eigvecs[:, idx])

# W-归一化：让 phi_i^T W phi_i = 1
Phi_from_eig = np.zeros_like(eigvecs)
for i in range(eigvecs.shape[1]):
    v = eigvecs[:, i:i+1]
    norm_w = np.sqrt((v.T @ W @ v).item())
    Phi_from_eig[:, i:i+1] = v / norm_w

print("\n=== 矩阵 M = Y Y^T W ===")
print(M)

print("\n=== 特征值 lambda ===")
print(eigvals)

print("\n=== 由特征值问题得到的 POD 模态（列向量）Phi_from_eig ===")
print(Phi_from_eig)

print("\n=== 检查 W-正交归一：Phi^T W Phi ===")
print(Phi_from_eig.T @ W @ Phi_from_eig)

# ============================================================
# 5) 计算模态系数
#
#    a_ij = <y_j, phi_i>_W = phi_i^T W y_j
#    矩阵写法：A = Phi^T W Y
# ============================================================
A_from_eig = Phi_from_eig.T @ W @ Y

print("\n=== 模态系数矩阵 A = Phi^T W Y ===")
print(A_from_eig)

# ============================================================
# 6) 用全部模态重构中心化快照，再加回均值
# ============================================================
Y_rec_full = Phi_from_eig @ A_from_eig
X_rec_full = x_bar + Y_rec_full

print("\n=== 用全部模态重构后的中心化快照 Y_rec_full ===")
print(Y_rec_full)

print("\n=== 用全部模态重构后的原始快照 X_rec_full ===")
print(X_rec_full)

# ============================================================
# 7) 只保留第一模态做截断重构
# ============================================================
Phi1 = Phi_from_eig[:, :1]
A1 = A_from_eig[:1, :]

Y_rec_r1 = Phi1 @ A1
X_rec_r1 = x_bar + Y_rec_r1

print("\n=== 只保留第1模态时的中心化重构 Y_rec_r1 ===")
print(Y_rec_r1)

print("\n=== 只保留第1模态时的原始快照重构 X_rec_r1 ===")
print(X_rec_r1)

# ============================================================
# 8) 计算总加权误差（验证手算结果）
#
#    误差定义：
#    sum_j || y_j - yhat_j ||_W^2
# ============================================================
def weighted_total_error(Y_true, Y_hat, W):
    err = 0.0
    for j in range(Y_true.shape[1]):
        d = (Y_true[:, j:j+1] - Y_hat[:, j:j+1])
        err += (d.T @ W @ d).item()
    return err

err_full = weighted_total_error(Y, Y_rec_full, W)
err_r1 = weighted_total_error(Y, Y_rec_r1, W)

print("\n=== 全模态重构的总加权误差 ===")
print(err_full)

print("\n=== 只保留第1模态时的总加权误差 ===")
print(err_r1)

# ============================================================
# 9) 方法二：用 SVD 再验证一遍
#
#    带权 POD 可转化为对 Y_tilde = W^(1/2) Y 做普通 SVD
#    若 Y_tilde = U Sigma V^T
#    则 Phi = W^(-1/2) U
# ============================================================
sqrtW = np.diag(np.sqrt(np.diag(W)))
inv_sqrtW = np.diag(1.0 / np.sqrt(np.diag(W)))

Y_tilde = sqrtW @ Y
U, s, Vt = np.linalg.svd(Y_tilde, full_matrices=False)

Phi_from_svd = inv_sqrtW @ U
A_from_svd = Phi_from_svd.T @ W @ Y

print("\n=== 加权后的快照矩阵 Y_tilde = W^(1/2) Y ===")
print(Y_tilde)

print("\n=== SVD 的左奇异向量 U ===")
print(U)

print("\n=== 奇异值 s ===")
print(s)

print("\n=== 由 SVD 恢复的 POD 模态 Phi_from_svd ===")
print(Phi_from_svd)

print("\n=== 检查 SVD 恢复出的模态是否 W-正交归一 ===")
print(Phi_from_svd.T @ W @ Phi_from_svd)

print("\n=== 由 SVD 得到的模态系数 A_from_svd ===")
print(A_from_svd)

# ============================================================
# 10) 能量与累计能量占比
#
#     特征值 lambda_i = s_i^2
# ============================================================
lambda_from_svd = s**2
energy_ratio = lambda_from_svd / np.sum(lambda_from_svd)
cum_energy_ratio = np.cumsum(energy_ratio)

print("\n=== 由 SVD 得到的特征值 lambda = s^2 ===")
print(lambda_from_svd)

print("\n=== 每个模态的能量占比 ===")
print(energy_ratio)

print("\n=== 累计能量占比 ===")
print(cum_energy_ratio)
```

结果：

```python
'''
=== 原始快照矩阵 X ===
[[2. 0. 1. 1.]
 [1. 1. 2. 0.]]

=== 均值快照 x_bar ===
[[1.]
 [1.]]

=== 中心化后的快照矩阵 Y = X - x_bar ===
[[ 1. -1.  0.  0.]
 [ 0.  0.  1. -1.]]

=== 矩阵 M = Y Y^T W ===
[[4. 0.]
 [0. 2.]]

=== 特征值 lambda ===
[4. 2.]

=== 由特征值问题得到的 POD 模态（列向量）Phi_from_eig ===
[[0.707107 0.      ]
 [0.       1.      ]]

=== 检查 W-正交归一：Phi^T W Phi ===
[[1. 0.]
 [0. 1.]]

=== 模态系数矩阵 A = Phi^T W Y ===
[[ 1.414214 -1.414214  0.        0.      ]
 [ 0.        0.        1.       -1.      ]]

=== 用全部模态重构后的中心化快照 Y_rec_full ===
[[ 1. -1.  0.  0.]
 [ 0.  0.  1. -1.]]

=== 用全部模态重构后的原始快照 X_rec_full ===
[[2. 0. 1. 1.]
 [1. 1. 2. 0.]]

=== 只保留第1模态时的中心化重构 Y_rec_r1 ===
[[ 1. -1.  0.  0.]
 [ 0.  0.  0.  0.]]

=== 只保留第1模态时的原始快照重构 X_rec_r1 ===
[[2. 0. 1. 1.]
 [1. 1. 1. 1.]]

=== 全模态重构的总加权误差 ===
1.9721522630525295e-31

=== 只保留第1模态时的总加权误差 ===
2.0

=== 加权后的快照矩阵 Y_tilde = W^(1/2) Y ===
[[ 1.414214 -1.414214  0.        0.      ]
 [ 0.        0.        1.       -1.      ]]

=== SVD 的左奇异向量 U ===
[[1. 0.]
 [0. 1.]]

=== 奇异值 s ===
[2.       1.414214]

=== 由 SVD 恢复的 POD 模态 Phi_from_svd ===
[[0.707107 0.      ]
 [0.       1.      ]]

=== 检查 SVD 恢复出的模态是否 W-正交归一 ===
[[1. 0.]
 [0. 1.]]

=== 由 SVD 得到的模态系数 A_from_svd ===
[[ 1.414214 -1.414214  0.        0.      ]
 [ 0.        0.        1.       -1.      ]]

=== 由 SVD 得到的特征值 lambda = s^2 ===
[4. 2.]

=== 每个模态的能量占比 ===
[0.666667 0.333333]

=== 累计能量占比 ===
[0.666667 1.      ]
'''
```
