# 手算 POD

> 参考资料：
> 
> [康斯坦茨大学数学与统计系 POD 讲义](https://www.math.uni-konstanz.de/numerik/personen/volkwein/teaching/POD-Book.pdf)
> 
> [Stanford University CME345 课程讲义](https://web.stanford.edu/group/frg/course_work/CME345/CA-AA216-CME345-Ch4.pdf)

## 1. 先造一个最小快照矩阵

设我们有两个快照，每个快照都是三维向量：

$$
y_1=
\begin{bmatrix}
2 \\
0 \\
0
\end{bmatrix},
\qquad
y_2=
\begin{bmatrix}
0 \\
1 \\
1
\end{bmatrix}.
$$

把它们排成快照矩阵：

$$
Y=[y_1,y_2]=
\begin{bmatrix}
2 & 0 \\
0 & 1 \\
0 & 1
\end{bmatrix}.
$$

这里：

- 行数 $m=3$，表示每个快照有 3 个自由度；
- 列数 $n=2$，表示一共只有 2 个快照。

为了先专注在 POD 主体上，这里先**不减均值**，也先使用最简单的 Euclidean 内积：

$$
\langle a,b\rangle = a^\top b.
$$

---

## 2. POD 的问题到底在求什么

先只求**第一模态** $\psi\in\mathbb{R}^3$，并要求它是单位向量：

$$
\|\psi\|=1.
$$

POD 要找的是：让所有快照在这个方向上的总投影能量最大的方向，也就是求

$$
\max_{\|\psi\|=1}\sum_{j=1}^2 \langle y_j,\psi\rangle^2.
$$

把这两个快照代进去，就是

$$
\max_{\|\psi\|=1}
\left(
\langle y_1,\psi\rangle^2+\langle y_2,\psi\rangle^2
\right).
$$

---

## 3. 先把目标函数写开

设

$$
\psi=
\begin{bmatrix}
a \\
b \\
c
\end{bmatrix},
\qquad
a^2+b^2+c^2=1.
$$

那么

$$
\langle y_1,\psi\rangle
=
\begin{bmatrix}
2 & 0 & 0
\end{bmatrix}
\begin{bmatrix}
a \\
b \\
c
\end{bmatrix}
=2a,
$$

$$
\langle y_2,\psi\rangle
=
\begin{bmatrix}
0 & 1 & 1
\end{bmatrix}
\begin{bmatrix}
a \\
b \\
c
\end{bmatrix}
=b+c.
$$

所以目标函数变成

$$
J(a,b,c)=(2a)^2+(b+c)^2=4a^2+(b+c)^2.
$$

约束是

$$
a^2+b^2+c^2=1.
$$

所以第一模态就是解下面这个约束优化问题：

$$
\max_{a^2+b^2+c^2=1}\left[4a^2+(b+c)^2\right].
$$

---

## 4. 直接看出最大值在哪

因为

$$
(b+c)^2\le 2(b^2+c^2),
$$

所以

$$
J(a,b,c)=4a^2+(b+c)^2
\le 4a^2+2(b^2+c^2).
$$

再用约束 $a^2+b^2+c^2=1$，得到

$$
J(a,b,c)\le 4a^2+2(1-a^2)=2+2a^2\le 4.
$$

要达到最大值 $4$，必须同时满足：

1. $a^2=1$；
2. $b=c=0$。

所以最大值在

$$
\psi_1=
\begin{bmatrix}
1 \\
0 \\
0
\end{bmatrix}
\quad\text{或}\quad
\begin{bmatrix}
-1 \\
0 \\
0
\end{bmatrix}
$$

处取得。模态的正负号等价，这里取

$$
\boxed{
\psi_1=
\begin{bmatrix}
1 \\
0 \\
0
\end{bmatrix}
}
$$

这就是**第一 POD 模态**。

它对应的最大投影能量是

$$
\lambda_1=4.
$$

---

## 5. 第二模态怎么求

第二模态要满足两个条件：

1. $\|\psi_2\|=1$；
2. $\psi_2\perp\psi_1$。

因为

$$
\psi_1=
\begin{bmatrix}
1 \\
0 \\
0
\end{bmatrix},
$$

所以正交条件说明第二模态一定形如

$$
\psi_2=
\begin{bmatrix}
0 \\
b \\
c
\end{bmatrix},
\qquad
b^2+c^2=1.
$$

这时目标函数变成

$$
J(0,b,c)=4\cdot 0^2+(b+c)^2=(b+c)^2.
$$

又因为

$$
(b+c)^2\le 2(b^2+c^2)=2,
$$

等号成立当且仅当 $b=c$。

再结合 $b^2+c^2=1$，得到

$$
b=c=\frac{1}{\sqrt{2}}.
$$

所以第二模态是

$$
\boxed{
\psi_2=
\frac{1}{\sqrt{2}}
\begin{bmatrix}
0 \\
1 \\
1
\end{bmatrix}
}
$$

对应能量是

$$
\lambda_2=2.
$$

---

## 6. 第三模态也能顺手得到

因为空间是三维的，前两个模态已经正交归一，所以第三模态就是同时与它们正交的单位向量。

设

$$
\psi_3=
\begin{bmatrix}
0 \\
x \\
y
\end{bmatrix},
$$

并要求与 $\psi_2$ 正交：

$$
\frac{1}{\sqrt{2}}(x+y)=0
\quad\Rightarrow\quad
x+y=0.
$$

再加上单位长度条件 $x^2+y^2=1$，可取

$$
x=\frac{1}{\sqrt{2}},
\qquad
y=-\frac{1}{\sqrt{2}}.
$$

因此

$$
\boxed{
\psi_3=
\frac{1}{\sqrt{2}}
\begin{bmatrix}
0 \\
1 \\
-1
\end{bmatrix}
}
$$

对应能量是

$$
\lambda_3=0.
$$

因为所有快照都落在 $\mathrm{span}(y_1,y_2)$ 张成的二维子空间里，所以第三个方向上没有任何能量。

---

## 7. 用矩阵方法重做一遍

上面是直接从优化问题手算出来的。现在用 POD 的标准矩阵形式再做一次。

POD 的特征值问题是

$$
YY^\top\psi=\lambda\psi.
$$

先算

$$
YY^\top=
\begin{bmatrix}
2 & 0 \\
0 & 1 \\
0 & 1
\end{bmatrix}
\begin{bmatrix}
2 & 0 & 0 \\
0 & 1 & 1
\end{bmatrix}
=
\begin{bmatrix}
4 & 0 & 0 \\
0 & 1 & 1 \\
0 & 1 & 1
\end{bmatrix}.
$$

所以我们要解

$$
\begin{bmatrix}
4 & 0 & 0 \\
0 & 1 & 1 \\
0 & 1 & 1
\end{bmatrix}
\psi
=
\lambda\psi.
$$

现在可以直接验证：

### 第一模态

$$
\psi_1=
\begin{bmatrix}
1 \\
0 \\
0
\end{bmatrix}
$$

时，

$$
YY^\top\psi_1=
\begin{bmatrix}
4 \\
0 \\
0
\end{bmatrix}
=4\psi_1,
$$

所以 $\lambda_1=4$。

### 第二模态

$$
\psi_2=
\frac{1}{\sqrt{2}}
\begin{bmatrix}
0 \\
1 \\
1
\end{bmatrix}
$$

时，

$$
YY^\top\psi_2=
\frac{1}{\sqrt{2}}
\begin{bmatrix}
0 \\
2 \\
2
\end{bmatrix}
=2\psi_2,
$$

所以 $\lambda_2=2$。

### 第三模态

$$
\psi_3=
\frac{1}{\sqrt{2}}
\begin{bmatrix}
0 \\
1 \\
-1
\end{bmatrix}
$$

时，

$$
YY^\top\psi_3=
\frac{1}{\sqrt{2}}
\begin{bmatrix}
0 \\
0 \\
0
\end{bmatrix}
=0\psi_3,
$$

所以 $\lambda_3=0$。

这和上面从优化问题直接算出来的结果完全一致。

---

## 8. 现在把它和 SVD 连起来

现在看最关键的一步。

因为

$$
Y=U\Sigma V^\top,
$$

所以

$$
YY^\top=U\Sigma^2U^\top.
$$

这说明：

- $YY^\top$ 的特征向量，就是 $Y$ 的左奇异向量；
- 特征值，就是奇异值的平方。

我们已经求出了特征值：

$$
\lambda_1=4,\qquad
\lambda_2=2,\qquad
\lambda_3=0.
$$

所以奇异值就是

$$
\sigma_1=\sqrt{4}=2,\qquad
\sigma_2=\sqrt{2},\qquad
\sigma_3=0.
$$

因为 $Y$ 是一个 $3\times 2$ 矩阵，真正非零的奇异值只有两个：

$$
\boxed{
\sigma_1=2,\qquad
\sigma_2=\sqrt{2}
}
$$

对应的左奇异向量就是我们刚刚得到的两个 POD 模态：

$$
u_1=\psi_1=
\begin{bmatrix}
1 \\
0 \\
0
\end{bmatrix},
\qquad
u_2=\psi_2=
\frac{1}{\sqrt{2}}
\begin{bmatrix}
0 \\
1 \\
1
\end{bmatrix}.
$$

所以这一步一定要记住：

$$
\boxed{
\text{POD 模态就是快照矩阵 } Y \text{ 的左奇异向量}
}
$$

---

## 9. 再把右奇异向量也手算出来

为了把 SVD 彻底打通，我们再求 $V$。

先算

$$
Y^\top Y=
\begin{bmatrix}
2 & 0 & 0 \\
0 & 1 & 1
\end{bmatrix}
\begin{bmatrix}
2 & 0 \\
0 & 1 \\
0 & 1
\end{bmatrix}
=
\begin{bmatrix}
4 & 0 \\
0 & 2
\end{bmatrix}.
$$

所以 $Y^\top Y$ 的特征值也是 $4,2$，对应的单位特征向量很简单：

$$
v_1=
\begin{bmatrix}
1 \\
0
\end{bmatrix},
\qquad
v_2=
\begin{bmatrix}
0 \\
1
\end{bmatrix}.
$$

因此

$$
V=
\begin{bmatrix}
1 & 0 \\
0 & 1
\end{bmatrix}
=I.
$$

于是这个例子的薄 SVD 就是

$$
Y=U\Sigma V^\top,
$$

其中

$$
U=
\begin{bmatrix}
1 & 0 \\
0 & \frac{1}{\sqrt{2}} \\
0 & \frac{1}{\sqrt{2}}
\end{bmatrix},
\qquad
\Sigma=
\begin{bmatrix}
2 & 0 \\
0 & \sqrt{2}
\end{bmatrix},
\qquad
V=I.
$$

乘开验证：

$$
U\Sigma=
\begin{bmatrix}
2 & 0 \\
0 & 1 \\
0 & 1
\end{bmatrix}
=Y.
$$

完全正确。

---

## 10. 模态系数怎么手算

POD 展开里，每个快照都写成

$$
y_j=\sum_{i=1}^r a_{ij}\psi_i.
$$

系数由投影给出：

$$
a_{ij}=\langle y_j,\psi_i\rangle.
$$

下面来逐项计算。

### 对快照 $y_1=[2,0,0]^\top$

第一模态系数：

$$
a_{11}
=
\langle y_1,\psi_1\rangle
=
\begin{bmatrix}
2 & 0 & 0
\end{bmatrix}
\begin{bmatrix}
1 \\
0 \\
0
\end{bmatrix}
=2.
$$

第二模态系数：

$$
a_{21}
=
\langle y_1,\psi_2\rangle
=
\begin{bmatrix}
2 & 0 & 0
\end{bmatrix}
\frac{1}{\sqrt{2}}
\begin{bmatrix}
0 \\
1 \\
1
\end{bmatrix}
=0.
$$

所以

$$
y_1=2\psi_1+0\psi_2.
$$

### 对快照 $y_2=[0,1,1]^\top$

第一模态系数：

$$
a_{12}
=
\langle y_2,\psi_1\rangle
=
\begin{bmatrix}
0 & 1 & 1
\end{bmatrix}
\begin{bmatrix}
1 \\
0 \\
0
\end{bmatrix}
=0.
$$

第二模态系数：

$$
a_{22}
=
\langle y_2,\psi_2\rangle
=
\begin{bmatrix}
0 & 1 & 1
\end{bmatrix}
\frac{1}{\sqrt{2}}
\begin{bmatrix}
0 \\
1 \\
1
\end{bmatrix}
=
\frac{2}{\sqrt{2}}
=
\sqrt{2}.
$$

所以

$$
y_2=0\psi_1+\sqrt{2}\,\psi_2.
$$

于是系数矩阵就是

$$
A=
\begin{bmatrix}
2 & 0 \\
0 & \sqrt{2}
\end{bmatrix}.
$$

你会发现它正好就是

$$
\Sigma V^\top=
\begin{bmatrix}
2 & 0 \\
0 & \sqrt{2}
\end{bmatrix}.
$$

这不是巧合，而是一般结论：

$$
\boxed{
A=\Sigma V^\top
}
$$

---

## 11. 看看只保留一个模态时会发生什么

现在做截断 POD，只保留第一模态 $\psi_1$。

### 对 $y_1$

一模态重构：

$$
\hat y_1
=
\langle y_1,\psi_1\rangle\psi_1
=
2\psi_1
=
\begin{bmatrix}
2 \\
0 \\
0
\end{bmatrix}.
$$

它被完全重构。

### 对 $y_2$

一模态重构：

$$
\hat y_2
=
\langle y_2,\psi_1\rangle\psi_1
=0.
$$

所以第二个快照完全丢失了。

---

## 12. 为什么第一模态能量占比是 $4/(4+2)=2/3$

POD 中每个模态的能量由特征值给出：

$$
\lambda_1=4,\qquad
\lambda_2=2.
$$

总能量为

$$
\lambda_1+\lambda_2=6.
$$

所以第一模态能量占比为

$$
\frac{4}{6}=\frac{2}{3}.
$$

第二模态占比为

$$
\frac{2}{6}=\frac{1}{3}.
$$

这表示：

- 第一模态解释了全部快照波动中的 $66.7\%$；
- 前两个模态解释了 $100\%$。

因为这个例子本来就是秩 2。

---

## 13. 为什么这就是“最优低秩逼近”

如果只允许用一个单位向量 $\psi$ 来逼近所有快照，那么总重构误差是

$$
\sum_{j=1}^2 \|y_j-\langle y_j,\psi\rangle\psi\|^2.
$$

POD 选出的 $\psi_1$ 会使这个量最小。

在这个例子里，若取

$$
\psi_1=
\begin{bmatrix}
1 \\
0 \\
0
\end{bmatrix},
$$

那么：

- 对 $y_1$，误差为 0；
- 对 $y_2$，误差为

$$
\|y_2\|^2=0^2+1^2+1^2=2.
$$

总误差为

$$
0+2=2.
$$

这恰好等于被丢弃的特征值之和：

$$
\lambda_2=2.
$$

这就是一般公式在这个例子里的具体体现：

$$
\boxed{
\text{截断到 } r \text{ 个模态后的总平方误差}
=
\sum_{i=r+1}^{m}\lambda_i
}
$$

---

## 14. 如果先减均值，会改哪一步

这一步很重要，因为真实流体里经常要这么做。

当前两个快照的均值是

$$
\bar y=\frac{1}{2}(y_1+y_2)
=
\frac{1}{2}
\begin{bmatrix}
2 \\
1 \\
1
\end{bmatrix}
=
\begin{bmatrix}
1 \\
0.5 \\
0.5
\end{bmatrix}.
$$

中心化后的快照变成

$$
y_1'=y_1-\bar y=
\begin{bmatrix}
1 \\
-0.5 \\
-0.5
\end{bmatrix},
\qquad
y_2'=y_2-\bar y=
\begin{bmatrix}
-1 \\
0.5 \\
0.5
\end{bmatrix}.
$$

然后你就不再对

$$
Y=[y_1,y_2]
$$

做 POD，而是对

$$
Y'=[y_1',y_2']
$$

做 POD。

也就是说，**推导完全一样，唯一变化只是把快照矩阵从 $Y$ 换成中心化后的 $Y'$**。

---

## 15. 这一套推导的逻辑链

这一遍手算里，你已经看到了一条完整链条：

### 第一步：定义优化问题

找单位向量 $\psi$，使

$$
\sum_j \langle y_j,\psi\rangle^2
$$

最大。

### 第二步：转成矩阵形式

这个目标等于

$$
\psi^\top YY^\top\psi.
$$

### 第三步：解特征值问题

于是得到

$$
YY^\top\psi=\lambda\psi.
$$

### 第四步：按特征值从大到小取模态

得到 $\psi_1,\psi_2,\dots$。

### 第五步：和 SVD 对应

因为

$$
YY^\top=U\Sigma^2U^\top,
$$

所以 POD 模态就是左奇异向量。

---

## 16. 用一句话总结这个手算例子

对

$$
Y=
\begin{bmatrix}
2 & 0 \\
0 & 1 \\
0 & 1
\end{bmatrix}
$$

做 POD，得到

$$
\psi_1=
\begin{bmatrix}
1 \\
0 \\
0
\end{bmatrix},
\qquad
\lambda_1=4,
$$

$$
\psi_2=
\frac{1}{\sqrt{2}}
\begin{bmatrix}
0 \\
1 \\
1
\end{bmatrix},
\qquad
\lambda_2=2,
$$

$$
\psi_3=
\frac{1}{\sqrt{2}}
\begin{bmatrix}
0 \\
1 \\
-1
\end{bmatrix},
\qquad
\lambda_3=0.
$$

而这正好就是快照矩阵 $Y$ 的左奇异向量与奇异值平方之间的对应关系。
