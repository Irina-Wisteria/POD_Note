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

=== 由特征值问题得到的 POD 模态(列向量)Phi_from_eig ===
[[0.707107 0.      ]
 [0.       1.      ]]

=== 检查 W-正交归一: Phi^T W Phi ===
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