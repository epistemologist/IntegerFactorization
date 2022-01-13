from utility import ModN

def add_h(A, B, C, X1, Z1, X2, Z2, Xminus, Zminus):
    if X1 == Z1 == 0:
        return (X2, Z2)
    if X2 == Z2 == 0:
        return (X1, Z1)
    tmp = X1 * X2 - A * Z1 * Z2
    Xplus = Zminus * (tmp * tmp - 4 * B * Z1 * Z2 * (X1 * Z2 + X2 * Z1 + C * Z1 * Z2))
    tmp2 = X1 * Z2 - X2 * Z1
    Zplus = Xminus * tmp2 * tmp2
    return (Xplus, Zplus)


def double_h(A, B, C, X1, Z1):
    tmp = X1 * X1 - A * Z1 * Z1
    Xplus = tmp * tmp - 4 * B * Z1 * Z1 * Z1 * (2 * X1 + C * Z1)
    Zplus = 4 * Z1 * (X1 * X1 * X1 + C * X1 * X1 * Z1 + A * X1 * Z1 * Z1 + B * Z1 * Z1 * Z1)
    return (Xplus, Zplus)


def point_multiplication(A, B, C, X, Z, n):
    if n == 0:
        return [0, 0]
    if n == 1:
        return [X, Z]
    if n == 2:
        return double_h(A, B, C, X, Z)
    U, V = X, Z
    T, W = double_h(A, B, C, X, Z)
    bit_length = int(n).bit_length()
    for j in range(bit_length - 2, 0, -1):
        if (n >> j) & 1:
            U, V = add_h(A, B, C, T, W, U, V, X, Z)
            T, W = double_h(A, B, C, T, W)
        else:
            T, W = add_h(A, B, C, U, V, T, W, X, Z)
            U, V = double_h(A, B, C, U, V)
    if n & 1:
        return add_h(A, B, C, U, V, T, W, X, Z)
    return double_h(A, B, C, U, V)

# y^2 = x^3 + 4x^2 + 2x + 13, P = (20 : 25 : 1)
n = 479
A = ModN(2, n)
B = ModN(13, n)
C = ModN(4, n)
X = ModN(17, n)
Z = ModN(1, n)

print(point_multiplication(A,B,C,X,Z,37))

