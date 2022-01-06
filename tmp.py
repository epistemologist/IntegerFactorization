def pow(x,n):
    y = 1
    while n > 1:
        if n%2 == 0:
            x = x*x
            n //= 2
        else:
            y = x*y
            x = x*x
            n = (n-1) // 2
    return x*y
