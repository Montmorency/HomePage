#x_1 = (Hu_0-a_0u_0)
#b_1 = <x_1|x_1>
#u_1 = x_1*(1.0/b_1)
#The earliest known use of a three-term recurrence was by the 
#Ancient Greeks to find the lowest common divisor of two 
#natural numbers (positive integers). Reorganizing the arithmetic 
#to emphasize its similarity to the recursion method, we replace 
#H by a real number, r, between nought and one. The natural numbers 
#U_1 < U_2 < U_3 < ... < U_n < ... replace the {u_{n} } and the recursion 
#proceeds by finding the smallest natural number, 
#An, such that, 
#\begin{equation}
#U_{n+1} = A_{n} U_{n} + U{n-1} 
#\end{equation}
#and $rU_{n+1}$  is no further from a natural number than $rU_n$. 
#In fact the sequence {rU_n} comes closer and closer to being a 
#sequence of natural numbers as n increases. If r is a 
#rational number, then the recursion terminates when $rU_{n}$
#is a natural number. There is a deep relationship 
#between the recurrence for the rational approximants to a real number 
#and the rational approximants to a resolvent.

a = 28
b = 21
#a = 35
#b = 22
a = 65
b = 26
#a = 1389
#b = 1003
u_nm1 = a
u_n = b
u_np1 = a
while (u_nm1 > u_n):
#find largest integral quotient this is an operation
#equivalent to a hamiltonian that is an integral operator.
    q = 0
    #while (u_nm1 > u_n): 
    while (u_nm1-q*u_n)>0: 
    #    u_nm1 -= u_n
        u_np1 = u_nm1- q*u_n
        q += 1
#with parallel assignment
    #u_nm1, u_n = u_n, u_nm1
    u_nm1, u_n = u_n, u_np1
#with a temp variable
    print u_nm1, u_n, q

