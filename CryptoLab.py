#!/usr/bin/env python
# coding: utf-8




import math
import numpy as np
W = np.array([[1, 1], [1, -1]],int)
p = (15, 11, 4, 1, 3, 15, 5, 9, 0, 10, 14, 7, 6, 8, 2, 15)





# Матрица Адамара-Сильвестра
def hadamard_matrix(dimension):
    if dimension == 1:
        return W
    else:
        return np.hstack((np.vstack((hadamard_matrix(dimension-1),hadamard_matrix(dimension-1))),np.vstack((hadamard_matrix(dimension-1),-hadamard_matrix(dimension-1)))))





# Получение массива значений функций
def function_value(function):
    values = [ [] for _ in range(len(bin(max(function))[2:]))]
    for value in function:
        for function_number in range(len(bin(max(function))[2:])):
            values[function_number].append(value >> function_number & 1)
    values.reverse()
    return values





# Генерирует все мономы в лексикографическом порядке
def generate_monoms(n):
    monoms = ['1']
    x = 1
    for i in range(2 ** n - 1):
        mon = ''
        for j in range(n):
            if x >> j & 1:
                mon += 'x' + str(4 - j)
        monoms.append(mon)
        x += 1
    return(monoms)





# Вывод значений функции в нормальном виде
def builds_function(function):
    answer = ''
    num_of_bits = len(bin(max(function))[2:])
    monoms = generate_monoms(num_of_bits)
    print_yi = "y{}"
    counter = 0
    for i in range(num_of_bits):
        answer += print_yi.format(counter) + ' = '
        coefficients = function_value(anf(function))[i]
        for j in range(len(function)):
            if coefficients[j]:
                answer += monoms[j] + ' ⊕ '
        counter += 1
        answer = answer[0:len(answer) - 3:]
        answer += '\n'
    return(answer)





# Многочлен Жегалкина
def anf(function):
    buff = []
    iterations = len(function)
    ANF = [ function[0] ]
    for _ in range(iterations-1):
        for i in range(len(function)-1):
            buff.append(function[i]^function[i+1])
        ANF.append(buff[0])
        function = buff
        buff = []
    return ANF





# Коэффцициенты Фурье
def fourier_spectrum(function):
    values = function_value(function)
    spectrum = [] 
    M = hadamard_matrix(int(math.log2(len(function))))
    for functions in range(len(bin(max(function))[2:])):
        spectrum.append(np.matmul(np.array(values[functions]).T,M).tolist())
    return spectrum





# Коэффициенты Уолша-Адамара
def walsch_spectrum(function):
    values = function_value(function)
    unit_array = np.array([1 for _ in range(len(function))])
    M = hadamard_matrix(int(math.log2(len(function))))
    spectrum = []
    for functions in range(len(bin(max(function))[2:])):
        spectrum.append(np.matmul((unit_array - 2 * np.array(values[functions])).T,M).tolist())
    return spectrum





# Задание 5 - строгий лавинный критерий
def strict_avalanche_test(function):
    values = function_value(function)
    var = [[0] for _ in range(int(math.log2(len(function))))]
    for k in range(int(math.log2(len(function)))): 
        for i in range(int(math.log2(len(function)))):
            num = 0
            for j in range(len(function)):
                 if values[k][j] != values[k][j^(2**i)]:
                    num=num+1
            var[k][0]=var[k][0]+num
    return var





# Задание 6
def bit_independence_criterion(function):
    values = function_value(function)
    g = {1:[0,1],2:[0,2],3:[1,2],4:[0,3],5:[1,3],6:[2,3]}
    var = [[0] for _ in range(6)]
    for k in g: 
        for i in range(int(math.log2(len(function)))):
            num = 0
            for j in range(len(function)):
                 if values[g[k][0]][j] != values[g[k][0]][j^(2**i)] and values[g[k][1]][j] != values[g[k][1]][j^(2**i)]:
                    num=num+1
            var[k-1][0]=var[k-1][0]+num
    return var,g





# Задание 7 - лавинный критерий
def avalanche_test(function):
    values = function_value(function)
    var = [[0] for _ in range(5)]
    for i in range(int(math.log2(len(function)))):
        for j in range(len(function)):
            num = 0
            for k in range(int(math.log2(len(function)))):
                if values[k][j] != values[k][j^(2**i)]:
                    num=num+1
            var[num][0]=var[num][0]+1
    return var





print("Значения разрядных функций:")
print(np.array(function_value(p)), end="\n\n")
print("Коэффициенты Жегалкина разрядных функций:")
print(np.array(function_value(anf(p))), end="\n\n")
print("Аналитический вид разрядных функций:")
print(builds_function(p), end="\n\n")
print("Коэффициенты Фурье разряжных функций:")
print(np.array(fourier_spectrum(p)), end="\n\n")
print("Коэффициенты Уолша-Адамара разрядных функций:")
print(np.array(walsch_spectrum(p)), end="\n\n")
print("Строгий лавинный критерий:")
print(np.array(strict_avalanche_test(p)), end="\n\n")
print("Bit independence criterion:")
print(np.array(bit_independence_criterion(p)), end="\n\n")
print("Лавинный критерий:")
print(np.array(avalanche_test(p)), end="\n\n")







