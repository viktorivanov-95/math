import numpy
import time
import tracemalloc
import math
from decimal import Decimal, getcontext, ROUND_DOWN, ROUND_HALF_UP

# Метод ортогональных базисов
def check_with_orthogonal_bases(bases, residues):
    print("\nМетод ортогональных базисов")
    bases_count = len(bases)
    full_diapazon = numpy.prod(bases.astype(object))
    print(f"Основания системы: {bases}")
    print(f"Остатки числа: {residues}")
    print(f"Полный диапазон (P): {full_diapazon}")
    
    # 1. Вычисление Pi = P/pi
    print("\n1. Вычисление Pi = P/pi для каждого основания")
    P = numpy.zeros(bases_count, dtype=object)
    for i in range(bases_count):
        P[i] = full_diapazon // bases[i]
        print(f"P[{i}] = {full_diapazon} // {bases[i]} = {P[i]}")
    
    # 2. Вычисление βi = Pi mod pi
    print("\n2. Вычисление βi = Pi mod pi")
    beta = numpy.zeros(bases_count, dtype=object)
    for i in range(bases_count):
        beta[i] = P[i] % bases[i]
        print(f"β[{i}] = {P[i]} mod {bases[i]} = {beta[i]}")
    
    # 3. Вычисление mi (обратные к βi по модулю pi)
    print("\n3. Вычисление весов mi (обратные к βi)")
    m = numpy.zeros(bases_count, dtype=object)
    for i in range(bases_count):
        m[i] = pow(int(beta[i]), -1, int(bases[i]))
        print(f"m[{i}] = inv({beta[i]}) mod {bases[i]} = {m[i]}")
    
    # 4. Вычисление ортогональных базисов Bi = mi * Pi
    print("\n4. Вычисление ортогональных базисов Bi")
    B = numpy.zeros(bases_count, dtype=object)
    for i in range(bases_count):
        B[i] = m[i] * P[i]
        print(f"B[{i}] = {m[i]} * {P[i]} = {B[i]}")
    
    # 5. Вычисление числа A
    print("\n5. Вычисление числа A")
    A = 0
    for i in range(bases_count):
        if i == 0:
            print(f"A[{i}] = 0 + {int(residues[i])} * {int(B[i])} = {int(residues[i]*B[i])}")
        else:
            print(f"A[{i}] = {int(A)} + {int(residues[i])} * {int(B[i])} = {int(A + residues[i]*B[i])}")
        A = int(A + residues[i] * B[i])
    A = numpy.mod(A, full_diapazon)
    print(f"\nA = {A} mod {full_diapazon} = {A}")
    
    return A

# Метод перевода числа в ОПСС
def check_with_opss(bases, residues):
    print("\nМетод перевода числа в ОПСС")
    bases_count = len(bases)
    print(f"Основания системы: {bases}")
    print(f"Остатки числа: {residues}")
    
    # 1. Вычисление констант c_ij (обратные элементы)
    print("\n1. Вычисление констант c_ij (обратные элементы)")
    c = numpy.zeros((bases_count, bases_count), dtype=object)
    for i in range(bases_count):
        for j in range(i + 1, bases_count):
            c[i][j] = pow(int(bases[i]), -1, int(bases[j]))
            print(f"c[{i}][{j}] = inv({bases[i]}) mod {bases[j]} = {c[i][j]}")
    
    # 2. Перевод в ОПСС
    print("\n2. Перевод числа в ОПСС")
    opss = numpy.zeros(bases_count, dtype=object)
    opss[0] = residues[0] % bases[0]
    print(f"opss[0] = {residues[0]} mod {bases[0]} = {opss[0]}")
    
    # Последующие разряды
    for i in range(1, bases_count):
        temp = residues[i]
        for j in range(i):
            temp = (temp - opss[j]) * c[j][i] % bases[i]
        opss[i] = temp % bases[i]
        print(f"opss[{i}] = {temp} mod {bases[i]} = {opss[i]}")
    print(f"\nЧисло в ОПСС: {opss}")
    
    # 3. Перевод в десятичную систему
    print("\n3. Перевод из ОПСС в десятичную систему")
    A_opss = 0
    product = 1
    for i in range(bases_count):
        term = opss[i] * product
        print(f"{opss[i]} * {product} = {term}")
        A_opss += term
        if i < bases_count - 1:
            product *= bases[i]
    print(f"\nИтоговое число: {A_opss}")
    return A_opss

# Метод относительных величин с поддержкой больших чисел
def relative_values_method(bases, residues, precision=10):
    print(f"\nМетод относительных величин (точность: {precision} знаков)")
    bases_count = len(bases)
    
    # Используем Decimal для вычисления диапазонов
    full_diapazon = Decimal(1)
    for base in bases:
        full_diapazon *= Decimal(str(base))
    
    print(f"Основания системы: {bases}")
    print(f"Остатки числа: {residues}")
    print(f"Полный диапазон (P): {full_diapazon}")
    
    # Настройка точности Decimal
    getcontext().prec = precision + 50
    getcontext().rounding = ROUND_HALF_UP
    
    # 1. Вычисление Pi = P/pi
    print("\n1. Вычисление Pi = P/pi для каждого основания")
    P = numpy.zeros(bases_count, dtype=object)
    for i in range(bases_count):
        base_decimal = Decimal(str(bases[i]))
        P_decimal = full_diapazon / base_decimal
        # Ограничиваем количество знаков после запятой
        P_decimal = P_decimal.quantize(Decimal(f'1.{"0" * precision}'), rounding=ROUND_HALF_UP)
        P[i] = P_decimal
        print(f"P[{i}] = {full_diapazon} / {base_decimal} = {P[i]}")
    
    # 2. Вычисление коэффициентов k[i]
    print("\n2. Вычисление коэффициентов k[i] = (Pi^(-1) mod pi) / pi")
    k = numpy.zeros(bases_count, dtype=object)
    for i in range(bases_count):
        # Вычисляем обратный элемент с целыми числами
        Pi_int = int(Decimal(P[i]).to_integral_value(rounding=ROUND_HALF_UP))
        base_int = int(bases[i])
        inv_Pi = pow(Pi_int, -1, base_int)
        # Вычисляем k[i] с ограниченной точностью
        k_decimal = Decimal(str(inv_Pi)) / Decimal(str(base_int))
        k_decimal = k_decimal.quantize(Decimal(f'1.{"0" * precision}'), rounding=ROUND_HALF_UP)
        k[i] = k_decimal
        print(f"k[{i}] = inv({Pi_int}) mod {base_int} / {base_int} = {k[i]}")
    
    # 3. Вычисление A/P
    print("\n3. Вычисление A/P = Σ(residues[i] * k[i])")
    A_relative = Decimal('0')
    for i in range(bases_count):
        term = Decimal(str(residues[i])) * k[i]
        term = term.quantize(Decimal(f'1.{"0" * precision}'), rounding=ROUND_HALF_UP)
        if i == 0:
            print(f"A/P[{i}] = 0 + {residues[i]} * {k[i]} = {term}")
        else:
            print(f"A/P[{i}] = {A_relative} + {residues[i]} * {k[i]} = {A_relative + term}")
        A_relative += term
        A_relative = A_relative.quantize(Decimal(f'1.{"0" * precision}'), rounding=ROUND_HALF_UP)
    
    # 4. Вычисление A
    print(f"\n4. Вычисление A = (A/P mod 1) * P")
    fractional_part = A_relative % Decimal('1')
    fractional_part = fractional_part.quantize(Decimal(f'1.{"0" * precision}'), rounding=ROUND_HALF_UP)
    print(f"Дробная часть A/P = {fractional_part}")
    A_relative = fractional_part * full_diapazon
    A_relative = A_relative.quantize(Decimal(f'1.{"0" * precision}'), rounding=ROUND_HALF_UP)
    print(f"A = {fractional_part} * {full_diapazon} = {A_relative}")
    
    return A_relative

# Начало выполнения программы
start_time = time.time()
tracemalloc.start()

print("Введите основания системы (все основания должны быть попарно взаимно простыми):")

# Ввод данных
bases_count = int(input('Введите количество оснований: '))
bases = numpy.zeros(bases_count, dtype=object)

for i in range(bases_count):
    bases[i] = int(input(f'Введите {i+1} основание: '))

residues = numpy.zeros(bases_count, dtype=object)
for i in range(bases_count):
    residues[i] = int(input(f'Введите остаток по основанию {bases[i]}: '))

# Настройка точности метода относительных величин
print("\nНастройка точности метода относительных величин")
while True:
    try:
        precision = int(input('Введите количество знаков после запятой для точности (1-50): '))
        if 1 <= precision <= 50:
            break
        else:
            print("Ошибка: введите число от 1 до 50")
    except ValueError:
        print("Ошибка: введите целое число")

# Вычисление полного диапазона
full_diapazon = numpy.prod(bases.astype(object))
full_diapazon_dec = Decimal(1)
for base in bases:
    full_diapazon_dec *= Decimal(str(base))

# Вывод введенных данных
print("\nВведенные данные")
print('Основания:', bases)
print('Остатки числа:', residues)
print('Полный диапазон P:', full_diapazon_dec)

# Метод относительных величин
A_relative = relative_values_method(bases, residues, precision)

# Метод ортогональных базисов
A_orth = check_with_orthogonal_bases(bases, residues)

# Метод ОПСС
A_opss = check_with_opss(bases, residues)

# Сравнение результатов
print(f"Метод относительных величин: A = {A_relative}")
print(f"Метод ортогональных базисов: A = {A_orth}")
print(f"Метод ОПСС: A = {A_opss}")

# Преобразуем точные результаты в Decimal для сравнения
A_orth_dec = Decimal(str(A_orth))
A_opss_dec = Decimal(str(A_opss))

print(f"\nРазница между относительными величинами и ортогональными базисами: {abs(A_relative - A_orth_dec)}")
print(f"Разница между относительными величинами и ОПСС: {abs(A_relative - A_opss_dec)}")
print(f"Разница между ортогональными базисами и ОПСС: {abs(A_orth_dec - A_opss_dec)}")

if A_orth == A_opss:
    print("\n Методы ортогональных базисов и ОПСС дали одинаковый результат")
else:
    print("\n Методы ортогональных базисов и ОПСС дали разные результаты!")

# Замер ресурсов
current, peak = tracemalloc.get_traced_memory()
tracemalloc.stop()
execution_time = time.time() - start_time

print(f'Пиковое использование памяти: {peak / 1024:.2f} KB')
print(f'Время выполнения: {execution_time:.6f} секунд')
print(f'Точность вычислений: {precision} знаков после запятой')
