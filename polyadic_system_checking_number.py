import numpy
import time
import tracemalloc

def opss_conversion(bases, residues):
    """
    Перевод числа из модулярной системы в десятичную с помощью ОПСС
    """
    bases_count = len(bases)
    full_diapazon = numpy.prod(bases.astype(object))
    
    print("\n1. Вычисление полного диапазона:")
    print(f"Полный диапазон P = {full_diapazon}")
    
    # 2. Поиск констант c_ij (обратных элементов)
    print("\n2. Поиск констант c_ij (обратных элементов):")
    c = numpy.zeros((bases_count-1, bases_count), dtype=object)
    
    for j in range(bases_count-1):
        for i in range(j+1):
            try:
                c[j][i] = pow(int(bases[i]), -1, int(bases[j+1]))
                print(f'c[{j}][{i}] = {c[j][i]} (поскольку {c[j][i]} * {bases[i]} ≡ 1 mod {bases[j+1]})')
            except ValueError:
                print(f'Невозможно найти c[{j}][{i}] - обратный элемент для {bases[i]} mod {bases[j+1]}')
                return None, None
    
    print("\nМатрица констант c_ij:")
    print(c[:,:bases_count-1])
    
    # 3. Расчет числа в ОПСС
    print("\n3. Расчет числа в ОПСС:")
    opss = numpy.zeros(bases_count, dtype=object)
    
    opss[0] = residues[0] % bases[0]
    print(f'opss[0] = {residues[0]} mod {bases[0]} = {opss[0]}')
    
    if bases_count > 1:
        temp = (residues[1] - opss[0]) * c[0][0]
        opss[1] = temp % bases[1]
        print(f'opss[1] = ({residues[1]} - {opss[0]}) * {c[0][0]} = {temp} mod {bases[1]} = {opss[1]}')
    
    for i in range(2, bases_count):
        # Первый шаг вычисления
        opss[i] = (residues[i] - opss[0]) * c[i-1][0]
        print(f'\nНачальное opss[{i}] = ({residues[i]} - {opss[0]}) * {c[i-1][0]} = {opss[i]}')
        
        # Коррекция с использованием предыдущих значений
        for j in range(1, i):
            opss[i] = (opss[i] - opss[j]) * c[i-1][j]
            print(f'Шаг {j}: opss[{i}] = (пред.результат - {opss[j]}) * {c[i-1][j]} = {opss[i]}')
        
        opss[i] %= bases[i]
        print(f'Финальное opss[{i}] = {opss[i]} mod {bases[i]}')
    
    print(f'\nЧисло в ОПСС: {opss}')
    
    # 4. Перевод из ОПСС в десятичную систему
    print("\n4. Перевод из ОПСС в десятичную систему:")
    decimal_number = 0
    for i in range(bases_count):
        product = numpy.prod(bases[:i].astype(object)) if i > 0 else 1
        term = opss[i] * product
        print(f'{opss[i]} * {product} = {term}')
        decimal_number += term
    
    print(f'\nИтоговое десятичное число: {decimal_number}')
    
    return opss, decimal_number

def check_with_orthogonal_bases(bases, residues):
    """
    Проверка методом ортогональных базисов (алгоритм из скрипта orthogonal_bases)
    """
    print("\n1. Исходные данные для проверки:")
    print(f"Основания системы: {bases}")
    print(f"Остатки числа: {residues}")
    
    # Вычисление полного диапазона
    print("\n2. Вычисление полного диапазона:")
    full_diapazon = numpy.prod(bases.astype(object))
    print(f"Полный диапазон P (произведение всех оснований) = {full_diapazon}")
    
    # Вычисление ортогональных базисов
    print("\n3. Вычисление ортогональных базисов:")
    bases_count = len(bases)
    
    # 3.1 Поиск величины P/pi = Pi для каждого основания
    print("\n3.1 Поиск величины P/pi = Pi для каждого основания:")
    P = numpy.zeros(bases_count, dtype=object)
    for i in range(bases_count):
        P[i] = full_diapazon // bases[i]
        print(f"P[{i}] = {full_diapazon} // {bases[i]} = {P[i]}")
    print(f"Величина P = {P}")
    
    # 3.2 Поиск величины βi = Pi(mod pi)
    print("\n3.2 Поиск величины βi = Pi(mod pi):")
    beta = numpy.zeros(bases_count, dtype=object)
    for i in range(bases_count):
        beta[i] = P[i] % bases[i]
        print(f"β[{i}] = {P[i]} mod {bases[i]} = {beta[i]}")
    print(f"Величина β = {beta}")
    
    # 3.3 Нахождение веса mi базисов (обратный элемент βi по модулю pi)
    print("\n3.3 Нахождение весов mi базисов (обратный элемент βi по модулю pi):")
    m = numpy.zeros(bases_count, dtype=object)
    for i in range(bases_count):
        try:
            m[i] = pow(int(beta[i]), -1, int(bases[i]))
            print(f"m[{i}] = inv({beta[i]}) mod {bases[i]} = {m[i]}")
        except ValueError:
            print(f"Ошибка: невозможно найти обратный элемент для beta[{i}] = {beta[i]} по модулю {bases[i]}")
            return None
    print(f"Веса базисов m = {m}")
    
    # 3.4 Вычисление ортогональных базисов системы Bi = mi * Pi
    print("\n3.4 Вычисление ортогональных базисов системы Bi = mi * Pi:")
    B = numpy.zeros(bases_count, dtype=object)
    for i in range(bases_count):
        B[i] = m[i] * P[i]
        print(f"B[{i}] = {m[i]} * {P[i]} = {B[i]}")
    print(f"Ортогональные базисы системы B = {B}")
    
    # 3.5 Нахождение числа A в десятичной системе
    print("\n3.5 Нахождение числа A в десятичной системе:")
    A = 0
    for i in range(bases_count):
        old_A = A
        A += residues[i] * B[i]
        print(f"A[{i}] = {old_A} + {residues[i]} * {B[i]} = {A}")
    
    A = A % full_diapazon
    print(f"\nA = {A} (после взятия по модулю {full_diapazon})")
    
    return A

# Начало программы
start_time = time.time()
tracemalloc.start()

print("Введите основания системы (все основания должны быть попарно взаимно простыми):")
bases_count = int(input('Введите количество оснований: '))
bases = numpy.zeros(bases_count, dtype=object)

for i in range(bases_count):
    bases[i] = int(input(f'Введите {i+1} основание: '))

print('\nВведите остатки числа по каждому основанию:')
residues = numpy.zeros(bases_count, dtype=object)
for i in range(bases_count):
    residues[i] = int(input(f'Введите остаток по основанию {bases[i]}: '))

full_diapazon = numpy.prod(bases.astype(object))

print('\nВведенные данные:')
print(f'Основания: {bases}')
print(f'Остатки числа: {residues}')
print(f'Полный диапазон P = {full_diapazon}')

# Перевод числа из СОК в десятичную систему через ОПСС
print('Перевод числа из модулярной системы счисления в позиционную (десятичную) методом перевода в ОПСС')
opss, decimal_number = opss_conversion(bases, residues)

if opss is not None:
    # Проверка методом ортогональных базисов
    print('Проверка методом ортогональных базисов')
    A_orth = check_with_orthogonal_bases(bases, residues)
    
    # Сравнение результатов
    if A_orth is not None:
        print('Сравнение результатов')
        print(f'Число из ОПСС: {decimal_number}')
        print(f'Число из ортогональных базисов: {A_orth}')
        
        if decimal_number == A_orth:
            print('Результаты совпадают!')
        else:
            print('Результаты не совпадают!')

# Замер памяти и времени
current, peak = tracemalloc.get_traced_memory()
tracemalloc.stop()
execution_time = time.time() - start_time

print(f'Использовано памяти: {current/1024:.2f} KB')
print(f'Пиковое использование памяти: {peak/1024:.2f} KB')
print(f'Время выполнения: {execution_time * 1e6:.2f} микросекунд')
print(f'Время выполнения: {execution_time:.6f} секунд')
