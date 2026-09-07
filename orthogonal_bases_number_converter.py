import numpy
import time
import tracemalloc

def modular_to_decimal(bases, residues):
    # Подсчет затрат памяти и времени
    start_time = time.time()
    tracemalloc.start()
    
    print("\n1. Исходные данные:")
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
    
    # Замер памяти и времени выполнения
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    execution_time = time.time() - start_time
    print("\nМетрики выполнения:")
    print(f"Пиковое использование памяти: {peak / 1024:.2f} KB")
    print(f"Время выполнения: {execution_time:.6f} секунд")
    
    return A

# Ввод данных
print("Введите основания системы (все основания должны быть попарно взаимно простыми):")
bases_count = int(input('Введите количество оснований: '))
bases = numpy.zeros(bases_count, dtype=object)
for i in range(bases_count):
    bases[i] = int(input(f'Введите {i+1} основание: '))

# Ввод остатков числа
print("\nВведите остатки числа для каждого основания:")
residues = numpy.zeros(bases_count, dtype=object)
for i in range(bases_count):
    residues[i] = int(input(f'Введите {i+1} остаток: '))

# Перевод числа из модулярной системы в десятичную
A = modular_to_decimal(bases, residues)
if A is not None:
    print("\nИтоговый результат:")
    print(f"Десятичное число A = {A}")
